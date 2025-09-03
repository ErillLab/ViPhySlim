import preprocessing
from Bio import SeqIO, Phylo
from Bio.Phylo.Consensus import *
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
import matplotlib.pyplot as plt
from mpi4py import MPI
import numpy as np
import shutil
import csv
import os
import time
import random

RUN_MODE = preprocessing.execution_type

# Initialize MPI communication. Each process gets a unique rank (ID) and 
# the total number of processes (size). It is used to distribute pairwise 
# tasks.
if RUN_MODE == 'parallel':
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
else:
    comm = None
    rank = None
    size = 1 


'''
    Returns True if serial run (rank is None) or if parallel run and rank is 0.
'''
def i_am_main_process():
    return not rank


"""
    Calculates the distance between two genomes using a selected distance formula.
    It takes as input:
    - distance_formula: the formula to use ("d0", "d4", or "d6"),
    - h_total: the total number of covered positions across both genomes,
    - i_xy: the average sequence identity between the two genomes,
    - lambda_xy: the total length of both genomes combined.

    Returns a normalized distance value between 0 and 1.
"""
def distance_calculator(distance_formula, h_total, i_xy, lambda_xy):
    if distance_formula == "d0":
        distance = 1 - (h_total / lambda_xy)
    elif distance_formula == "d4":
        distance = 1 - ((2 * i_xy) / h_total)
    elif distance_formula == "d6":
        distance = 1 - ((2 * i_xy) / lambda_xy)
    return distance


"""
    Generates a bootstrap sample from the given vector.

    This function creates a new list by randomly sampling the input vector. 
    It is used to generate resampled coverage vectors for bootstrapping support 
    in phylogenetic trees.

    Parameters:
    - vector (list[float]): A list of values to sample from.

    Returns:
    - list[float]: A new list of the same length as the input, where each element
    is randomly selected from the original vector.
"""
def bootstrap_vector(vector):
    return [random.choice(vector) for _ in vector]


"""
    Main execution function for pairwise genome comparison and phylogenetic tree construction.
    
    This function performs the following steps:

    1. Loads genome sequences from a specified directory.
    2. Initializes a distance matrix and distributes pairwise genome comparison tasks
        across multiple MPI processes.
    3. Each process performs reciprocal BLAST between assigned genome pairs:
        - Parses BLAST results to extract alignments and sequence identities.
        - Builds per-genome coverage vectors representing position-wise identity coverage.
        - Aggregates coverage to calculate sequence identity and alignment metrics.
        - Computes a pairwise distance using a selectable formula (d0, d4, or d6).
    4. The root process gathers all pairwise distances to assemble a symmetric distance matrix.
    5. Outputs:
        - Saves the distance matrix as a CSV file.
        - Optionally constructs a Neighbor Joining tree and saves it in Newick format.
        - Optionally generates a PNG image of the phylogenetic tree.
        - If bootstrapping is enabled:
            - Stores the coverage vectors for each genome pair.
            - Generates multiple bootstrap replicate trees.
            - Calculates bootstrap support values and annotates the original tree.
            - Saves the bootstrapped tree in both Newick format and as a PNG image.

    Parallelization using MPI allows efficient processing of large genome datasets
    by distributing the computational load of pairwise comparisons.
"""
def algorithm():
    genomes = []
    for g in os.listdir(preprocessing.translated_path):
        if g.endswith(".fasta"):
            genomes.append(os.path.join(preprocessing.translated_path, g))

    n = len(genomes)
    distance_matrix = np.zeros((n, n))
    names = [os.path.basename(g).split(".fasta")[0] for g in genomes]
    
    tree_names = names if preprocessing.accession_name else preprocessing.accession_names

    if i_am_main_process():
        if RUN_MODE == "parallel":
            # Create genome pairs (i, j).
            # In parallel mode, the task list is split into chunks, one per worker process.
            tasks = [(i, j) for i in range(n) for j in range(i + 1, n)]
            chunks = [tasks[i::size - 1] for i in range(size - 1)]

            for m in range(1, size):
                # Master process sends each worker its subset of genome pairs to process.
                comm.send((genomes, chunks[m - 1]), dest=m)

            for m in range(1, size):
                # Master process collects computed distances from each worker
                # and fills the distance matrix.
                results = comm.recv(source=m)
                for i, j, distance in results:
                    distance_matrix[i][j] = distance
                    distance_matrix[j][i] = distance
        
        else:
            tasks = [(i, j) for i in range(n) for j in range(i + 1, n)]
            for i, j in tasks:
                gi = genomes[i]
                gj = genomes[j]
                gi_name = os.path.basename(gi).split(".fasta")[0]
                gj_name = os.path.basename(gj).split(".fasta")[0]
                gi_db = os.path.join(preprocessing.db_path, gi_name)
                gj_db = os.path.join(preprocessing.db_path, gj_name)

                # Run BLASTP in one direction (gi vs gj).
                # Output is saved in XML format for later parsing.
                first_blast_out = os.path.join(preprocessing.db_path, f"{gi_name}_{gj_name}.xml")
                first_blast = NcbiblastpCommandline(query=gi, db=gj_db, evalue=preprocessing.e_value, word_size=preprocessing.word_size, outfmt=5, out=first_blast_out)
                first_blast()
                time.sleep(0.5)

                # Run BLASTP in other direction (gj vs gi).
                second_blast_out = os.path.join(preprocessing.db_path, f"{gj_name}_{gi_name}.xml")
                second_blast = NcbiblastpCommandline(query=gj, db=gi_db, evalue=preprocessing.e_value, word_size=preprocessing.word_size, outfmt=5, out=second_blast_out)
                second_blast()
                time.sleep(0.5)

                with open(gi, "r") as handle_gi:
                    gi_length = sum(len(record.seq) for record in SeqIO.parse(handle_gi, "fasta"))
                with open(gj, "r") as handle_gj:
                    gj_length = sum(len(record.seq) for record in SeqIO.parse(handle_gj, "fasta"))

                # Coverage vectors store identity values for each position in the genome sequence. 
                # They represent how well each base is covered by reciprocal BLAST alignments.
                coverage_vector_gi = [0] * gi_length
                coverage_vector_gj = [0] * gj_length

                # Variables for distance calculation:
                    # - len_x / len_y: alignment lengths for each direction.
                    # - v_x / v_y: cumulative coverage values across aligned positions.
                    # - h_xy / h_yx: number of covered positions (initialized at 1 to avoid division by zero).
                    # - i_xy: average identity score.
                identities = 0
                len_x = len_y = 0
                h_xy = h_yx = 1
                i_xy = 0
                v_x = v_y = 0

                # Parse BLAST XML output.
                # Iterate through all alignments and high-scoring pairs (HSPs).
                with open(first_blast_out, "r") as handle:
                    blast_records = NCBIXML.parse(handle)
                    for blast_record in blast_records:
                        for alignment in blast_record.alignments:
                            for hsp in alignment.hsps:
                                len_x = (hsp.query_end - hsp.query_start) + 1
                                identity = (hsp.identities + ((hsp.positives - hsp.identities) * preprocessing.positives)) / len_x
                                identities += identity
                                for r in range(hsp.query_start - 1, hsp.query_end - 1):
                                    if coverage_vector_gi[r] < identity:
                                        coverage_vector_gi[r] = identity

                with open(second_blast_out, "r") as handle:
                    blast_records = NCBIXML.parse(handle)
                    for blast_record in blast_records:
                        for alignment in blast_record.alignments:
                            for hsp in alignment.hsps:
                                len_y = (hsp.query_end - hsp.query_start) + 1
                                identity = (hsp.identities + ((hsp.positives - hsp.identities) * preprocessing.positives)) / len_y
                                identities += identity
                                for r in range(hsp.query_start - 1, hsp.query_end - 1):
                                    if coverage_vector_gj[r] < identity:
                                        coverage_vector_gj[r] = identity

                for l in range(len(coverage_vector_gi)):
                    if coverage_vector_gi[l] != 0:
                        v_x += coverage_vector_gi[l]
                        h_xy += 1

                for l in range(len(coverage_vector_gj)):
                    if coverage_vector_gj[l] != 0:
                        v_y += coverage_vector_gj[l]
                        h_yx += 1

                if preprocessing.bootstrapping > 0:
                    # Save coverage vectors for this genome pair.
                    # These will be re-sampled later during bootstrapping to build replicate trees.
                    os.makedirs('Output/CoverageVectors', exist_ok=True)
                    np.save(f'Output/CoverageVectors/{gi_name}_{gj_name}.npy', np.array(coverage_vector_gi, dtype=np.float32))
                    np.save(f'Output/CoverageVectors/{gj_name}_{gi_name}.npy', np.array(coverage_vector_gj, dtype=np.float32))

                h_total = h_xy + h_yx
                if len_x != 0 or len_y != 0:
                    i_xy = (v_x + v_y) / 2

                lambda_xy = gi_length + gj_length
                distance = distance_calculator(preprocessing.distance_formula, h_total, i_xy, lambda_xy)
                distance_matrix[i][j] = distance
                distance_matrix[j][i] = distance

    else:
        genomes, my_tasks = comm.recv(source=0)
        results = []

        for i, j in my_tasks:
            gi = genomes[i]
            gj = genomes[j]
            gi_name = os.path.basename(gi).split(".fasta")[0]
            gj_name = os.path.basename(gj).split(".fasta")[0]
            gi_db = os.path.join(preprocessing.db_path, gi_name)
            gj_db = os.path.join(preprocessing.db_path, gj_name)

            # Run BLASTP in one direction (gi vs gj).
            # Output is saved in XML format for later parsing.
            first_blast_out = os.path.join(preprocessing.db_path, f"{gi_name}_{gj_name}.xml")
            first_blast = NcbiblastpCommandline(query=gi, db=gj_db, evalue=preprocessing.e_value, word_size=preprocessing.word_size, outfmt=5, out=first_blast_out)
            first_blast()
            time.sleep(0.5)

            # Run BLASTP in other direction (gj vs gi).
            second_blast_out = os.path.join(preprocessing.db_path, f"{gj_name}_{gi_name}.xml")
            second_blast = NcbiblastpCommandline(query=gj, db=gi_db, evalue=preprocessing.e_value, word_size=preprocessing.word_size, outfmt=5, out=second_blast_out)
            second_blast()
            time.sleep(0.5)

            with open(gi, "r") as handle_gi:
                gi_length = sum(len(record.seq) for record in SeqIO.parse(handle_gi, "fasta"))

            with open(gj, "r") as handle_gj:
                gj_length = sum(len(record.seq) for record in SeqIO.parse(handle_gj, "fasta"))

            # Coverage vectors store identity values for each position in the genome sequence. 
            # They represent how well each base is covered by reciprocal BLAST alignments.
            coverage_vector_gi = [0] * gi_length
            coverage_vector_gj = [0] * gj_length

            # Variables for distance calculation:
                # - len_x / len_y: alignment lengths for each direction.
                # - v_x / v_y: cumulative coverage values across aligned positions.
                # - h_xy / h_yx: number of covered positions (initialized at 1 to avoid division by zero).
                # - i_xy: average identity score.
            identities = 0
            len_x = len_y = 0
            h_xy = h_yx = 1
            i_xy = 0
            v_x = v_y = 0

            # Parse BLAST XML output.
            # Iterate through all alignments and high-scoring pairs (HSPs).
            with open(first_blast_out, "r") as handle:
                blast_records = NCBIXML.parse(handle)
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            len_x = (hsp.query_end - hsp.query_start) + 1
                            identity = (hsp.identities + ((hsp.positives - hsp.identities) * preprocessing.positives)) / len_x
                            identities += identity
                            for r in range(hsp.query_start - 1, hsp.query_end - 1):
                                if coverage_vector_gi[r] < identity:
                                    coverage_vector_gi[r] = identity

            with open(second_blast_out, "r") as handle:
                blast_records = NCBIXML.parse(handle)
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            len_y = (hsp.query_end - hsp.query_start) + 1
                            identity = (hsp.identities + ((hsp.positives - hsp.identities) * preprocessing.positives)) / len_y
                            identities += identity
                            for r in range(hsp.query_start - 1, hsp.query_end - 1):
                                if coverage_vector_gj[r] < identity:
                                    coverage_vector_gj[r] = identity

            for l in range(len(coverage_vector_gi)):
                if coverage_vector_gi[l] != 0:
                    v_x += coverage_vector_gi[l]
                    h_xy += 1

            for l in range(len(coverage_vector_gj)):
                if coverage_vector_gj[l] != 0:
                    v_y += coverage_vector_gj[l]
                    h_yx += 1

            if preprocessing.bootstrapping > 0:
                # Save coverage vectors for this genome pair.
                # These will be re-sampled later during bootstrapping to build replicate trees. 
                os.makedirs('Output/CoverageVectors', exist_ok=True)
                np.save(f'Output/CoverageVectors/{gi_name}_{gj_name}.npy', np.array(coverage_vector_gi, dtype=np.float32))
                np.save(f'Output/CoverageVectors/{gj_name}_{gi_name}.npy', np.array(coverage_vector_gj, dtype=np.float32))

            h_total = h_xy + h_yx
            if len_x != 0 or len_y != 0:
                i_xy = (v_x + v_y) / 2

            lambda_xy = gi_length + gj_length
            distance_formula = preprocessing.distance_formula
            distance = distance_calculator(distance_formula, h_total, i_xy, lambda_xy)

            results.append((i, j, distance))

        comm.send(results, dest=0)

    if i_am_main_process():
        os.makedirs("Output", exist_ok=True)

        with open(f'Output/{preprocessing.output_file_prefix}.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow([f'{str(n)}'])
            for i in range(n):
                row = [names[i]]
                for j in range(n):
                    row.append(f'{distance_matrix[i][j]:.8f}')
                writer.writerow(row)

        if preprocessing.generate_tree:
            triangle = []
            for i in range(n):
                triangle.append([float(distance_matrix[i][j]) for j in range(i + 1)])

            dm = DistanceMatrix(tree_names, triangle)
            constructor = DistanceTreeConstructor()
            tree = constructor.nj(dm)
            tree.rooted = preprocessing.root_tree

            for node in tree.get_nonterminals():
                node.name = None

            fig_width = max(10, 0.12 * max(len(name) for name in tree_names))
            fig_height = max(5, 0.3 * len(tree_names))
            fig = plt.figure(figsize=(fig_width, fig_height))
            ax = fig.add_subplot(1, 1, 1)

            Phylo.write(tree, file=f'Output/{preprocessing.output_file_prefix}_tree.nhx', format="newick")

            if preprocessing.generate_image:
                Phylo.draw(tree, do_show=False, axes=ax)
                plt.savefig(f"Output/{preprocessing.output_file_prefix}_tree.png")

    if preprocessing.bootstrapping > 0:
        vector_pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]

        if i_am_main_process():
            if RUN_MODE == "parallel":
                num_replicas = preprocessing.bootstrapping
                num_workers = size - 1
                replicas_por_worker = [num_replicas // num_workers] * num_workers
                for i in range(num_replicas % num_workers):
                    replicas_por_worker[i] += 1

                for worker_rank in range(1, size):
                    # Master sends all information required for bootstrapping to each worker:
                    # genome names, vector pairs, number of replicates, and distance formula.
                    comm.send(names, dest=worker_rank, tag=0)
                    comm.send(tree_names, dest=worker_rank, tag=99)
                    comm.send(vector_pairs, dest=worker_rank, tag=1)
                    comm.send(replicas_por_worker[worker_rank - 1], dest=worker_rank, tag=2)
                    comm.send(preprocessing.distance_formula, dest=worker_rank, tag=3)

                list_of_trees = []
                for worker_rank in range(1, size):
                    trees = comm.recv(source=worker_rank, tag=4)
                    list_of_trees.extend(trees)

                Phylo.write(list_of_trees, "Output/bootstrapped_trees.nhx", "newick")

            else:  
                list_of_trees = []
                # Generates bootstrap replicates by resampling coverage vectors,
                # recalculating pairwise distances, and building replicate Neighbor Joining trees.
                for _ in range(preprocessing.bootstrapping):
                    bs_distance_matrix = np.zeros((n, n))

                    for i, j in vector_pairs:
                        gi_name = names[i]
                        gj_name = names[j]

                        vec_gi = bootstrap_vector(np.load(f'Output/CoverageVectors/{gi_name}_{gj_name}.npy'))
                        vec_gj = bootstrap_vector(np.load(f'Output/CoverageVectors/{gj_name}_{gi_name}.npy'))

                        v_x = sum(vec_gi[l] for l in range(len(vec_gi)) if vec_gi[l] != 0)
                        v_y = sum(vec_gj[l] for l in range(len(vec_gj)) if vec_gj[l] != 0)
                        h_xy = sum(1 for l in vec_gi if l != 0) + 1
                        h_yx = sum(1 for l in vec_gj if l != 0) + 1
                        h_total = h_xy + h_yx
                        i_xy = (v_x + v_y) / 2
                        lambda_xy = len(vec_gi) + len(vec_gj)

                        distance = distance_calculator(preprocessing.distance_formula, h_total, i_xy, lambda_xy)

                        bs_distance_matrix[i][j] = distance
                        bs_distance_matrix[j][i] = distance

                    triangle = [[float(bs_distance_matrix[i][j]) for j in range(i + 1)] for i in range(n)]
                    tree = DistanceTreeConstructor().nj(DistanceMatrix(tree_names, triangle))

                    for node in tree.get_nonterminals():
                        node.name = None
                    list_of_trees.append(tree)

                Phylo.write(list_of_trees, "Output/bootstrapped_trees.nhx", "newick")

        else:
            names = comm.recv(source=0, tag=0)
            tree_names = comm.recv(source=0, tag=99)
            vector_pairs = comm.recv(source=0, tag=1)
            num_replicas = comm.recv(source=0, tag=2)
            distance_formula = comm.recv(source=0, tag=3)
            n = len(names)
            list_of_trees = []

            for _ in range(num_replicas):
                bs_distance_matrix = np.zeros((n, n))

                for i, j in vector_pairs:
                    gi_name = names[i]
                    gj_name = names[j]

                    vec_gi = bootstrap_vector(np.load(f'Output/CoverageVectors/{gi_name}_{gj_name}.npy'))
                    vec_gj = bootstrap_vector(np.load(f'Output/CoverageVectors/{gj_name}_{gi_name}.npy'))

                    v_x = v_y = 0
                    h_xy = h_yx = 1

                    for l in range(len(vec_gi)):
                        if vec_gi[l] != 0:
                            v_x += vec_gi[l]
                            h_xy += 1

                    for l in range(len(vec_gj)):
                        if vec_gj[l] != 0:
                            v_y += vec_gj[l]
                            h_yx += 1

                    h_total = h_xy + h_yx
                    i_xy = (v_x + v_y) / 2
                    lambda_xy = len(vec_gi) + len(vec_gj)
                    distance = distance_calculator(distance_formula, h_total, i_xy, lambda_xy)

                    bs_distance_matrix[i][j] = distance
                    bs_distance_matrix[j][i] = distance

                triangle = [[float(bs_distance_matrix[i][j]) for j in range(i + 1)] for i in range(n)]
                dm = DistanceMatrix(tree_names, triangle)
                tree = DistanceTreeConstructor().nj(dm)

                for node in tree.get_nonterminals():
                    node.name = None

                list_of_trees.append(tree)

            comm.send(list_of_trees, dest=0, tag=4)

        if RUN_MODE == "parallel":
            comm.Barrier()

        if i_am_main_process():
            trees = list(Phylo.parse("Output/bootstrapped_trees.nhx", "newick"))
            original_tree = next(Phylo.parse(f"Output/{preprocessing.output_file_prefix}_tree.nhx", "newick"))
            # Combine original tree with bootstrap replicates to compute branch support values.
            tree_with_support = get_support(original_tree, list_of_trees)

            Phylo.write(tree_with_support, f'Output/{preprocessing.output_file_prefix}_bootstrapped_tree.nhx', 'newick')

            if preprocessing.generate_image:
                fig = plt.figure(figsize=(fig_width, fig_height))
                ax = fig.add_subplot(1, 1, 1)
                Phylo.draw(tree_with_support, do_show=False, axes=ax)
                plt.savefig(f"Output/{preprocessing.output_file_prefix}_bootstrapped_tree.png")

            if os.path.exists("Output/CoverageVectors"):
                shutil.rmtree("Output/CoverageVectors")