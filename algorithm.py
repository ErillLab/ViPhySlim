import preprocessing
from Bio import SeqIO, Phylo
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
import matplotlib.pyplot as plt
from mpi4py import MPI
import numpy as np
import csv
import os
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

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
    Main execution function that performs the pairwise comparison of all genomes in a directory.
    It does the following:
    1. Loads genome file paths and initializes a distance matrix.
    2. Distributes genome pair comparison tasks across multiple MPI processes.
    3. For each genome pair (gi, gj), performs:
    - Reciprocal BLAST (gi -> gj and gj -> gi),
    - Parsing of BLAST outputs to extract alignments and identities,
    - Construction of coverage vectors for both genomes,
    - Calculation of total identities and covered positions,
    - Computation of a pairwise distance using distance_calculator.
    4. Collects results in the root process and builds a symmetric distance matrix.
    5. Optionally:
    - Writes the distance matrix to a CSV file,
    - Constructs a phylogenetic tree using the Neighbor Joining method,
    - Saves the tree in Newick format and optionally as a PNG image.

    This function is parallelized with MPI to efficiently handle large genome datasets.
"""
def algorithm():
    genomes = []
    for g in os.listdir(preprocessing.translated_path):
        if g.endswith(".fasta"):
            genomes.append(os.path.join(preprocessing.translated_path, g))

    n = len(genomes)
    distance_matrix = np.zeros((n, n))
    names = [os.path.basename(g).split(".fasta")[0] for g in genomes]

    if rank == 0:
        tasks = [(i, j) for i in range(n) for j in range(i + 1, n)]
        chunks = [tasks[i::size - 1] for i in range(size - 1)]

        for m in range(1, size):
            comm.send((genomes, chunks[m - 1]), dest=m)

        for m in range(1, size):
            results = comm.recv(source=m)
            for i, j, distance in results:
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

            first_blast_out = os.path.join(preprocessing.db_path, f"{gi_name}_{gj_name}.xml")
            first_blast = NcbiblastpCommandline(query=gi, db=gj_db, evalue=preprocessing.e_value, word_size=preprocessing.word_size, outfmt=5, out=first_blast_out)
            first_blast()
            time.sleep(0.5)

            second_blast_out = os.path.join(preprocessing.db_path, f"{gj_name}_{gi_name}.xml")
            second_blast = NcbiblastpCommandline(query=gj, db=gi_db, evalue=preprocessing.e_value, word_size=preprocessing.word_size, outfmt=5, out=second_blast_out)
            second_blast()
            time.sleep(0.5)

            with open(gi, "r") as handle_gi:
                gi_length = sum(len(record.seq) for record in SeqIO.parse(handle_gi, "fasta"))

            with open(gj, "r") as handle_gj:
                gj_length = sum(len(record.seq) for record in SeqIO.parse(handle_gj, "fasta"))

            coverage_vector_gi = [0] * gi_length
            coverage_vector_gj = [0] * gj_length

            identities = 0
            len_x = 0
            h_xy = 1
            h_yx = 1
            len_y = 0
            i_xy = 0
            v_x = 0
            v_y = 0

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

            h_total = h_xy + h_yx
            if len_x != 0 or len_y != 0:
                i_xy = (v_x + v_y) / 2

            lambda_xy = gi_length + gj_length
            distance_formula = preprocessing.distance_formula
            distance = distance_calculator(distance_formula, h_total, i_xy, lambda_xy)

            results.append((i, j, distance))

        comm.send(results, dest=0)

    if rank == 0:
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
            if not preprocessing.accession_name:
                names = preprocessing.accession_names

            triangle = []
            for i in range(n):
                triangle.append([float(distance_matrix[i][j]) for j in range(i + 1)])
                
            dm = DistanceMatrix(names, triangle)
            constructor = DistanceTreeConstructor()
            tree = constructor.nj(dm)
            tree.rooted = preprocessing.root_tree

            fig_width = max(10, 0.12 * max(len(name) for name in names))
            fig_height = max(5, 0.3 * len(names))
            fig = plt.figure(figsize=(fig_width, fig_height))
            ax = fig.add_subplot(1, 1, 1)

            Phylo.write(tree, file=f'Output/{preprocessing.output_file_prefix}_tree.nhx', format="newick")

            if preprocessing.generate_image:
                Phylo.draw(tree, do_show=False, axes=ax)
                plt.savefig(f"Output/{preprocessing.output_file_prefix}_tree.png")