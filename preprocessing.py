import os
from Bio import Entrez, SeqIO
import yaml
import shutil
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbimakeblastdbCommandline
import time

'''
    Load configuration from YAML and create
    working directories
'''
yaml_file = "configuration.yml"
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

with open(yaml_file, "r") as file:
    document = yaml.safe_load(file)
    email = document["entrez_parameters"]["user_email"]
    genomes_path = os.path.join(BASE_DIR, document["input_parameters"]["genomes_path"])
    working_path = os.path.join(BASE_DIR, document["operational_parameters"]["working_path"])
    output_file_prefix = document["output_parameters"]["output_file_name"]
    erase_working_folder = document["operational_parameters"]["erase_working_folder"]
    accessions = document["input_parameters"]["accession_numbers"]
    
    distance_formula = document['algorithm_parameters']['distance_formula']
    bootstrapping = document['algorithm_parameters']['bootstrap_replicates']
    
    word_size = document["BLAST_parameters"]["word_size"]
    positives = document["BLAST_parameters"]["positives"]
    e_value = document["BLAST_parameters"]["e_value"]
    
    generate_tree = document["tree_parameters"]["generate_tree"]
    root_tree = document["tree_parameters"]["root_tree"]
    accession_name = document["tree_parameters"]["accession_name"]
    generate_image = document["tree_parameters"]["generate_image"]
    
Entrez.email = email

translated_path = os.path.join(working_path, "Translated")
db_path = os.path.join(working_path, "DB")
os.makedirs(translated_path, exist_ok=True)
os.makedirs(db_path, exist_ok=True)

if not accession_name:
    accession_names = []

'''
    Translate the given sequence record into 
    six reading frames
'''
def translate_complementary(record):
    forward1 = record.seq.translate()
    forward2 = record.seq[1:].translate()
    forward3 = record.seq[2:].translate()

    reverse = record.seq.reverse_complement()
    reverse1 = reverse.translate()
    reverse2 = reverse[1:].translate()
    reverse3 = reverse[2:].translate()

    sequence = str(forward1) + str(forward2) + str(forward3) + str(reverse1) + str(reverse2) + str(reverse3)

    return SeqRecord(Seq(sequence), id=record.id, description=f'{record.id} 6-frames translation')


'''
    Reads the list of accession numbers and if the
    corresponding genome file does not exist, it 
    downloads them from the NCBI database
'''
def check_and_download_genomes():
    for accession in accessions:
        fasta_path = os.path.join(genomes_path, f"{accession}.fasta")
        genbank_path = os.path.join(genomes_path, f"{accession}.gb")

        if os.path.exists(fasta_path) or os.path.exists(genbank_path):
            continue

        with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
            record = SeqIO.read(handle, "genbank")

        with open(fasta_path, "w") as fasta_file:
            SeqIO.write(record, fasta_file, "fasta")
        
        time.sleep(1)
        

'''
    Loads each genome file and translates it
    using the "translate_six_frames" function.
    The translated sequence is saved in the
    "Translated" folder inside the working directory. 
    Finally, it creates an individual BLAST database
    for each genome and stores it in the "DB" folder 
    in the working directory
'''
def process_and_translate_genomes():
    for filename in os.listdir(genomes_path):
        if filename.endswith(".fasta") or filename.endswith(".gb"):
            filepath = os.path.join(genomes_path, filename)
            accession_number = os.path.splitext(filename)[0]

            if not accession_name:
                if filename.endswith(".fasta"):
                    with open(filepath, "r") as file:
                        header = file.readline().strip() 
                        virus_name = header.split(None, 1)[1].split(',')[0]
                elif filename.endswith(".gb"):
                    with open(filepath, "r") as file:
                        for line in file:
                            if line.startswith("DEFINITION"):
                                virus_name = line.split("  ", 1)[1].split(',')[0]
                                break
                accession_names.append(virus_name)

            record = SeqIO.read(filepath, "fasta" if filename.endswith(".fasta") else "genbank")
            translated_record = translate_complementary(record)
            
            translated_file_path = os.path.join(translated_path, f"{accession_number}.fasta")
            with open(translated_file_path, "w") as result_file:
                SeqIO.write(translated_record, result_file, "fasta")

            individual_db_path = os.path.join(db_path, accession_number)
            create_db_cmd = NcbimakeblastdbCommandline(dbtype="prot", input_file=translated_file_path, out=individual_db_path)
            create_db_cmd()
            time.sleep(0.5)

            
'''
    Runs the preprocessing stage:
        - Erases both folders in WorkingDirectory path
        - Checks and downloads missing genome files
        - Translates genome sequences into six frames and saves them
        - Creates individual BLAST databases for each translated genome
'''
def preprocessing():
    if erase_working_folder:
        shutil.rmtree(db_path)
        shutil.rmtree(translated_path)
        
    if not os.path.exists(db_path):
        os.mkdir(db_path)
    if not os.path.exists(translated_path):
        os.mkdir(translated_path)

    check_and_download_genomes()
    process_and_translate_genomes()