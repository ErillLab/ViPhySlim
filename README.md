# ViPhySlim

**Phylogenomic inference for viruses using whole-genome distance metrics.**

[![GitHub](https://img.shields.io/badge/GitHub-ErillLab%2FViPhySlim-blue?logo=github)](https://github.com/ErillLab/ViPhySlim)

## Overview

**ViPhySlim** is an open-source Python package for viral phylogenetic inference based on whole-genome comparisons. It is inspired by the [VICTOR](https://ggdc.dsmz.de/victor.php) tool, providing a lightweight, scalable, and locally deployable alternative. ViPhySlim is designed to handle a large number of viral genomes, leveraging parallel computation via MPI.

### Why ViPhySlim?

Inferring the evolutionary relationships of viruses is a challenging task. Unlike cellular organisms, viral evolution is heavily influenced by horizontal gene transfer, making traditional gene-based phylogenetics less reliable. Genome-wide approaches offer better accuracy, but existing tools have limitations:

- Closed-source  
- Web-only interface  
- Limited to 100 genomes per analysis  

**ViPhySlim** addresses these challenges by offering:

- Local and scalable phylogenetic analysis  
- Parallel processing with MPI  
- Open-source Python implementation  
- Conda-packaged for easy installation and reproducibility  

## Features

- Pairwise genome distance computation  
- MPI-based parallel processing for high-throughput scalability  
- Output compatible with phylogenetic tree building tools  
- Easily deployable via `venv` or Conda (BioConda channel coming soon)

## Distance Metrics

ViPhySlim implements **whole-genome distance estimation** using the **Genome BLAST Distance Phylogeny (GBDP)** approach. After computing BLASTp alignments between translated viral genomes, distances are calculated using one of three formulas:

- **d0**: Proportion of the genome covered by high-scoring segment pairs (HSPs).  
- **d4**: Fraction of identical amino acid pairs within HSPs relative to the total aligned length. This measure is more robust when working with incomplete genomes.  
- **d6**: Fraction of identical amino acid pairs relative to the entire genome size. This metric preserves the highest amount of evolutionary information.  

The choice of distance depends on the dataset characteristics and the research focus.  
By default, **d6** is recommended for complete genomes, while **d4** is more reliable for fragmented or incomplete assemblies.

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/ErillLab/ViPhySlim.git
cd ViPhySlim
```

### 2. Create and activate the Conda environment
```bash
conda env create -f viphyslim_environment.yml
conda activate viphyslim
```

## Usage
Type of execution: "parallel":
```bash
mpiexec -np <num_processes> python main.py
```
### Example
```bash
mpiexec -np 5 python main.py
```

Type of execution: "serial":
```bash
python main.py
```

**Note:** When running the code, make sure the execution mode (parallel or serial) corresponds to the setting chosen in the configuration file.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Citation
If you use **ViPhySlim** in your research, please cite this repository.
