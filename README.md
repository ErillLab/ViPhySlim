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

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/ErillLab/ViPhySlim.git
cd ViPhySlim
```

### 2. Create and activate a virtual environment
```bash
python3 -m venv viphyslim-env
source viphyslim-env/bin/activate
```

### 3. Install dependencies 
```bash
pip install -r requirements.txt
```

## Usage 
```bash
mpiexec -np <num_processes> python main.py
```

## Example
```bash
mpiexec -np 5 python main.py
```

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Citation
If you use **ViPhySlim** in your research, please cite this repository.
