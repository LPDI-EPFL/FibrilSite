# FibrilSite
**FibrilSite** is a computational pipeline to define fibril sites (i.e., pockets) and identify shared features across a database of defined sites. The results of applying this framework to 61 sites from 19 fibril structures is reported in
> Mapping the Structural Landscape of Amyloid Fibrils to Guide Polymorph-Specific Therapeutics \
Ahmed Sadek, Bruno E. Correia#, Hilal A. Lashuel# \
https://doi.org/10.1101/2025.05.08.652887

Fibril surface sites are extracted from fibril surface as point clouds featurized with surface properties including Poisson–Boltzmann continuum electrostatics, hydrophobicity, shape index and hydrogen bond donors/acceptor patterns. 

### Fibril structures preparation
Fibril structures were pre-processed and surface features were computed using the **MaSIF site** tool within the MaSIF (Molecular Surface Interaction Fingerprinting) framework available at (https://github.com/LPDI-EPFL/masif.git), as described in 
> Gainza, P., Sverrisson, F., Monti, F. et al. Deciphering interaction fingerprints from protein molecular surfaces using geometric deep learning. Nat Methods 17, 184–192 (2020). https://doi.org/10.1038/s41592-019-0666-6

### Fibril sites definition
Fibril site definition is peformed using the *fibril_grooves_extractor.ipynb* notebook.

To run the notebook you need to install the code package containing all the required funtions as follows: 

1- Create a conda environment with python version 3.6

    conda create -n fibrilsite python=3.6

2- Install the The StructuralBioInformatics Library

    pip install StrBioInfo==0.9a0.dev1


3- Install the code package in the **fibrilsite** environment 

    pip install git+https://github.com/LPDI-EPFL/FibrilSite.git

4- Install Pandas v.0.25.3
    
    pip install pandas==0.25.3

5- Install **PyMesh** in the **fibrilsite** environment following the instruction on https://pymesh.readthedocs.io/en/latest/installation.html 

6- Export the **fibrilsite** environment to Jupyter as follows:

    conda install -c anaconda ipykernel -y
    python -m ipykernel install --user --name=fibrilsite

### Fibril sites alignment 
The defined sites are aligned in an all vs all manner using Open3d, then analysed based on the shared surface fraction and the similarity of their surface geometrical and physicochemical features using the provided scripts and instructions in the site_alignment folder

### Data for the work
All data are freely accessible on Zenodo (https://doi.org/10.5281/zenodo.15192320).

