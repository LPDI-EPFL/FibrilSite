# FibrilSite  – Site Alignment
To run the defined fibril sites alignment and alignment analysis:

## Provided scripts:
- **01_run_fibril_site_surface_feature_mapping.sh :** \
    executable script to run *[./scripts/fibril_site_input_feats_mapper.py]* that will extract the surface features for the defined fibril site points following MaSIF processing of the fibril pdb files.
    
    **params**: \
        - site_src_dir    : Path to folder containing the defined fibril sites \
        - input_feats_src : Path to folder contaning MaSIF calculated surface features *[/data_preparation/04b-precomputation_12A/precomputation]* \
        - output_folder   : Path to root folder for exporting

    **output**: \
        - {date}_all_sites_parsed.csv      : CSV containing information of all defined sites \
        - {date}_all_sites_input_feats.csv : CSV containing information of all defined sites including the calculated surface features
    
    **exported information**: \
        - fibril               : fibril PDB id \
        - pocket_id            : pocket (site) name \
        - isolation            : level of pocket definition \
        - MaSIF_index          : pocket (site) surface points MaSIF index -> needed for mapping \
        - atom_type            : surface point atom type – if applicable \
        - chain                : fibril chain \
        - coords               : surface point atomic coordinates – if applicable \
        - point_direction      : dot product of surface normals of one anchor point and its neighbour -> needed for pocket expansion \
        - resid                : residue id \
        - resname              : residue name \
        - sasa                 : per residue solvent-accessible surface area \
        - surf_charge          : ply-file parsed surface electrostatics value \
        - surf_coords          : ply-file parsed surface coordinates \
        - surf_hbond           : ply-file parsed hydrogen-bond donor/acceptor potential \
        - surf_hphob           : ply-file parsed hydropathy value \
        - surf_norm_fibril_dot : dot product of point surface normal and fibril elongation axis \
        - surf_normals         : surface point normal \
        - input_si             : surface point calculated shape index \
        - input_charge         : surface point calculated electrostatics value \
        - input_hphob          : surface point calculated hydropathy value \
        - input_hbonds         : surface point calculated hydrogen-bond donor/acceptor potential 

- **02_run_fibril_site_registeration.sh :** \
    executable script to run *[./scripts/fibril_site_registeration.py]* that will align the defined fibril sites based on surface features. 

    **params**: \
        - info_file     : CSV containing information of all defined sites including the calculated surface features \
        - sites_folder  : Path to folder containing the defined fibril sites \
        - output_folder : Path to root folder for exporting
     
     **output**: \
        - o3d_objects_npy                    : Folder containing the numpy site alignment results for site-pairs \
        - {date}_all_sites_alignment_results : CSV containing information of fibril sites alignment
   
    **exported information**: \
        - source_pocket | target_pocket : aligned pockets (sites) names\
        - ransac_rmse                   : RMSE (Root Mean Square Error) computed on inlier correspondences following RANSAC (Random Sample Consensus) \
        - ransac_fitness                : alignment quality following RANSAC \
        - ransac_nb_corres              : number of correspondences following RANSAC \
        - icp_rmse                      : RMSE computed on inlier correspondences following ICP (Iterative Closest Point) \
        - icp_fitness                   : alignment quality following ICP \
        - icp_nb_corres                 : number of correspondences following ICP \
        - size_source | size_target     : number of points of the size and target pockets (sites)

- **03_run_fibril_site_alignment_analysis.sh :**
    executable script to run *[./scripts/fibril_site_alignment_analysis.py]* that will analyse the sites alignment results. \

    **params**: \
        - site_src_dir     : Path to folder containing the defined fibril sites \
        - fibrils_pdb_src  : Path to folder containing fibrils pdb files \
        - reg_results_src  : Path to folder containing alignment results npy files \
        - sites_info_csv   : Path to CSV file containing fibril sites info \
        - sites_align_csv  : Path to CSV file containing fibril sites alignment results \
        - SSmax            : Threshold for Site Surface overlap between site matches, recommended >= 0.5 \
        - Fdiff            : Threshold for Surface feature difference between site matches, recommended <= 0.6 \
        - output_folder    : Path to root folder for exporting

     **output**: \
        - all_vs_all   : Folder containing CSV files with the top 5 matches for each pocket based on surface feature similarity only  \
        - aSyn_to_aSyn : Folder containing two CSV files for site matches among ex vivo aSyn fibrils \
        - aSyn_to_invitro_aSyn : Folder containing two CSV files for site matches between ex vivo and in vitro aSyn fibrils \
        - aSyn_to_other_amyloids : Folder containing two CSV files for site matches between ex vivo fibrils of aSyn and other amyloids \
        - invitro_aSyn_to_other_amyloids : Folder containing two CSV files for site matches between in vitro aSyn fibrils and ex vivo fibrils of other amyloids \
             [ 1-"all_{}_matches.csv" : all site matches | 2-"sel_{}_matches.csv" : site matches satisfying the site similarity condition ] \
        - identified_matches_alignments : Folder containing fibril and site alignments for identified site matches
      
    **exported information**: \
        - pocket_matches       : a tuple of pocket (site) matches \
        - pocket_pairs         : a list of pocket (site) matches \
        - source_pocket | target_pocket : aligned pockets (sites) names\
        - icp_rmse             : RMSE computed on inlier correspondences following ICP (Iterative Closest Point) \
        - icp_fitness_source   : alignment quality following ICP, calculated based on source pocket \
        - icp_fitness_target   : alignment quality following ICP, calculated based on target pocket \
        - icp_nb_corres        : number of correspondences following ICP \
        - size_source | size_target : number of points of the size and target pockets (sites) \
        - src2target_size_ratio : ratio of the pocket match sizes, calculated as number of points of source pocket / number of points of target pocket \
        - source_pocket_fibril | target_pocket_fibril    : name of respective fibril for the defined pockets (sites) \
        - combined_fitness_score (**SSmax**): best alignment between two pockets (sites), SSmax = max(icp_fitness_source, icp_fitness_target) \
        - icp_mean_input_diff (**Fdiff**): surface features similarity between the pocket (site) matches, calculated as the Euclidean distance between the surface features computed on inlier correspondences following ICP 
    
The *[./scripts/fibril_site_alignment_analysis.py]* script is tailored for the current study. For more customisation in case of site database expansion with sites from other fibrils, you would need to use the provided notebook **fibril_site_alignment_analysis.ipynb**. 

> To reproduce the results from this paper, use the provided dataset at Zenodo (https://doi.org/10.5281/zenodo.15192320) and run the scripts 01, 02, 03 sequentially (~20 minutes).

## Usage
1- Create a conda environment with python version 3.8
    
    conda create -n alignsites python=3.8

2- Install the code package in the **fibrilsite** environment 

    pip install git+https://github.com/A-Sadek/FibrilSite.git

Or

    pip install fibrilsite

3- Install Openbabel

    conda install conda-forge::openbabel
    
4- Export the **alignsites** environment to Jupyter as follows:

    conda install -c anaconda ipykernel -y
    python -m ipykernel install --user --name=alignsites
 
5- Activate the **alignsites** environment

6- Map the defined sites to their calculted surface features using 

    source 01_run_fibril_site_surface_feature_mapping.sh

7- Run the all-vs-all alignment 
    
    source 02_run_fibril_site_registeration.sh

8- Run the alignment analysis 
    
    source 03_run_fibril_site_alignment_analysis.sh

or through the **fibril_site_alignment_analysis.ipynb** notebook

9- Validate the matches through visualization 

