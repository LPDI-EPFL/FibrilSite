# Libraries
import os, glob, prody, datetime

import open3d as o3d
from openbabel import pybel

import numpy as np
import pandas as pd
from tqdm import tqdm

from Bio.PDB import PDBParser

# Alignmnet analysis functions

def get_input_desc_diff(df_all_info:pd.DataFrame, df:pd.DataFrame, idx:int, input_feats_npy:list) -> pd.DataFrame:
    """ This function calculates the Input features difference """ 
    
    # define the open3d results file path 
    path_ = [p for p in input_feats_npy if df.at[idx, 'source_pocket'] in p
            if df.at[idx, 'target_pocket'] in p]
    assert len(path_) == 1, 'Path finding ERROR'
    path_ = path_[0]
    
    # load the results file
    reg_results = np.load(path_, allow_pickle=True).item()['icp_test']
    
    # get the correspondences indicies
    source_corres_set = reg_results['corres_set'][0][:,0]
    target_corres_set = reg_results['corres_set'][0][:,1]
    assert len(source_corres_set) == len(target_corres_set)    
    
    # get the corres set properties
    ## for source pocket
    df_src = df_all_info[df_all_info.pocket_id == df.at[idx, 'source_pocket']].reset_index(drop=True)
    source_input_si    = df_src['input_si'].to_numpy()[source_corres_set]
    source_surf_charge = df_src['input_charge'].to_numpy()[source_corres_set]
    source_surf_hphob  = df_src['input_hphob'].to_numpy()[source_corres_set]
    source_surf_hbonds = df_src['input_hbonds'].to_numpy()[source_corres_set]
    
    feats_source = np.stack(( source_surf_charge, source_surf_hphob, source_surf_hbonds, source_input_si), axis=1)
    
    ## for target pocket
    df_target = df_all_info[df_all_info.pocket_id == df.at[idx, 'target_pocket']].reset_index(drop=True)
    target_input_si    = df_target['input_si'].to_numpy()[target_corres_set]
    target_surf_charge = df_target['input_charge'].to_numpy()[target_corres_set]
    target_surf_hphob  = df_target['input_hphob'].to_numpy()[target_corres_set]
    target_surf_hbonds = df_target['input_hbonds'].to_numpy()[target_corres_set]
    
    feats_target = np.stack(( target_surf_charge, target_surf_hphob, target_surf_hbonds, target_input_si), axis=1)
    
    assert feats_source.shape == feats_target.shape, 'Feats parsing failure'
    
    input_desc_diff = np.mean(np.linalg.norm(feats_target - feats_source, axis = 1))

    return input_desc_diff

def calc_input_feat_diff(df:pd.DataFrame, df_all_info:pd.DataFrame, input_feats_npy:list, output:str, export:bool=False, ident=None ) -> pd.DataFrame :
    """ This function is to calculate the input feat differece, raw and weighted """

    # calculate the input feature difference
    input_feat_diff_vessel = []

    for idx01 in df.index:
        input_feat_diff_v = get_input_desc_diff(df_all_info=df_all_info, df=df, idx=idx01, input_feats_npy=input_feats_npy)

        # append    
        input_feat_diff_vessel.append(input_feat_diff_v)

    # add to the dataframe
    df.insert(9, 'icp_mean_input_diff', input_feat_diff_vessel)

    assert len(input_feat_diff_vessel) == df.shape[0]

    if export and ident != None:
        df.to_csv(os.path.join(output, str(ident) + '_matches.csv'))
    
    return df

def fix_xyz(name:str, file:str, output:str):
    """ """
    # read the file
    lines = [f.strip() for f in open(file, 'r')]

    # new file name
    out_path = f'{output}/{name}_fixed_coords.xyz'

    if len(lines[0].split()) == 4 :
        with open(out_path, 'w') as w:
            w.write(str(len(lines)) + '\n')
            w.write('mod' + '\n')
            for i in range(0,len(lines),1):
                w.write(lines[i][0:])
                if i < len(lines)-1:
                    w.write('\n')
    else:
        with open(out_path, 'w') as w:
            w.write(str(len(lines)) + '\n')
            w.write('mod' + '\n')
            for i in range(0,len(lines),1):
                w.write('C ' + lines[i])
                if i < len(lines)-1:
                    w.write('\n')

    return out_path

def xyz2pdb(file:str, output:str):
    """ Covnvert xyz files to pdb files using obabel """
    
    out_path = f"{output}/{os.path.basename(file).replace('xyz','pdb')}"
    for molecule in pybel.readfile("xyz", file):
        molecule.write("pdb", out_path)
    return out_path

def get_trans_matrix(mobile:str, target:str):
    """ This function is to calculate the transformation matrix required to get the pockets from their original place to the aligned one """
    struct_start = prody.parsePDB(mobile)
    struct_end   = prody.parsePDB(target)
    return prody.superpose(mobile=struct_start, target=struct_end)[1].getMatrix()

def transform_fibrils(file:str, pocket:str, mat:np.array, output:str):
    """ apply calculated transformation matrix to fibril's originating pockets"""
    struct = prody.parsePDB(file)
    prody.moveAtoms(atoms=struct, by=mat)
    prody.writePDB(filename=f"{output}/{os.path.basename(file).split('_')[0]}_{pocket}_transformed.pdb", atoms=struct)
    return None

def write_xyz(pocket:str, coords:np.array, output:str):
    """ This function is to write xyz files for the correspondece sets from their np.array """
    with open(f"{output}/{pocket}_corres_set.xyz", 'w') as w:
        for item in coords:
            w.write(f"C {item[0]} {item[1]} {item[2]}" + '\n')
    return None   

def align_site_pcd(df:pd.DataFrame, ply_files:list, input_feats_npy:list, output:str):
    """ 
    This function is to align the sites point clouds from the pocket matches
    
    df               : dataframe containing the site matches 
    ply_files        : list containing paths to sites convex hull ply files
    input_feats_npy  : list containing paths to the alignment outputs for registering pocket point cloud pairs
    output           : path to export folder root
    """
    
    # align all defined site matches
    for idx00 in tqdm(df.index, desc="aligning sites"):
        src_pocket    = df.at[idx00, "source_pocket"]
        target_pocket = df.at[idx00, "target_pocket"]
        
        # make each match output folder
        root_out_temp = os.path.join(output, src_pocket + "_" + target_pocket)
        os.makedirs(root_out_temp, exist_ok=False)
        
        out_temp = os.path.join(root_out_temp, "aligned_sites")
        os.makedirs(out_temp, exist_ok=False)
        
        # get the pockets ply file paths
        src_ply    = [p for p  in ply_files if src_pocket.split("_")[-1]    in os.path.basename(p).split("_convex")[0].split("_")][0]
        target_ply = [p for p  in ply_files if target_pocket.split("_")[-1] in os.path.basename(p).split("_convex")[0].split("_")][0]
        
        # get the pockets registeration output py file path
        reg_npy = [
            p for p in input_feats_npy 
            if src_pocket in p
            if target_pocket in p
        ]
        
        #print(reg_npy)
        assert len(reg_npy) == 1, f"Many registeration output have been found for {src_pocket} and {target_pocket}"
        reg_npy = reg_npy[0] # extract the path in case all good
        
        
        # load the icp registeration results
        reg_results = np.load(reg_npy, allow_pickle=True).item()['icp_test']
            
        # load the point clouds
        source = o3d.io.read_point_cloud(src_ply, format='ply')
        target = o3d.io.read_point_cloud(target_ply, format='ply')
            
        # prepare source pc
        source_corr = o3d.geometry.PointCloud()
        source_corr.points  = o3d.utility.Vector3dVector(np.asarray(source.points) - np.mean(np.asarray(source.points), axis=0))        ## remove center of mass
        source_corr.transform(reg_results['transformation'][0])
        o3d.io.write_point_cloud(f"{out_temp}/{src_pocket.split('_')[-1]}_transformed.xyz", source_corr)
            
        # prepare target pc
        target_corr = o3d.geometry.PointCloud()
        target_corr.points  = o3d.utility.Vector3dVector(np.asarray(target.points) - np.mean(np.asarray(target.points), axis=0))
        o3d.io.write_point_cloud(f"{out_temp}/{target_pocket.split('_')[-1]}_aligned.xyz", target_corr)
        
    return None

def align_site_fibrils(align_out:str, src_pockets_paths:list, ply_files:list, fibrils_paths:list, fibril_pocket_map:dict, reg_results_paths:list):
    """ 
    This function is to align the fibril where shared features were identified
    
    align_out          : path to aligned sites point clouds
    src_pockets_paths  : list of paths to the defined fibril sites
    ply_files          : list of paths to the sites ply files
    fibrils_paths      : list of paths to the fibrils pdb files
    fibril_pocket_map  : dictionary mapping the fibrils to the defined sites
    reg_results_paths  : list of paths to the alignmnets output
    """

    for file in tqdm(sorted([os.path.abspath(p).strip() for p in glob.iglob(os.path.join(align_out, "*"))]), desc="aligning fibrils"):
        
        # make output dirs
        temp_out01 = os.path.join(file, "aligned_fibrils")
        os.makedirs(temp_out01, exist_ok=False)
    
        temp_out02 = os.path.join(temp_out01, "_temp")
        os.makedirs(temp_out02, exist_ok=False)
        
        # get the pockets xyz files
        pocket01_end_xyz = [p.strip() for p in glob.iglob(os.path.join(file, "aligned_sites" ,'*.xyz')) if 'transformed' in os.path.basename(p)][0]
        pocket02_end_xyz = [p.strip() for p in glob.iglob(os.path.join(file, "aligned_sites" ,'*.xyz')) if 'aligned' in os.path.basename(p)][0]
    
        pocket01_name = pocket01_end_xyz.split('/')[-1].split('_transformed')[0]
        pocket02_name = pocket02_end_xyz.split('/')[-1].split('_aligned')[0]

        pocket01_start_xyz = [p.strip() for p in src_pockets_paths if pocket01_name in p.split("/")[-2].split("_")][0]
        pocket02_start_xyz = [p.strip() for p in src_pockets_paths if pocket02_name in p.split("/")[-2].split("_")][0]
        
        # fix the xyz files to match the obabel format
        pocket01_start_xyz = fix_xyz(name=pocket01_name+'_start', file=pocket01_start_xyz, output=temp_out02)
        pocket01_end_xyz   = fix_xyz(name=pocket01_name+'_end',   file=pocket01_end_xyz,   output=temp_out02)
        pocket02_start_xyz = fix_xyz(name=pocket02_name+'_start', file=pocket02_start_xyz, output=temp_out02)
        pocket02_end_xyz   = fix_xyz(name=pocket02_name+'_end',   file=pocket02_end_xyz,   output=temp_out02)
        
        # generate pdb files from xyz files
        pocket01_start_pdb = xyz2pdb(file=pocket01_start_xyz, output=temp_out02)
        pocket01_end_pdb   = xyz2pdb(file=pocket01_end_xyz,   output=temp_out02)
        pocket02_start_pdb = xyz2pdb(file=pocket02_start_xyz, output=temp_out02)
        pocket02_end_pdb   = xyz2pdb(file=pocket02_end_xyz,   output=temp_out02)
        
        # get the transformation matrix between pocket start and end states
        pocket01_mat = get_trans_matrix(mobile=pocket01_start_pdb, target=pocket01_end_pdb).T
        pocket02_mat = get_trans_matrix(mobile=pocket02_start_pdb, target=pocket02_end_pdb).T
    
        # get the respective fibrils paths
        pocket01_fibril = [path_ for path_ in fibrils_paths if [k for k,v in fibril_pocket_map.items() for item in v if item == pocket01_name.split('_')[-1]][0] == os.path.basename(path_).split("_")[0]][0]
        pocket02_fibril = [path_ for path_ in fibrils_paths if [k for k,v in fibril_pocket_map.items() for item in v if item == pocket02_name.split('_')[-1]][0] == os.path.basename(path_).split("_")[0]][0]
    
        # transform fibrils
        transform_fibrils(file=pocket01_fibril, pocket=pocket01_name, mat=pocket01_mat, output=temp_out01)
        transform_fibrils(file=pocket02_fibril, pocket=pocket02_name, mat=pocket02_mat, output=temp_out01)
    
        # load the pocket point clouds
        pocket01_pcd_path = [p for p in ply_files if pocket01_name == os.path.basename(p).split("_convex")[0].split("_")[-1]][0]
        pocket02_pcd_path = [p for p in ply_files if pocket02_name == os.path.basename(p).split("_convex")[0].split("_")[-1]][0]
    
        pocket01_pcd = o3d.io.read_point_cloud(pocket01_pcd_path, format='ply')
        pocket02_pcd = o3d.io.read_point_cloud(pocket02_pcd_path, format='ply')
    
        reg_npy = [p for p in reg_results_paths 
                   if pocket01_name in os.path.basename(p).split("_input-feats")[0].split("_")
                   if pocket02_name in os.path.basename(p).split("_input-feats")[0].split("_")][0]
        
        reg_results = np.load(reg_npy, allow_pickle=True).item()['icp_test']
    
        # prepare pocket01
        pocket01_pcd_corr = o3d.geometry.PointCloud()
        pocket01_pcd_corr.points  = o3d.utility.Vector3dVector(np.asarray(pocket01_pcd.points) - np.mean(np.asarray(pocket01_pcd.points), axis=0))
        pocket01_pcd_corr.transform(reg_results['transformation'][0])
    
        # prepare pocket02
        pocket02_pcd_corr = o3d.geometry.PointCloud()
        pocket02_pcd_corr.points  = o3d.utility.Vector3dVector(np.asarray(pocket02_pcd.points) - np.mean(np.asarray(pocket02_pcd.points), axis=0))
    
        # get correspondence sets for each pocket
        pocket01_corres_set = reg_results['corres_set'][0][:,0]
        pocket02_corres_set = reg_results['corres_set'][0][:,1]
        pocket01_pcd_corr_corres_set = np.asarray(pocket01_pcd_corr.points)[pocket01_corres_set]
        pocket02_pcd_corr_corres_set = np.asarray(pocket02_pcd_corr.points)[pocket02_corres_set]
    
        # export correspondence sets as xyz
        write_xyz(pocket=pocket01_name, coords=pocket01_pcd_corr_corres_set, output=temp_out01)
        write_xyz(pocket=pocket02_name, coords=pocket02_pcd_corr_corres_set, output=temp_out01)
    
    return None
