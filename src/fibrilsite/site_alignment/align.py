# Libraries
import os, glob, shutil, copy, prody, datetime

import open3d as o3d

import numpy as np
import pandas as pd
from shutil import copyfile
from tqdm import tqdm

from Bio.PDB import PDBParser

# Functions
# Alignmnet functions
def execute_global_registration(source_pcd, target_pcd, source_feats,
                                target_feats, distance_threshold = 1.5, verbose:bool = False):
    """
    This function is to use Ransac to align the pockets' point clouds
    
    source_pcd           : source pcd that will be the target to align on
    target_pcd           : target pcd that will be aligned to source pcd
    source_feats         : source pcd features that will be used in the alignment
    target_feats         : target pcd features that will be used in the alignment
    distance threshold   : Ransac distance threshold to align two pcd, 1.5 was adopted from MaSIF seed search  
    verbose                    : Parameter controlling verbosity
    """
    
    if verbose:
        print(":: RANSAC registration on downsampled point clouds.")
        print("   we use a liberal distance threshold %.3f." % distance_threshold)
    
    result = o3d.pipelines.registration.registration_ransac_based_on_feature_matching(
        source_pcd, target_pcd, source_feats, target_feats, True,
        distance_threshold,
        o3d.pipelines.registration.TransformationEstimationPointToPoint(False),
        3, 
        [
            o3d.pipelines.registration.CorrespondenceCheckerBasedOnEdgeLength(0.9),
            o3d.pipelines.registration.CorrespondenceCheckerBasedOnDistance(distance_threshold),
            o3d.pipelines.registration.CorrespondenceCheckerBasedOnNormal(np.pi/2),
        ], 
        o3d.pipelines.registration.RANSACConvergenceCriteria(20000, 0.999))
    return result

def refine_registration(source_pcd, target_pcd, result_ransac, distance_threshold = 1.0, verbose:bool = False, generalized: bool = False):
    
    """
    This function is to use ICP to refine the aligned the pockets' point clouds
    
    source_pcd           : source pcd that will be aligned to the target pcd (apply the transformation matrix to this pcd)
    target_pcd           : target pcd that the source pcd will aligned to
    source_feats         : source pcd features that will be used in the alignment
    target_feats         : target pcd features that will be used in the alignment
    distance threshold   : Ransac distance threshold to align two pcd, 1.0 was adopted from MaSIF seed search  
    """
    
    if verbose:
        print(":: Point-to-plane ICP registration is applied on original point")
        print("   clouds to refine the alignment. This time we use a strict")
        print("   distance threshold %.3f." % distance_threshold)
    
    if generalized:
        result = o3d.pipelines.registration.registration_generalized_icp(
            source_pcd, target_pcd, distance_threshold, result_ransac.transformation,
        )
    
    else: 
        result = o3d.pipelines.registration.registration_icp(
            source_pcd, target_pcd, distance_threshold, result_ransac.transformation,
            o3d.pipelines.registration.TransformationEstimationPointToPlane(),
            o3d.pipelines.registration.ICPConvergenceCriteria(max_iteration=100)
        )
    return result

def ply_parser_hull(ply_file:str, out_path=None ) -> pd.DataFrame:
    """
    This function is to parse convex hull MaSIF ply files and return a dataframe with the surface point coords and surface points normals.
    
    out_path : File exporting desitination
    ply_file : ply file to be parsed
    """
    # info container
    dic = {'surf_coords' :[],
           'surf_normals':[],
          }
    # Parse the ply file
    with open(ply_file, 'r') as f:
        lines = [line.strip() for line in f.readlines() if not line[0].isalpha() and len(line.split()) == 9]
    
    for line in lines:
        # get surface point coords
        surf_x, surf_y, surf_z = float(line.split()[0]), float(line.split()[1]), float(line.split()[2])
        surf_coords = np.array([surf_x, surf_y, surf_z])
                
        # get surface normals
        surf_nx, surf_ny, surf_nz = float(line.split()[6]), float(line.split()[7]), float(line.split()[8])
        surf_norm_coords = np.array([surf_nx, surf_ny, surf_nz])

        # add to dict       
        dic['surf_coords'].append(surf_coords)
        dic['surf_normals'].append(surf_norm_coords)
    
    df_surf = pd.DataFrame.from_dict(dic)
    df_surf.reset_index(drop=False, inplace=True)
    df_surf.rename(columns={'index':'MaSIF_index'}, inplace=True)
    if out_path != None:
        df_surf.to_csv(f'{out_path}/{datetime.date.today()}_{os.path.basename(ply_file).split("_")[0]}_masif_params_from_ply.csv')
    
    return df_surf

def best_registration(source_pcd:o3d.geometry.PointCloud, target_pcd:o3d.geometry.PointCloud, source_feats:o3d.pipelines.registration.Feature, target_feats:o3d.pipelines.registration.Feature, n:int = 3, verbose:bool = False):
    """Takes the source and target point cloud, and performs n instances of global registration. In cases where the n 
    attempts fail to find any correspondence, the process is repeated. Outputs the registration result that yields the biggest
    number of correspondences
    
        Args:
        source_pcd (o3d.geometry.PointCloud): Source point cloud
        target_pcd (o3d.geometry.PointCloud): Target point cloud
        source_feats (o3d.pipelines.registration.Feature): Source features
        target_feats (o3d.pipelines.registration.Feature): Target Features
        n (int, optional): Number of global registrations to perform. Defaults to 3.
        verbose (bool, optional): Whether or not to print the performance of the registration. Defaults to False

    Returns:
        (o3d.pypeline.registration.registration_results): _description_ [0] for ransac results, [1] for icp results
    """
    # containers
    results_ransac  = []
    results_icp     = []
    
    n_corres_ransac = []
    n_corres_icp    = []
    
    measures_icp    = []

    for i in range(n):
        result_ransac = execute_global_registration(
            source_pcd=source_pcd,
            target_pcd=target_pcd,
            source_feats=source_feats,
            target_feats=target_feats,
            distance_threshold= 1.5,
        )
        result_icp = refine_registration(
            source_pcd=source_pcd,
            target_pcd=target_pcd,
            result_ransac=result_ransac,
            distance_threshold=1
        )
        
        results_ransac.append(result_ransac)
        results_icp.append(result_icp)
        
        n_corres_ransac.append(len(np.asarray(result_ransac.correspondence_set)))
        n_corres_icp.append(len(np.asarray(result_icp.correspondence_set)))
        
        measures_icp.append([result_icp.inlier_rmse, result_icp.fitness])   # for the verbose function
        
    if verbose:
        print('Chosen alignment measures', measures_icp[np.argmax(n_corres_icp)])
        
    return results_ransac[np.argmax(n_corres_icp)], results_icp[np.argmax(n_corres_icp)]

def global_reg_pipeline(source_pocket:str, target_pocket:str, path_dic:dict, df_pockets:pd.DataFrame, output:str, n_regs:int = 5 ):
    """This function takes the name of two pockets (source and target) and performs global registration (RANSAC followed by ICP)
        It outputs different metrics on the registration

    Args:
        source_pocket (string): Name of the source pocket
        target_pocket (string): Name of the target pocket
        n_regs (int, optional): Number of time to perform the global registration for each pair of pockets. Defaults to 5
        matching_on_corres (bool, optional): What type of point matching to use. Defaults to True

    Returns:
        pd.DataFrame: Contains valuable metrics on the registration process

    Modifications:
        Removed "matching on all" doesn't make sense the way it was coded, you need a distance cutoff otherwise the matching can be very off, besides one target point can be matched to multiple source points.
    """

    # format paths
    source_path = f"{path_dic[source_pocket]}"
    target_path = f"{path_dic[target_pocket]}"

    # parse ply files
    df_source = ply_parser_hull(ply_file=source_path)
    df_target = ply_parser_hull(ply_file=target_path)

    target_surf_normals = df_target.surf_normals.to_numpy()
    source_surf_normals = df_source.surf_normals.to_numpy()

    # fetching the input feats from the pocket DataFrame
    source_input_si     = df_pockets[df_pockets['pocket'] == source_pocket]['input_si'].to_numpy()
    source_surf_charge  = df_pockets[df_pockets['pocket'] == source_pocket]['input_charge'].to_numpy()
    source_surf_hphob   = df_pockets[df_pockets['pocket'] == source_pocket]['input_hphob'].to_numpy()
    source_surf_hbonds  = df_pockets[df_pockets['pocket'] == source_pocket]['input_hbonds'].to_numpy()

    target_input_si     = df_pockets[df_pockets['pocket'] == target_pocket]['input_si'].to_numpy()       
    target_surf_charge  = df_pockets[df_pockets['pocket'] == target_pocket]['input_charge'].to_numpy()       
    target_surf_hphob   = df_pockets[df_pockets['pocket'] == target_pocket]['input_hphob'].to_numpy()       
    target_surf_hbonds  = df_pockets[df_pockets['pocket'] == target_pocket]['input_hbonds'].to_numpy()         

    feats_source = np.stack(( source_surf_charge, source_surf_hphob, source_surf_hbonds, source_input_si), axis=1)
    feats_target = np.stack(( target_surf_charge, target_surf_hphob, target_surf_hbonds, target_input_si), axis=1)

    # get the features for the alignment, not the calculated ones
    source_feats       = o3d.pipelines.registration.Feature()
    source_feats.data  = feats_source.T

    target_feats       = o3d.pipelines.registration.Feature()
    target_feats.data  = feats_target.T

    ### --- preprocessing --- 
    
    # substract the center of mass of each PC from itself
    source = o3d.io.read_point_cloud(source_path, format='ply')
    target = o3d.io.read_point_cloud(target_path, format='ply')

    # remove center of mass
    source_corr = o3d.geometry.PointCloud()
    source_corr.points  = o3d.utility.Vector3dVector(np.asarray(source.points) - np.mean(np.asarray(source.points), axis=0))
    source_corr.normals = o3d.utility.Vector3dVector(source_surf_normals)

    target_corr = o3d.geometry.PointCloud()
    target_corr.points  = o3d.utility.Vector3dVector(np.asarray(target.points) - np.mean(np.asarray(target.points), axis=0))
    target_corr.normals = o3d.utility.Vector3dVector(target_surf_normals)

    # storing the array form of the source and target
    source_points = np.asarray(source_corr.points)
    target_points = np.asarray(target_corr.points)
    size_source = len(source_points)
    size_target = len(target_points)

    ### --- running the registration ---

    result_ransac, result_icp = best_registration(source_pcd=source_corr,
            target_pcd=target_corr,
            source_feats=source_feats,
            target_feats=target_feats,
            n = n_regs
    )

    # storing the performances and attributes of the results
    ## from the ransac run
    ransac_rmse = result_ransac.inlier_rmse
    ransac_fitness = result_ransac.fitness

    ## from the icp run
    icp_rmse = result_icp.inlier_rmse
    icp_fitness = result_icp.fitness
    
    # export the o3d registration object
    dic = {}
    dic['ransac_test'] = {'corres_set':[np.asarray(result_ransac.correspondence_set)],
                          'transformation':[result_ransac.transformation],
                          'fitness': [result_ransac.fitness],
                          'rmse': [result_ransac.inlier_rmse]}

    dic['icp_test'] = {'corres_set':[np.asarray(result_icp.correspondence_set)],
                        'transformation':[result_icp.transformation],
                        'fitness': [result_icp.fitness],
                        'rmse': [result_icp.inlier_rmse]}
    
    ## a dir for the npy files
    npy_out = f'{output}/o3d_objects_npy/'
    if not os.path.exists(npy_out):
        os.makedirs(npy_out)
    
    ## save the objects
    np.save(f"{npy_out}/{source_pocket}_{target_pocket}_input-feats.npy", dic, allow_pickle=True)


    # creating a dataframe that includes all informations we want to analyse
    output_df = pd.DataFrame({
        'source_pocket':[source_pocket], 
        'target_pocket':[target_pocket], 
        'ransac_rmse':[ransac_rmse], 
        'ransac_fitness':[ransac_fitness],
        'ransac_nb_corres':[len(np.asarray(result_ransac.correspondence_set))],
        'icp_rmse':[icp_rmse], 
        'icp_fitness':[icp_fitness], 
        'icp_nb_corres':[len(np.asarray(result_icp.correspondence_set))], 
        'size_source':[size_source], 
        'size_target':[size_target], 
    })

    return output_df

def registrate_all_pockets(n_regs:int, path_dic:dict, df_pockets:pd.DataFrame, output:str) -> pd.DataFrame:
    """Registrates all pockets

    Args:
        n_regs (int): See doc global_reg_pipeline()
    """
    # container
    vessel = []

    # loop over every pocket in the training set (about 8min)
    for n in tqdm(range(len(path_dic.keys())), desc="aligning fibril sites"):
        for m in range(n+1, len(list(path_dic.keys()))):

            source = list(path_dic.keys())[n]
            target = list(path_dic.keys())[m]
           
            result = global_reg_pipeline(source_pocket=source, target_pocket=target, path_dic=path_dic, df_pockets=df_pockets, n_regs = n_regs, output=output)
            vessel.append(result)
    
    # export a dataframe with all results
    df = pd.concat(vessel).reset_index(drop=True)
    df.to_csv(os.path.join(output, str(datetime.date.today())+ "_all_sites_alignment_results.csv"))

    return None