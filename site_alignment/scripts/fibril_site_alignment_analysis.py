# This script is to analyse the output of the fibril site alignments

# libraries
import os, glob, datetime, argparse
import pandas as pd
from fibrilsite.site_alignment import ply_parser_hull
from fibrilsite.site_alignment_analysis import *

# functions
def create_parser():
    """
    Create a CLI parser.
    :return: the parser object.
    """
    parse = argparse.ArgumentParser()
    parse.add_argument("--site_src_dir",     type=str,   nargs=1, help="Path to folder containing the defined fibril sites.")
    parse.add_argument("--fibrils_pdb_src",  type=str,   nargs=1, help="Path to folder containing fibrils pdb files.")
    parse.add_argument("--reg_results_src",  type=str,   nargs=1, help="Path to folder containing alignment results npy files.")
    parse.add_argument("--sites_info_csv",   type=str,   nargs=1, help="Path to CSV file containing fibril sites info.")
    parse.add_argument("--sites_align_csv",  type=str,   nargs=1, help="Path to CSV file containing fibril sites alignment results.")
    parse.add_argument("--SSmax",            type=float, nargs=1, help="Threshold for Site Surface overlap between site matches, recommended >= 0.5.")
    parse.add_argument("--Fdiff",            type=float, nargs=1, help="Threshold for Surface feature difference between site matches, recommended <= 0.6.")
    parse.add_argument("--output_folder",    type=str,   nargs=1, help="Path to root folder for exporting.")
    return parse

os.path.abspath('./2025-04-08_registration_outputs/2025-04-08_all_sites_alignment_results.csv')

def parse_args(parser):
    """
    Parse the arguments of a parser object.
    :param parser: the parser object.
    :return: the specified arguments
    """
    args = parser.parse_args()
    return args

def fibril_name_str(x):
    """ This function is to add the fibril name for the respective pocket """
    if x.split('_')[0] in ['PHF','CTEII']:
        return 'Tau'
    elif x.split('_')[0] in ['A53T', 'Pol', 'LF', 'G51D', 'MSA', 'E46K', 'H50Q', 'pY39']:
        return 'aSyn'
    else:
        return x.split('_')[0]

def add_size_metrics(df:pd.DataFrame) -> pd.DataFrame:
    """  """
    df['src2target_size_ratio']      = round(df['size_source'] / df['size_target'],2)
    df['target_matched_percent_icp'] = round(df['icp_nb_corres']/df['size_target'],2)
    return df

def add_pocket_matches(df:pd.DataFrame) -> pd.DataFrame:
    """  """
    vessel = []
    for item00 in df.pocket_pairs:
        sorted_items = sorted([item00[0],item00[2]])
        vessel.append((sorted_items[0], sorted_items[1]))
    
    # add to Dataframe
    df.insert(0, 'pocket_matches', vessel)

    # check if there are duplicated matches
    assert len(df.pocket_matches.to_list()) == len(list(set(df.pocket_matches.to_list()))), "Duplicated Matches Found"
    print("Unique Matches")

    return df

def main():
    """ Main execution point """

     # Parse arguments
    args = parse_args(create_parser())
    site_src         = os.path.abspath(args.site_src_dir[0])
    fibrils_pdb_src  = os.path.abspath(args.fibrils_pdb_src[0])
    reg_results_src  = os.path.abspath(args.reg_results_src[0])
    sites_info_csv   = os.path.abspath(args.sites_info_csv[0])
    sites_align_csv  = os.path.abspath(args.sites_align_csv[0])
    ssmax_thresh     = float(args.SSmax[0])
    fdiff_thresh     = float(args.Fdiff[0])
    output_root      = os.path.abspath(args.output_folder[0])

    # make main output folder
    main_output = os.path.join(output_root, str(datetime.date.today()) + "_site_alignment_analysis")
    os.makedirs(main_output, exist_ok=False)

    # make output for all vs all comparison
    allvsall_out = os.path.join(main_output, "all_vs_all")
    os.makedirs(allvsall_out, exist_ok=False)

    # make output for ex vivo aSyn to aSyn fibril sites comparison
    asyn2asyn_out = os.path.join(main_output, "aSyn_to_aSyn")
    os.makedirs(asyn2asyn_out, exist_ok=False)

    # make output for ex vivo aSyn to other amyloid fibrils sites comparison
    asyn2others_out = os.path.join(main_output, "aSyn_to_other_amyloids")
    os.makedirs(asyn2others_out, exist_ok=False)

    # make output for ex vivo to in vitro aSyn fibril sites comparison
    asyn2invitro_out = os.path.join(main_output, "aSyn_to_invitro_aSyn")
    os.makedirs(asyn2invitro_out, exist_ok=False)

    # make output for in vitro aSyn to other amyloid fibrils sites comparison
    asyninvitro2others_out = os.path.join(main_output, "invitro_aSyn_to_other_amyloids")
    os.makedirs(asyninvitro2others_out, exist_ok=False)

    # make output for aligned sites
    aligned_sites_output = os.path.join(main_output, "identified_matches_alignments")
    os.makedirs(aligned_sites_output, exist_ok=False)

    # ---------------------------------------------------------------------- #

    # Parse file paths
    # parse the output npy files after site registration
    reg_results_paths = [l.strip() for l in glob.iglob(os.path.join(reg_results_src , "o3d_objects_npy", "*.npy"))]
    
    # parse the fibril pdb file paths
    fibrils_paths = [p.strip() for p in glob.iglob(os.path.join(fibrils_pdb_src, "*.pdb"))]

    # parse the paths for the defined fibril sites ply files 
    ply_files = [p.strip() for p in glob.iglob(os.path.join(site_src, "*", "*.ply")) if "convex" in p]

    # ---------------------------------------------------------------------- #
    
    # get the paths to the defined site (pocket) files
    ## we need to get the refined site version if available, other wise to load the isolated site version

    # container
    src_pockets_paths = []

    for folder00 in glob.iglob(os.path.join(site_src, '*')):
        try:
            [p.strip() for p in glob.iglob(os.path.join(folder00, '*.xyz')) if 'refined' in p][0]
        except IndexError:
            src_pockets_paths.append([p.strip() for p in glob.iglob(os.path.join(folder00, '*.xyz')) if 'isolated' in p][0])
        else:
            src_pockets_paths.append([p.strip() for p in glob.iglob(os.path.join(folder00, '*.xyz')) if 'refined' in p][0])

    assert len(src_pockets_paths) == len(ply_files), "Defined sites xyz and ply files count mismatch"

    # ---------------------------------------------------------------------- #

    # Map fibril sites to respective fibril pdbs
    done = []
    fibril_pocket_map = {}

    for path01 in glob.iglob(os.path.join(site_src, '*')):
        p_name = os.path.basename(path01).split("_")
        if p_name[1] not in done :
            done.append(p_name[1])
            fibril_pocket_map[p_name[1]] = []
            fibril_pocket_map[p_name[1]].append(p_name[-1])
        else:
            fibril_pocket_map[p_name[1]].append(p_name[-1])

    assert len(fibril_pocket_map.keys()) == len(fibrils_paths), "Parsed fibrils and provided fibril pdb paths count mismatch"
    print("Mapped fibril sites: ", fibril_pocket_map)

    # load all fibril site points with their features
    df_all_info = pd.read_csv(sites_info_csv, index_col=0)

    # load all fibril site alignment results  
    df_input_all = pd.read_csv(sites_align_csv, index_col=0)

    # define the loaded site names
    pockets = list(set(df_all_info.pocket_id))
    assert len(pockets) == len(src_pockets_paths), "Parsed sites xyz files and site info mismatch"

    # add fibril source for source and target pockets
    df_input_all['source_pocket_fibril'] = df_input_all['source_pocket'].apply(lambda x: fibril_name_str(x))
    df_input_all['target_pocket_fibril'] = df_input_all['target_pocket'].apply(lambda x: fibril_name_str(x))

    # add the alternate metrics from the target pocket side
    df_input_all = add_size_metrics(df=df_input_all)

    # rename the fitness score columns
    df_input_all.rename(columns={"icp_fitness":"icp_fitness_source", "target_matched_percent_icp":"icp_fitness_target"}, inplace=True)

    # add the combined fitness score
    comb_fit_sc_vessel = []
    for idx00 in df_input_all.index:
        comb_fit_sc_vessel.append(round(max(df_input_all.at[idx00, "icp_fitness_source"], df_input_all.at[idx00, "icp_fitness_target"]),2))

    df_input_all.insert(df_input_all.shape[1], "combined_fitness_score", comb_fit_sc_vessel)

    # get rid of the ransac columns
    df_input_all = df_input_all[[
        'source_pocket', 'target_pocket',
        'icp_rmse', 'icp_fitness_source', 'icp_fitness_target', 'combined_fitness_score',
        'icp_nb_corres', 'size_source', 'size_target', 'src2target_size_ratio',
        'source_pocket_fibril', 'target_pocket_fibril']]

    # get the pocket pairs
    df_input_all.insert(0, 'pocket_pairs', [[s,sf,t,tf] for s,sf,t,tf in zip(df_input_all.source_pocket, df_input_all.source_pocket_fibril, df_input_all.target_pocket, df_input_all.target_pocket_fibril)])
    df_input_all.sort_values(by='pocket_pairs', inplace=True)
    df_input_all.reset_index(drop=True, inplace=True)

    # calculate the input feat diff
    df_input_all = calc_input_feat_diff(df=df_input_all, df_all_info=df_all_info, input_feats_npy=reg_results_paths, output=main_output, export=False)

    # add the pocket matches
    df_input_all = add_pocket_matches(df=df_input_all)

    # get the approximations
    df_input_all["combined_fitness_score"] = df_input_all["combined_fitness_score"].apply(lambda x: round(x, 2))
    df_input_all["icp_mean_input_diff"] = df_input_all["icp_mean_input_diff"].apply(lambda x: round(x, 2))

    # export 
    df_input_all.to_csv(os.path.join(main_output, str(datetime.date.today())+'_all_input_reg_pockets.csv'))

    # ---------------------------------------------------------------------- #
    
    # Find the all vs all matches
    # base the selections on the ICP mean input feats (Fdiff)
    # container
    allvsall_vessel = []

    for poc in tqdm(pockets):
        df_temp = df_input_all.copy()[(df_input_all.source_pocket == poc) | (df_input_all.target_pocket == poc)].sort_values(by="icp_mean_input_diff", ascending=True).reset_index(drop=True).head(5)
        allvsall_vessel.append(df_temp)
        df_temp.to_csv(os.path.join(allvsall_out, poc + "_top5_nghs.csv"))

    # put the selected nghs into a df
    df_allvsall_sel_ngh = pd.concat(allvsall_vessel).reset_index(drop=True)
    df_allvsall_sel_ngh.to_csv(os.path.join(allvsall_out, "allvsall_top5_nghs.csv"))

    # ---------------------------------------------------------------------- #

    # Compare sites among ex vivo aSyn fibrils
    # get the pocket names in the brain derived structures
    asyn_brain_pockets = [p for p in pockets if p.split("_")[0] in ["MSA", "LF"]]

    # get the information for the brain derived pockets
    df_input_b2b = df_input_all.copy()[(df_input_all.source_pocket.isin(asyn_brain_pockets)) & (df_input_all.target_pocket.isin(asyn_brain_pockets))].sort_values(by='combined_fitness_score', ascending=False).reset_index(drop=1)

    # export all aSyn brain pockets matches
    df_input_b2b.to_csv(os.path.join(asyn2asyn_out, 'all_asyn_brain_to_asyn_brain_matches.csv'))

    # get the matches that pass the simirity threshold
    # SSmax : combined_fitness_score >= 0.5
    # Fdiff : icp_mean_input_diff <= 0.6

    df_b2b_sel = df_input_b2b.copy()[(df_input_b2b.combined_fitness_score >= ssmax_thresh) & (df_input_b2b.icp_mean_input_diff <= fdiff_thresh)].sort_values(by="combined_fitness_score", ascending=False).reset_index(drop=True)
    print(f"Identified {df_b2b_sel.shape[0]} site matched among ex vivo aSyn fibril sites")

    # export
    df_b2b_sel.to_csv(os.path.join(asyn2asyn_out, "sel_asyn_brain_to_asyn_brain_matches.csv"))

    # align site point clouds
    align_site_pcd(df=df_b2b_sel, ply_files=ply_files, input_feats_npy=reg_results_paths, output=aligned_sites_output)

    # ---------------------------------------------------------------------- #

    # Compare sites between ex vivo fibrils of aSyn and other amyloid proteins
    # get the pocket names in the brain derived structures
    asyn_brain_and_other_amyloids_pockets = [p for p in pockets if p.split("_")[0] not in ['A53T', 'Pol', 'G51D', 'E46K', 'H50Q', 'pY39']]

    # get the information for the desired pockets
    df_input_b2others = df_input_all.copy()[
        (df_input_all.source_pocket.isin(asyn_brain_and_other_amyloids_pockets)) & 
        (df_input_all.target_pocket.isin(asyn_brain_and_other_amyloids_pockets))
        ].sort_values(by='combined_fitness_score', ascending=False).reset_index(drop=1)

    # export all aSyn brain pockets matches
    df_input_b2others.to_csv(os.path.join(asyn2others_out, 'all_asyn_brain_to_other_amyloids_matches.csv'))
    
    # get the matches that pass the simirity threshold
    df_b2o_sel = df_input_b2others.copy()[(df_input_b2others.combined_fitness_score >= ssmax_thresh) & (df_input_b2others.icp_mean_input_diff <= fdiff_thresh)].sort_values(by="combined_fitness_score", ascending=False).reset_index(drop=True)

    # make sure that the matches are not among ex vivo aSyn fibrils 
    df_b2o_sel["fibril_src_target_match"] = df_b2o_sel["source_pocket_fibril"] == df_b2o_sel["target_pocket_fibril"]
    df_b2o_sel = df_b2o_sel[df_b2o_sel.fibril_src_target_match == False]
    print(f"Identified {df_b2o_sel.shape[0]} site matched between sites from ex vivo fibrils of aSyn and other amyloids proteins")

    # export
    df_b2o_sel.to_csv(os.path.join(asyn2others_out, 'sel_asyn_brain_to_other_amyloids_matches.csv'))

    # ---------------------------------------------------------------------- #

    ### Compare sites between ex vivo and in vitro fibrils of aSyn

    # get the pocket names
    asyn_brain_and_other_invitro_pockets = [p for p in pockets if p.split("_")[0] in ['MSA', 'LF', 'A53T', 'Pol', 'G51D', 'E46K', 'H50Q', 'pY39']]

    # get the information for the desired pockets
    df_input_b2invitro = df_input_all.copy()[
        (df_input_all.source_pocket.isin(asyn_brain_and_other_invitro_pockets)) & 
        (df_input_all.target_pocket.isin(asyn_brain_and_other_invitro_pockets))
        ].sort_values(by='combined_fitness_score', ascending=False).reset_index(drop=1)

    # export all aSyn brain pockets matches
    df_input_b2invitro.to_csv(os.path.join(asyn2invitro_out, 'all_asyn_brain_to_invitro_matches.csv'))

    # get the matches that pass the simirity threshold
    df_b2invitro_sel = df_input_b2invitro.copy()[(df_input_b2invitro.combined_fitness_score >= ssmax_thresh) & (df_input_b2invitro.icp_mean_input_diff <= fdiff_thresh)].sort_values(by="combined_fitness_score", ascending=False).reset_index(drop=True)

    # ensure that the matches are with the ex vivo aSyn structures
    df_b2invitro_sel = df_b2invitro_sel[(df_b2invitro_sel.source_pocket.isin(asyn_brain_pockets)) | (df_b2invitro_sel.target_pocket.isin(asyn_brain_pockets))].reset_index(drop=True)

    # make sure that the matches are not among ex vivo aSyn fibrils 
    df_b2invitro_sel_rej = df_b2invitro_sel[(df_b2invitro_sel.source_pocket.isin(asyn_brain_pockets)) & (df_b2invitro_sel.target_pocket.isin(asyn_brain_pockets))].reset_index(drop=True)
    df_b2invitro_sel = df_b2invitro_sel[~df_b2invitro_sel.pocket_matches.isin(df_b2invitro_sel_rej.pocket_matches.to_list())].reset_index(drop=True)
    print(f"Identified {df_b2invitro_sel.shape[0]} site matched between sites from ex vivo and in vitro fibrils of aSyn")

    # export
    df_b2invitro_sel.to_csv(os.path.join(asyn2invitro_out, 'sel_asyn_brain_to_invitro_matches.csv'))

    # align site point clouds
    align_site_pcd(df=df_b2invitro_sel, ply_files=ply_files, input_feats_npy=reg_results_paths, output=aligned_sites_output)

    # ---------------------------------------------------------------------- #

    # Compare sites between in vitro fibrils of aSyn and ex vivo fibrils of other amyloids
    # get the pocket names
    asyn_invitro_and_other_exvivo_pockets = [p for p in pockets if p.split("_")[0] not in ['MSA', 'LF']]

    # get the information for the desired pockets
    df_input_invitro2others = df_input_all.copy()[
        (df_input_all.source_pocket.isin(asyn_invitro_and_other_exvivo_pockets)) & 
        (df_input_all.target_pocket.isin(asyn_invitro_and_other_exvivo_pockets))
        ].sort_values(by='combined_fitness_score', ascending=False).reset_index(drop=1)
    
    # make sure that the matches are not among ex vivo aSyn fibrils 
    df_input_invitro2others["fibril_src_target_match"] = df_input_invitro2others["source_pocket_fibril"] == df_input_invitro2others["target_pocket_fibril"]
    df_input_invitro2others = df_input_invitro2others[df_input_invitro2others.fibril_src_target_match == False].reset_index(drop=True)

    # export
    df_input_invitro2others.to_csv(os.path.join(asyninvitro2others_out, 'all_asyn_invitro_to_other_amyloids_matches.csv'))

    # get the matches that pass the simirity threshold
    df_invitro2others_sel = df_input_invitro2others.copy()[(df_input_invitro2others.combined_fitness_score >= ssmax_thresh) & (df_input_invitro2others.icp_mean_input_diff <= fdiff_thresh)].sort_values(by="combined_fitness_score", ascending=False).reset_index(drop=True)
    print(f"Identified {df_invitro2others_sel.shape[0]} site matched between sites from in vitro fibrils of aSyn and ex vivo fibrils of other amyloid proteins")

    # export
    df_invitro2others_sel.to_csv(os.path.join(asyninvitro2others_out, 'sel_asyn_invitro_to_other_amyloids_matches.csv'))

    # ---------------------------------------------------------------------- #

    ## Align Fibrils based on shared features for matched sites
    # align site source fibrils
    align_site_fibrils(
        align_out=aligned_sites_output, 
        src_pockets_paths=src_pockets_paths, 
        ply_files=ply_files, 
        fibrils_paths=fibrils_paths, 
        fibril_pocket_map=fibril_pocket_map, 
        reg_results_paths=reg_results_paths
        )
    
    return None

if __name__ == "__main__":
    main()






    









