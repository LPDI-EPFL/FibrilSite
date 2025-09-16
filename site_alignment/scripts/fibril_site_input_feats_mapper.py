# this script is map the defined sites points to their calculated surface features 

# libraries
import os, glob, datetime, argparse
import numpy as np
import pandas as pd
from tqdm import tqdm

# functions
def create_parser():
    """
    Create a CLI parser.
    :return: the parser object.
    """
    parse = argparse.ArgumentParser()
    parse.add_argument("--site_src_dir",     type=str, nargs=1, help="Path to folder containing the defined fibril sites.")
    parse.add_argument("--input_feats_src",  type=str, nargs=1, help="Path to folder contaning MaSIF calculated surface features.")
    parse.add_argument("--output_folder",    type=str, nargs=1, help="Path to root folder for exporting.")
    return parse

def parse_args(parser):
    """
    Parse the arguments of a parser object.
    :param parser: the parser object.
    :return: the specified arguments
    """
    args = parser.parse_args()
    return args

## Parse the isolated pockets
def site_parser(site_src_dir:str) -> pd.DataFrame:
    """
    This function is to parse the defined sites point clouds

    site_src_dir : path to the folder containing the defined fibril sites
    """
    # container for storing the parsed sites
    vessel01 = [] 

    for pocket in tqdm(sorted(glob.iglob(os.path.join(site_src_dir, "*", "*.csv"))), desc='parsing defined sites'):
        if 'refined' in pocket or 'isolate' in pocket:
            df_temp01 = pd.read_csv(pocket, index_col=0)
            df_temp01.insert(0, 'fibril', os.path.basename(pocket).replace('.csv','').split('_')[1])
            df_temp01.insert(1, 'pocket_id', "_".join(os.path.basename(pocket).replace('.csv','').split('_')[2:-1]))
            df_temp01.insert(2, 'isolation', os.path.basename(pocket).replace('.csv','').split('_')[-1])
            vessel01.append(df_temp01)

    # add all these parsed pockets to a dataframe
    df_pockets_crude = pd.concat(vessel01).reset_index(drop=1)
    
    return df_pockets_crude

## Remove the isolate points for the pockets that underwent a round of refinement 
def pocket_cleaner(df_pockets_crude:pd.DataFrame, output:str) -> pd.DataFrame:
    """
    This function is to remove the isolate points for the pockets that underwent a round of refinement 

    df_pockets_crude : dataframe containing the the parsed sites info
    output           : path to export folder
    """

    # container to contain the filtered points info
    vessel02 = [] 

    for g01 in tqdm(df_pockets_crude.groupby(by='pocket_id'), desc='cleaning'):
        if len(set(g01[1].isolation)) < 2 :
            vessel02.append(g01[1])
        else:
            vessel02.append(g01[1][g01[1].isolation == 'refined'])

    df_pockets = pd.concat(vessel02).reset_index(drop=1)
    df_pockets.to_csv(os.path.join(output, str(datetime.date.today())+"_all_sites_parsed.csv"))

    return df_pockets

## Map surface features to site points
def feats_mapper(df_pockets:pd.DataFrame, input_feats_paths:list, output:str):
    """
    This function is to parse the MaSIF calculated input features and map them to the defined points

    *surface feature parsing* 

    MaSIF files are organzed in the manner of MaSIF point index that is the center of the patch, 
    the 200 points within the patch and their computed surface features. 
    We map the defined site points to their corresponding features based on the MaSIF indexing.
    
    df_pockets        : dataframe containing the info for the parsed site points
    input_feats_paths : list containing paths to the MaSIF calculated surface features
    """
    # container
    vessel03 = [] 

    for g02 in tqdm(df_pockets.groupby(by='fibril'), desc='getting feats'):
    
        # load the input feat files
        input_feats = np.load([f for f in input_feats_paths if g02[0] in f if os.path.basename(f) == "p1_input_feat.npy"][0])[:,0,:]

        # add the values
        df_temp02 = g02[1]
        df_temp02.insert(df_temp02.shape[1], 'input_si',     [input_feats[x][0] for x in df_temp02.MaSIF_index])
        df_temp02.insert(df_temp02.shape[1], 'input_charge', [input_feats[x][3] for x in df_temp02.MaSIF_index])
        df_temp02.insert(df_temp02.shape[1], 'input_hphob',  [input_feats[x][4] for x in df_temp02.MaSIF_index])
        df_temp02.insert(df_temp02.shape[1], 'input_hbonds', [input_feats[x][2] for x in df_temp02.MaSIF_index])

        vessel03.append(df_temp02)

    df_all_feats = pd.concat(vessel03).reset_index(drop=1)
    df_all_feats.to_csv(os.path.join(output, str(datetime.date.today())+"_all_sites_input_feats.csv"))

    return None

def main():
    """ Main execution point """
    
    # Parse arguments
    args = parse_args(create_parser())
    site_src_dir    = os.path.abspath(args.site_src_dir[0])
    input_feats_src = os.path.abspath(args.input_feats_src[0])
    output_root     = os.path.abspath(args.output_folder[0])

    # make output folder
    output = os.path.join(output_root, str(datetime.date.today())+"_sites_parsed_info/")
    os.makedirs(output, exist_ok=0)

    # paths to the input feats folder
    input_feats_paths = [f.strip() for f in glob.iglob(os.path.join(input_feats_src, "*", "*"))]

    # Parse the isolated pockets
    df_pockets_crude = site_parser(site_src_dir=site_src_dir)

    # Remove the isolate points for the pockets that underwent a round of refinement 
    df_pockets = pocket_cleaner(df_pockets_crude=df_pockets_crude, output=output)

    # Map the surface features to the identified site points
    feats_mapper(df_pockets=df_pockets, input_feats_paths=input_feats_paths, output=output)

    return None

if __name__ == "__main__":
    main()