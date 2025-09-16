#!/bin/bash

# Activate "alignsites" environment

python ./scripts/fibril_site_alignment_analysis.py --site_src_dir "../sel_fibril_sites" \
                                                   --fibrils_pdb_src "../pdbs" \
                                                   --reg_results_src "./2025-04-08_registration_outputs" \
                                                   --sites_info_csv "../run_alignment/2025-04-08_sites_parsed_info/2025-04-08_all_sites_input_feats.csv" \
                                                   --sites_align_csv "./2025-04-08_registration_outputs/2025-04-08_all_sites_alignment_results.csv" \
                                                   --SSmax 0.5 \
                                                   --Fdiff 0.6 \
                                                   --output_folder "."
