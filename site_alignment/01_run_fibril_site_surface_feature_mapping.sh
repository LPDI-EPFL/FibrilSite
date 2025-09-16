#!/bin/bash

# Activate "fibrilsite" environment

python ./scripts/fibril_site_input_feats_mapper.py --site_src_dir "../sel_fibril_sites" \
                                                   --input_feats_src "../data_preparation/04b-precomputation_12A/precomputation"  \
                                                   --output_folder "."