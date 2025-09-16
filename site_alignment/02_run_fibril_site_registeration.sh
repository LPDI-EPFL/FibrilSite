#!/bin/bash

# Activate "alignsites" environment

python ./scripts/fibril_site_registeration.py --info_file "./2025-04-08_sites_parsed_info/2025-04-08_all_sites_input_feats.csv" \
                                              --sites_folder "../sel_fibril_sites/"  \
                                              --output_folder "."