#!/bin/bash

module load R/3.6.3
export SOURCETRACKER_PATH=/users/0/gould209/.local/bin/sourcetracker
Rscript sourcetracker_for_qiime.r -i study.txt -m study.txt -o ./output_dir
Rscript sourcetracker_for_qiime.r -i study.txt -m study.txt -o ./output_dir
