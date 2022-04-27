#!/bin/bash
#SBATCH --account=
#SBATCH --cpus-per-task=10

#load modules
ml python/3.7-2019.10
source /users/PAS1475/zhenyuwu/virtualenv/caper/bin/activate 
ml java/11.0.8
ml singularity/current

cd /HOME/$working_directory
#encode syntax
caper run \
/users/PAS1475/zhenyuwu/Software/chip-seq-pipeline2/chip.wdl \
-i  1_ENCODE_input json.json \
--singularity

find . -name "*metadata.json" | xargs croo 