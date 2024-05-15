#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=8gb
#SBATCH --partition=Centaurus

module purge
module load anaconda3

### Define a shortcut to the annotation file
annot=/users/jaileru/final/Xenopus_tropicalis_annotation.gtf

### Activate the htseq environment
source activate htseq

### Run HTSeq-count to get the read counts
htseq-count -f bam -t gene -i gene_id  control_1.bam  $annot > ./counts-files/control_1.htseq.out
