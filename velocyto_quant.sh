#!/bin/sh
# LSF directivesy
#BSUB -J "velocytoquanti"		# name of job
#BSUB -e loom24229.err                  # error file path, where % will be replaced with job ID
#BSUB -o loom24229.out                  # output file path, where % will be replaced with job ID
#BSUB -q verylong                     # queue 
##BSUB -q night
##BSUB -B                             # start of the job email
#BSUB -N                              # end of the job enmail
#BSUB -u julia.ruehle@dkfz-heidelberg.de   # user email
#BSUB -W 180:00                       # requested time
#BSUB -n 4                            # one CPU
#BSUB -M 40GB                        # memory limit
#BSUB -R "rusage[mem=40GB]"
#BSUB -P genanno


module load anaconda3/2019.07
eval "$(conda shell.bash hook)" 
conda activate velocyto

RUN_PATH='/dkfz/groups/OE0540/users/ruehle/rnavelocity/velocyto'
cd $RUN_PATH

velocyto run-smartseq2 -o /dkfz/groups/OE0540/users/ruehle/rnavelocity/velocyto/looms -e /icgc/dkfzlsdf/analysis/B260/projects/HipSci/openAccess/endoderm_differentation/data_raw/run_24229/star/*/*2pass.Aligned.sortedByCoord.dedup.bam  /dkfz/groups/OE0540/users/ruehle/rnavelocity/pre_kallisto/velocity_index/Homo_sapiens.GRCh37.75.gtf  


