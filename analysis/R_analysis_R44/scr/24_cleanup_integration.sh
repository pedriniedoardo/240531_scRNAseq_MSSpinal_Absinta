#!/bin/bash
#SBATCH --job-name=scProcess
#SBATCH --account pedrini.edoardo
#SBATCH --mem=64GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=6  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="24_cleanup_integration.err"
#SBATCH --output="24_cleanup_integration.out"

echo "my job strart now" > 24_cleanup_integration.log;

date >> 24_cleanup_integration.log;

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_R44;

Rscript /home/pedrini.edoardo/scr/project_edoardo/240531_scRNAseq_MSSpinal_Absinta/analysis/R_analysis_R44/scr/24_cleanup_integration.R

date >> 24_cleanup_integration.log;
echo "all done!!" >> 24_cleanup_integration.log