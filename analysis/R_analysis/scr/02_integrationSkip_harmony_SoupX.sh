#!/bin/bash
#SBATCH --job-name=scProcess
#SBATCH --account pedrini.edoardo
#SBATCH --mem=64GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=6  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="02_integrationSkip_harmony_SoupX.err"
#SBATCH --output="02_integrationSkip_harmony_SoupX.out"

echo "my job strart now" > 02_integrationSkip_harmony_SoupX.log;

date >> 02_integrationSkip_harmony_SoupX.log;

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_R4.2_recover2;

Rscript /home/pedrini.edoardo/scr/project_edoardo/240531_scRNAseq_MSSpinal_Absinta/analysis/R_analysis/scr/02_integrationSkip_harmony_SoupX.R

date >> 02_integrationSkip_harmony_SoupX.log;
echo "all done!!" >> 02_integrationSkip_harmony_SoupX.log