#!/bin/bash
#SBATCH --job-name=scProcess
#SBATCH --account pedrini.edoardo
#SBATCH --mem=128GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=8  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="01_apply_fixed_thresholds_00500_07000_05_SoupX.err"
#SBATCH --output="01_apply_fixed_thresholds_00500_07000_05_SoupX.out"

echo "my job strart now" > 01_apply_fixed_thresholds_00500_07000_05_SoupX.log;

date >> 01_apply_fixed_thresholds_00500_07000_05_SoupX.log;

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_R4.2_recover2;

Rscript /beegfs/scratch/ric.cosr/pedrini.edoardo/project_edoardo/240531_scRNAseq_MSSpinal_Absinta/analysis/R_analysis/scr/01_apply_fixed_thresholds_00500_07000_05_SoupX.R

date >> 01_apply_fixed_thresholds_00500_07000_05_SoupX.log;
echo "all done!!" >> 01_apply_fixed_thresholds_00500_07000_05_SoupX.log