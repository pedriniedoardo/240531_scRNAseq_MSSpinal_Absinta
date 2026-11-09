#!/bin/bash
#SBATCH --job-name=miloR_adj
#SBATCH --account pedrini.edoardo
#SBATCH --mem=64GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=2  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="08c_miloR_DA_adjusted_fullObject.err"
#SBATCH --output="08c_miloR_DA_adjusted_fullObject.out"

PROJECT_ROOT="/beegfs/scratch/ric.cosr/pedrini.edoardo/project_edoardo/240531_scRNAseq_MSSpinal_Absinta/analysis/R_analysis_R45_pixi"

echo "my job strart now" > 08c_miloR_DA_adjusted_fullObject.log;

date >> 08c_miloR_DA_adjusted_fullObject.log;

# this project uses pixi (not conda) to select the R executable/system deps;
# renv (activated from .Rprofile) then manages the package library from
# inside the R session — pixi must select R first, so cd into PROJECT_ROOT
# (where .Rprofile/renv.lock/pixi.toml live) before invoking Rscript via
# `pixi run`, regardless of where this job was submitted from.
# the relative paths in the R script (../../out/...) are resolved from PROJECT_ROOT.
cd "$PROJECT_ROOT" || exit 1

# requires the milo object of 08_miloR_DA_fullObject.R (out/object/analysis_R45_pixi/08_miloR/); only the test is rerun, so the job is short
# the R script creates its own output folders (08_miloR_adjusted)
pixi run Rscript "$PROJECT_ROOT/scr/08c_miloR_DA_adjusted_fullObject.R"
# exit code of the R script: without it a failed R script still ends with "all done!!" and a COMPLETED job
status=$?

date >> 08c_miloR_DA_adjusted_fullObject.log;
if [ $status -eq 0 ]; then
  echo "all done!!" >> 08c_miloR_DA_adjusted_fullObject.log
else
  echo "FAILED: the R script exited with status $status" >> 08c_miloR_DA_adjusted_fullObject.log
fi
exit $status
