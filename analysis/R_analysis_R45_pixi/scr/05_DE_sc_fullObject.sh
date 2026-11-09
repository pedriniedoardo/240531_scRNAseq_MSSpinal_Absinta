#!/bin/bash
#SBATCH --job-name=scDEsc
#SBATCH --account pedrini.edoardo
#SBATCH --mem=128GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=6  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="05_DE_sc_fullObject.err"
#SBATCH --output="05_DE_sc_fullObject.out"

PROJECT_ROOT="/beegfs/scratch/ric.cosr/pedrini.edoardo/project_edoardo/240531_scRNAseq_MSSpinal_Absinta/analysis/R_analysis_R45_pixi"

echo "my job strart now" > 05_DE_sc_fullObject.log;

date >> 05_DE_sc_fullObject.log;

# this project uses pixi (not conda) to select the R executable/system deps;
# renv (activated from .Rprofile) then manages the package library from
# inside the R session — pixi must select R first, so cd into PROJECT_ROOT
# (where .Rprofile/renv.lock/pixi.toml live) before invoking Rscript via
# `pixi run`, regardless of where this job was submitted from.
# the relative paths in the R script (../../out/...) are resolved from PROJECT_ROOT.
cd "$PROJECT_ROOT" || exit 1

# the R script does not create the output folders
mkdir -p ../../out/table/analysis_R45_pixi/05_DE_sc/
mkdir -p ../../out/plot/analysis_R45_pixi/05_DE_sc/

pixi run Rscript "$PROJECT_ROOT/scr/05_DE_sc_fullObject.R"

date >> 05_DE_sc_fullObject.log;
echo "all done!!" >> 05_DE_sc_fullObject.log
