#!/bin/bash
#SBATCH --time=5:00
#SBATCH --partition=cpar
#SBATCH --constraint=k20

nvprof ./bin/cuda < inputdata.txt