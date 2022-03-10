#65;6203;1c#!/bin/bash
#SBATCH --job-name=MHIRX1
##Select the desired partition:
#SBATCH --partition=serial-top
##Define the màximum number of cores you need to use:
##Keep in mind that most of matlab functions use internal multithreading.
##You probably need to send your matlab jobs to parallel queues if you use
##a loot of parallelizable functions.
#SBATCH -c 1
##Define the màximun memory you will use:
#SBATCH --mem=120G
##Define the maximum time your job will run.
##With format Days-Hours (1day-0hours)
#SBATCH --time=1-00
##Tell the queue manager to not consider the number of cpus in use from
##other jobs.
#SBATCH --oversubscribe
echo "Start:" `date`
echo "I ran on:"
cd $SLURM_SUBMIT_DIR
echo $SLURM_NODELIST
##load the desired matlab version
module load matlab/R2017b
##The matlab script should be called without the .m extension
matlab -nodisplay -r "mainMHIReducedX(4,3,1000,[10:40:1000 1200:200:5000],12,1e-2,0.2,103); quit"
echo "End:" `date`
