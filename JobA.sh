#!/bin/bash
#SBATCH --job-name=C1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4000m
#SBATCH --time=0:05:00
##SBATCH --licenses=abaqus@slurmdb:25
#SBATCH --account=awaas1

myhost=`hostname`
echo "Hostname: $myhost"

if [[ $SLURM_JOB_NODELIST ]] ; then
 echo "SLURM_NPROCS=$SLURM_NPROCS"
 echo "SLURM_JOB_NUM_NODES=$SLURM_JOB_NUM_NODES"
fi

##module load abaqus/2019

##module load intel/19.1

module load matlab

cd /home/ekheng/clustertest/

#unset SLURM_GTIDS
#rm -rfv abq1.*

##abaqus job=cube-purematrix input=cube-purematrix	scratch=/gpfs/accounts/awaas_root/awaas1/ekheng/final \
##	interactive user=RJD-VUMAT-master-tension.f double=both

##abaqus job=FE3LMS input=/gpfs/accounts/awaas_root/awaas1/ekheng/RoyanTest/FE3LMS interactive \
##user=/gpfs/accounts/awaas_root/awaas1/ekheng/RoyanTest/RJD-VUMAT-master-tension3.f scratch=/gpfs/accounts/awaas_root/awaas1/ekheng/RoyanTest \
##double=both cpus=30 parallel=domain domains=30

##abaqus cae noGUI="/home/ekheng/ondemand/data/sys/myjobs/projects/default/33/abaqusMacros.py"

matlab -nodisplay -r record

sbatch --dependency=afterany:$SLURM_JOB_ID  /home/ekheng/clustertest/Job2.sh