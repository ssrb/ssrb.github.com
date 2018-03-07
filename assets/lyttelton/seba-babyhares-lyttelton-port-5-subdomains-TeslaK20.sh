#!/bin/bash
################################################################################
#################################### Microway Cluster Management Software (MCMS)
################################################################################
# Sample submission script for batch jobs. (SLURM)
#
# With the exception of interactive shell sessions, all jobs you run must be
# through a batch script such as this.
#
# YOU MUST EDIT THIS FILE TO CONTAIN YOUR COMMANDS.
#
# Once your script is complete, submit it with the sbatch utility.
# For example:
#  sbatch example-custom-job.sh
#
#----------------------------------------------------------------
#
# Using a job scheduler allows all users to share the compute cluster without
# slowing down each other's jobs. Without it, it's not possible to be certain
# two programs won't try to use the same CPUs, GPUs or coprocessors.
#
# A SLURM batch script is just a *nix shell script with additional
# information describing the desired cluster resources.
#
# The complete documentation is available here:
# http://slurm.schedmd.com/sbatch.html
#
################################################################################


# The following sections (starting with #SBATCH) inform the scheduler what
# capabilities and resources your job will require. All other lines beginning
# with '#' characters are comments.
#
# Add comments so you don't forget what your script is doing.


# Request the number of CPUs and compute nodes needed by your job.
# You must instruct the number of nodes and the number of cores.
#
# Request one processor core:
#   --nodes=1 --ntasks-per-node=1
#
# Request two compute nodes, each with ten processor cores:
#   --nodes=2 --ntasks-per-node=10
#
# Request all twenty processor cores on one compute node:
#   --nodes=1 --ntasks-per-node=20
#
#SBATCH --ntasks=5 --cpus-per-task=4


# Request number of GPUs and the type of GPU
#
# Request nodes with one Tesla K40 GPU
#   --gres=gpu:1
#   --constraint=K40
#
# Request nodes with two Tesla K20 GPUs
#   --gres=gpu:2
#   --constraint=K20
#
#SBATCH --gres=gpu:3
#SBATCH --constraint=K20


# Estimate how much memory you expect your job to use (in megabytes).
# Common values:
#     4GB    4096
#     8GB    8192
#    16GB   16384
#    32GB   32768
#    64GB   65536
#
#SBATCH --mem=32768


# Specify how long you expect your job to run. By default SLURM will kill jobs
# that over-run their reservation, so you need to make a realistic estimate.
# PLAY NICE!
#
# Request 1 day:
#   --time=24:00:00
#
# Request 1 hour:
#   --time=1:00:00
#
# Request 15 minutes:
#   --time=15:00
#
#SBATCH --time=10:00


# Set the job name
#SBATCH --job-name=LytteltonWaveSimulation


# If desired, you may set the filenames for your program's output.
# The %j variable will be replaced with the Job ID # of the job.
#SBATCH --error=LytteltonWaveSimulation.%j.output.errors
#SBATCH --output=LytteltonWaveSimulation.%j.output.log


# Edit this to set number of GPUs per node
gpus_per_node=2



# Print helpful start-up and shutdown messages:
#   * Start and Stop times
#   * Node counts and names
#   * CPU / Core counts and information
#   * GPU counts and information
#   * Memory information
#   * Hard Drive / SSD information
#
# Also, change to the working directory
source /mcms/core/slurm/scripts/libexec/slurm.jobstart_messages.sh



################################################################################
# THIS IS WHERE YOU SPECIFY THE JOBS/TASKS TO RUN
################################################################################

# customize MPI version and additional modules if necessary
module load intel
module load openmpi
module load cuda
module list 2>&1
echo; echo;


################
### Option 1 ###
################
# Run a parallel MPI command.
mpirun -c 5 /home/seba/work/ssrb.github.com/assets/lyttelton/BabyHares/babyhares_k20 /home/seba/work/ssrb.github.com/assets/lyttelton/BabyHares/lyttelton.mesh /home/seba/work/ssrb.github.com/assets/lyttelton/BabyHares/interface2.vids

################
### Option 2 ###
################
# If you're running a single command on a single compute node.
#/path/to/your/program --options

################################################################################
################################################################################
