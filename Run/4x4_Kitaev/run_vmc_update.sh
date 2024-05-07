#!/bin/bash 
#SBATCH --job-name=Kitaev-vmc-update-Jx=1-Jy=1e-2-Jz=1
##SBATCH --mem=20gb 
#SBATCH --partition=hebhcnormal01 
#SBATCH --time=99-00:00:00 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=64 
#SBATCH --cpus-per-task=1
#SBATCH --error=/public/home/xyliu/Kitaev/Run/4x4_Kitaev/log/error.err
#SBATCH --output=/public/home/xyliu/Kitaev/Run/4x4_Kitaev/log/Kitaev-vmc-update-Jx=1-Jy=1e-2-Jz=1-D=8.log
mpirun -np 64 /public/home/xyliu/Kitaev/build/Kitaev_vmc_update /public/home/xyliu/Kitaev/Run/4x4_Kitaev/vmc_update_params.json
