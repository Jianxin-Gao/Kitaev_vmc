#!/bin/bash 
#SBATCH --job-name=Kitaev-simple-update-Jx=1-Jy=1-Jz=1
#SBATCH --mem=1gb 
#SBATCH --partition=hebhcnormal01 
#SBATCH --time=99-00:00:00 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=1 
#SBATCH --output=/public/home/xyliu/Kitaev/Run/4x4_Kitaev/log/Kitaev-simple-update-Jx=1-Jy=1-Jz=1-TruncErr=1e-7-Dmin=4-Dmax=4-Tau=0.2-Step=100.log 
/public/home/xyliu/Kitaev/build/Kitaev_simple_update /public/home/xyliu/Kitaev/Run/4x4_Kitaev/simple_update_params.json
