#!/bin/bash

#SBATCH --job-name=theCurrentJob
#SBATCH --qos=short
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=30
#SBATCH --output=theCurrentJob.out
#SBATCH --error=theCurrentJob.err
#SBATCH --time=200

export config=        # path to DRUSvar/DRUSvar/configs/imagenet_256_1c.yml
export modelPath=     # vitro1c || both1c (trained on both vitro and vivo datasets) || image1c (trained only on ImageNet)
export MATLAB_PATH=   # path to DRUSvar/MATLABfiles/ 

# activate micromamba
source ~/.bashrc
micromamba activate ddrm

# debug an virtual environment library error
echo $LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/micromamba/yzhang2018@ec-nantes.fr/envs/ddrm/lib/
echo $LD_LIBRARY_PATH

# launch the python script
python -u pathToDRUSvar/DRUSvar/main.py 
--ni 
--config $config  
--doc $modelPath  
--matlab_path $MATLAB_PATH
--ckpt   # the name of the selected ckpt in the modelPath  
--deg    # DRUS || DENO || NumericalCysts || NumericalScatterers


