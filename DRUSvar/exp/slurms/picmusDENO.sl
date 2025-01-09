#!/bin/bash

#SBATCH --job-name=picmusDENO
#SBATCH --qos=short
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=30
#SBATCH --output=picmusDENO.out
#SBATCH --error=picmusDENO.err
#SBATCH --time=200

export CONFIG=/home/yzhang2018@ec-nantes.fr/DRUSvar/DRUSvar/configs/imagenet_256_1c.yml
export MODEL_PATH=both1c
export MATLAB_PATH=/home/yzhang2018@ec-nantes.fr/DRUSvar/MATLABfiles/ 
echo CONFIG
echo MODEL_PATH
echo MATLAB_PATH

# activate micromamba
source ~/.bashrc
micromamba activate ddrm


# launch python script
# python -c "from PIL import Image; print('ok')" #no problem
echo $LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/micromamba/yzhang2018@ec-nantes.fr/envs/ddrm/lib/
echo $LD_LIBRARY_PATH

python -u /home/yzhang2018@ec-nantes.fr/DRUSvar/DRUSvar/main.py --ni --config $CONFIG  --doc $MODEL_PATH  --ckpt model004000.pt --matlab_path $MATLAB_PATH --deg DENO

