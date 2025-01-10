#!/bin/bash

#SBATCH --job-name=picmusDENOvitro
#SBATCH --qos=short
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=30
#SBATCH --output=picmusDENOvitro.out
#SBATCH --error=picmusDENOvitro.err
#SBATCH --time=60

export CONFIG=/home/yzhang2018@ec-nantes.fr/DRUSvar/DRUSvar/configs/imagenet_256_3c.yml
export MODEL_PATH=vitro
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

python -u /home/yzhang2018@ec-nantes.fr/DRUSvar/DRUSvar/main.py --ni --config $CONFIG  --doc $MODEL_PATH  --ckpt model006000.pt --matlab_path $MATLAB_PATH --deg DENOvitro

