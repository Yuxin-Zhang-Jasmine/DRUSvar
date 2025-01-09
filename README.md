# Diffusion Reconstruction of Ultrasound Images with Informative Uncertainty

The structure of the repository is:
```bash
├── <DRUS-v2>  
│   ├── <configs> 
│   ├── <exp>  
│   │   ├──<image_samples>  
│   │   ├──<logs>             # the fined-tuned diffusion models, can be downloaded from:
│   │   │   ├──<vitro>        # https://uncloud.univ-nantes.fr/index.php/s/BRbDYsq2CNjJgoR
│   │   │   ├──<vivo>         # https://uncloud.univ-nantes.fr/index.php/s/FJAtmeN6QXaj6yD
│   │   │   ├──<CAROTIDcross> # https://uncloud.univ-nantes.fr/index.php/s/e53P2LboafkzimH
│   │   ├──<slurms>           # the job example/template if use HPC
│   ├── <functions>
│   ├── <guided_diffusion>
│   ├── <runners>  
│   ├── main.py
│   │ 
├── <Observation_SVDresults> 
│   ├── <Observation> 
│   │   ├──<01_simulation>    # the simulated fetus measurements (for discussion)
│   │   ├──<02_picmus>        # PICMUS measurements (for experiments)
│   ├── <SVD>                 # can be downloaded from: 
│   │   ├──<01_simulation>    # https://drive.google.com/drive/folders/1mCMaL6OR9mLvoRaQxmeWGLnUADYU32KF?usp=sharing
│   │   ├──<02_picmus>        # https://drive.google.com/drive/folders/10KwoH5G-s8Gk_aCj7WxTZ_L3596u44dI?usp=sharing
│   ├──compute_discussionSimu_Hty.m  # the script for computing the fetus measurements
│   ├──compute_mainDRUS_svdBH.m      # the script for computing the PICMUS DRUS SVD(BH)
│   ├──compute_mainDRUS_By.m         # the script for computing the PICMUS DRUS measurements
│   ├──compute_mainDeno_By.m         # the script for computing the PICMUS Deno measurements
├── environment.yml 
```
### Fine-tuning
The ultrasound datasets and the code for fine-tuning can be found in [this repository](https://gitlab.univ-nantes.fr/zhang-y-7/guided-diffusion-us).

### Results in the paper
All the results and the scripts for showing the results is [here](https://uncloud.univ-nantes.fr/index.php/s/PkWaC3NDF7ocE95).

### Reproduction steps
1. compute the measurements and the SVD results of the inverse problem model matrix using the scripts in `<Observation_SVDresults>`
2. edit `DRUS-v2/runners/diffusion.py` for a specific task
3. do restoration on HPC using slurm or as follows:
```
python main.py --ni --config {CONFIG.yml} --doc {MODELFOLDER} --ckpt {CKPT} --timesteps {STEPS} --deg {DEG} --image_folder {SAVEDIR} --matlab_path {MATLABPATH}
```
where
- `CONFIG` is the name of the config file (see `DRUS-v2/configs/`), including hyperparameters such as batch size and network architectures.
- `MODELFOLDER` is the name of the folder saving the diffusion model checkpoints 
- `CKPT` is the name of the selected diffusion model checkpoint
- `STEPS` controls how many timesteps (less than 1000) used in the process. (e.g. 50) 
- `DEG` is the degradation type, such as DRUS or Deno 
- `SAVEDIR` is name of the directory saving restored images
- `MATLABPATH` is the path of the folder `<Observation_SVDresults/>`. 

For example
```
python main.py --ni --config imagenet_256_3c.yml --doc vitro  --ckpt model006000.pt   --timesteps 50 --deg DRUS --image_folder vitroDRUS --matlab_path /home/.../Observation_SVDresults/
```

## References and Acknowledgements
```
@misc{zhang2023diffusion_v2,
    title={Diffusion Reconstruction of Ultrasound Images with Informative Uncertainty}, 
    author={Yuxin Zhang and Clément Huneau and Jérôme Idier and Diana Mateus},
    year={2023},
    eprint={2310.20618},
    archivePrefix={arXiv},
}
@misc{zhang2023diffusion_v1,
    title={Ultrasound Image Reconstruction with Denoising Diffusion Restoration Models}, 
    author={Yuxin Zhang and Clément Huneau and Jérôme Idier and Diana Mateus},
    year={2023},
    eprint={2307.15990},
    archivePrefix={arXiv},
}
@inproceedings{kawar2022denoising,
    title={Denoising Diffusion Restoration Models},
    author={Bahjat Kawar and Michael Elad and Stefano Ermon and Jiaming Song},
    booktitle={Advances in Neural Information Processing Systems},
    year={2022}
}
```

This implementation is based on / inspired by:
- [https://ddrm-ml.github.io/](https://ddrm-ml.github.io/) (the DDRM repo)
