# Ultrasound Imaging Based on the Variance of a Diffusion Restoration Model

The structure of the repository is:
```bash
├── <DRUSvar>  
│   ├── <configs> 
│   ├── <exp>  
│   │   ├──<logs>             # The fined-tuned diffusion models, can be downloaded from: 
│   │   │   ├──<image1c>      # https://uncloud.univ-nantes.fr/index.php/s/SWamKLe3W5JTbSo   
│   │   │   ├──<vitro1c>         
│   │   │   ├──<both1c> 
│   │   ├──<slurms>           # the job examples &template if use HPC
│   ├── <functions>
│   ├── <guided_diffusion>
│   ├── <runners>  
│   ├── main.py
│   │ 
├── <MATLABfiles> 
│   ├── <numerical>           # for numerical experiments
│   │   ├──<SimulatedData>    # simulated data
│   │   ├──<Test_cysts>       # results
│   │   ├──<Test_scatterers>  # results
│   ├── <picmus>              # for real-data test 
│   │   ├──<DAS>              # the baseline
│   │   ├──<Observation>      # BH (measurements) for diffusion models
│   │   ├──<SVD>              # https://drive.google.com/drive/folders/10KwoH5G-s8Gk_aCj7WxTZ_L3596u44dI?usp=sharing
│   │   ├──<Test_picmus>      # results
│   ├── <src>                 # help resources
├── environment.yml 
```
### Fine-tuning
The ultrasound datasets and the code for fine-tuning can be found in [this repository](https://github.com/Yuxin-Zhang-Jasmine/guided-diffusion-ultrasound).


### Reproduction steps
on HPC using slurm, or as follows:
```
python main.py --ni --config {CONFIG.yml} --doc {MODELFOLDER} --ckpt {CKPT} --deg {DEG} --matlab_path {MATLABPATH}
```
where
- `CONFIG` is the name of the config file (see `DRUSvar/configs/`), including hyperparameters such as batch size and network architectures.
- `MODELFOLDER` is the name of the folder saving the diffusion model checkpoints 
- `CKPT` is the name of the selected diffusion model checkpoint
- `DEG` is the degradation type, such as DRUS, DENO, NumericalCysts, or NumericalScatterers
- `MATLABPATH` is the path of the folder `<MATLABfiles/>`. 

For example
```
python /home/.../DRUSvar/main.py --ni --config /home/.../DRUSvar/configs/imagenet_256_1c.yml --doc vitro1c  --ckpt model004000.pt  --deg DRUS --matlab_path /home/.../MATLABfiles/
```

## References and Acknowledgements
```
@inproceedings{DRUSvar,
    title={Ultrasound Imaging based on the Variance of a Diffusion Restoration Model},
    author={Zhang, Yuxin and Huneau, Cl{\'e}ment and Idier, J{\'e}r{\^o}me and Mateus, Diana},
    booktitle={EUSIPCO},
    year={2024}
}
```

This implementation is based on / inspired by:
- [https://ddrm-ml.github.io/](https://ddrm-ml.github.io/) (the DDRM repo)

