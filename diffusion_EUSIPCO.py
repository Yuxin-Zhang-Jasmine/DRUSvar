import os
from math import sqrt, pi, ceil
import numpy as np
import torch
from functions.ckpt_util import download
from functions.denoising import efficient_generalized_steps
from guided_diffusion.script_util import create_model
import mat73
from scipy.io import savemat


def get_beta_schedule(beta_schedule, *, beta_start, beta_end, num_diffusion_timesteps):
    def sigmoid(x):
        return 1 / (np.exp(-x) + 1)

    if beta_schedule == "quad":
        betas = (
                np.linspace(
                    beta_start ** 0.5,
                    beta_end ** 0.5,
                    num_diffusion_timesteps,

                    dtype=np.float64,
                )
                ** 2
        )
    elif beta_schedule == "linear":
        betas = np.linspace(
            beta_start, beta_end, num_diffusion_timesteps, dtype=np.float64
        )
    elif beta_schedule == "const":
        betas = beta_end * np.ones(num_diffusion_timesteps, dtype=np.float64)
    elif beta_schedule == "jsd":  # 1/T, 1/(T-1), 1/(T-2), ..., 1
        betas = 1.0 / np.linspace(
            num_diffusion_timesteps, 1, num_diffusion_timesteps, dtype=np.float64
        )
    elif beta_schedule == "sigmoid":
        betas = np.linspace(-6, 6, num_diffusion_timesteps)
        betas = sigmoid(betas) * (beta_end - beta_start) + beta_start
    else:
        raise NotImplementedError(beta_schedule)
    assert betas.shape == (num_diffusion_timesteps,)
    return betas


class Diffusion(object):
    def __init__(self, args, config, device=None):
        self.args = args
        self.config = config
        if device is None:
            device = (
                torch.device("cuda")
                if torch.cuda.is_available()
                else torch.device("cpu")
            )
        self.device = device

        self.model_var_type = config.model.var_type
        betas = get_beta_schedule(
            beta_schedule=config.diffusion.beta_schedule,
            beta_start=config.diffusion.beta_start,
            beta_end=config.diffusion.beta_end,
            num_diffusion_timesteps=config.diffusion.num_diffusion_timesteps,
        )
        betas = self.betas = torch.from_numpy(betas).float().to(self.device)
        self.num_timesteps = betas.shape[0]

        alphas = 1.0 - betas
        alphas_cumprod = alphas.cumprod(dim=0)
        alphas_cumprod_prev = torch.cat(
            [torch.ones(1).to(device), alphas_cumprod[:-1]], dim=0
        )
        self.alphas_cumprod_prev = alphas_cumprod_prev
        posterior_variance = (
                betas * (1.0 - alphas_cumprod_prev) / (1.0 - alphas_cumprod)
        )
        if self.model_var_type == "fixedlarge":
            self.logvar = betas.log()
            # torch.cat(
            # [posterior_variance[1:2], betas[1:]], dim=0).log()
        elif self.model_var_type == "fixedsmall":
            self.logvar = posterior_variance.clamp(min=1e-20).log()

    def sample(self):
        cls_fn = None

        config_dict = vars(self.config.model)
        model = create_model(**config_dict)
        if self.config.model.use_fp16:
            model.convert_to_fp16()
        ckpt = os.path.join(self.args.log_path, self.args.ckpt)
        print('Path of the current ckpt: ' + ckpt)
        if not os.path.exists(ckpt):
            print('The model does not exist, downloading an Imagenet 3c ckpt...')
            download(
                'https://openaipublic.blob.core.windows.net/diffusion/jul-2021/%dx%d_diffusion_uncond.pt' % (
                    self.config.data.image_size, self.config.data.image_size), ckpt)
        model.load_state_dict(torch.load(ckpt, map_location=self.device))
        model.to(self.device)
        model.eval()
        model = torch.nn.DataParallel(model)
        self.sample_sequence(model, cls_fn)

    def sample_sequence(self, model, cls_fn=None):
        args, config = self.args, self.config
        print("data channels : " + str(self.config.data.channels))
        print("model in_channels : " + str(self.config.model.in_channels))
        print('The corresponding MATLAB path: ' + self.args.matlab_path)

        # ** define the dasLst and the gammaLst **
        # ---EUSIPCO synthetic Test (gaussianMultiNoise)-----------------------------------------------
        speckletypes = ['m5']
        gammaLevels = [0.02, 0.05, 0.08, 0.11, 0.14, 0.17, 0.2, 0.23, 0.26, 0.29, 0.32, 0.35]  # contrast_speckle
        # gammaLevels = [0,0.003,0.006,0.009,0.012,0.015,0.018, 0.04, 0.06, 0.08, 0.1]  # resolution
        repeat = 20
        timestepsLst = [50] * len(speckletypes) * len(gammaLevels) * repeat  # 1 img, 20 repeat
        gammaLst = np.zeros((len(gammaLevels), repeat*len(speckletypes)))
        for i in range(repeat*len(speckletypes)):
            gammaLst[:, i] = gammaLevels
        gammaLst = list(gammaLst.reshape(len(gammaLevels)*repeat*len(speckletypes)))
        dasLst = []
        for s in range(len(speckletypes)):
            dasLst += ['signed_' + speckletypes[s] + '_' + str(self.args.imgIdx)] * repeat
        dasLst = dasLst * len(gammaLevels)

        print('len(dasLst) = ', len(dasLst))
        print('len(gammaLst) = ', len(gammaLst))
        print('len(timestepsLst) = ', len(timestepsLst))
        # # -------------------------------------------------------------------------------------------

        # ** get SVD results of the model matrix **
        print(f'Loading the degradation function/matrix  and it\'s svd (' + self.args.deg + ')')
        if self.args.deg == "us_VarianceAnalysis":
            from functions.svd_replacement import ultrasound0
            svdPath = self.args.matlab_path + 'SVD/02_picmus/'
            Up = torch.from_numpy(mat73.loadmat(svdPath + 'Ud.mat')['Ud'])
            lbdp = torch.from_numpy(mat73.loadmat(svdPath + 'Sigma.mat')['Sigma'])
            Vp = torch.from_numpy(mat73.loadmat(svdPath + 'Vd.mat')['Vd'])
            H_funcs = ultrasound0(config.data.channels, Up, lbdp, Vp, self.device)
        elif self.args.deg == "deno_VarianceAnalysis":
            from functions.svd_replacement import Denoising
            H_funcs = Denoising(config.data.channels, self.config.data.image_size, self.device)
        elif self.args.deg == 'deblur_aniso_VarianceAnalysis':
            from functions.svd_replacement import Deblurring2D
            sigma = self.args.sigma  # 1.2
            kerHalfLen = ceil(2 * sigma)
            pdf = lambda x: torch.exp(torch.Tensor([-0.5 * (x / sigma) ** 2]))
            kernel2 = torch.Tensor(
                [pdf(i) for i in range(-kerHalfLen, 0, 1)] + [pdf(i) for i in range(kerHalfLen + 1)]).to(self.device)
            # kernel2 = torch.Tensor([pdf(-3), pdf(-2), pdf(-1), pdf(0), pdf(1), pdf(2), pdf(3)]).to(self.device)
            sigma = self.args.sigma  # 1.2
            kerHalfLen = ceil(2 * sigma)
            pdf = lambda x: torch.cos(torch.Tensor([2*pi*0.5*x])) * torch.exp(torch.Tensor([-0.5 * (x/sigma)**2]))
            kernel1 = torch.Tensor([pdf(i) for i in range(-kerHalfLen, 0, 1)] + [pdf(i) for i in range(kerHalfLen+1)]).to(self.device)
            H_funcs = Deblurring2D(kernel1 / kernel1.sum(), kernel2 / kernel2.sum(), config.data.channels, self.config.data.image_size, self.device)
        else:
            print("ERROR: problem_model (--deg) type not supported")
            quit()


        idx_so_far = 0
        print(f'Start restoration')
        for _ in range(len(dasLst)):
            timesteps = timestepsLst[idx_so_far]
            dasSaveName = dasLst[idx_so_far] + '.mat'
            gamma = gammaLst[idx_so_far]

            # ** load the ground Truth x_orig **
            x_orig = torch.from_numpy(
                mat73.loadmat(self.args.matlab_path + 'Observation/01_simulation/Deno/' + dasSaveName)['x'])
            x_orig = (x_orig.view(1, 3, 256, 256)).to(self.device)
            if self.config.model.in_channels == 1:
                x_orig = torch.mean(x_orig, 1, True)
            # ** apply the inverse problem model (y_0 becomes with shape=(1,xxx))**
            y_0 = H_funcs.H(x_orig)
            amp = y_0.abs().max()
            if idx_so_far == 0:
                print('y_0.shape = ', y_0.shape)
            # **  load and truncate the additive noise to the same size as y_0 **
            additiveNoise = torch.from_numpy(
                mat73.loadmat(self.args.matlab_path+'Observation/01_simulation/Deno/additiveNoise_'+str(gamma)+'.mat')['n'])
            additiveNoise = (additiveNoise.view(1, 3*65536)).to(self.device)
            additiveNoise = additiveNoise[:, :y_0.size(1)]
            # ** add the additive noise **
            y_0 = y_0 + additiveNoise * amp

            # y_0 = y_0 + gamma * amp * torch.randn_like(y_0)

            # ** save the degraded image as .mat **
            savemat(os.path.join(self.args.image_folder, f"{idx_so_far + 1}.mat"), {'y_0': y_0[0].cpu().numpy()})
            
            # ===========================================================================================================
            # ** Begin DDRM **
            x = torch.randn(
                y_0.shape[0],
                config.data.channels,
                config.data.image_size,
                config.data.image_size,
                device=self.device,
            )
            with torch.no_grad():
                x, _ = self.sample_image(x, model, H_funcs, y_0, float(gamma*amp), last=False, cls_fn=cls_fn, timesteps=timesteps)
            # ** save the DDRM restored image as .mat **
            savemat(os.path.join(self.args.image_folder, f"{idx_so_far + 1}_{-1}.mat"),
                    {'x': x[-1][0].detach().cpu().numpy()})

            idx_so_far += y_0.shape[0]  # iterate multiple images
            print(f'Finish {idx_so_far}')

    def sample_image(self, x, model, H_funcs, y_0, gamma, last=False, cls_fn=None, classes=None, timesteps=50):
        skip = self.num_timesteps // timesteps
        seq = range(0, self.num_timesteps, skip)
        x = efficient_generalized_steps(x, seq, model, self.betas, H_funcs, y_0, gamma, etaB=self.args.etaB,
                                        etaA=self.args.eta, etaC=self.args.eta, cls_fn=cls_fn,
                                        classes=classes)
        if last:
            x = x[0][-1]
        return x
