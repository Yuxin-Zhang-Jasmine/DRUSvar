import os
from math import sqrt, pi, ceil
import numpy as np
import torch
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
        model.load_state_dict(torch.load(ckpt, map_location=self.device))
        model.to(self.device)
        model.eval()
        model = torch.nn.DataParallel(model)
        self.sample_sequence(model, cls_fn)

    def sample_sequence(self, model, cls_fn=None):
        args, config = self.args, self.config
        print("data channels : " + str(config.data.channels))
        print("model in_channels : " + str(config.model.in_channels))
        print('The corresponding MATLAB path: ' + args.matlab_path)

        # ** define the dasLst, the gammaLst, and the savepath
        if args.deg in ["DRUS", "DENO"]:
            repeat = 20
            dasLst = ['simu_reso', 'simu_cont', 'expe_reso', 'expe_cont', 'expe_cross', 'expe_long'] * repeat
            if args.deg == "DENO":
                gammaLst = [0.11, 0.12, 0.16, 0.025, 0.02, 0.05] * repeat
            elif args.deg == "DRUS":
                gammaLst = [30, 65, 90, 55, 12, 1.5] * repeat
            resultsPath = args.matlab_path + 'picmus/Test_picmus/results/' + args.deg + '/'
        elif args.deg in ["NumericalCysts", "NumericalScatterers"]:
            repeat = 20
            if args.deg == "NumericalCysts":
                numeType = 'cysts'
                gammaLevels = [0.02, 0.05, 0.08, 0.11, 0.14, 0.17, 0.2, 0.23, 0.26, 0.29, 0.32, 0.35]  # Cysts
            elif args.deg == "NumericalScatterers":
                numeType = 'scatterers'
                gammaLevels = [0.0, 0.003,0.006,0.009,0.012,0.015,0.018, 0.04, 0.06, 0.08, 0.1]  # Scatterers
            gammaLst = np.zeros((len(gammaLevels), repeat))
            for i in range(repeat):
                gammaLst[:, i] = gammaLevels
            gammaLst = list(gammaLst.reshape(len(gammaLevels) * repeat))
            dasLst = [numeType + '_mp_' + str(args.mIdx)] * repeat * len(gammaLevels)
            resultsPath = args.matlab_path + 'numerical/Test_' + numeType + '/results/sigma1.2/m' + str(args.mIdx) + '/'
        # to validate that the basic info. is correct
        print('len(dasLst) = ', len(dasLst))
        print('len(gammaLst) = ', len(gammaLst))
        print('Results will be saved to: ' + resultsPath)
        os.makedirs(resultsPath, exist_ok=True)

        # ** get SVD results of the model matrix **
        print(f'Loading the SVD of the degradation matrix (' + args.deg + ')')
        if args.deg == "DRUS":  # apply on the PICMUS dataset with the full SVD of model matrix BH
            from functions.svd_replacement import ultrasound0
            svdPath = args.matlab_path + 'picmus/SVD/'
            Up = torch.from_numpy(mat73.loadmat(svdPath + 'Ud.mat')['Ud'])
            lbdp = torch.from_numpy(mat73.loadmat(svdPath + 'Sigma.mat')['Sigma'])
            Vp = torch.from_numpy(mat73.loadmat(svdPath + 'Vd.mat')['Vd'])
            H_funcs = ultrasound0(config.data.channels, Up, lbdp, Vp, self.device)
        elif args.deg == "DENO":  # apply on the PICMUS dataset with simple denoising (BH \approx Identity)
            from functions.svd_replacement import Denoising
            H_funcs = Denoising(config.data.channels, config.data.image_size, self.device)
        elif args.deg in ["NumericalCysts", "NumericalScatterers"]:  # apply on the numerical dataset (Cysts and Scatterers) degraded with a 2D blurring kernel
            from functions.svd_replacement import Deblurring2D
            sigma = 1.2
            kerHalfLen = ceil(2 * sigma)
            pdf = lambda x: torch.exp(torch.Tensor([-0.5 * (x / sigma) ** 2]))
            kernel2 = torch.Tensor(
                [pdf(i) for i in range(-kerHalfLen, 0, 1)] + [pdf(i) for i in range(kerHalfLen + 1)]).to(self.device)
            pdf = lambda x: torch.cos(torch.Tensor([2*pi*0.5*x])) * torch.exp(torch.Tensor([-0.5 * (x/sigma)**2]))
            kernel1 = torch.Tensor([pdf(i) for i in range(-kerHalfLen, 0, 1)] + [pdf(i) for i in range(kerHalfLen+1)]).to(self.device)
            H_funcs = Deblurring2D(kernel1 / kernel1.sum(), kernel2 / kernel2.sum(), config.data.channels, config.data.image_size, self.device)
        else:
            print("ERROR: the degradation model type (--deg) is not supported")
            quit()


        idx_so_far = 0
        print(f'Start restoration')
        for _ in range(len(dasLst)):
            dasSaveName = dasLst[idx_so_far] + '.mat'
            gamma = gammaLst[idx_so_far]
            # load/simulate the observation y_0  & define the save path
            if args.deg in ["DRUS", "DENO"]:
                y_0 = torch.from_numpy(mat73.loadmat(args.matlab_path + 'picmus/Observation/'+args.deg+'/' + dasSaveName)['By'])
                y_0 = (y_0.view(1, -1)).repeat(1, config.data.channels).to(self.device)
                gamma = gamma * sqrt(3)
            elif args.deg in ["NumericalCysts", "NumericalScatterers"]:
                # load the ground Truth x_orig, (dasSaveName has format signed_m5_XX.mat)
                x_orig = torch.from_numpy(
                    mat73.loadmat(args.matlab_path + 'numerical/SimulatedData/'+numeType+'/' + dasSaveName)['x'])
                x_orig = (x_orig.view(1, 1, 256, 256)).to(self.device)
                # apply the inverse problem model (y_0 becomes with shape=(1,xxx))
                y_0 = H_funcs.H(x_orig)
                amp = y_0.abs().max()
                if idx_so_far == 0:
                    print('y_0.shape = ', y_0.shape)
                # load and truncate the additive noise to the same size as y_0
                additiveNoise = torch.from_numpy(mat73.loadmat(args.matlab_path+'numerical/SimulatedData/additiveNoises/additiveNoise_' + str(gamma) + '.mat')['n'])
                additiveNoise = (additiveNoise.view(1, 65536)).to(self.device)
                # add the additive noise
                y_0 = y_0 + additiveNoise * amp
                gamma = float(gamma * amp)
                # save the degraded image as .mat
                savemat(os.path.join(resultsPath, f"{idx_so_far + 1}.mat"), {'y_0': y_0[0].cpu().numpy()})

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
                x, _ = self.sample_image(x, model, H_funcs, y_0, gamma, last=False, cls_fn=cls_fn, timesteps=args.timesteps)
            # save the DDRM restored image as .mat
            savemat(os.path.join(resultsPath, f"{idx_so_far + 1}_{-1}.mat"),{'x': x[-1][0].detach().cpu().numpy()})
            idx_so_far += y_0.shape[0]  # iterate multiple images
            print(f'Finish {idx_so_far}')

    def sample_image(self, x, model, H_funcs, y_0, gamma, last=False, cls_fn=None, classes=None, timesteps=50):
        skip = self.num_timesteps // timesteps
        seq = range(0, self.num_timesteps, skip)
        x = efficient_generalized_steps(x, seq, model, self.betas, H_funcs, y_0, gamma, etaB=self.args.etaB,
                                        etaA=self.args.eta, etaC=self.args.eta, cls_fn=cls_fn, classes=classes)
        if last:
            x = x[0][-1]
        return x
