saveto = '/numerical/SimulatedData/additiveNoises/';
n1 = single(randn(256*256,1));

gammaLevelsCysts = [0.02, 0.05, 0.08, 0.11, 0.14, 0.17, 0.2, 0.23, 0.26, 0.29, 0.32, 0.35];
for gamma = gammaLevelsCysts
    n = single(n1 .* gamma);
    save([pwd saveto 'additiveNoise_' num2str(gamma) '.mat'], 'n')
end

gammaLevelsScatterers = [0,0.003,0.006,0.009,0.012,0.015,0.018, 0.04, 0.06, 0.08, 0.1];
for gamma = gammaLevelsScatterers
    n = single(n1 .* gamma);
    save([pwd saveto 'additiveNoise_' num2str(gamma) '.mat'], 'n')
end

