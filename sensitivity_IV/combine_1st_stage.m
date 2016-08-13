% combine the 4 first-stage estimators and export to a csv file.
load('DGP_Han_BB_n_240_m_160_Rep_500_C1_0.5_seed_204.mat')
load('GMM_DGP_Han_B_n_240_m_160_Rep_500_C1_0.5_seed_204.mat')

BB = dataset(B1, B2, B_GMM0, B_GMM1);

title = ['combine_1st_stage_DGP_', DGP, '_B_n_', num2str(n), '_m_', num2str(m), '_Rep_', num2str(Rep), '_C1_', ...
    num2str(C1), '_seed_', num2str(seed) ];
export(BB, 'File', [title, '.csv'], 'delimiter', ',');
