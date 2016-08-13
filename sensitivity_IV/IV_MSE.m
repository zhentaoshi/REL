MSE = zeros(4, 8);
b_true = 1;



load('DGP_linearIV_BB_n_120_m_80_Rep_500_C1_0.25.mat')
MSE(1, 1:2) = [120, 80];
A = B2(:,1) - b_true;
MSE(1, 3:4) = [std(A), mean(A)];
load('DGP_linearIV_BB_n_120_m_80_Rep_500_C1_0.5.mat');
A = B2(:,1) - b_true;
MSE(1, 5:6) = [std(A), mean(A)];
load('DGP_linearIV_BB_n_120_m_80_Rep_500_C1_0.75.mat');
A = B2(:,1) - b_true;
MSE(1, 7:8) = [std(A), mean(A)];
load('DGP_linearIV_BB_n_120_m_80_Rep_500_C1_1.mat');
A = B2(:,1) - b_true;
MSE(1, 9:10) = [std(A), mean(A)];



load('DGP_linearIV_BB_n_120_m_160_Rep_500_C1_0.25.mat')
MSE(2, 1:2) = [120, 160];
A = B2(:,1) - b_true;
MSE(2, 3:4) = [std(A), mean(A)];
load('DGP_linearIV_BB_n_120_m_160_Rep_500_C1_0.5.mat');
A = B2(:,1) - b_true;
MSE(2, 5:6) = [std(A), mean(A)];
load('DGP_linearIV_BB_n_120_m_160_Rep_500_C1_0.75.mat');
A = B2(:,1) - b_true;
MSE(2, 7:8) = [std(A), mean(A)];
load('DGP_linearIV_BB_n_120_m_160_Rep_500_C1_1.mat');
A = B2(:,1) - b_true;
MSE(2, 9:10) = [std(A), mean(A)];



load('DGP_linearIV_BB_n_240_m_80_Rep_500_C1_0.25.mat')
MSE(3, 1:2) = [240, 80];
A = B2(:,1) - b_true;
MSE(3, 3:4) = [std(A), mean(A)];
load('DGP_linearIV_BB_n_240_m_80_Rep_500_C1_0.5.mat');
A = B2(:,1) - b_true;
MSE(3, 5:6) = [std(A), mean(A)];
load('DGP_linearIV_BB_n_240_m_80_Rep_500_C1_0.75.mat');
A = B2(:,1) - b_true;
MSE(3, 7:8) = [std(A), mean(A)];
load('DGP_linearIV_BB_n_240_m_80_Rep_500_C1_1.mat');
A = B2(:,1) - b_true;
MSE(3, 9:10) = [std(A), mean(A)];




load('DGP_linearIV_BB_n_240_m_160_Rep_500_C1_0.25.mat')
MSE(4, 1:2) = [240, 160];
A = B2(:,1) - b_true;
MSE(4, 3:4) = [std(A), mean(A)];
load('DGP_linearIV_BB_n_240_m_160_Rep_500_C1_0.5.mat');
A = B2(:,1) - b_true;
MSE(4, 5:6) = [std(A), mean(A)];
load('DGP_linearIV_BB_n_240_m_160_Rep_500_C1_0.75.mat');
A = B2(:,1) - b_true;
MSE(4, 7:8) = [std(A), mean(A)];
load('DGP_linearIV_BB_n_240_m_160_Rep_500_C1_1.mat');
A = B2(:,1) - b_true;
MSE(4, 9:10) = [std(A), mean(A)];

export( dataset(MSE), 'file', 'IV_REL_MSE.csv', 'Delimiter', ',')

