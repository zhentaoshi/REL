# to run this file, type in linux: 
# bash: file_name.bat



matlab -nosplash -nodesktop -r " load('DGP_linearIV_BB_n_120_m_80_Rep_500_C1_0.5.mat');master_boost; exit;"  &
matlab -nosplash -nodesktop -r "load('DGP_linearIV_BB_n_120_m_160_Rep_500_C1_0.5.mat');master_boost; exit;"  &
matlab -nosplash -nodesktop -r " load('DGP_linearIV_BB_n_240_m_80_Rep_500_C1_0.5.mat');master_boost; exit;"  &
matlab -nosplash -nodesktop -r "load('DGP_linearIV_BB_n_240_m_160_Rep_500_C1_0.5.mat');master_boost; exit;"  &

matlab -nosplash -nodesktop -r " load('DGP_linearIV_BB_n_120_m_80_Rep_500_C1_0.25.mat');master_boost; exit;"  &
matlab -nosplash -nodesktop -r "load('DGP_linearIV_BB_n_120_m_160_Rep_500_C1_0.25.mat');master_boost; exit;"  &
matlab -nosplash -nodesktop -r " load('DGP_linearIV_BB_n_240_m_80_Rep_500_C1_0.25.mat');master_boost; exit;"  &
matlab -nosplash -nodesktop -r "load('DGP_linearIV_BB_n_240_m_160_Rep_500_C1_0.25.mat');master_boost; exit;"  &

matlab -nosplash -nodesktop -r " load('DGP_linearIV_BB_n_120_m_80_Rep_500_C1_0.75.mat');master_boost; exit;"  &
matlab -nosplash -nodesktop -r "load('DGP_linearIV_BB_n_120_m_160_Rep_500_C1_0.75.mat');master_boost; exit;"  &
matlab -nosplash -nodesktop -r " load('DGP_linearIV_BB_n_240_m_80_Rep_500_C1_0.75.mat');master_boost; exit;"  &
matlab -nosplash -nodesktop -r "load('DGP_linearIV_BB_n_240_m_160_Rep_500_C1_0.75.mat');master_boost; exit;"  &

matlab -nosplash -nodesktop -r " load('DGP_linearIV_BB_n_120_m_80_Rep_500_C1_1.mat');master_boost; exit;"  &
matlab -nosplash -nodesktop -r "load('DGP_linearIV_BB_n_120_m_160_Rep_500_C1_1.mat');master_boost; exit;"  &
matlab -nosplash -nodesktop -r " load('DGP_linearIV_BB_n_240_m_80_Rep_500_C1_1.mat');master_boost; exit;"  &
matlab -nosplash -nodesktop -r "load('DGP_linearIV_BB_n_240_m_160_Rep_500_C1_1.mat');master_boost; exit;"  &
