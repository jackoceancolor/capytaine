% function to compare Capytaine output with another BEM code.
%
% Adam Keester, 11/5/2020
clc; clear all; close all;
hydro = struct();

cpt_dir = 'C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases/';

%% Compare sphere to disordered capt output
hydro = Read_capytaine_v1(hydro, [cpt_dir '/sphere_test.nc']);
% adjust cg, cb, vo because hydrostatics aren't dont yet
hydro(1).cg = [0; 0; -2];
hydro(1).cb = [0; 0; -1.874];
hydro(1).Vo = 261.3639;

hydro = Read_capytaine_v1(hydro, [cpt_dir '/sphere_test_disordered.nc']);
% adjust cg, cb, vo because hydrostatics aren't dont yet
hydro(2).cg = [0; 0; -2];
hydro(2).cb = [0; 0; -1.874];
hydro(2).Vo = 261.3639;

hydro = Combine_BEM(hydro); % Combine hydro structs to plot comparison together
hydro = Radiation_IRF(hydro,60,[],[],[],[]);
hydro = Radiation_IRF_SS(hydro,[],[]);
hydro = Excitation_IRF(hydro,160,[],[],[],[]);
Write_H5(hydro)
Plot_BEMIO(hydro)
hydro_tests = hydro;
clear hydro

%% Compare sphere to wamit, nemoh
hydro = struct();
hydro = Read_capytaine_v1(hydro, [cpt_dir '/sphere_1.nc']);
% adjust cg, cb, vo because hydrostatics aren't dont yet
hydro(1).cg = [0; 0; -2];
hydro(1).cb = [0; 0; -1.874];
hydro(1).Vo = 261.3639;

hydro = Read_WAMIT(hydro,'..\WAMIT\Sphere\sphere.out',[]); % Compare to WAMIT
hydro = Read_NEMOH(hydro,'..\NEMOH\Sphere\'); % Compare to NEMOH
hydro = Combine_BEM(hydro); % Combine hydro structs to plot comparison together
hydro = Radiation_IRF(hydro,60,[],[],[],[]);
hydro = Radiation_IRF_SS(hydro,[],[]);
hydro = Excitation_IRF(hydro,160,[],[],[],[]);
Write_H5(hydro)
Plot_BEMIO(hydro)

%%
% % for RM3
% hydro.cg = [0,0; 0,0; -0.72,-21.29];
% hydro.cb = [0,0; 0,0; -1.293,-15.604];
% hydro.Vo = [725.834, 886.691];
% file = [cpt_dir '/rm3.nc']; % short (w=3x1) RM3 run, two bodies, no b2b
% % file = [cpt_dir '/rm3_combo.nc']; % short RM3 run, two bodies, attempted b2b (doesn't work?)
% % file = [cpt_dir '/rm3_combo_long.nc']; % long RM3 run, two bodies, attempted b2b (doesn't work?)


%% individual plots comparing hydro structs
clear;clc;
load('hydro_sphere_test_capy.mat');
hydro_c = hydro;

load('hydro_sphere_wamit.mat');
hydro_w = hydro;

load('hydro_sphere_nemoh.mat');
hydro_n = hydro;

for i=1:length(hydro_c.w)
    ind_w(i) = find(hydro_w.w==hydro_c.w(i));
end

figure()
plot(hydro_c.w, hydro_c.A(1,1,:),...
    hydro_c.w, hydro_w.A(1,1,ind_w),...
    hydro_c.w, hydro_n.A(1,1,ind_w));
legend('capy','wamit','nemoh');
title('Surge added mass');

figure()
plot(hydro_c.w, hydro_c.A(1,1,:),...
    hydro_c.w, hydro_w.A(1,1,ind_w),...
    hydro_c.w, hydro_n.A(1,1,ind_w));
legend('capy','wamit','nemoh');
title('Heave added mass');

figure()
plot(hydro_c.w, hydro_c.A(1,1,:),...
    hydro_c.w, hydro_w.A(1,1,ind_w),...
    hydro_c.w, hydro_n.A(1,1,ind_w));
legend('capy','wamit','nemoh');
title('Pitch added mass');

figure()
plot(hydro_c.w, hydro_c.ex_ma(1,:,:),...
    hydro_c.w, hydro_w.ex_ma(1,:,ind_w),...
    hydro_c.w, hydro_n.ex_ma(1,:,ind_w));
legend('capy','wamit','nemoh');
title('Surge F_exc magnitude');

figure()
plot(hydro_c.w, hydro_c.ex_ma(3,:,:),...
    hydro_c.w, hydro_w.ex_ma(3,:,ind_w),...
    hydro_c.w, hydro_n.ex_ma(3,:,ind_w));
legend('capy','wamit','nemoh');
title('Heave F_exc magnitude');

figure()
plot(hydro_c.w, hydro_c.ex_ma(5,:,:),...
    hydro_c.w, hydro_w.ex_ma(5,:,ind_w),...
    hydro_c.w, hydro_n.ex_ma(5,:,ind_w));
legend('capy','wamit','nemoh');
title('Pitch F_exc magnitude');



%%
% % wamit
% hydro = Read_WAMIT(hydro,'rm3.out',[]);
% hydro = Radiation_IRF(hydro,20,[],[],[],[]);
% hydro = Radiation_IRF_SS(hydro,[],[]);
% hydro = Excitation_IRF(hydro,20,[],[],[],[]);
% Write_H5(hydro)
% Plot_BEMIO(hydro)

% % nemoh
% hydro = Read_NEMOH(hydro,'..\RM3\');
% % hydro = Read_WAMIT(hydro,'..\..\WAMIT\RM3\rm3.out',[]);
% % hydro = Combine_BEM(hydro); % Compare WAMIT
% hydro = Radiation_IRF(hydro,60,[],[],[],1.9);
% hydro = Radiation_IRF_SS(hydro,[],[]);
% hydro = Excitation_IRF(hydro,157,[],[],[],1.9);
% Write_H5(hydro)
% Plot_BEMIO(hydro)

% % aqwa
% hydro = Read_AQWA(hydro,'aqwa_example_data.AH1','aqwa_example_data.LIS');
% hydro = Radiation_IRF(hydro,15,[],[],[],[]);
% hydro = Radiation_IRF_SS(hydro,[],[]);
% hydro = Excitation_IRF(hydro,15,[],[],[],[]);
% Write_H5(hydro)
% Plot_BEMIO(hydro)

% these from addditional functions called in bemio.m files?
% hydro_nemoh.Ainf, % from Normalize() - simply gets A(max frequency)
% hydro_nemoh.ra_K, % radiation IRF (K), from Radiation_IRF()
% hydro_nemoh.ra_t, % radiation IRF time, from Radiation_IRF()
% hydro_nemoh.ra_w, % radiation IRF frequency, from Radiation_IRF()
% hydro_nemoh.ss_A, % from Radiation_IRF_SS()
% hydro_nemoh.ss_B, % from Radiation_IRF_SS()
% hydro_nemoh.ss_C, % from Radiation_IRF_SS()
% hydro_nemoh.ss_D, % from Radiation_IRF_SS()
% hydro_nemoh.ss_K, % from Radiation_IRF_SS()
% hydro_nemoh.ss_conv, % from Radiation_IRF_SS()
% hydro_nemoh.ss_R2, % from Radiation_IRF_SS()
% hydro_nemoh.ss_O, % from Radiation_IRF_SS()
% hydro_nemoh.ex_K, % from Excitation_IRF()
% hydro_nemoh.ex_t, % from Excitation_IRF()
% hydro_nemoh.ex_w % from Excitation_IRF()
