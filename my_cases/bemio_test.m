% clc; clear all; close all;
hydro = struct();

% hydro = Read_CAPYTAINE_b2b_v4(hydro,'.\RM3\temp.nc'); % test new hydrostatics implementation
% hydro(1).body = {'float_cpt_b2b','spar_cpt_b2b'};

% hydro = Read_CAPYTAINE_b2b_v4(hydro,'.\RM3\test.nc');
% hydro(1).body = {'float_cpt_b2b','spar_cpt_b2b'};
% 
hydro = Read_CAPYTAINE_b2b_v4(hydro,'.\RM3\rm3_b2b_test.nc');
hydro(end).body = {'float_cpt_old_b2b','spar_cpt_old_b2b'};

hydro = Combine_BEM(hydro); % Compare to NEMOH and WAMIT
hydro = Radiation_IRF(hydro,20,[],[],[],[]);
hydro = Radiation_IRF_SS(hydro,[],[]);
hydro = Excitation_IRF(hydro,200,[],[],[],[]);
% % Write_H5(hydro)
Plot_BEMIO(hydro)

