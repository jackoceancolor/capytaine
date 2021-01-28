% clc; clear all; close all;
clear
hydro = struct();

load('../WEC-Sim/hydro_rm3_wamit.mat');
hydro.body = {'float_wamit','spar_wamit'};

hydro = Read_CAPYTAINE_b2b_v3(hydro,'.\rm3_b2b_full.nc'); % 2 bodies, run all 12 dofs for each but not combined
hydro(end).body = {'float_cpt_oldB2B','spar_cpt_oldB2B'};

% hydro = Read_CAPYTAINE_b2b(hydro,'.\rm3_full.nc'); % 2 bodies, full run, run separately and problems appended
% hydro(end).body = {'float_cpt','spar_cpt'};

hydro = Read_CAPYTAINE_b2b_v4(hydro,'.\rm3_b2b_new_full.nc'); % runs all 12 dofs for the combined body
hydro(end).body = {'float_cpt_newB2B','spar_cpt_newB2B'};

% hydro = Read_NEMOH(hydro,'..\..\NEMOH\RM3\');
% hydro(end+1) = Read_WAMIT(hydrow,'..\..\WAMIT\RM3\rm3.out',[]);
% hydrow = load('../WEC-Sim/hydro_rm3_wamit.mat');
% hydro(end) = hydrow.hydro;
% hydro(end).body = {'float_wamit','spar_wamit'};

hydro = Combine_BEM(hydro); % Compare to NEMOH and WAMIT
hydro = Radiation_IRF(hydro,20,[],[],[],[]);
hydro = Radiation_IRF_SS(hydro,[],[]);
hydro = Excitation_IRF(hydro,200,[],[],[],[]);
% % Write_H5(hydro)
Plot_BEMIO(hydro)

