clc; clear all; close all;
hydro = struct();

hydro = Read_CAPYTAINE_b2b_v2(hydro,'.\cubes_b2b_full.nc');
hydro(1).body = {'r_cube_cpt_b2b','t_cube_cpt_b2b'};

hydro = Read_CAPYTAINE(hydro,'.\cubes_full.nc');
hydro(2).body = {'r_cube_cpt','t_cube_cpt'};

% hydro = Read_NEMOH(hydro,'..\..\NEMOH\Cubes\');
% hydro = Read_WAMIT(hydro,'..\..\WAMIT\Cubes\cubes.out',[]);

hydro = Combine_BEM(hydro); % Compare to NEMOH and WAMIT
hydro = Radiation_IRF(hydro,20,[],[],[],[]);
hydro = Radiation_IRF_SS(hydro,[],[]);
hydro = Excitation_IRF(hydro,200,[],[],[],[]);
% Write_H5(hydro)
Plot_BEMIO(hydro)

