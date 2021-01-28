% function to compare Capytaine output with another BEM code.
%
% Adam Keester, 11/5/2020
clc; clear all; close all;
hydro = struct();

%% Compare sphere to disordered capt output 
% The 'disordered' capytaine output refers to an .nc file that contains 
% jumbled direction dimensions (e.g. [roll, surge, sway, pitch, yaw, heav])
% Correctly ordered dimensions is not guaranteed in capytaine output so
% the correction in read_capytaine is tested here.

hydro = Read_CAPYTAINE(hydro, './sphere_short.nc');
% adjust cg, cb, vo because hydrostatics aren't done yet
% values for sphere:
hydro(1).cg = [0; 0; -2];
hydro(1).cb = [0; 0; -1.874];
hydro(1).Vo = 261.3639;

hydro = Read_CAPYTAINE(hydro, './sphere_short_disordered.nc');
% adjust cg, cb, vo because hydrostatics aren't done yet
% values for sphere:
hydro(2).cg = [0; 0; -2];
hydro(2).cb = [0; 0; -1.874];
hydro(2).Vo = 261.3639;
hydro(2).body = 'sphere_capytaine_disordered';

hydro = Read_CAPYTAINE(hydro, './sphere_test.nc');
% adjust cg, cb, vo because hydrostatics aren't done yet
% values for sphere:
hydro(3).cg = [0; 0; -2];
hydro(3).cb = [0; 0; -1.874];
hydro(3).Vo = 261.3639;
hydro(3).body = 'sphere_test';

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
hydro = Read_CAPYTAINE(hydro, './sphere.nc'); % cpt rotation about cog, nemoh mesh
% adjust cg, cb, vo because hydrostatics aren't dont yet
hydro(1).cg = [0; 0; -2];
hydro(1).cb = [0; 0; -1.874];
hydro(1).Vo = 261.3639;

hydro = Read_WAMIT(hydro,'..\..\WAMIT\Sphere\sphere.out',[]); % Compare to WAMIT
hydro(2).body = 'sphere_wamit'; % change name for plotting clarity

hydro = Read_NEMOH(hydro,'..\..\NEMOH\Sphere\'); % Compare to NEMOH
hydro(3).body = 'sphere_nemoh'; % change name for plotting clarity

hydro = Combine_BEM(hydro); % Combine hydro structs to plot comparison together
hydro = Radiation_IRF(hydro,60,[],[],[],[]);
hydro = Radiation_IRF_SS(hydro,[],[]);
hydro = Excitation_IRF(hydro,160,[],[],[],[]);
Write_H5(hydro)
Plot_BEMIO(hydro)
