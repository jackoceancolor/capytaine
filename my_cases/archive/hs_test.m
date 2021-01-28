

fn = 'sphere_hs_test.nc';
load('hydro_sphere_nemoh.mat');
hydro_n = hydro;
load('hydro_sphere_wamit.mat');
hydro_w = hydro;
clear hydro

info = ncinfo(fn);
vo = ncread(fn,'displaced_volume');
k_hs = ncread(fn,'hydrostatic_stiffness');
cg = ncread(fn,'center_of_mass');
cb = ncread(fn,'center_of_buoyancy');

[hydro_n.cb hydro_w.cb cb]
[hydro_n.cg hydro_w.cg cg]
[hydro_n.Vo hydro_w.Vo vo]
[hydro_n.C -9.99e6*ones(6,1) hydro_w.C -9.99e6*ones(6,1) k_hs]


%%
clear;clc;close all;

hydro = struct();
hydro = Read_CAPYTAINE(hydro,'sphere_hs_test.nc');




