% test capytaine b2b output

%% RM3
load('./WEC-Sim/hydro_rm3_nemoh.mat');
hydro_n = hydro;
load('./WEC-Sim/hydro_rm3_wamit.mat');
hydro_w = hydro;
clear hydro

% fn = './RM3/rm3_b2b_full.nc';
fn = './RM3/temp.nc';
info = ncinfo(fn);
hydro = struct();
hydro = Read_CAPYTAINE_b2b_v4(hydro,fn);
hydro_c = hydro;
clear hydro

% check if b2b interactions are even present
am = [abs(hydro_c.A(1:6,7:12,:)) abs(hydro_c.A(7:12,1:6,:))];
rd = [abs(hydro_c.B(1:6,7:12,:)) abs(hydro_c.B(7:12,1:6,:))];
err_a = (hydro_c.A(:,:,1)-hydro_w.A(:,:,1))./hydro_w.A(:,:,1)*100;
err_b = (hydro_c.B(:,:,1)-hydro_w.B(:,:,1))./hydro_w.B(:,:,1)*100;
max(abs(err_a),[],'all')
max(abs(err_b),[],'all')

%% Cubes
load('./WEC-Sim/hydro_cubes_nemoh.mat');
hydro_n = hydro;
load('./WEC-Sim/hydro_cubes_wamit.mat');
hydro_w = hydro;
clear hydro

fn = './Cubes/cubes_b2b_full.nc';
info = ncinfo(fn);
hydro = struct();
hydro = Read_CAPYTAINE_b2b_v4(hydro,fn);
hydro_c = hydro;
clear hydro

% check if b2b interactions are even present
err_a = max( [abs(hydro_c.A(1:6,7:12,:)) abs(hydro_c.A(7:12,1:6,:))],[],'all')
err_b = max( [abs(hydro_c.B(1:6,7:12,:)) abs(hydro_c.B(7:12,1:6,:))],[],'all')


%% Sphere (single body)
fn = './Sphere/sphere_full.nc';
load('./WEC-Sim/hydro_sphere_wamit.mat');
hydro_w = hydro;
clear hydro

info = ncinfo(fn);
hydro = struct();
hydro = Read_CAPYTAINE_b2b_v4(hydro,fn);
hydro = Read_CAPYTAINE(hydro,fn);
hydro = Read_WAMIT(hydro,'../../WEC-Sim/examples/BEMIO/WAMIT/Sphere/sphere.out',[]);
hydro = Combine_BEM(hydro); % Compare to NEMOH and WAMIT
hydro = Radiation_IRF(hydro,20,[],[],[],[]);
hydro = Radiation_IRF_SS(hydro,[],[]);
hydro = Excitation_IRF(hydro,200,[],[],[],[]);
% Write_H5(hydro)
Plot_BEMIO(hydro)


%% compare reading directly from nc file to hydro output from read_capytaine
% w_c = ncread(fn,'omega');
% 
% am_c = ncread(fn,'added_mass')/1000;
% max(abs(sum(am_c,4)/1000-hydro_c.A),[],'all')
% err_amc = max( [abs(am_c(1:6,7:12,:,:)) abs(am_c(7:12,1:6,:,:))],[],'all')
% 
% tmpw = zeros(1,1,length(w_c),1);
% tmpw(1,1,:,1) = w_c;
% rd_c = ncread(fn,'radiation_damping')./(1000*tmpw);
% max(abs(sum(rd_c,4)-hydro_c.B),[],'all')
% err_rdc = max( [abs(rd_c(1:6,7:12,:,:)) abs(rd_c(7:12,1:6,:,:))],[],'all')



