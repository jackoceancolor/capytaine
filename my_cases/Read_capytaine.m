% function to start testing the reading of Capytaine output
% assuming netcdf format

% import WAMIT RM3 simulation for comparison
hydro_wamit = load('hydro_rm3_wamit.mat');
hydro_wamit = hydro_wamit.hydro;
hydro_nemoh = load('hydro_rm3_nemoh.mat');
hydro_nemoh = hydro_nemoh.hydro;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in NETCDF file\
% matlab functions:
% nccreate	Create variable in NetCDF file
% ncdisp	Display contents of NetCDF data source in Command Window
% ncinfo	Return information about NetCDF data source
% ncread	Read data from variable in NetCDF data source
% ncreadatt	Read attribute value from NetCDF data source
% ncwrite	Write data to NetCDF file
% ncwriteatt	Write attribute to NetCDF file
% ncwriteschema	Add NetCDF schema definitions to NetCDF file

% choose file to read
% filename = '"C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases/capytaine1.nc"';
% fn2 = eval(filename);

% filename = '"C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases/rm3.nc"';
% fn2 = eval(filename);

filename = '"C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases/rm3_combo.nc"';
fn2 = eval(filename);

% filename = '"C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases/rm3_combo_long.nc"';
% fn2 = eval(filename);

hydro = struct();

BEGIN TESTING SCRIPT USING THE RM3_COMBO.NC AND RM3_COMBO_LONG.NC OUTPUT

[a,b] = size(hydro);  % Check on what is already there
if b==1
    if isfield(hydro(b),'Nb')==0;  F = 1;
    else;  F = 2;
    end
elseif b>1;  F = b+1;
end

%%
p = waitbar(0,'Reading AQWA output file...'); %Progress bar

hydro(F).code = 'CAPYTAINE';
[filepath,name,ext] = fileparts(filename);
hydro(F).file = name;  % Base name

% get names of all variables and dimensions in Capytaine output
info = ncinfo(fn2);
cpt_vars = cell(length(info.Variables),1);
cpt_dims = cell(length(info.Dimensions),1);
for ii=1:length(info.Variables)
    cpt_vars{ii} = info.Variables(ii).Name;
end
for ii=1:length(info.Dimensions)
    cpt_dims{ii} = info.Dimensions(ii).Name;
end

% list of variables to get from netcdf file
req_vars = {
    'complex',...
    'body_name',...
    'g',...
    'rho',...
    'water_depth',...
    'influenced_dof'...
    'radiating_dof',...
    'theta',...
    'omega',...
    'wave_direction',...
    'added_mass',...
    'radiation_damping',...
    'diffraction_force',...
    'Froude_Krylov_force',...
    'kochin_diffraction',...
    'kochin',...     % variables
    };
% {'wavenumber';'wavelength';'radiation_damping';} % other variables in .nc file not currently used

% check that all required Capytaine output is present
for i=1:length(req_vars)
    if ~contains(lower(req_vars{i}), lower(cpt_vars))
        error('capytaine output does not contain: %s',req_vars{i});
    end
end


% begin parsing netcdf file to hydro struct
cpt = struct();

% cg, 
% cb, 
% Vo
hydro(F).rho = ncread(fn2,'rho');
hydro(F).g = ncread(fn2,'g');
hydro(F).h = ncread(fn2,'water_depth');
hydro(F).Nb = info.Dimensions(getInd(info.Dimensions,'body_name')).Length;
hydro(F).Nf = info.Dimensions(getInd(info.Dimensions,'omega')).Length;
hydro(F).Nh = info.Dimensions(getInd(info.Dimensions,'wave_direction')).Length;

tmp = ncread(fn2,'body_name')';
for i=1:hydro(F).Nb
    hydro(F).body{i} = tmp(i,:);
end
hydro(F).w = ncread(fn2,'omega')';
hydro(F).T = 1./hydro(F).w;
hydro(F).beta = ncread(fn2,'wave_direction');

dof_i = info.Dimensions(getInd(info.Dimensions,'influenced_dof')).Length;
dof_r = info.Dimensions(getInd(info.Dimensions,'radiating_dof')).Length;
hydro(F).dof = [dof_i, dof_r];
waitbar(1/8);


%%
hydro(F).cg = zeros(3,hydro(F).Nb);
hydro(F).cb = zeros(3,hydro(F).Nb);
hydro(F).cb = zeros(1,hydro(F).Nb);
% Vo = Displacement volume
waitbar(2/8);

%%
% hydro(F).C(i,:,m) = tmp{1,1}(1:6);  % Linear restoring stiffness
waitbar(3/8);

%%
% hydro(F).A(i,:,k) = tmp{1,1}(2:2:end);  % Added mass
% hydro(F).B(i,:,k) = tmp{1,1}(3:2:end);  % Radiation damping
waitbar(4/8);

%%
tmp = ncread(fn2,'complex')';
if tmp(1,:) == "re" && tmp(2,:) == "im"
    i_re = 1;
    i_im = 2;
elseif tmp(1,:) == "im" && tmp(2,:) == "re"
    i_im = 1;
    i_re = 2;
else
    error('check complex dimension indices');
end

%% Excitation Force [6*Nb,Nh,Nf];
% CHECK if the kochin_diffraction is actually the excitation force
% i_var = getInd(info.Variables,'kochin_diffraction');
% dim = info.Variables(i_var).Dimensions;
% if (dim(1).Name~="influenced_dof" || dim(2).Name~="wave_direction" || ...
%         dim(3).Name~="omega" || dim(4).Name~="body_name" || dim(5).Name~="complex")
%     error("capytain dimensions for diffraction_force incorrect");
% end
% 
% tmp = ncread(fn2,'kochin_diffraction');
% for n=1:hydro(F).Nb
%     hydro(F).sc_re(6*(n-1)+1:6*n,:,:) = tmp(:,:,:,n,1);         % Real part of diffraction force
%     hydro(F).sc_im(6*(n-1)+1:6*n,:,:) = tmp(:,:,:,n,2);         % Imaginary part of diffraction force
% end
% hydro(F).sc_ma = (hydro(F).sc_re.^2 + hydro(F).sc_im.^2).^0.5;  % Magnitude of diffraction force
% hydro(F).sc_ph = atan(hydro(F).sc_im./hydro(F).sc_re);          % Phase of diffraction force
% waitbar(5/8);

%% Diffraction Force (scattering) [6*Nb,Nh,Nf];
i_var = getInd(info.Variables,'diffraction_force');
dim = info.Variables(i_var).Dimensions;
if (dim(1).Name~="influenced_dof" || dim(2).Name~="wave_direction" || ...
        dim(3).Name~="omega" || dim(4).Name~="body_name" || dim(5).Name~="complex")
    error("capytain dimensions for diffraction_force incorrect");
end

tmp = ncread(fn2,'diffraction_force');
for n=1:hydro(F).Nb
    hydro(F).sc_re(6*(n-1)+1:6*n,:,:) = tmp(:,:,:,n,1);         % Real part of diffraction force
    hydro(F).sc_im(6*(n-1)+1:6*n,:,:) = tmp(:,:,:,n,2);         % Imaginary part of diffraction force
end
hydro(F).sc_ma = (hydro(F).sc_re.^2 + hydro(F).sc_im.^2).^0.5;  % Magnitude of diffraction force
hydro(F).sc_ph = atan(hydro(F).sc_im./hydro(F).sc_re);          % Phase of diffraction force
waitbar(6/8);

%% Froude-Krylov force file
i_var = getInd(info.Variables,'Froude_Krylov_force');
dim = info.Variables(i_var).Dimensions;
if (dim(1).Name~="influenced_dof" || dim(2).Name~="wave_direction" || ...
        dim(3).Name~="omega" || dim(4).Name~="body_name" || dim(5).Name~="complex")
    error("capytain dimensions for Froude_Krylov_force incorrect");
end

tmp = ncread(fn2,'Froude_Krylov_force');
for n=1:hydro(F).Nb
    hydro(F).fk_re(6*(n-1)+1:6*n,:,:) = tmp(:,:,:,n,i_re);      % Real part of diffraction force
    hydro(F).fk_im(6*(n-1)+1:6*n,:,:) = tmp(:,:,:,n,i_im);      % Imaginary part of diffraction force
end
hydro(F).fk_ma = (hydro(F).fk_re.^2 + hydro(F).fk_im.^2).^0.5;  % Magnitude of diffraction force
hydro(F).fk_ph = atan(hydro(F).fk_im./hydro(F).fk_re);          % Phase of diffraction force
waitbar(7/8);

%% ================= READING KOCHIN FILES ===================%
theta(ntheta)= Kochin(3*(ntheta-1)+1); % theta
Kochin_BVP(ntheta,1,x)= Kochin(3*(ntheta-1)+2); % magnitude
Kochin_BVP(ntheta,2,x)= Kochin(3*(ntheta-1)+3); % phase

% these from addditional functions called in bemio.m files?
% hydro_nemoh.Ainf,
% hydro_nemoh.ra_K, 
% hydro_nemoh.ra_t, 
% hydro_nemoh.ra_w, 
% hydro_nemoh.ss_A, 
% hydro_nemoh.ss_B, 
% hydro_nemoh.ss_C, 
% hydro_nemoh.ss_D, 
% hydro_nemoh.ss_K, 
% hydro_nemoh.ss_conv, 
% hydro_nemoh.ss_R2, 
% hydro_nemoh.ss_O, 
% hydro_nemoh.ex_K, 
% hydro_nemoh.ex_t, 
% hydro_nemoh.ex_w

waitbar(8/8);

%%
% get all required variables from netcdf file
% for i=1:length(req_vars)
%     current_var = req_vars{i};
% %     f = strcat('hydro.', current_var, ' = ncread(', fn2, ', "', current_var, '");');
%     f = ['cpt.', current_var, ' = ncread(', filename, ', "', current_var, '");'];
% %     fprintf([f '\n']);
%     eval(f);
% end; clear i;
% cpt.radiating_dof = cpt.radiating_dof';      % transpose char strings to correct orientation
% cpt.influenced_dof = cpt.influenced_dof';    % transpose char strings to correct orientation

% hydro(F).Nb = info.Dimensions(getInd(info.Dimensions,'body_name')).Length;
% hydro(F).Nf = info.Dimensions(getInd(info.Dimensions,'body_name')).Length;

% bodies = 


function ind = getInd(dimStruct, str2find)
    j = 0;
    for j=1:length(dimStruct)
        if string(dimStruct(j).Name) == str2find
            ind = j;
        end
    end
end

