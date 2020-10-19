% function to start testing the reading of Capytaine output
% assuming netcdf format

% import WAMIT RM3 simulation for comparison
hydro_wamit = load('hydro_rm3_wamit.mat');
hydro_wamit = hydro_wamit.hydro;

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

filename = '"C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases/rm3.nc"';
fn2 = eval(filename);

hydro = struct();

[a,b] = size(hydro);  % Check on what is already there
if b==1
    if isfield(hydro(b),'Nb')==0;  F = 1;
    else;  F = 2;
    end
elseif b>1;  F = b+1;
end

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
    'omega',...
    'radiating_dof',...   % dimensions
    'theta',...
    'wave_direction',...
    'influenced_dof'...        % dimensions
    'added_mass',...
    'radiation_damping',...
    'diffraction_force',...    % variables
    'Froude_Krylov_force',...
    'kochin_diffraction',...
    'kochin',...     % variables
    };

% check that all required Capytaine output is present
for i=1:length(req_vars)
    if ~contains(lower(req_vars{i}), lower(cpt_vars))
        error('capytaine output does not contain: %s',req_vars{i});
    end
end

% being parsing netcdf file to hydro struct
cpt.struct();

% get all required variables from netcdf file
for i=1:length(req_vars)
    current_var = req_vars{i};
%     f = strcat('hydro.', current_var, ' = ncread(', fn2, ', "', current_var, '");');
    f = ['cpt.', current_var, ' = ncread(', filename, ', "', current_var, '");'];
%     fprintf([f '\n']);
    eval(f);
end; clear i;
cpt.radiating_dof = cpt.radiating_dof';      % transpose char strings to correct orientation
cpt.influenced_dof = cpt.influenced_dof';    % transpose char strings to correct orientation

% hydro(F).Nb = info.Dimensions(getInd(info.Dimensions,'body_name')).Length;
% hydro(F).Nf = info.Dimensions(getInd(info.Dimensions,'body_name')).Length;

bodies = 


function ind = getInd(dimStruct, str2find)
    j = 0;
    for j=1:length(dimStruct)
        if dimStruct(j).Name == str2find
            ind = j;
        end
    end
end

