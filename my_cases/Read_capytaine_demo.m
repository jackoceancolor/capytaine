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





% list of variables to get from netcdf file
hydro = struct();
var_to_read = {'complex',...
    'omega',...
    'radiating_dof',...   % dimensions
    'theta',...
    'wave_direction',...
    'influenced_dof'...        % dimensions
    ...
    'added_mass','radiation_damping','diffraction_force',...    % variables
    'Froude_Krylov_force','kochin_diffraction','kochin',...     % variables
    };

% choose file to read
% filename = '"C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases/capytaine1.nc"';
% fn2 = eval(filename);

filename = '"C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases/rm3.nc"';
fn2 = eval(filename);

for i=1:length(var_to_read)
    current_var = var_to_read{i};
%     f = strcat('hydro.', current_var, ' = ncread(', fn2, ', "', current_var, '");');
    f = ['hydro.', current_var, ' = ncread(', filename, ', "', current_var, '");'];
%     fprintf([f '\n']);
    eval(f);
end; clear i;
        
% cpt_omega = ncread(filename, var_to_read);



