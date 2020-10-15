
% function to start testing the reading of Capytaine output
% assuming netcdf format
hydro_wamit = load('wamit_sphere_hydro_struct.mat');
hydro_wamit = hydro_wamit.hydro;
hydro = struct;

var_to_read = {'complex','omega', 'radiating_dof',...   % dimensions
    'theta','wave_direction','influenced_dof'...        % dimensions
    ...
    'added_mass','radiation_damping','diffraction_force',...    % variables
    'Froude_Krylov_force','kochin_diffraction','kochin',...     % variables
    };
filename = '"C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases/capytaine1.nc"';
fn2 = "C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases/capytaine1.nc";

for i=1:length(var_to_read)
    current_var = var_to_read{i};
%     f = strcat('cpt_', current_var, ' = ncread(', filename, ', "', current_var, '");');
    f = ['cpt_', current_var, ' = ncread(', filename, ', "', current_var, '");'];
%     fprintf([f '\n']);
    eval(f);
end
        
% cpt_omega = ncread(filename, var_to_read);



