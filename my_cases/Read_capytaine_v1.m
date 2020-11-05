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

file_bf = '"C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases/rm3.nc"';
fbf = eval(file_bf);
info_bf = ncinfo(fbf);

% file_cs = '"C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases/rm3_combo.nc"';
% fcs = eval(file_cs);
% info_cs = ncinfo(fcs);
% 
% file_cl = '"C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases/rm3_combo_long.nc"';
% fcl = eval(file_cl);
% info_cl = ncinfo(fcl);


fn2 = fbf; % fbf, fcs, cl
filename = file_bf; % file_bf, file_cs, file_cl
hydro = struct();

%%
[a,b] = size(hydro);  % Check on what is already there
if b==1 && ~isfield(hydro(b),'Nb')
    F = 1;
elseif b>=1
    F = b+1;
end

%
% p = waitbar(0,'Reading AQWA output file...'); %Progress bar

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
        error('Capytaine output does not contain: %s',req_vars{i});
    end
end

% test dofs reordering function
test = ["Surge","Heave","Yaw","Sway","Roll","Pitch"];
var = ([1:6]'*[1 3 6 2 4 5])';
nv = reset_dofs(var,test,1);


% begin parsing netcdf file to hydro struct

% Read center of gravity, center of buoyancy, displaced volume
% cg, 
% cb, 
% Vo

% Read density, gravity and water depth
hydro(F).rho = ncread(fn2,'rho');
hydro(F).g = ncread(fn2,'g');
hydro(F).h = ncread(fn2,'water_depth');

% Read number of bodies and body names
tmp = getInd(info.Dimensions,'body_name');
if tmp==0
    hydro(F).Nb = 1;
else
    hydro(F).Nb = info.Dimensions(tmp).Length;
end
tmp = ncread(fn2,'body_name')';
for i=1:hydro(F).Nb
    hydro(F).body{i} = tmp(i,:);
end

% Read number of frequencies and wave headings
hydro(F).Nf = info.Dimensions(getInd(info.Dimensions,'omega')).Length;
hydro(F).Nh = info.Dimensions(getInd(info.Dimensions,'wave_direction')).Length;

% Read frequency array, wave direction and calculate period from frequency
hydro(F).w = ncread(fn2,'omega')';
hydro(F).T = 2*pi./hydro(F).w;
hydro(F).beta = ncread(fn2,'wave_direction');

% Read number of dofs (should be 6x6)
dof_i = info.Dimensions(getInd(info.Dimensions,'influenced_dof')).Length;
dof_r = info.Dimensions(getInd(info.Dimensions,'radiating_dof')).Length;
hydro(F).dof = [dof_i, dof_r];

%% Reordering parameters
% check the ordering of the 'complex' dimension
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

% check the ordering of the 'radiating_dof' dimension
rdofs = lower(string(ncread(fn2,'radiating_dof')'));
rdofs = erase(rdofs, char(0));
if strcmp(rdofs(1), "surge") && ... 
        strcmp(rdofs(2), "sway") && ...
        strcmp(rdofs(3), "heave") && ...
        strcmp(rdofs(4), "roll") && ...
        strcmp(rdofs(5), "pitch") && ...
        strcmp(rdofs(6), "yaw")
    reset_rdofs = false;
else
    reset_rdofs = true;
end

% check the ordering of the 'influenced_dof' dimension
idofs = lower(string(ncread(fn2,'influenced_dof')'));
idofs = erase(idofs,char(0));
if strcmp(idofs(1), "surge") && ... 
        strcmp(idofs(2), "sway") && ...
        strcmp(idofs(3), "heave") && ...
        strcmp(idofs(4), "roll") && ...
        strcmp(idofs(5), "pitch") && ...
        strcmp(idofs(6), "yaw")
    reset_idofs = false;
else
    reset_idofs = true;
end

% Throw error if dofs in wrong order
% todo - create function that reorders dofs for all necessary variables
if reset_rdofs || reset_idofs
    tmp = sprintf(['Radiating or influenced dofs not in the correct order. '...
        'Rerun Capytaine with the following command before writing to Netcdf file: \n'...
        '>> sorted_dofs = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]\n'...
        '>> data = data.sel(radiating_dof=sorted_dofs, influenced_dof=sorted_dofs))']);
    error(tmp)
end
% waitbar(1/8);


%%
% Capytaine doesn't currently have hydrostatics data
% hydro(F).cg = zeros(3,hydro(F).Nb);
% hydro(F).cb = zeros(3,hydro(F).Nb);
% Vo = Displacement volume,
%waitbar(2/8);

%% Linear restoring stiffness [6, 6, Nb]
% todo - hydrostatics data from Capytaine (not currently output)
% calculate from cpt mesh?
hydro(F).C = zeros(6,6);  % Linear restoring stiffness
%waitbar(3/8);

%% Radiation added mass [6*Nb, 6*Nb, Nf]
% Get index of variable
i_var = getInd(info.Variables,'added_mass');

% get dimensions of the variable
dim = info.Variables(i_var).Dimensions;
i_infdof = getInd(dim,'influenced_dof');
i_raddof = getInd(dim,'radiating_dof');
i_w = getInd(dim,'omega');
i_bod = getInd(dim,'body_name');

% read variable
tmp = ncread(fn2,'added_mass');

% permute variable to correct dimensions if incorrect
if hydro(F).Nb == 1 || i_bod == 0
    tmp = permute(tmp,[i_infdof, i_raddof, i_w]);
else
    tmp = permute(tmp,[i_infdof, i_raddof, i_w, i_bod]);
end

% permute the influenced dof direction is not output by Capytaine correctly
if reset_idofs
    tmp = reset_dofs(tmp,idofs,1);
end
if reset_rdofs
    tmp = reset_dofs(tmp,rdofs,2);
end

hydro(F).A = size(6*hydro(F).Nb, 6*hydro(F).Nb, hydro(F).Nf);
for n=1:hydro(F).Nb
    hydro(F).A(6*(n-1)+1:6*n,6*(n-1)+1:6*n,:) = tmp(:,:,:,n); % Radiation added mass matrix
end

%% Radiation damping [6*Nb, 6*Nb, Nf]
% Get index of variable
i_var = getInd(info.Variables,'radiation_damping');

% get dimensions of the variable
dim = info.Variables(i_var).Dimensions;
i_infdof = getInd(dim,'influenced_dof');
i_raddof = getInd(dim,'radiating_dof');
i_w = getInd(dim,'omega');
i_bod = getInd(dim,'body_name');

% read variable
tmp = ncread(fn2,'radiation_damping');

% permute variable to correct dimensions if incorrect
if hydro(F).Nb == 1 || i_bod == 0
    tmp = permute(tmp,[i_infdof, i_raddof, i_w]);
else
    tmp = permute(tmp,[i_infdof, i_raddof, i_w, i_bod]);
end

% permute the influenced dof direction is not output by Capytaine correctly
if reset_idofs
    tmp = reset_dofs(tmp,idofs,1);
end
if reset_rdofs
    tmp = reset_dofs(tmp,rdofs,2);
end

hydro(F).A = size(6*hydro(F).Nb, 6*hydro(F).Nb, hydro(F).Nf);
for n=1:hydro(F).Nb
    hydro(F).A(6*(n-1)+1:6*n,6*(n-1)+1:6*n,:) = tmp(:,:,:,n); % Radiation added mass matrix
end
%waitbar(4/8);

%% Excitation Force [6*Nb,Nh,Nf];
% CHECK if the kochin_diffraction is actually the excitation force??
% i_var = getInd(info.Variables,'kochin_diffraction');
% waitbar(5/8);

%% Diffraction Force (scattering) [6*Nb,Nh,Nf];
% Get index of variable
i_var = getInd(info.Variables,'diffraction_force');

% get dimensions of the variable
dim = info.Variables(i_var).Dimensions;
i_infdof = getInd(dim,'influenced_dof');
i_dir = getInd(dim,'wave_direction');
i_w = getInd(dim,'omega');
i_bod = getInd(dim,'body_name');
i_comp = getInd(dim,'complex');

% read variable
tmp = ncread(fn2,'diffraction_force');

% permute variable to correct dimensions if incorrect
if hydro(F).Nb == 1 || i_bod == 0
    tmp = permute(tmp,[i_infdof, i_dir, i_w, i_comp]);
else
    tmp = permute(tmp,[i_infdof, i_dir, i_w, i_comp, i_bod]);
end

% permute the influenced dof direction is not output by Capytaine correctly
if reset_idofs
    tmp = reset_dofs(tmp,idofs,1); 
end

% Set real and imaginary components of variable. Calculate magnitude and
% phase from components
for n=1:hydro(F).Nb
    hydro(F).sc_re(6*(n-1)+1:6*n,:,:) = tmp(:,:,:,i_re,n);      % Real part of diffraction force
    hydro(F).sc_im(6*(n-1)+1:6*n,:,:) = tmp(:,:,:,i_im,n);      % Imaginary part of diffraction force
end
hydro(F).sc_ma = (hydro(F).sc_re.^2 + hydro(F).sc_im.^2).^0.5;  % Magnitude of diffraction force
hydro(F).sc_ph = atan(hydro(F).sc_im./hydro(F).sc_re);          % Phase of diffraction force
%waitbar(6/8);

%% Froude-Krylov force file [6*Nb,Nh,Nf];
% Get index of variable
i_var = getInd(info.Variables,'Froude_Krylov_force');

% get dimensions of the variable
dim = info.Variables(i_var).Dimensions;
i_infdof = getInd(dim,'influenced_dof');
i_dir = getInd(dim,'wave_direction');
i_w = getInd(dim,'omega');
i_bod = getInd(dim,'body_name');
i_comp = getInd(dim,'complex');

% read variable
tmp = ncread(fn2,'diffraction_force');

% permute variable to correct dimensions if incorrect
if hydro(F).Nb == 1 || i_bod == 0
    tmp = permute(tmp,[i_infdof, i_dir, i_w, i_comp]);
else
    tmp = permute(tmp,[i_infdof, i_dir, i_w, i_comp, i_bod]);
end

% permute the influenced dof direction is not output by Capytaine correctly
if reset_idofs
    tmp = reset_dofs(tmp,idofs,1);
end

% Set real and imaginary components of variable. Calculate magnitude and
% phase from components
for n=1:hydro(F).Nb
    hydro(F).fk_re(6*(n-1)+1:6*n,:,:) = tmp(:,:,:,i_re,n);      % Real part of Froude Krylov force
    hydro(F).fk_im(6*(n-1)+1:6*n,:,:) = tmp(:,:,:,i_im,n);      % Imaginary part of Froude Krylov force
end
hydro(F).fk_ma = (hydro(F).fk_re.^2 + hydro(F).fk_im.^2).^0.5;  % Magnitude of Froude Krylov force
hydro(F).fk_ph = atan(hydro(F).fk_im./hydro(F).fk_re);          % Phase of Froude Krylov force
%waitbar(7/8);


%% ================= READING KOCHIN FILES ===================%
% theta(ntheta)= Kochin(3*(ntheta-1)+1); % theta
% Kochin_BVP(ntheta,1,x)= Kochin(3*(ntheta-1)+2); % magnitude
% Kochin_BVP(ntheta,2,x)= Kochin(3*(ntheta-1)+3); % phase

% CHECK if the kochin_diffraction is actually the excitation force??

%waitbar(8/8);

hydro = Normalize(hydro);  % Normalize the data according the WAMIT convention

%% Next function calls for WAMIT, NEMOH, AQWA after call to READ_''.m
hydro = Radiation_IRF(hydro,60,[],[],[],[]);
hydro = Radiation_IRF_SS(hydro,[],[]);
hydro = Excitation_IRF(hydro,160,[],[],[],[]);
Write_H5(hydro)
Plot_BEMIO(hydro)

% % wamit
% hydro = Read_WAMIT(hydro,'rm3.out',[]);
% hydro = Radiation_IRF(hydro,20,[],[],[],[]);
% hydro = Radiation_IRF_SS(hydro,[],[]);
% hydro = Excitation_IRF(hydro,20,[],[],[],[]);
% Write_H5(hydro)
% Plot_BEMIO(hydro)

% % nemoh
% hydro = Read_NEMOH(hydro,'..\RM3\');
% % hydro = Read_WAMIT(hydro,'..\..\WAMIT\RM3\rm3.out',[]);
% % hydro = Combine_BEM(hydro); % Compare WAMIT
% hydro = Radiation_IRF(hydro,60,[],[],[],1.9);
% hydro = Radiation_IRF_SS(hydro,[],[]);
% hydro = Excitation_IRF(hydro,157,[],[],[],1.9);
% Write_H5(hydro)
% Plot_BEMIO(hydro)

% % aqwa
% hydro = Read_AQWA(hydro,'aqwa_example_data.AH1','aqwa_example_data.LIS');
% hydro = Radiation_IRF(hydro,15,[],[],[],[]);
% hydro = Radiation_IRF_SS(hydro,[],[]);
% hydro = Excitation_IRF(hydro,15,[],[],[],[]);
% Write_H5(hydro)
% Plot_BEMIO(hydro)

% these from addditional functions called in bemio.m files?
% hydro_nemoh.Ainf, % from Normalize() - simply gets A(max frequency)
% hydro_nemoh.ra_K, % radiation IRF (K), from Radiation_IRF()
% hydro_nemoh.ra_t, % radiation IRF time, from Radiation_IRF()
% hydro_nemoh.ra_w, % radiation IRF frequency, from Radiation_IRF()
% hydro_nemoh.ss_A, % from Radiation_IRF_SS()
% hydro_nemoh.ss_B, % from Radiation_IRF_SS()
% hydro_nemoh.ss_C, % from Radiation_IRF_SS()
% hydro_nemoh.ss_D, % from Radiation_IRF_SS()
% hydro_nemoh.ss_K, % from Radiation_IRF_SS()
% hydro_nemoh.ss_conv, % from Radiation_IRF_SS()
% hydro_nemoh.ss_R2, % from Radiation_IRF_SS()
% hydro_nemoh.ss_O, % from Radiation_IRF_SS()
% hydro_nemoh.ex_K, % from Excitation_IRF()
% hydro_nemoh.ex_t, % from Excitation_IRF()
% hydro_nemoh.ex_w % from Excitation_IRF()


function new_var = reset_dofs(variable, old_dofs, dim)
% this function will rearrange dimension dim of a variable if the dofs 
% are not in the correct order: surge, sway, heave, roll, pitch, yaw
new_dofs = ["surge", "sway", "heave", "roll", "pitch", "yaw"];
new_inds = zeros(6,1);
for i=1:6
    for j=1:6
        if lower(old_dofs(j)) == new_dofs(i)
            new_inds(i) = j;
            continue
        end
    end
end

new_var = zeros(size(variable));
str = [repmat(':,',[1,dim-1]) 'new_inds,' repmat(':,',[1,7]) ':'];
eval(['new_var = variable(' str ');']);

% new_var = variable(new_inds,:,:,:,:,:,:,:);

end



