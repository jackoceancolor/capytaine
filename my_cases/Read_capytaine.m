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

file_cs = '"C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases/rm3_combo.nc"';
fcs = eval(file_cs);
info_cs = ncinfo(fcs);

file_cl = '"C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases/rm3_combo_long.nc"';
fcl = eval(file_cl);
info_cl = ncinfo(fcl);


fn2 = fcs; % fbf, fcs, cl
filename = file_cs; % file_bf, file_cs, file_cl
hydro = struct();

% BEGIN TESTING SCRIPT USING THE RM3_COMBO.NC AND RM3_COMBO_LONG.NC OUTPUT

[a,b] = size(hydro);  % Check on what is already there
if b==1
    if isfield(hydro(b),'Nb')==0;  F = 1;
    else;  F = 2;
    end
elseif b>1;  F = b+1;
end

%%
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
        error('capytaine output does not contain: %s',req_vars{i});
    end
end


test = ["Surge","Heave","Yaw","Sway","Roll","Pitch"];
var = ([1:6]'*[1 3 6 2 4 5])';
nv = reset_dofs(var,test,1);


% begin parsing netcdf file to hydro struct
cpt = struct();

% cg, 
% cb, 
% Vo
hydro(F).rho = ncread(fn2,'rho');
hydro(F).g = ncread(fn2,'g');
hydro(F).h = ncread(fn2,'water_depth');
tmp = getInd(info.Dimensions,'body_name');
if tmp==0
    hydro(F).Nb = 1;
else
    hydro(F).Nb = info.Dimensions(tmp).Length;
end
hydro(F).Nf = info.Dimensions(getInd(info.Dimensions,'omega')).Length;
hydro(F).Nh = info.Dimensions(getInd(info.Dimensions,'wave_direction')).Length;

tmp = ncread(fn2,'body_name')';
for i=1:hydro(F).Nb
    hydro(F).body{i} = tmp(i,:);
end
hydro(F).w = ncread(fn2,'omega')';
hydro(F).T = 2*pi./hydro(F).w;
hydro(F).beta = ncread(fn2,'wave_direction');

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
rdofs = lower(string(ncread(fn2,'radiating_dofs')'));
if rdofs(1)=="surge" && rdofs(2)=="sway" && rdofs(3)=="heave" && ...
        rdofs(4)=="roll" && rdofs(5)=="pitch" && rdofs(6)=="yaw"
    reset_rdofs = False;
else
    reset_rdofs = True;
end

% check the ordering of the 'influenced_dof' dimension
idofs = lower(string(ncread(fn2,'influenced_dofs')'));
if idofs(1)=="surge" && idofs(2)=="sway" && idofs(3)=="heave" && ...
        idofs(4)=="roll" && idofs(5)=="pitch" && idofs(6)=="yaw"
    reset_idofs = False;
else
    reset_idofs = True;
end

%waitbar(1/8);


%%
% Capytaine doesn't currently have hydrostatics data
% hydro(F).cg = zeros(3,hydro(F).Nb);
% hydro(F).cb = zeros(3,hydro(F).Nb);
% Vo = Displacement volume,
%waitbar(2/8);

%%
% Capytaine doesn't currently have hydrostatics data
% hydro(F).C(i,:,m) = tmp{1,1}(1:6);  % Linear restoring stiffness
%waitbar(3/8);

%%
% hydro(F).A(i,:,k) = tmp{1,1}(2:2:end);  % Added mass
% hydro(F).B(i,:,k) = tmp{1,1}(3:2:end);  % Radiation damping
%waitbar(4/8);


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
% if reset_rdof; tmp = reset_dofs(tmp,rdofs,1); end
% if reset_idof; tmp = reset_dofs(tmp,idofs,2); end
% for n=1:hydro(F).Nb
%     hydro(F).sc_re(6*(n-1)+1:6*n,:,:) = tmp(:,:,:,n,i_re);         % Real part of diffraction force
%     hydro(F).sc_im(6*(n-1)+1:6*n,:,:) = tmp(:,:,:,n,i_im);         % Imaginary part of diffraction force
% end
% hydro(F).sc_ma = (hydro(F).sc_re.^2 + hydro(F).sc_im.^2).^0.5;  % Magnitude of diffraction force
% hydro(F).sc_ph = atan(hydro(F).sc_im./hydro(F).sc_re);          % Phase of diffraction force
% %waitbar(5/8);

%% Diffraction Force (scattering) [6*Nb,Nh,Nf];
i_var = getInd(info.Variables,'diffraction_force');

dim = info.Variables(i_var).Dimensions;
i_infdof = getInd(dim,'influenced_dof');
i_dir = getInd(dim,'wave_direction');
i_w = getInd(dim,'omega');
i_bod = getInd(dim,'body_name');
i_comp = getInd(dim,'complex');

tmp = ncread(fn2,'diffraction_force');
if hydro(F).Nb == 1 || i_bod == 0
    tmp = permute(tmp,[i_infdof, i_dir, i_w, i_comp]);
else
    tmp = permute(tmp,[i_infdof, i_dir, i_w, i_comp, i_bod]);
end
if reset_idof
    tmp = reset_dofs(tmp,idofs,1); 
end

for n=1:hydro(F).Nb
    hydro(F).sc_re(6*(n-1)+1:6*n,:,:) = tmp(:,:,:,i_re,n);      % Real part of diffraction force
    hydro(F).sc_im(6*(n-1)+1:6*n,:,:) = tmp(:,:,:,i_im,n);      % Imaginary part of diffraction force
end
hydro(F).sc_ma = (hydro(F).sc_re.^2 + hydro(F).sc_im.^2).^0.5;  % Magnitude of diffraction force
hydro(F).sc_ph = atan(hydro(F).sc_im./hydro(F).sc_re);          % Phase of diffraction force
%waitbar(6/8);

%% Froude-Krylov force file [6*Nb,Nh,Nf];
i_var = getInd(info.Variables,'Froude_Krylov_force');

dim = info.Variables(i_var).Dimensions;
i_infdof = getInd(dim,'influenced_dof');
i_dir = getInd(dim,'wave_direction');
i_w = getInd(dim,'omega');
i_bod = getInd(dim,'body_name');
i_comp = getInd(dim,'complex');

tmp = ncread(fn2,'diffraction_force');
if hydro(F).Nb == 1 || i_bod == 0
    tmp = permute(tmp,[i_infdof, i_dir, i_w, i_comp]);
else
    tmp = permute(tmp,[i_infdof, i_dir, i_w, i_comp, i_bod]);
end
if reset_idof
    tmp = reset_dofs(tmp,idofs,1); 
end

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



