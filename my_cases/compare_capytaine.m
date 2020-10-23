% function to compare Capytaine output with another BEM code. Run after
% successfully creating a hydro struct with the Read_Capytaine.m script.
%
% Adam Keester, 10/22/2020
hydro_wamit = load('hydro_rm3_wamit.mat');
hydro_wamit = hydro_wamit.hydro;
hydro_nemoh = load('hydro_rm3_nemoh.mat');
hydro_nemoh = hydro_nemoh.hydro;

file_cs = '"C:/Users/akeeste/Documents/Software/GitHub/capytaine/my_cases/rm3_combo.nc"';
fcs = eval(file_cs);
info_cs = ncinfo(fcs);


a_w = hydro_wamit.A(:,:,[1 130 260]);
a_n = hydro_nemoh.A(:,:,[1 130 260]);

i_dof = ncread(fcs,'influenced_dof')';
r_dof = ncread(fcs,'radiating_dof')';
% 
% a_c = ncread(fcs,'added_mass');
% tmp = a_c;
% tmp(1:6,1:6,:) = a_c(7:12,7:12,:);
% tmp(7:12,7:12,:) = a_c(1:6,1:6,:);
% tmp(1:6,7:12,:) = a_c(7:12,1:6,:);
% tmp(7:12,1:6,:) = a_c(1:6,7:12,:);
% a_c = tmp;

tmp = ncread(fbf,'added_mass');
tmp_b = tmp(:,:,:,1);
tmp_p = tmp(:,:,:,2);
z = zeros([12,12,3]);
a_c = z;
a_c(1:6,1:6,:) = tmp_b;
a_c(7:12,7:12,:) = tmp_p;

for iw = 1:3
    figure(iw)
    subplot(1,3,1);
    heatmap(a_w(:,:,iw));
    caxis([0 1e4]);
    title('wamit');

    subplot(1,3,2);
    heatmap(a_n(:,:,iw));
    caxis([0 1e4]);
    title('nemoh');

    subplot(1,3,3);
    heatmap(a_c(:,:,iw));
    caxis([0 1e4]);
    title('capytaine');
end

