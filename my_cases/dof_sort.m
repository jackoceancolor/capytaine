% script to test dof sorting for GBM and B2B and Nbodies
% FIRST check that capytaine b2b output is actually right and add that into
% read_capytaine.
% SECOND check that GBM stand alone is right.
% THIRD check that its even possible to output dofs in a disordered way 
%     i.e. not (b1_surge, b1_sway, ..., b1_yaw, b1_gbm1, b1_gbm2, b2_surge, ..., b2_yaw)

rdofs00 = ["surge"; "sway"; "heave"; "roll"; "pitch"; "yaw"]; % 1 2 3 4 5 6 7
rdofs0 = ["surge"; "sway"; "heave"; "roll"; "pitch"; "yaw"; "gbm1"]; % 1 2 3 4 5 6 7
rdofs1 = ["sway"; "surge"; "pitch"; "roll"; "gbm1"; "heave"; "yaw"]; % 2 1 5 4 7 3 6

rdofs2 = ["flap__sway"; "flap__surge"; "flap__pitch";...
    "flap__roll"; "flap__heave"; "flap__yaw";...
    "base__sway"; "base__surge"; "base__pitch";...
    "base__roll"; "base__heave"; "base__yaw"]; % 2 1 5; 4 3 6; 8 7 11; 10 9 12

rdofs3 = ["flap__sway"; "base__gbm1"; "flap__surge"; "flap__pitch";...
    "base__sway"; "base__surge"; "base__pitch";...
    "flap__roll"; "flap__heave"; "flap__yaw";...
    "base__roll"; "base__heave"; "base__yaw"; "flap__gbm1"]; % 2 1 5; 8 7 11; 4 3 6; 10 9 12

rdofs4 = ["flap__surge"; "flap__sway"; "flap__heave"; "flap__roll"; "flap__pitch"; "flap__yaw";...
    "base__surge"; "base__sway"; "base__heave"; "base__roll"; "base__pitch"; "base__yaw"]; % 2 1 5; 8 7 11; 4 3 6; 10 9 12

body01 = {'flap'};
body234 = {'flap','base'};

std_dofs = ["surge", "sway", "heave", "roll", "pitch", "yaw"];

[tmp00,i00,reset00] = sorted_dof_list(rdofs00,body01);
[tmp0,i0,reset0] = sorted_dof_list(rdofs0,body01);
[tmp1,i1,reset1] = sorted_dof_list(rdofs1,body01);
[tmp2,i2,reset2] = sorted_dof_list(rdofs2,body234);
[tmp3,i3,reset3] = sorted_dof_list(rdofs3,body234);
[tmp4,i4,reset4] = sorted_dof_list(rdofs4,body234);




function [sorted_dofs,inds,reset_tf] = sorted_dof_list(old_dofs, body_names)
% 1. sort by body name: 'bodyName{i}__dofName'
% 2. sort each body's dofs by std + gbm: surge, sway, heave, roll, pitch, yaw, gbm1, gbm2, ...

% this function will rearrange dimension 'dim' of a 'variable' if the dofs
% ('old_dofs') are not in the correct order: 
%    [surge, sway, heave, roll, pitch, yaw, gbm, ...]

% list of standard dofs
std_dofs = ["surge", "sway", "heave", "roll", "pitch", "yaw"];
tmp = '__';

if length(body_names)==1
    body_names = {''};
    tmp = '';
end

sorted_dofs = [];
for k=1:length(body_names)
    body_dofs = old_dofs(contains(old_dofs,body_names{k})); % all dofs for body k
    
    std_body_dofs = strcat(body_names{k},tmp,std_dofs); % standard 6 dofs for body k
    gbm_dofs = body_dofs(~contains(body_dofs,std_body_dofs)); % any gbm dofs for body k (i.e. not in std list)
    if isempty(gbm_dofs); gbm_dofs=[]; end
    sorted_dofs = [sorted_dofs std_body_dofs gbm_dofs]; % concatenate [std(k) gbm(k) std(k+1) gbm(k+1)...]
end

% set the indices that sort the old dofs/variables into the correct order
for j=1:length(old_dofs)
    for i=1:length(sorted_dofs)
        if lower(old_dofs(j)) == sorted_dofs(i)
            inds(i) = j;
            continue
        end
    end
end

% check that inds is setup correctly. test should match sorted_dofs
% test = old_dofs(inds); 
reset_tf = any(inds~=1:length(sorted_dofs));

end




