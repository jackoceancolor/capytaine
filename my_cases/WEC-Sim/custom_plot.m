% custom plot bemio script
set(groot,'defaultLineLineWidth',1.5)
set(groot,'defaultFigureWindowStyle','docked')

% load variables
load hydro_sphere_nemoh.mat
hydro_n = hydro;

load hydro_sphere_wamit.mat
hydro_w = hydro;

load hydro_sphere_cpt_owm.mat
hydro_c_owm = hydro;

load hydro_sphere_cpt_cgwm.mat
hydro_c_cgwm = hydro;

load hydro_sphere_cpt_cgnm.mat
hydro_c_cgnm = hydro;
clear hydro;

vars = {'hydro_n','hydro_w','hydro_c_cgnm','hydro_c_cgwm','hydro_c_owm'};
leg = {'nemoh', 'wamit', 'capy w/ cg, nemoh mesh', 'capy w/ cg, wamit mesh', 'capy w/ org, wamit mesh'};
sty = {'b--','k--','r--','g--','m--'};

%% plot
% Added Mass
for a = 1:1:3
    figure()
    hold on
    for i = 1:length(vars)
        eval(['plot(' vars{i} '.w,squeeze(' vars{i} '.A(a,a,:)),sty{i});']);
    end
    hold off
    xlabel('Frequency (1/s)');
    ylabel('Added mass');
    legend(leg{1:i});
    title(['Added mass in ' int2str(a) ' direction']);
end

%%
% Excitation force magnitude
for a = 1:2:5
    figure()
    hold on
    for i = 1:3
        eval(['plot(' vars{i} '.w,squeeze(' vars{i} '.ex_ma(a,1,:)),sty{i});']);
    end
    hold off
    xlabel('Frequency (1/s)');
    ylabel('excitation force magnitude');
    legend(leg{1:i});
    title(['Excitation force magnitude in ' int2str(a) ' direction']);
    hold on
    plot(hydro_c_cgnm.w,squeeze(hydro_c_cgnm.fk_ma(a,1,:)),'g--');
    plot(hydro_c_cgnm.w,squeeze(hydro_c_cgnm.sc_ma(a,1,:)),'m--');
end



%% Added Mass
clear X Y Legends
Fig1 = figure('Position',[50,500,975,521]);
Title = ['Normalized Added Mass: $$\bar{A}_{i,j}(\omega) = {\frac{A_{i,j}(\omega)}{\rho}}$$'];
Subtitles = {'Surge','Heave','Pitch'};
XLables = {'$$\omega (rad/s)$$','$$\omega (rad/s)$$','$$\omega (rad/s)$$'};
YLables = {'$$\bar{A}_{1,2}(\omega)$$','$$\bar{A}_{2,1}(\omega)$$','$$\bar{A}_{5,5}(\omega)$$'};
X = hydro.w;
a = 0;
for i = 1:hydro.Nb    
    m = hydro.dof(i);
    Y(1,i,:) = squeeze(hydro.A(a+1,a+2,:));
    Legends{1,i} = [hydro.body{i}];
    Y(2,i,:) = squeeze(hydro.A(a+2,a+1,:));
    Legends{2,i} = [hydro.body{i}];
    Y(3,i,:) = squeeze(hydro.A(a+5,a+5,:));
    Legends{3,i} = [hydro.body{i}];
    a = a + m;
end
Notes = {'Notes:',...
    ['$$\bullet$$ $$\bar{A}_{i,j}(\omega)$$ should tend towards a constant, ',...
    '$$A_{\infty}$$, within the specified $$\omega$$ range.'],...
    ['$$\bullet$$ Only $$\bar{A}_{i,j}(\omega)$$ for the surge, heave, and ',...
    'pitch DOFs are plotted here. If another DOF is significant to the system, ',...
    'that $$\bar{A}_{i,j}(\omega)$$ should also be plotted and verified before ',...
    'proceeding.']};
FormatPlot(Fig1,Title,Subtitles,XLables,YLables,X,Y,Legends,Notes)








%% Format
function FormatPlot(fig,heading,subtitle,x_lables,y_lables,X_data,Y_data,legends,notes)

axes1 = axes('Parent',fig,'Position',[0.0731 0.3645 0.2521 0.4720]);
hold(axes1,'on');
box(axes1,'on');
title(subtitle(1));
xlabel(x_lables(1),'Interpreter','latex');
ylabel(y_lables(1),'Interpreter','latex');

axes2 = axes('Parent',fig,'Position',[0.3983 0.3645 0.2521 0.4720]);
hold(axes2,'on');
box(axes2,'on');
title(subtitle(2));
xlabel(x_lables(2),'Interpreter','latex');
ylabel(y_lables(2),'Interpreter','latex');

axes3 = axes('Parent',fig,'Position',[0.7235 0.3645 0.2521 0.4720]);
hold(axes3,'on');
box(axes3,'on');
title(subtitle(3));
xlabel(x_lables(3),'Interpreter','latex');
ylabel(y_lables(3),'Interpreter','latex');

annotation(fig,'textbox',[0.0 0.9 1.0 0.1],...
    'String',heading,...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'FontWeight','bold',...
    'FontSize',12,...
    'EdgeColor','none');
annotation(fig,'textbox',[0.0 0.0 1.0 0.2628],...
    'String',notes,...
    'Interpreter','latex',...
    'FitBoxToText','off',...
    'EdgeColor','none');

[p,b,s]=size(Y_data);
for i = 1:b
    plot(X_data,squeeze(Y_data(1,i,:)),'LineWidth',1,'Parent',axes1);
    plot(X_data,squeeze(Y_data(2,i,:)),'LineWidth',1,'Parent',axes2);
    plot(X_data,squeeze(Y_data(3,i,:)),'LineWidth',1,'Parent',axes3);
end
legend(axes1,legends(1,:),'location','best','Box','off','Interpreter','none')
legend(axes2,legends(2,:),'location','best','Box','off','Interpreter','none')
legend(axes3,legends(3,:),'location','best','Box','off','Interpreter','none')

end