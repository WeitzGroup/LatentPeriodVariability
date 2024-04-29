% Code by Marian Dominguez-Mirazo, 2024
% The script ParameterInference/VirusHostParams/step2_diagnostic.jl 
% Should be run in advance as part of the parameter inference section. 
% t's output is required to plot this figure. 
clear all; close all; clc;

%% Define data ids and simulation id correspondence to CV
% List of dataset IDs, see Table 3 main text for parameter details
% See the Parameter Inference folder for details on data generation
% data1: E.coli and lambda
% data2: P.marinus and PHM2
% data3: E.hux and EhV
dataids = ["data1";"data2";"data3"];
% Simulation ids 
ids = [1,2,3,4,5,6,7];
% Corresponding CV values
cvs = [0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15];
cmap=parula(10);
%%
figure('Position',[10,10,800,800]);
tcl = tiledlayout(3,2,'TileSpacing','compact');
%% Retrieve diagnostics from files
% Loop through dataset
for j =1:size(dataids,1)
    % Current dataset ID
    dataid=dataids(j,:);
    % Read file
    file = strjoin(['../ParameterInference/VirusHostParams/step2_MCMC/round2/',dataid,'/diagnostic.csv'],'');
    tab = readtable(file, 'ReadVariableNames', false);
    tab = table2array(tab);
    rhat = tab(:,2:5);
    ESSratio = tab(:,6:9);

    %% Plot 
    
    % Rhat
    nexttile(2*(j-1)+1);
    for i = ids
        scatter(rhat(i,:),1:4,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:)); hold on;
    end
    xline(1,'LineStyle','--','Color','#808080','LineWidth',1.5);
    xlim([0.975,1.05]);
    yticks(1:4);
    yticklabels({'CV','\eta','\beta','\phi'});
    set(gca,'FontName','Latin Modern Roman');
    set(gca,'TickDir','out');
    set(gca,'TickLength',[0.025,0.025]);
    set(gca,'LineWidth',0.6);
    ax=gca;
    ax.FontSize=17;
    ylabel('Parameter','interpreter','latex','FontSize',20);

    % ESS ratio
    nexttile(2*(j-1)+2);
    for i = ids
        scatter(ESSratio(i,:),1:4,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:)); hold on;
    end
    xline(0.1,'LineStyle','--','Color','#808080','LineWidth',1.5); 
    xlim([0,0.2]);
    yticks(1:4);
    yticklabels({'CV','\eta','\beta','\phi'});
    set(gca,'FontName','Latin Modern Roman');
    set(gca,'TickDir','out');
    set(gca,'TickLength',[0.025,0.025]);
    set(gca,'LineWidth',0.6);
    ax=gca;
    ax.FontSize=17;

end

% Add x axis label to the bottom plots only
nexttile(1);
title('$\hat{R}$','interpreter','latex','FontSize',20);
nexttile(2);
title('ESS ratio','interpreter','latex','FontSize',20);

nexttile(4);
hL = legend('0.5','0.45','0.4','0.35','0.3','0.25','0.2');
title(hL,'CV');
hL.Layout.Tile = 'East';

nexttile(5);
xlabel('$\hat{R}$','interpreter','latex','FontSize',20);
nexttile(6);
xlabel('ESS ratio','interpreter','latex','FontSize',20);

%% Save figure
saveas(gcf,'../Figures/FigureS5.svg');