% Code by Marian Dominguez-Mirazo, 2024
clear all; close all; clc
%%
figure('Position',[10,10,800,700]);
tcl = tiledlayout(2,2);
%% Simulate lysogen induction

pars.mu = 1; % Irrelevant parameter
pars.K = 1;% Irrelevant parameter
pars.phi = 1; % Irrelevant parameter
pars.beta = 100; % Necessary to run but normalized afterwards
pars.eta = 1; % Necessary to run but normalized afterwards
pars.initS = 1e4; %Necessary to run but normalized afterwards
% Approximate CV (see methods)
cvs = 0.04:0.02:0.2;

options = odeset('AbsTol',1e-6,'RelTol',1e-6);
t = 0:0.001:3;

%% Run simulations
% Initialize saving arrays
firstburst=zeros(numel(cvs),1);
turbidity=zeros(numel(cvs),1);
% Loop through different various of cv
for i=1:numel(cvs)
    cv = cvs(i);
    pars.n = round(1./(cv.^2)-1);
    x0 = zeros(pars.n+3,1);
    x0(2) = pars.initS;
    [tsol,ysol] = ode45(@ODE_SEnIV,t,x0,options,pars);

    % Calculate turbidity cutoff
    lala = sum(ysol(:,1:end-1),2);
    yesidx=find(lala>1);
    lala=lala(yesidx);
    lala = lala/pars.initS;
    difflala = diff(lala);
    minlogorder = floor(log10(-min(difflala)));
    difflala(difflala<0) = difflala(difflala<0).*-1;
    % Get points with that order of change
    idx=find(log10(difflala) > minlogorder);
    turbidity(i) = tsol(idx(1));

    % Find timepoint when free virus starts to increase
    firstburst(i)=tsol(find((ysol(yesidx,end)+pars.initS) ./pars.initS -1 ...
        > pars.beta*0.0001,1));

    % Plot some runs 
    if mod(cv,0.04)==0
        nexttile(3);
        semilogy(tsol(yesidx),(ysol(yesidx,end) + pars.initS )./pars.initS,'Color','#808080'); hold on;
        scatter(firstburst(i),1,'filled','MarkerFaceColor','#D95319','Marker','diamond');

        nexttile(4);
        semilogy(tsol(yesidx),sum(ysol(yesidx,1:end-1),2),'Color','#808080'); hold on;
        scatter(turbidity(i),sum(ysol(idx(1),1:end-1),2),'filled','MarkerFaceColor','#7E2F8E','Marker','diamond');
    end
end

%% Load Dennehy and Wang data
file = '../Data/WangDennehy.csv';
virnoise = readtable(file, 'ReadVariableNames', true);
virnoise = table2array(virnoise(:,2:end));
popmean = virnoise(:,1);
popse = virnoise(:,2);
popN = virnoise(:,3);
singlemean = virnoise(:,4);
singlecv = virnoise(:,5);
singlesd = virnoise(:,6);
singleN = virnoise(:,7);
%% Sample the single cell measurement distribution assuming normality
% Number of strains
nstrains=numel(singleN);
% Number of samples 
nsamples=1000;
% Initialize storage
save_mean=zeros(nsamples,nstrains);
save_sd=zeros(nsamples,nstrains);
save_cv=zeros(nsamples,nstrains);
CI95_mean = zeros(2,nstrains);
CI95_sd = zeros(2,nstrains);
CI95_cv = zeros(2,nstrains);

%Loop by strain
for j = 1:nstrains
    % Multiple samples
    for i = 1:nsamples
        % Sample normal dist with mean, sd and number of observations
        lasample=normrnd(singlemean(j),singlesd(j),[1,singleN(j)]);
        % Save mean of sample
        save_mean(i,j)=mean(lasample);
        % Save sd of sample
        save_sd(i,j)=std(lasample);
        % Save cv of sample
        save_cv(i,j) = save_sd(i,j)/save_mean(i,j);
    end
    % Calculate 95% CI for strain 
    CI95_mean(:,j) = [quantile(save_mean(:,j),0.025); quantile(save_mean(:,j),0.975)];
    CI95_sd(:,j) = [quantile(save_sd(:,j),0.025) quantile(save_sd(:,j),0.975)];
    CI95_cv(:,j) = [quantile(save_cv(:,j),0.025) quantile(save_cv(:,j),0.975)];
end

%% Sample the turbidity measurement distribution assuming normality
% Initialize storage
save_decay =zeros(nsamples,nstrains);
CI95_decay = zeros(2,nstrains);

%Loop by strain
for j = 1:nstrains
    % Multiple samples
    for i = 1:nsamples
        % Sample normal dist with mean, sd and number of observations
        lasample=normrnd(popmean(j),popse(j)*sqrt(popN(j)-1),[1,popN(j)]);
        % Save mean of sample
        save_decay(i,j)=mean(lasample);
    end
    % Calculate 95% CI for strain 
    CI95_decay(:,j) = [quantile(save_decay(:,j),0.025); quantile(save_decay(:,j),0.975)];
end

%% Get ratio of pop and single cell estimates propagating error
% Initialize storage
ratio_95 = zeros(1,nstrains);
% Loop by strain
for j = 1:nstrains
    % SD of the population measurement
    this_popsd = std(save_decay(:,j)); %popse(j)*sqrt(popN(j)-1);
    % SD of the bootstrapped single-cell mean measurement
    this_meansd = std(save_mean(:,j));
    % SD of the ratio
    ratio_sd = sqrt((this_popsd/popmean(j))^2 + (this_meansd/singlemean(j))^2) * popmean(j)/singlemean(j);
    % Get 95% using twice the SD (cause normal distribution)
    ratio_95(j) = ratio_sd*2;
end


%% Plot theory
nexttile([1 2])

%Plot data
scatter(singlecv(1:5),popmean(1:5)./singlemean(1:5),'blue','filled','DisplayName','S107 mutants','SizeData',60,'Marker','o');  hold on;
scatter(singlecv(6:8),popmean(6:8)./singlemean(6:8),'blue','DisplayName','Other mutants','SizeData',60,'Marker','o','LineWidth',1.2);
scatter(singlecv(9),popmean(9)./singlemean(9),'red','DisplayName','Wildtype','SizeData',60,'Marker','o','LineWidth',1.2);
% Plot theory (simulations)
plot(cvs,firstburst,'Color','#D95319','DisplayName','PFU assay','LineWidth',1.5);
plot(cvs,turbidity,'Color','#7E2F8E','DisplayName','Turbidity assay','LineWidth',1.5);
% Plot erros
errorbar(singlecv,popmean./singlemean,singlecv' - CI95_cv(1,:), CI95_cv(2,:) - singlecv','horizontal','LineStyle','none','Color','k','DisplayName', '95% CI'); hold on;
errorbar(singlecv,popmean./singlemean,ratio_95 ,'LineStyle','none','Color','k','DisplayName','');

% Labels and aesthetics
legend('boxoff');
xlabel('Coefficient of Variation (CV)','FontSize',20);
ylabel({'Population to Cellular';'lysis time measurement ratio'},'Interpreter','tex','FontSize',20)
box off;
ax=gca;
ax.FontSize=17;
set(gca,'FontName','Latin Modern Roman');
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);

%%
nexttile(3);
title('PFU assay');
xlabel('Time, hr','FontSize',20);
ylabel ('PFU/ml','FontSize',20);
ylim([1,5e2]);
legend('Viral dynamics','First burst');
legend('boxoff');
box off;
ax=gca;
ax.FontSize=17;
set(gca,'FontName','Latin Modern Roman');
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);

nexttile(4);
title('Turbidity assay')
xlabel('Time, hr','FontSize',20);
ylabel ('CFU/ml','FontSize',20);
ylim([1e1 1e5]);
yticks([1e1 1e2 1e3 1e4]);
legend('Total host dynamics','Start of collapse');
legend('boxoff');
box off;
ax=gca;
ax.FontSize=17;
set(gca,'FontName','Latin Modern Roman');
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);

%%
saveas(gcf,'../Figures/FigureS2.svg');