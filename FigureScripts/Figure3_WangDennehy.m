% Code by Marian Dominguez-Mirazo, 2024
clear all; close all; clc
%%
figure('Position',[10,10,500,500]);

%% Load Dennehy and Wang data
file = '../Data/WangDennehy.csv';
virnoise = readtable(file, 'ReadVariableNames', true);
virnoise = table2array(virnoise(:,2:end));
% Store population-level measurements
popmean = virnoise(:,1);
popse = virnoise(:,2);
popN = virnoise(:,3);
% Store cellular-level measurements
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


for j = 1:nstrains
    % Prepare "true" means   
    single_mean_range = (0.9:0.001:1.1).*singlemean(j);
    save_mean_range=zeros(nsamples,numel(single_mean_range));
    CI95_mean_range = zeros(2,numel(single_mean_range));
    
    % Loop over a range of "true" means
    for k = 1:numel(single_mean_range)
        % Assign this experiment mean
        this_mean = single_mean_range(k);
        % Run multiple experiments
        for i = 1:nsamples
            lasample=normrnd(this_mean,singlesd(j),[1,singleN(j)]);%numsamps
            save_mean_range(i,k)=mean(lasample);
        end
        % Get the 95% CI of all experiments
        CI95_mean_range(:,k) = [quantile(save_mean_range(:,k),0.025); quantile(save_mean_range(:,k),0.975)];
    end
    
    idx = CI95_mean_range(1,:)<=singlemean(j) & CI95_mean_range(2,:)>=singlemean(j);
    single_mean_range_idx = single_mean_range(idx);
    CI95_mean(:,j) = [single_mean_range_idx(1), single_mean_range_idx(end)]; 

end
% And compare with the other methods see if there are substantital
% differences

% Prepare "true" sd   
single_sd_range = (0.9:0.001:1.1).*singlesd(j);
save_sd_range=zeros(nsamples,numel(single_sd_range));
CI95_sd_range = zeros(2,numel(single_sd_range));

% Loop over a range of "true" sds
for k = 1:numel(single_sd_range)
    % Assign this experiment sd
    this_sd = single_sd_range(k);
    % Run multiple experiments
    for i = 1:nsamples
        lasample=normrnd(singlemean(j),this_sd,[1,singleN(j)]);%numsamps
        single_sd_range(i,k)=std(lasample);
    end
    % Get the 95% CI of all experiments
    CI95_sd_range(:,k) = [quantile(single_sd_range(:,k),0.025); quantile(single_sd_range(:,k),0.975)];
end

idx = CI95_sd_range(1,:)<=singlesd(j) & CI95_sd_range(2,:)>=singlesd(j);
single_sd_range_idx = single_sd_range(idx);
CI95_sd(:,j) = [single_sd_range_idx(1), single_sd_range_idx(end)]; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ok, based on the mean CI looking so similar from each other we will just keep it as is to avoid complications....
%% Hopefully I'll still remember this some other time. 

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

%% Plot cellular to Pop level relationship
scatter(singlemean(1:5),popmean(1:5),'blue','filled','DisplayName','S107 mutants','SizeData',60,'Marker','o'); hold on;
scatter(singlemean(6:8),popmean(6:8),'blue','DisplayName','Other mutants','SizeData',60,'Marker','o','LineWidth',1.2);
scatter(singlemean(9),popmean(9),'red','DisplayName','Wildtype','SizeData',60,'Marker','o','LineWidth',1.2);
plot(25:85,25:85,'Color','#808080','LineStyle','--','LineWidth',1.7,'DisplayName', '1:1 ratio')
errorbar(singlemean,popmean,singlemean'-CI95_mean(1,:),CI95_mean(2,:)-singlemean','horizontal','LineStyle','none','Color','k','DisplayName','95% CI');
errorbar(singlemean,popmean,popmean'-CI95_decay(1,:),CI95_decay(2,:)-popmean','LineStyle','none','Color','k');


%labels and aesthetics
xlim([25,85]);
ylim([25,85]);
xlabel({'Cellular-level lysis time';'measurement (min)'},'interpreter','tex','FontSize',20);
ylabel({'Population-level lysis time';'measurement (min)'},'interpreter','tex','FontSize',20);
xticks(30:10:80);
yticks(30:10:80);
box off;
ax=gca;
ax.FontSize=19;
set(gca,'FontName','Latin Modern Roman');
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);
pbaspect([1,1,1]);
legend('Location','NorthWest');
legend('boxoff');
%%
saveas(gcf,'../Figures/Figure3.svg');