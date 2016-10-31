task = 'ListSQ';
imageX = 800;
imageY = 600;

skagg_prctile = [];
corr_half = [];
corr_half_prctile = [];
corr_even_odds = [];
corr_even_odds_prctile = [];
trial_count = [];

max_firing_rate = [];
which_monkey = [];

for monkey = 1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML[
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    
    for session = 1:length(session_data)
        
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{session},task);
        if isempty(task_file)
            continue
        end
        
        disp(task_file(1:8))
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials','stability_attribute','cfg','hdr','data');
        
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        
        if num_units == 0
            continue
        end
        
        if exist([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],'file')
            load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],...
                'spatial_info','nvr','peak_firing_rate')
        else
            continue
        end
        
        for unit = 1:num_units
            
            %unsable units
            if multiunit(unit) ==1
                continue
            end
            
            if isnan(spatial_info.rate(unit)) || spatial_info.rate(unit) == 0 %didn't run on unit since not stable
                continue
            end
            
            
            trial_count = [trial_count length(nvr{unit})];
            max_firing_rate = [max_firing_rate peak_firing_rate(3,unit)]; %all image firng rate
            which_monkey = [which_monkey monkey];
            
            
            skagg_prctile = [skagg_prctile spatial_info.shuffled_rate_prctile(unit)];
            corr_half = [corr_half spatial_info.spatialstability_halves(:,unit)];
            corr_half_prctile = [corr_half_prctile spatial_info.spatialstability_halves_prctile(:,unit)];
            corr_even_odds = [corr_even_odds spatial_info.spatialstability_even_odd(:,unit)];
            corr_even_odds_prctile = [corr_even_odds_prctile spatial_info.spatialstability_even_odd_prctile(:,unit)];
        end
    end
end

%%
%---Figure for Skaggs vs Number of Images and Peak Firing Rate---%
figure
subplot(2,3,1)
plot(skagg_prctile,max_firing_rate,'k.','markersize',7)
xlabel('Skagg Score %')
ylabel('Peak Firing Rate (Hz)')
title('Skagg Score % vs. Peak Firing Rate')

subplot(2,3,2)
plot(skagg_prctile,trial_count,'k.','markersize',7)
xlabel('Skagg Score %')
ylabel('# of Images')
title('Skagg Score % vs. # of Images')

%just see if there's a spurious relationship
subplot(2,3,3)
plot(max_firing_rate,trial_count,'k.','markersize',7)
xlabel('Peak Firing Rate (Hz)')
ylabel('# of Trials')
title('Peak Firing Rate vs. # of Trials')

subplot(2,3,4)
hist(skagg_prctile,0:2:100);
hold on
yl = ylim;
plot([95 95],[0 yl(2)],'k--')
xlabel('Skagg Score %')
ylabel('# of Neurons')
title('Distribution of Skagg Scores')
xlim([-1 101])

subplot(2,3,5)
hist(max_firing_rate,25);
xlabel('Peak Firing Rate (Hz)')
ylabel('# of Neurons')
title('Distribution of Peak Firing Rates')

subplot(2,3,6)
hist(trial_count,25);
xlabel('# of Images')
ylabel('# of Neurons')
title('Distribution of # of Images')


%---Figure for Different Correlations---%

figure
%Correlation values for first and 2nd half
subplot(2,3,1)
hold on
plot([-1 1],[-1 1],'k--')
plot(corr_half(1,:),corr_half(2,:),'b.','markersize',7)
hold off
xlabel('Spearman''s \rho')
ylabel('Kendall''s \tau')
title('Corr. 1st vs 2nd half')

subplot(2,3,4)
hold on
plot([-0 100],[0 100],'k--')
plot( corr_half_prctile(1,:), corr_half_prctile(2,:),'b.','markersize',7)
hold off
xlabel('Spearman''s \rho %')
ylabel('Kendall''s \tau %')
title('Corr. Percentile 1st vs 2nd half')

%Correlation values for even and odd trials
subplot(2,3,2)
hold on
plot([-1 1],[-1 1],'k--')
plot(corr_even_odds(1,:),corr_even_odds(2,:),'b.','markersize',7)
hold off
xlabel('Spearman''s \rho')
ylabel('Kendall''s \tau')
title('Corr. Even/Odd Trials')

subplot(2,3,5)
hold on
plot([-0 100],[0 100],'k--')
plot(corr_even_odds_prctile(1,:),corr_even_odds_prctile(2,:),'b.','markersize',7)
hold off
xlabel('Spearman''s \rho %')
ylabel('Kendall''s \tau %')
title('Corr. Percentile  Even/Odd Trials')

%Plots for 1st/2nd half vs Even/Odd trial correlations
subplot(2,3,3)
hold on
plot([-1 1],[-1 1],'k--')
p(1) = plot(corr_even_odds(1,:),corr_half(1,:),'.b');
p(2) = plot(corr_even_odds(2,:),corr_half(2,:),'.r');
hold off
legend(p,'Spearman','Kendal','Location','Best')
xlabel('Corr. Even/Odd Trials')
ylabel('Corr. 1st vs 2nd half')

subplot(2,3,6)
hold on
plot([-0 100],[0 100],'k--')
p(1) = plot(corr_even_odds_prctile(1,:),corr_half_prctile(1,:),'.b');
p(2) = plot(corr_even_odds_prctile(2,:),corr_half_prctile(2,:),'.r');
hold off
legend(p,'Spearman','Kendal','Location','Best')
xlabel('Corr. Percentile Even/Odd Trials')
ylabel('Corr. Percentile 1st vs 2nd half')

%---Figure for Correlations vs Skaggs---%

figure
subplot(2,2,1)
plot(skagg_prctile,corr_half_prctile(1,:),'.b');
xlabel('Skagg Percentile')
ylabel('Spearman''s \rho Percentile')
title('Skagg vs 1st/2nd Half Correlation')

subplot(2,2,3)
plot(skagg_prctile,corr_half_prctile(1,:),'.b');
xlabel('Skagg Percentile')
ylabel('Kendall''s \tau Percentile')

subplot(2,2,2)
plot(skagg_prctile,corr_even_odds_prctile(1,:),'.b');
xlabel('Skagg Percentile')
ylabel('Spearman''s \rho Percentile')
title('Skagg vs Even/Odd Correlation')

subplot(2,2,4)
plot(skagg_prctile,corr_even_odds_prctile(1,:),'.b');
xlabel('Skagg Percentile')
ylabel('Kendall''s \tau Percentile')

%%
%---Show Distribution of Correlations for each potential (~)place cell category--%

%spearmans stuff
cat = cell(1,9);
cat{1} = find((corr_half_prctile(1,:) > 95) & (skagg_prctile > 95) & (corr_even_odds_prctile(1,:) > 95));%ultra conservative
cat{2} = find((corr_half_prctile(1,:) > 95) & (corr_even_odds_prctile(1,:) > 95));%half + even/odds
cat{3} = find((corr_half_prctile(1,:) > 95) & (skagg_prctile > 95));%half + skaggs
cat{4} = find((skagg_prctile > 95) & (corr_even_odds_prctile(1,:) > 95));%even/odds +skaggs
cat{5} = find((corr_half_prctile(1,:) > 95));% %half
cat{6} = find((corr_even_odds_prctile(1,:) > 95));%even odds
cat{7} = find(skagg_prctile > 99); %skaggs 99%
cat{8} = find(skagg_prctile > 95);%skaggs 95%
cat{9} = find((corr_half_prctile(1,:) < 95) & (skagg_prctile < 95) & (corr_even_odds_prctile(1,:) < 95));%nothing useful

figure
subplot(2,2,1)
hold on
for c = 1:9
plot(c,skagg_prctile(cat{c}),'.k')
plot(c,nanmean(skagg_prctile(cat{c})),'ok','markersize',7)
end
set(gca,'Xtick',1:9)
set(gca,'XtickLabel',{'All','\rho_{1/2}+\rho_{e/o}','\rho_{1/2}+Skg_{95}','\rho_{e/o}+Skg_{95}',...
    '\rho_{1/2}','\rho_{e/o}','Skg_{99}','Skg_{95}','~sig'})
set(gca,'XTickLabelRotation',-45);
xlim([0.75 9.25])
ylabel('Skaggs Percentile')

subplot(2,2,2)
hold on
for c = 1:9
plot(c,corr_half(1,cat{c}),'.k')
plot(c,nanmean(corr_half(1,cat{c})),'ok','markersize',7)
end
set(gca,'Xtick',1:9)
set(gca,'XtickLabel',{'All','\rho_{1/2}+\rho_{e/o}','\rho_{1/2}+Skg_{95}','\rho_{e/o}+Skg_{95}',...
    '\rho_{1/2}','\rho_{e/o}','Skg_{99}','Skg_{95}','~sig'})
set(gca,'XTickLabelRotation',-45);
xlim([0.75 9.25])
ylabel('\rho 1st vs 2nd half')

subplot(2,2,3)
hold on
for c = 1:9
plot(c,corr_even_odds(1,cat{c}),'.k')
plot(c,nanmean(corr_even_odds(1,cat{c})),'ok','markersize',7)
end
set(gca,'Xtick',1:9)
set(gca,'XtickLabel',{'All','\rho_{1/2}+\rho_{e/o}','\rho_{1/2}+Skg_{95}','\rho_{e/o}+Skg_{95}',...
    '\rho_{1/2}','\rho_{e/o}','Skg_{99}','Skg_{95}','~sig'})
set(gca,'XTickLabelRotation',-45);
xlim([0.75 9.25])
ylabel('\rho Even/Odd Trials')

subplot(2,2,4)
hold on
for c = 1:9
plot(c,max_firing_rate(cat{c}),'.k')
plot(c,nanmean(max_firing_rate(cat{c})),'ok','markersize',7)
end
set(gca,'Xtick',1:9)
set(gca,'XtickLabel',{'All','\rho_{1/2}+\rho_{e/o}','\rho_{1/2}+Skg_{95}','\rho_{e/o}+Skg_{95}',...
    '\rho_{1/2}','\rho_{e/o}','Skg_{99}','Skg_{95}','~sig'})
set(gca,'XTickLabelRotation',-45);
xlim([0.75 9.25])
ylabel('Peak Firing Rate (Hz)')

subtitle('Spearman''s Categorical Spatial Results')

%Kendalls stuff
cat = cell(1,9);
cat{1} = find((corr_half_prctile(2,:) > 95) & (skagg_prctile > 95) & (corr_even_odds_prctile(2,:) > 95));%ultra conservative
cat{2} = find((corr_half_prctile(2,:) > 95) & (corr_even_odds_prctile(2,:) > 95));%half + even/odds
cat{3} = find((corr_half_prctile(2,:) > 95) & (skagg_prctile > 95));%half + skaggs
cat{4} = find((skagg_prctile > 95) & (corr_even_odds_prctile(2,:) > 95));%even/odds +skaggs
cat{5} = find((corr_half_prctile(2,:) > 95));% %half
cat{6} = find((corr_even_odds_prctile(2,:) > 95));%even odds
cat{7} = find(skagg_prctile > 99); %skaggs 99%
cat{8} = find(skagg_prctile > 95);%skaggs 95%
cat{9} = find((corr_half_prctile(2,:) < 95) & (skagg_prctile < 95) & (corr_even_odds_prctile(2,:) < 95));%nothing useful

figure
subplot(2,2,1)
hold on
for c = 1:9
plot(c,skagg_prctile(cat{c}),'.k')
plot(c,nanmean(skagg_prctile(cat{c})),'ok','markersize',7)
end
set(gca,'Xtick',1:9)
set(gca,'XtickLabel',{'All','\tau_{1/2}+\tau_{e/o}','\tau_{1/2}+Skg_{95}','\tau_{e/o}+Skg_{95}',...
    '\tau_{1/2}','\tau_{e/o}','Skg_{99}','Skg_{95}','~sig'})
xlim([0.75 9.25])
set(gca,'XTickLabelRotation',-45);
ylabel('Skaggs Percentile')

subplot(2,2,2)
hold on
for c = 1:9
plot(c,corr_half(2,cat{c}),'.k')
plot(c,nanmean(corr_half(2,cat{c})),'ok','markersize',7)
end
set(gca,'Xtick',1:9)
set(gca,'XtickLabel',{'All','\tau_{1/2}+\tau_{e/o}','\tau_{1/2}+Skg_{95}','\tau_{e/o}+Skg_{95}',...
    '\tau_{1/2}','\tau_{e/o}','Skg_{99}','Skg_{95}','~sig'})
set(gca,'XTickLabelRotation',-45);
xlim([0.75 9.25])
ylabel('\tau 1st vs 2nd half')

subplot(2,2,3)
hold on
for c = 1:9
plot(c,corr_even_odds(2,cat{c}),'.k')
plot(c,nanmean(corr_even_odds(2,cat{c})),'ok','markersize',7)
end
set(gca,'Xtick',1:9)
set(gca,'XtickLabel',{'All','\tau_{1/2}+\tau_{e/o}','\tau_{1/2}+Skg_{95}','\tau_{e/o}+Skg_{95}',...
    '\tau_{1/2}','\tau_{e/o}','Skg_{99}','Skg_{95}','~sig'})
set(gca,'XTickLabelRotation',-45);
xlim([0.75 9.25])
ylabel('\tau Even/Odd Trials')

subplot(2,2,4)
hold on
for c = 1:9
plot(c,max_firing_rate(cat{c}),'.k')
plot(c,nanmean(max_firing_rate(cat{c})),'ok','markersize',7)
end
set(gca,'Xtick',1:9)
set(gca,'XtickLabel',{'All','\tau_{1/2}+\tau_{e/o}','\tau_{1/2}+Skg_{95}','\tau_{e/o}+Skg_{95}',...
    '\tau_{1/2}','\tau_{e/o}','Skg_{99}','Skg_{95}','~sig'})
set(gca,'XTickLabelRotation',-45);
xlim([0.75 9.25])
ylabel('Peak Firing Rate (Hz)')

subtitle('Kendalls''s Categorical Spatial Results')
%%
%---Make a Table showing the summary data---%

%neuron counts
total_neurons = length(trial_count);
PW_total_neurons = sum(which_monkey == 1);
TO_total_neurons = sum(which_monkey == 2);

%skagg score stuff
skagg_95 = sum(skagg_prctile > 95);
PW_skagg_95 = sum((skagg_prctile > 95) & (which_monkey == 1));
TO_skagg_95 = sum((skagg_prctile > 95) & (which_monkey == 2));
skagg_99 = sum(skagg_prctile > 99);
PW_skagg_99 = sum((skagg_prctile > 99) & (which_monkey == 1));
TO_skagg_99 = sum((skagg_prctile > 99) & (which_monkey == 2));

%first vs second half correlations
Sprmn_half_95 = sum(corr_half_prctile(1,:) > 95);
PW_Sprmn_half_95 = sum((corr_half_prctile(1,:) > 95) & (which_monkey == 1));
TO_Sprmn_half_95 = sum((corr_half_prctile(1,:) > 95) & (which_monkey == 2));
Kendall_half_95 = sum(corr_half_prctile(2,:) > 95);
PW_Kendall_half_95 = sum((corr_half_prctile(2,:) > 95) & (which_monkey == 1));
TO_Kendall_half_95 = sum((corr_half_prctile(2,:) > 95) & (which_monkey == 2));

%even/odd trial correlations
Sprmn_even_odds_95 = sum(corr_even_odds_prctile(1,:) > 95);
PW_Sprmn_even_odds_95 = sum((corr_even_odds_prctile(1,:) > 95) & (which_monkey == 1));
TO_Sprmn_even_odds_95 = sum((corr_even_odds_prctile(1,:) > 95) & (which_monkey == 2));
Kendall_even_odds_95 = sum(corr_even_odds_prctile(2,:) > 95);
PW_Kendall_even_odds_95 = sum((corr_even_odds_prctile(2,:) > 95) & (which_monkey == 1));
TO_Kendall_even_odds_95 = sum((corr_even_odds_prctile(2,:) > 95) & (which_monkey == 2));

%correlations and skaggs
skagg_95_Sprmn_95_half = sum((corr_half_prctile(1,:) > 95) & (skagg_prctile > 95));
PW_skagg_95_Sprmn_95_half = sum((corr_half_prctile(1,:) > 95) & (skagg_prctile > 95) & (which_monkey == 1));
TO_skagg_95_Sprmn_95_half = sum((corr_half_prctile(1,:) > 95) & (skagg_prctile > 95) & (which_monkey == 2));
skagg_95_Kendall_95_half = sum((corr_half_prctile(2,:) > 95) & (skagg_prctile > 95));
PW_skagg_95_Kendall_95_half = sum((corr_half_prctile(2,:) > 95) & (skagg_prctile > 95) & (which_monkey == 1));
TO_skagg_95_Kendall_95_half = sum((corr_half_prctile(2,:) > 95) & (skagg_prctile > 95) & (which_monkey == 2));
skagg_95_Sprmn_95_even_odds = sum((corr_even_odds_prctile(1,:) > 95) & (skagg_prctile > 95));
PW_skagg_95_Sprmn_95_even_odds = sum((corr_even_odds_prctile(1,:) > 95) & (skagg_prctile > 95) & (which_monkey == 1));
TO_skagg_95_Sprmn_95_even_odds = sum((corr_even_odds_prctile(1,:) > 95) & (skagg_prctile > 95) & (which_monkey == 2));
skagg_95_Kendall_95_even_odds = sum((corr_even_odds_prctile(2,:) > 95) & (skagg_prctile > 95));
PW_skagg_95_Kendall_95_even_odds = sum((corr_even_odds_prctile(2,:) > 95) & (skagg_prctile > 95) & (which_monkey == 1));
TO_skagg_95_Kendall_95_even_odds = sum((corr_even_odds_prctile(2,:) > 95) & (skagg_prctile > 95) & (which_monkey == 2));

%ultra conservative all 3 measures significant
skagg_95_Sprmn_95_half_even_odds = sum((corr_half_prctile(1,:) > 95) & (skagg_prctile > 95) ...
     & (corr_even_odds_prctile(1,:) > 95));
PW_skagg_95_Sprmn_95_half_even_odds = sum((corr_half_prctile(1,:) > 95) & (skagg_prctile > 95) ...
     & (corr_even_odds_prctile(1,:) > 95) & (which_monkey == 1));
TO_skagg_95_Sprmn_95_half_even_odds = sum((corr_half_prctile(1,:) > 95) & (skagg_prctile > 95) ...
     & (corr_even_odds_prctile(1,:) > 95) & (which_monkey == 2));
skagg_95_Kendall_95_half_even_odds = sum((corr_half_prctile(2,:) > 95) & (skagg_prctile > 95) ...
     & (corr_even_odds_prctile(2,:) > 95));
PW_skagg_95_Kendall_95_half_even_odds = sum((corr_half_prctile(2,:) > 95) & (skagg_prctile > 95) ...
     & (corr_even_odds_prctile(2,:) > 95) & (which_monkey == 1));
TO_skagg_95_Kendall_95_half_even_odds = sum((corr_half_prctile(2,:) > 95) & (skagg_prctile > 95) ...
     & (corr_even_odds_prctile(2,:) > 95) & (which_monkey == 2));

%really strong correlations but no skaggs
no_skagg_Sprmn_99_half = sum((corr_half_prctile(1,:) > 99) & (skagg_prctile < 95));
no_skagg_Kendall_99_half = sum((corr_half_prctile(2,:) > 99) & (skagg_prctile < 95));
no_skagg_Sprmn_99_even_odds = sum((corr_even_odds_prctile(1,:) > 99) & (skagg_prctile < 95));
no_skagg_Kendall_99_even_odds = sum((corr_even_odds_prctile(2,:) > 99) & (skagg_prctile < 95));

%strong correlations but no skaggs
no_skagg_Sprmn_95 = sum((corr_half_prctile(1,:) > 95) & (corr_even_odds_prctile(1,:) > 95) & (skagg_prctile < 95));
no_skagg_Kendall_95 = sum((corr_half_prctile(2,:) > 95) & (corr_even_odds_prctile(2,:) > 95) & (skagg_prctile < 95));



clc
fprintf('%%------------------Spatial Analysis Results------------------%% \n')
%neuron counts
fprintf(['Total Units: \t' num2str(total_neurons) '\n'])
fprintf(['\t\t\t PW:' num2str(PW_total_neurons) '\t TO: ' num2str(TO_total_neurons) '\n'])
fprintf('\n')

%skagg only stuff
fprintf(['Skaggs 95%%: ' num2str(skagg_95) ' (' num2str(round(100*skagg_95/total_neurons),2) '%%)\n'])
fprintf(['\t\t\t PW:' num2str(PW_skagg_95) ' (' num2str(round(100*PW_skagg_95/PW_total_neurons),2) '%%)'...
    '\t\t TO:' num2str(TO_skagg_95) ' (' num2str(round(100*TO_skagg_95/TO_total_neurons),2) '%%)\n'])
fprintf(['Skaggs 99%%: ' num2str(skagg_99) ' (' num2str(round(100*skagg_99/total_neurons),2) '%%)\n'])
fprintf(['\t\t\t PW:' num2str(PW_skagg_99) ' (' num2str(round(100*PW_skagg_99/PW_total_neurons),2) '%%)'...
    '\t\t TO:' num2str(TO_skagg_99) ' (' num2str(round(100*TO_skagg_99/TO_total_neurons),2) '%%)\n'])
fprintf('\n')

%correlation only stuff
fprintf(['Corr. 1st/2nd Half Spearman: ' num2str(Sprmn_half_95) '(' ...
    num2str(round(100*Sprmn_half_95/total_neurons),2) '%%)\n'])
fprintf(['\t\t\t PW:' num2str(PW_Sprmn_half_95) ' (' num2str(round(100*PW_Sprmn_half_95/PW_total_neurons),2) '%%)'...
    '\t\t TO:' num2str(TO_Sprmn_half_95) ' (' num2str(round(100*TO_Sprmn_half_95/TO_total_neurons),2) '%%)\n'])
fprintf(['Corr. 1st/2nd Half Kendall:  ' num2str(Kendall_half_95) ...
    '('  num2str(round(100*Kendall_half_95/total_neurons),2) '%%)\n'])
fprintf(['\t\t\t PW:' num2str(PW_Kendall_half_95) ' (' num2str(round(100*PW_Kendall_half_95/PW_total_neurons),2) '%%)'...
    '\t\t TO:' num2str(TO_Kendall_half_95) ' (' num2str(round(100*TO_Kendall_half_95/TO_total_neurons),2) '%%)\n'])
fprintf('\n')
fprintf(['Corr. Even/Odd Trials Spearman: ' num2str(Sprmn_even_odds_95) '(' ...
    num2str(round(100*Sprmn_even_odds_95/total_neurons),2) '%%)\n'])
fprintf(['\t\t\t PW:' num2str(PW_Sprmn_even_odds_95) ' (' num2str(round(100*PW_Sprmn_even_odds_95/PW_total_neurons),2) '%%)'...
    '\t\t TO:' num2str(TO_Sprmn_even_odds_95) ' (' num2str(round(100*TO_Sprmn_even_odds_95/TO_total_neurons),2) '%%)\n'])
fprintf(['Corr. Even/Odd Trials Kendall:  ' num2str(Kendall_even_odds_95) ...
    '('  num2str(round(100*Kendall_even_odds_95/total_neurons),2) '%%)\n'])
fprintf(['\t\t\t PW:' num2str(PW_Kendall_even_odds_95) ' (' num2str(round(100*PW_Kendall_even_odds_95/PW_total_neurons),2) '%%)'...
    '\t\t TO:' num2str(TO_Kendall_even_odds_95) ' (' num2str(round(100*TO_Kendall_even_odds_95/TO_total_neurons),2) '%%)\n'])
fprintf('\n\n')

%correlations and skaggs
fprintf(['Spearman Corr 1st/2nd Half & Skaggs 95%%: ' num2str(skagg_95_Sprmn_95_half) '('... 
    num2str(round(100*skagg_95_Sprmn_95_half/total_neurons)) '%%)\n'])
fprintf(['\t\t\t PW:' num2str(PW_skagg_95_Sprmn_95_half) ' (' num2str(round(100*PW_skagg_95_Sprmn_95_half/PW_total_neurons),2) '%%)'...
    '\t\t TO:' num2str(TO_skagg_95_Sprmn_95_half) ' (' num2str(round(100*TO_skagg_95_Sprmn_95_half/TO_total_neurons),2) '%%)\n'])
fprintf(['Kendall Corr 1st/2nd Half & Skaggs 95%%: ' num2str(skagg_95_Kendall_95_half) '('... 
    num2str(round(100*skagg_95_Kendall_95_half/total_neurons)) '%%)\n'])
fprintf(['\t\t\t PW:' num2str(PW_skagg_95_Kendall_95_half) ' (' num2str(round(100*PW_skagg_95_Kendall_95_half/PW_total_neurons),2) '%%)'...
    '\t\t TO:' num2str(TO_skagg_95_Kendall_95_half) ' (' num2str(round(100*TO_skagg_95_Kendall_95_half/TO_total_neurons),2) '%%)\n'])
fprintf(['Spearman Corr Even/Odd Trials & Skaggs 95%%: ' num2str(skagg_95_Sprmn_95_even_odds) '(' ...
    num2str(round(100*skagg_95_Sprmn_95_even_odds/total_neurons)) '%%)\n'])
fprintf(['\t\t\t PW:' num2str(PW_skagg_95_Sprmn_95_even_odds) ' (' num2str(round(100*PW_skagg_95_Sprmn_95_even_odds/PW_total_neurons),2) '%%)'...
    '\t\t TO:' num2str(TO_skagg_95_Sprmn_95_even_odds) ' (' num2str(round(100*TO_skagg_95_Sprmn_95_even_odds/TO_total_neurons),2) '%%)\n'])
fprintf(['Kendall Corr Even/Odd Trials & Skaggs 95%%: ' num2str(skagg_95_Kendall_95_even_odds) '('...
    num2str(round(100*skagg_95_Kendall_95_even_odds/total_neurons)) '%%)\n'])
fprintf(['\t\t\t PW:' num2str(PW_skagg_95_Kendall_95_even_odds) ' (' num2str(round(100*PW_skagg_95_Kendall_95_even_odds/PW_total_neurons),2) '%%)'...
    '\t\t TO:' num2str(TO_skagg_95_Kendall_95_even_odds) ' (' num2str(round(100*TO_skagg_95_Kendall_95_even_odds/TO_total_neurons),2) '%%)\n'])
fprintf('\n\n')

%really strong correlations without skaggs
fprintf(['Spearman Corr 1st/2nd Half 99%% & ~Skaggs: ' num2str(no_skagg_Sprmn_99_half) '('... 
    num2str(round(100*no_skagg_Sprmn_99_half/total_neurons)) '%%)\n'])
fprintf(['Kendall Corr 1st/2nd Half 99%% & ~Skaggs: ' num2str(no_skagg_Kendall_99_half) '('... 
    num2str(round(100*no_skagg_Kendall_99_half/total_neurons)) '%%)\n'])
fprintf(['Spearman Corr Even/Odd Trials 99%% & ~Skaggs: ' num2str(no_skagg_Sprmn_99_even_odds) '(' ...
    num2str(round(100*no_skagg_Sprmn_99_even_odds/total_neurons)) '%%)\n'])
fprintf(['Kendall Corr Even/Odd Trials 99%% & ~Skaggs: ' num2str(no_skagg_Kendall_99_even_odds) '('...
    num2str(round(100*no_skagg_Kendall_99_even_odds/total_neurons)) '%%)\n'])
%half and even/odd correlations without skaggs
fprintf(['Spearman Corr 1st/2nd Half & Even/Odds 95%% ~Skaggs: ' num2str(no_skagg_Sprmn_95) ...
    '(' num2str(round(100*no_skagg_Sprmn_95/total_neurons)) '%%)\n']);
fprintf(['Kendall Corr 1st/2nd Half & Even/Odds 95%% ~Skaggs: ' num2str(no_skagg_Kendall_95) ...
    '(' num2str(round(100*no_skagg_Kendall_95/total_neurons)) '%%)\n']);
fprintf('\n\n')

%ultra conservative passing all 3 measures
fprintf(['Spearman Corr Half & Even/Odd, & Skagg 95%%: ' num2str(skagg_95_Sprmn_95_half_even_odds) ...
    '(' num2str(round(100*skagg_95_Sprmn_95_half_even_odds/total_neurons)) '%%) \n'])
fprintf(['\t\t\t PW:' num2str(PW_skagg_95_Sprmn_95_half_even_odds) ...
    ' (' num2str(round(100*PW_skagg_95_Sprmn_95_half_even_odds/PW_total_neurons),2) '%%)'...
    '\t\t TO:' num2str(TO_skagg_95_Sprmn_95_half_even_odds)...
    ' (' num2str(round(100*TO_skagg_95_Sprmn_95_half_even_odds/TO_total_neurons),2) '%%)\n'])
fprintf(['Kendall Corr Half & Even/Odd, & Skagg 95%%: ' num2str(skagg_95_Kendall_95_half_even_odds) ...
    '(' num2str(round(100*skagg_95_Kendall_95_half_even_odds/total_neurons)) '%%) \n'])
fprintf(['\t\t\t PW:' num2str(PW_skagg_95_Kendall_95_half_even_odds) ...
    ' (' num2str(round(100*PW_skagg_95_Kendall_95_half_even_odds/PW_total_neurons),2) '%%)'...
    '\t\t TO:' num2str(TO_skagg_95_Kendall_95_half_even_odds)...
    ' (' num2str(round(100*TO_skagg_95_Kendall_95_half_even_odds/TO_total_neurons),2) '%%)\n'])





