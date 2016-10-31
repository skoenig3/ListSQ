%get recording locations for view cells
clar
monkeys = {'Vivian','Tobii'};

task = 'ListSQ';
Fs = 1000;




spatial_type = [];%1 >95 skagg, 2 > 99 skagg, 3 > 95 spatial corr, 4 >95 spatial corr and skagg, 5 normal
corr_half = [];
corr_even_odds = [];
for monk = 2
    monkey = monkeys{monk};
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';
        
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working on importing the data
    end
    
    for sess = 1:length(session_data)
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{sess},task);
        if isempty(task_file)
            continue
        end
        
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials','stability_attribute','cfg','hdr','data');
        
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        
        if num_units == 0
            continue
        end
        
        load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'])
        disp(task_file(1:8))
        
        filter_size = filter_width*10;
        H = fspecial('gaussian',filter_size,filter_width);
        
        for unit = 1:size(valid_trials,2)
            
            %unsable units
            if multiunit(unit) ==1
                continue
            end
            
            if isnan(spatial_info.rate(unit)) || spatial_info.rate(unit) == 0 %didn't run on unit since not stable
                continue
            end
            
            
            ntrials = floor(size(spike_times{unit},1)/2);
            %---get rate map for first half vs second half---%
            %first half
            filtered_time1 = filter_time(eyepos{unit}(1:ntrials*2,:),imageX,imageY,Fs,binsize,H);
            filtered_time1(filtered_time1 <= min_bin_dur) = NaN;
            filtered_space1 = filter_space(eyepos{unit}(1:ntrials*2,:),...
                spike_times{unit}(1:ntrials,:),imageX,imageY,binsize,H);
            ratemap1 = filtered_space1./filtered_time1;
            %             ratemap1(isnan(ratemap1)) = 0;
            
            %second half
            filtered_time2 = filter_time(eyepos{unit}(ntrials*2+1:ntrials*2*2,:),...
                imageX,imageY,Fs,binsize,H);
            filtered_time2(filtered_time2 <= min_bin_dur) = NaN;
            filtered_space2 = filter_space(eyepos{unit}(ntrials*2+1:ntrials*2*2,:),...
                spike_times{unit}(ntrials+1:ntrials*2,:),imageX,imageY,binsize,H);
            ratemap2 = filtered_space2./filtered_time2;
            %             ratemap2(isnan(ratemap2)) = 0;
            
            %---Get rate map for even and odd trials---%
            x = eyepos{unit}(1:2:end,:);
            y = eyepos{unit}(2:2:end,:);
            
            x_odd = x(1:2:end,:);
            x_even = x(2:2:end,:);
            y_odd = y(1:2:end,:);
            y_even = y(2:2:end,:);
            
            eye_even = NaN(2*size(x_even,1),size(x_even,2));
            eye_even(1:2:end,:) = x_even;
            eye_even(2:2:end,:) = y_even;
            
            eye_odd = NaN(2*size(x_odd,1),size(x_odd,2));
            eye_odd(1:2:end,:) = x_odd;
            eye_odd(2:2:end,:) = y_odd;
            
            
            filtered_time_even = filter_time(eye_even,imageX,imageY,Fs,binsize,H);
            filtered_time_even(filtered_time_even <= min_bin_dur) = NaN;
            filtered_space_even = filter_space(eye_even,...
                spike_times{unit}(2:2:end,:),imageX,imageY,binsize,H);
            ratemap_even = filtered_space_even./filtered_time_even;
            
            filtered_time_odd = filter_time(eye_odd,imageX,imageY,Fs,binsize,H);
            filtered_time_odd(filtered_time_odd <= min_bin_dur) = NaN;
            filtered_space_odd = filter_space(eye_odd,...
                spike_times{unit}(1:2:end,:),imageX,imageY,binsize,H);
            ratemap_odd = filtered_space_odd./filtered_time_odd;
            
            max_fr = max([ratemap1(:)' ratemap2(:)' ratemap_even(:)' ratemap_odd(:)']);
            min_fr = min([ratemap1(:)' ratemap2(:)' ratemap_even(:)' ratemap_odd(:)']);
            
            r = corr(ratemap1(:),ratemap2(:),'row','pairwise','type','Pearson');
            rho = corr(ratemap1(:),ratemap2(:),'row','pairwise','type','Spearman');
            tau = corr(ratemap1(:),ratemap2(:),'row','pairwise','type','Kendall');
            
            corr_half = [corr_half [r rho tau]'];
            
            r = corr(ratemap_even(:),ratemap_odd(:),'row','pairwise','type','Pearson');
            rho = corr(ratemap_even(:),ratemap_odd(:),'row','pairwise','type','Spearman');
            tau = corr(ratemap_even(:),ratemap_odd(:),'row','pairwise','type','Kendall');
            
            corr_even_odds = [corr_even_odds [r rho tau]'];
            
            if spatial_info.shuffled_rate_prctile(unit) > 95
                if spatial_info.shuffled_spatialstability_prctile(unit) > 95
                    spatial_type = [spatial_type 1];
                elseif spatial_info.shuffled_rate_prctile(unit) > 99
                    spatial_type = [spatial_type 2];
                else
                    spatial_type = [spatial_type 3];
                end
            elseif spatial_info.shuffled_spatialstability_prctile(unit) > 95
                spatial_type = [spatial_type 4];
            else
                spatial_type = [spatial_type 5];%normal
            end  
            
            %             figure
            %             subplot(2,2,1)
            %             h = imagesc(ratemap1);
            %             set(h,'alphadata',~isnan(ratemap1));
            %             axis off
            %             axis equal
            %             %     caxis([min_fr max_fr])
            %             colormap('jet')
            %             title('First 1/2')
            %
            %             subplot(2,2,3)
            %             h = imagesc(ratemap2);
            %             set(h,'alphadata',~isnan(ratemap2));
            %             axis off
            %             axis equal
            %             %     caxis([min_fr max_fr])
            %             colormap('jet')
            %             title('Second Half')
            %
            %             subplot(2,2,2)
            %             h = imagesc(ratemap_odd);
            %             set(h,'alphadata',~isnan(ratemap_odd));
            %             axis off
            %             axis equal
            %             %     caxis([min_fr max_fr])
            %             colormap('jet')
            %             title('Odd trials')
            %
            %             subplot(2,2,4)
            %             h = imagesc(ratemap_even);
            %             set(h,'alphadata',~isnan(ratemap_even));
            %             axis off
            %             axis equal
            %             %     caxis([min_fr max_fr])
            %             colormap('jet')
            %             title('Even Trials')
            
        end
    end
end
%%
figure
plot(corr_even_odds,corr_half,'.')
xlabel('Correlation Even Odds'),ylabel('Correlation 1st 2nd halfs')
legend('r','rho','tau')
%% Even Odd trials
figure
subplot(2,2,1)
hold on
for type = 1:5
   plot( type ,corr_even_odds(1,spatial_type == type),'.k')
   plot(type,nanmean(corr_even_odds(1,spatial_type == type)),'ko','markersize',10)
end
hold off
set(gca,'Xtick',[0:6]);
set(gca,'XtickLabel',{'','Corr/Skagg 95','Skagg 99','Skagg 95','Corr 95','~sig'})
xlim([0 6])
ylabel('Pearsons r')

subplot(2,2,2)
hold on
for type = 1:5
   plot( type ,corr_even_odds(2,spatial_type == type),'.k')
   plot(type,nanmean(corr_even_odds(2,spatial_type == type)),'ko','markersize',10)
end
hold off
set(gca,'Xtick',[0:6]);
set(gca,'XtickLabel',{'','Corr/Skagg 95','Skagg 99','Skagg 95','Corr 95','~sig'})
xlim([0 6])
ylabel('Spearmans rho')

subplot(2,2,3)
hold on
for type = 1:5
   plot( type ,corr_even_odds(3,spatial_type == type),'.k')
   plot(type,nanmean(corr_even_odds(3,spatial_type == type)),'ko','markersize',10)
end
hold off
set(gca,'Xtick',[0:6]);
set(gca,'XtickLabel',{'','Corr/Skagg 95','Skagg 99','Skagg 95','Corr 95','~sig'})
xlim([0 6])
ylabel('Kendalls tau')

subtitle('Even Odd Trials')

%% First Half Second Half
figure
subplot(2,2,1)
hold on
for type = 1:5
   plot( type ,corr_half(1,spatial_type == type),'.k')
   plot(type,nanmean(corr_half(1,spatial_type == type)),'ko','markersize',10)
end
hold off
set(gca,'Xtick',[0:6]);
set(gca,'XtickLabel',{'','Corr/Skagg 95','Skagg 99','Skagg 95','Corr 95','~sig'})
xlim([0 6])
ylabel('Pearsons r')

subplot(2,2,2)
hold on
for type = 1:5
   plot( type ,corr_half(2,spatial_type == type),'.k')
   plot(type,nanmean(corr_half(2,spatial_type == type)),'ko','markersize',10)
end
hold off
set(gca,'Xtick',[0:6]);
set(gca,'XtickLabel',{'','Corr/Skagg 95','Skagg 99','Skagg 95','Corr 95','~sig'})
xlim([0 6])
ylabel('Spearmans rho')

subplot(2,2,3)
hold on
for type = 1:5
   plot( type ,corr_half(3,spatial_type == type),'.k')
   plot(type,nanmean(corr_half(3,spatial_type == type)),'ko','markersize',10)
end
hold off
set(gca,'Xtick',[0:6]);
set(gca,'XtickLabel',{'','Corr/Skagg 95','Skagg 99','Skagg 95','Corr 95','~sig'})
xlim([0 6])
ylabel('Kendalls tau')

subtitle('1st Half vs 2nd Half')
%%
figure
subplot(2,2,1)
plot(corr_even_odds(1,:),corr_even_odds(2,:),'.')
hold on
plot([-1 1],[-1 1],'k--')
hold off
xlabel('Pearsons r')
ylabel('Spearmans rho')
axis square

subplot(2,2,2)
plot(corr_even_odds(1,:),corr_even_odds(3,:),'.')
hold on
plot([-1 1],[-1 1],'k--')
hold off
xlabel('Pearsons r')
ylabel('Kendalls tau')
axis square

subplot(2,2,3)
plot(corr_even_odds(2,:),corr_even_odds(3,:),'.')
hold on
plot([-1 1],[-1 1],'k--')
hold off
xlabel('Spearmans rho')
ylabel('Kendalls tau')
axis square

subtitle('Even vs Odd Trials')


