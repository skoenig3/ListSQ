%get recording locations for view cells
clar
monkeys = {'Vivian','Tobii'};

task = 'ListSQ';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)


%---Firing Rate Curves for List---%
all_firing_rates = []; %all fixations
all_in_rates = [];%fixations out-> in
all_out_rates = []; %fixaiton out-> out
all_in_out_rates = [];%fixaions out-> in
all_peaks = [];%peak firing time of place cells for out-> in fixations

%---Firing Rates Between Sequence and List Tasks---%
all_context_gain = []; %change in firing rate (list_fr-seq_fr)/list_fr
all_context_gain2 = [];%gain list_fr/seq_fr
sig_p_list = []; %whether firing rates were significantly different for in vs out for list
sig_p_seq = []; %whether firing rates were significantly different for in vs out for seq
all_peak_list = []; %peak firing rate during list
all_peak_seq = []; %peak firing rate during seq
all_seq_peak_time = [];%peak firing time seq
all_list_peak_time = [];%peak firing time list

other_sig_p = []; %for units that had corr or skagg sig were they also significant, type 1: skagg, 2: corr 1/2, 3: corr e/o, 4: either corr+skgg
%---place field properties---%
all_areas = []; %area of place field
place_coverage = {}; %place field locations
coverage = {}; %eye data coverage for place cells

%---Other Parameters---%
monkey_count = zeros(2,2);%row 1 place row 2 non-place
all_unit_names = {}; %place cell unit names
all_monkeys = []; %1s and 2s
AP_location = []; %AP location of recorded place cell


Fs = 1000;
imageX = 800;
imageY = 600;

for monk = 2:-1:1
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
        
        num_trials = length(cfg.trl);
        valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
        if all(isnan(valid_trials(1,:)))
            continue
        end
        
        disp(task_file(1:8))
        load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],...
            'spatial_info','spike_times','eyepos','binsize','filter_width')
        
        load([data_dir task_file(1:8) '-Place_Cell_Analysis.mat'],...
            'twin1','twin2','smval','list_fixation_locked_firing','in_out',...
            'p_in_out','sliding_step','sliding_window','p_thresh',...
            'task_file','area','contextual_gain','unit_stats','sequence_fixation_locked_firing',...
            'all_place_field_matrix','seq_p_in_out')
        
        smval = 40;%smval/2;
        for unit = 1:num_units
            if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                    && (spatial_info.spatialstability_even_odd_prctile(unit) > 95) % ... %spatial consistency
                %                     && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                
                if contextual_gain(2,unit) < 1 %firing rate aligned in time is less than 1 Hze
                    continue
                end
                %---Other Parameters---%
                monkey_count(1,monk) = monkey_count(1,monk)+1;
                all_unit_names = [all_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}];
                all_monkeys = [all_monkeys monk]; %1s and 2s
                AP_location = [AP_location chamber_zero(1)+ session_data{sess}.location(1)]; %AP location of recorded place cell
                
                
                %---place field properties---%
                H = define_spatial_filter(filter_width);
                trial_data{1} = eyepos{unit};
                trial_data{2} = spike_times{unit};
                [~,timemaps] = get_firing_rate_map(trial_data,imageX,imageY,binsize,H,Fs,'all');
                coverage = [coverage {timemaps}];  %eye data coverage for place cells
                
                place_coverage = [place_coverage all_place_field_matrix(unit)]; %place field locations
                
                %---Firing Rates Between Sequence and List Tasks---%
                all_context_gain = [all_context_gain (contextual_gain(2,unit)-contextual_gain(4,unit))/contextual_gain(2,unit)];%change in firing rate (list_fr-seq_fr)/list_fr
                all_context_gain2 = [all_context_gain2 contextual_gain(2,unit)/contextual_gain(4,unit)];%gain list_fr/seq_fr
                all_peak_list = [all_peak_list contextual_gain(2,unit)]; %peak firing rate during list
                all_peak_seq = [all_peak_seq contextual_gain(4,unit)]; %peak firing rate during seq
                all_seq_peak_time = [all_seq_peak_time contextual_gain(3,unit)];%peak firing time seq
                all_list_peak_time = [all_list_peak_time contextual_gain(1,unit)];%peak firing time list
                
                
                if p_in_out(unit,1) < 0.05
                    sig_p_list = [sig_p_list 1]; %whether firing rates were significantly different for in vs out for list
                elseif  p_in_out(unit,2) < 0.05
                     sig_p_list = [sig_p_list 2]; %whether firing rates were significantly different for in vs out for list
                else
                    sig_p_list = [sig_p_list 0]; %whether firing rates were significantly different for in vs out for list
                end
                
                if ~isnan(contextual_gain(1,unit)) %if had at least 1 item in the field
                    if seq_p_in_out(unit) < 0.05
                        sig_p_seq = [sig_p_seq 1]; %whether firing rates were significantly different for in vs out for seq
                    else
                        sig_p_seq = [sig_p_seq 0]; %whether firing rates were significantly different for in vs out for seq
                    end
                else
                    sig_p_seq = [sig_p_seq NaN]; %whether firing rates were significantly different for in vs out for seq
                end
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %---Firing Rate Curves for List---%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %1) first fixation in: out-> in
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 1,:);
                [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                if isnan(all_peak_list(end))
                   all_peak_list(end) = max(firing_rate(twin1:end)); 
                end
                firing_rate = firing_rate-nanmean(firing_rate(1:twin1));%normalize to firing rate before fixation inside field
                firing_rate = firing_rate/max(firing_rate);
                all_in_rates = [all_in_rates; firing_rate];
                
                %find peaks
                [PKS,LOCS]= findpeaks(firing_rate,'MinPeakWidth',50);
                
                %remove peaks less than 1/2 the max
                LOCS(PKS < 0.5) = [];
                PKS(PKS < 0.5) = [];
                
                if isempty(LOCS)
                    LOCS = NaN;
                end
                
                
                %take the first peak
                all_peaks = [all_peaks LOCS(1)];
                
                %all fixations
                firing_rate = list_fixation_locked_firing{unit};
                firing_rate(sum(firing_rate,2) == 0,:) = [];
                [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                firing_rate = firing_rate-nanmean(firing_rate);%normalize to firing rate before fixation inside field
                firing_rate = firing_rate/max(firing_rate);
                all_firing_rates = [all_firing_rates; firing_rate];
                
                all_areas = [all_areas area(unit)];
                
                
                %3) fixations in -> out
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 3,:);
                [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                firing_rate = firing_rate-nanmean(firing_rate(1:twin1));
                firing_rate = firing_rate/max(firing_rate);
                all_in_out_rates = [all_in_out_rates; firing_rate];
                
                
                %4) fixation out of field but not first: out-> out
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 4,:);
                [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                firing_rate = firing_rate-nanmean(firing_rate);
                firing_rate = firing_rate/max(firing_rate);
                all_out_rates = [all_out_rates; firing_rate];
            elseif ~isnan(spatial_info.shuffled_rate_prctile(unit))
                monkey_count(2,monk) = monkey_count(2,monk)+1;
                
                if p_in_out(unit,1) < 0.05
                    p = 1;
                elseif  p_in_out(unit,2) < 0.05
                    p = 2;
                else
                    p = 0;
                end
                
                %for units that had corr or skagg sig were they also significant, type 1: skagg, 2: corr 1/2, 3: corr e/o, 4: either corr+skgg, 5: both corr
                if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                        && ~(spatial_info.spatialstability_even_odd_prctile(unit) > 95) ... %spatial consistency
                        && ~(spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                    other_sig_p = [other_sig_p [p;1]];
                elseif (spatial_info.spatialstability_even_odd_prctile(unit) > 95) ... %spatial consistency
                        && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                    other_sig_p = [other_sig_p [p;5]];
                elseif (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                        && ((spatial_info.spatialstability_even_odd_prctile(unit) > 95) ... %spatial consistency
                        || (spatial_info.spatialstability_halves_prctile(unit) > 95)) %spatial stability
                    other_sig_p = [other_sig_p [p;4]];
                elseif (spatial_info.spatialstability_even_odd_prctile(unit) > 95)
                    other_sig_p = [other_sig_p [p;2]];
                elseif (spatial_info.spatialstability_halves_prctile(unit) > 95)
                    other_sig_p = [other_sig_p [p;3]];
                end
            end
        end
    end
end
%% Plot Population Firing Rate Curves
figure
subplot(2,2,1)
[m,i] = max(all_in_rates(:,twin1:end),[],2);
[mm,ii] = sort(all_peaks);
[mm,ii] = sort(i);
imagesc([-twin1:twin2-1],[1:size(all_in_rates,1)],all_in_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_in_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('out -> in')

subplot(2,2,2)
% [m,i] = max(all_out_rates,[],2);
% [mm,ii] = sort(i);
imagesc([-twin1:twin2-1],[1:size(all_out_rates,1)],all_out_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_out_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('out -> out')

subplot(2,2,3)
% [m,i] = max(all_in_out_rates,[],2);
% [mm,ii] = sort(i);
imagesc([-twin1:twin2-1],[1:size(all_in_out_rates,1)],all_in_out_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_in_out_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('in->out')

subplot(2,2,4)
% [m,i] = max(all_firing_rates,[],2);
% [mm,ii] = sort(i);
imagesc([-twin1:twin2-1],[1:size(all_firing_rates,1)],all_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_firing_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('All')

%% Distribution of Peak Times
figure
hist(all_peaks-twin1,22)
xlabel('Peak Delay from Fixation Start (ms)')
ylabel('Neruon Count')
hold on
plot([mean(all_peaks-twin1) nanmean(all_peaks-twin1)],[0 10],'k--')
title(['Mean Peak Delay From Fixation Start = ' num2str(round(nanmean(all_peaks-twin1))) ' ms'])
xlim([-twin1 twin2])
axis square

%add FWHM
%% Eye Coverage and Place Cell Coverage
all_coverage = zeros(size(coverage{1}));
coverage_count =  zeros(size(coverage{1}));
for c = 1:length(coverage)
    cv = coverage{c};
    coverage_count(~isnan(cv)) =  coverage_count(~isnan(cv))+1;
    cv(isnan(cv)) = 0;
    all_coverage = all_coverage+cv;
end
all_coverage = all_coverage./coverage_count;

all_place_coverage  = zeros(imageY,imageX);
% all_place_coverage2  = zeros(imageY,imageX);
coverage_count =  zeros(imageY,imageX);
for c = 1:length(place_coverage);
    cv = place_coverage{c};
    cov = coverage{c};
    cov = imresize(cov,[imageY,imageX],'method','nearest');%upsample to images size in pixel coordinates
     coverage_count(~isnan(cov)) =  coverage_count(~isnan(cov))+1;
%      all_place_coverage2 = all_place_coverage2+cv;
      cv(isnan(cov)) = 0;
    all_place_coverage = all_place_coverage+cv;
end
all_place_coverage  = all_place_coverage./coverage_count;

figure
subplot(1,2,1)
imagesc(all_coverage)
axis off
colormap('jet')
title('Place Cell Eye Data Coverage')

subplot(1,2,2)
imagesc(all_place_coverage)
axis off
colormap('jet')
title('Place Cell Place Field Coverage')
p99 = prctile(all_place_coverage(:),99);
clims = caxis;
caxis([clims(1) p99]);

%%
total_in_out = sum(~isnan(sig_p_seq));
firing_rate_not_sig = [nanmean(all_peak_seq(sig_p_seq == 0)) nanmedian(all_peak_seq(sig_p_seq == 0)) nanmin((sig_p_seq == 0)) nanmax(all_peak_seq(sig_p_seq == 0))];
    firing_rate_sig = [nanmean(all_peak_seq(sig_p_seq == 1)) nanmedian(all_peak_seq(sig_p_seq == 1)) nanmin((sig_p_seq == 1)) nanmax(all_peak_seq(sig_p_seq == 1))];
    
%%
x = all_list_peak_time(sig_p_seq == 1)-twin1;
y = all_seq_peak_time(sig_p_seq == 1)-twin2;
x(isnan(y)) = [];
y(isnan(y)) = [];
[r,p] = corrcoef(x,y);

figure
plot(x,y,'.k')
xlabel('List Delay to Peak from Fixation Start (ms)')
ylabel('Sequence Delay to Peak from Fixation Start (ms)')
title('Significant Response in Sequence and List')

%%
air = all_in_rates;
ap = all_peaks;
sgpl = sig_p_list;
air(isnan(ap),:) = [];
sgpl(isnan(ap)) = [];
ap(isnan(ap)) = [];
ap(sgpl == 0) = [];
air(sgpl == 0,:) = [];
figure
[m,i] = max(air(:,twin1:end),[],2);
[mm,ii] = sort(ap);
[mm,ii] = sort(i);
imagesc([-twin1:twin2-1],[1:size(air,1)],air(ii,:))
hold on,plot(mm,1:length(mm),'w*'),hold off
colormap('jet')
hold on
plot([0 0],[1 size(air,1)],'w--');
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('Neuron #')
% pct5 = prctile(air(:),5);
% caxis([pct5 1])
vals = air(:,1:twin1);
caxis([-std(vals(:)) 1])
xlim([-twin1 twin2])
%%

[m,i] = max(nanmean(air));
figure
plot([-twin1:twin2-1],nanmean(air))
hold on
plot(i-twin1,m,'r*')
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--');
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('Normalized Firing Rate')
title(['Average Place Cell Fixation Aligned Firing Rate, peak @ ' num2str(i-twin1)])
xlim([-twin1 twin2])

%%
FWHM = NaN(1,length(air));
for n = 1:size(air,1)
   cv = air(n,:); 
   peak = ap(n);
   tims = find(cv < 0.5*cv(peak));
   diffs = tims-peak;
   neg_diffs = diffs(diffs < 0);
   pos_diffs = diffs(diffs > 0);
   neg = min(abs(neg_diffs));
   pos = min(pos_diffs);
   
   if isempty(pos)
       pos = twin2;
        FWHM(n) = NaN;
   else
        FWHM(n) = neg+pos;
        
%         if FWHM(n) > 200
%                      figure
%                plot(cv)
%                hold on
%                plot(peak,cv(peak),'r*')
%                plot(peak-neg,cv(peak-neg),'k*')
%                plot(peak+pos,cv(peak+pos),'k*')
%                plot([peak-neg peak+pos],[cv(peak-neg) cv(peak+pos)],'k')
%                hold off
%                title(['FWHM = ' num2str(neg+pos)])
%         end
   end
   
end
figure
hist(FWHM,25)
xlabel('FWHM @ peak')
ylabel('Count')
title(['Median FWHM = ' num2str(nanmedian(FWHM),3)])
%%
not_sig = find(sig_p_list > 0);
peak_fr_not_sig = all_peak_list(not_sig);
peak_fr_seq = all_peak_seq(not_sig);