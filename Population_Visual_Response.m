% Population Visual Response
% written Seth Konig 9/1/16

clar
task = 'ListSQ';
Fs = 1000;%Hz
smval = 60;
imageX = 800;
imageY = 600;


cross_on_firing_rates = [];
cross_fixation_firing_rates = [];
image_on_firing_rates = [];
long_image_on_firing_rates = [];
image_off_firing_rates = [];
all_unit_names = {};
image_off_unit_names = {};
monkey_count = zeros(1,2);

figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Population Figures\Place Cell Eye Movements\';
for monkey = 1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        
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
        
        
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials','cfg','hdr','data','fixationstats');
        
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        
        if exist([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat'],'file') %want to remove later
            load([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat']);
        else
            continue
        end
        
        disp(['Importing ' task_file(1:8)])
        
        for unit = 1:num_units
            if ~isempty(time_lock_firing{unit,1})
                
                
                %for cross on
                if (epoch_data.rate_prctile(unit,1) > 95) && (epoch_data.temporalstability_prctile(unit,1) > 95)
                    firing_rate = time_lock_firing{unit,1};
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    cross_on_firing_rates = [cross_on_firing_rates ; firing_rate];
                end
                
                %for fixation on cross hair
                if (epoch_data.rate_prctile(unit,2) > 95) && (epoch_data.temporalstability_prctile(unit,2) > 95)
                    firing_rate = time_lock_firing{unit,2}(:,1:twin1+twin3);
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    cross_fixation_firing_rates = [cross_fixation_firing_rates; firing_rate];
                end
                
                
                %for image on
                if (epoch_data.rate_prctile(unit,3) > 95) && (epoch_data.temporalstability_prctile(unit,3) > 95)
                    firing_rate = time_lock_firing{unit,3};
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    image_on_firing_rates = [image_on_firing_rates ; firing_rate];
                end
                
                
                %for long image on
                if (epoch_data.rate_prctile(unit,5) > 95) && (epoch_data.temporalstability_prctile(unit,5) > 95)
                    firing_rate = time_lock_firing{unit,5};
                    [firing_rate,~]= nandens(firing_rate,smval2,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    long_image_on_firing_rates = [long_image_on_firing_rates ; firing_rate];
                    all_unit_names = [all_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}];
                    monkey_count(monkey) = monkey_count(monkey)+1;
                end
                
                %for image off
                if (epoch_data.rate_prctile(unit,4) > 95) && (epoch_data.temporalstability_prctile(unit,4) > 95)
                    firing_rate = time_lock_firing{unit,4};
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    image_off_firing_rates = [image_off_firing_rates ; firing_rate];
                    image_off_unit_names = [image_off_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}];
                end
            end
        end
    end
end
%% All Cells Aligned to Max
figure

subplot(2,2,1)
[m,i] = max(cross_on_firing_rates,[],2);
[mm,ii] = sort(i);
imagesc([-twin1:twin3-1],[1:size(cross_on_firing_rates,1)],cross_on_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(cross_on_firing_rates,1)],'w--');
hold off
xlabel('Time from Fixation  (ms)')
ylabel('Neuron #')
title('Crosshair Onset')



subplot(2,2,2)
[m,i] = max(cross_fixation_firing_rates,[],2);
[mm,ii] = sort(i);
imagesc([-twin1:twin3-1],[1:size(cross_fixation_firing_rates,1)],cross_fixation_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(cross_fixation_firing_rates,1)],'w--');
hold off
xlabel('Time from Fixation  (ms)')
ylabel('Neuron #')
title('Crosshair Fixation')


subplot(2,2,3)
[m,i] = max(image_on_firing_rates,[],2);
[mm,ii] = sort(i);
imagesc([-twin1:twin2-1],[1:size(image_on_firing_rates,1)],image_on_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(image_on_firing_rates,1)],'w--');
hold off
xlabel('Time from Image On (ms)')
ylabel('Neuron #')
title('Image Onset')

subplot(2,2,4)
[m,i] = max(image_off_firing_rates,[],2);
[mm,ii] = sort(i);
imagesc([-twin3:twin3-1],[1:size(image_off_firing_rates,1)],image_off_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(image_off_firing_rates,1)],'w--');
hold off
xlabel('Time from Image Off (ms)')
ylabel('Neuron #')
title('Image Offset')
%%

figure
[m,i] = max(long_image_on_firing_rates,[],2);
[mm,ii] = sort(i);
imagesc([-twin1:twin4-1],[1:size(long_image_on_firing_rates,1)],long_image_on_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(long_image_on_firing_rates,1)],'w--');
hold off
xlabel('Time from Image On (ms)')
ylabel('Neuron #')
title('Image Onset')


