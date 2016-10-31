% Population Visual Response and List/Sequence Eye Movements
% written Seth Konig 9/5/16

clar
task = 'ListSQ';
Fs = 1000;%Hz
smval = 60;
imageX = 800;
imageY = 600;


cross_on_image_on = cell(1,2); %aka visual response
cross_fixation_sequence_fixation = cell(1,2); %fixation response long period 500+ms so equavalent
cross_fixation_list_fixation = cell(1,2); %fixation response short period
image_on_sequence_fixation = cell(1,2); %longish time period
image_on_image_off = cell(1,2); %sort of a visual repsonse


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
            load([data_dir task_file(1:8) '-Fixation_locked_Sequence_results.mat']);
            seq_fix = fixation_locked_firing; %since same variable
            seq_info = fixation_info;
            seq_info2 = fixation_info2;
            load([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat']);
        else
            continue
        end
        
        disp(['Importing ' task_file(1:8)])
        
        for unit = 1:num_units
            if ~isempty(time_lock_firing{unit,1})
                
                
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %---For cross on and image on---%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (epoch_data.rate_prctile(unit,1) > 95) && (epoch_data.temporalstability_prctile(unit,1) > 95)...
                        && (epoch_data.rate_prctile(unit,3) > 95) && (epoch_data.temporalstability_prctile(unit,3) > 95)
                    
                    %---for cross on---%
                    firing_rate = time_lock_firing{unit,1};
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    cross_on_image_on{1} = [cross_on_image_on{1} ; firing_rate];
                    
                    %for image on
                    firing_rate = time_lock_firing{unit,3};
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    cross_on_image_on{2} = [cross_on_image_on{2} ; firing_rate];
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %---For Fixation on cross and Sequence Fixation---%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (epoch_data.rate_prctile(unit,2) > 95) && (epoch_data.temporalstability_prctile(unit,2) > 95) ...
                        && ( seq_info.rate_prctile(unit) > 95 ||  seq_info2.rate_prctile(unit) > 95)
                    
                    %--for fixation on cross hair--%
                    firing_rate = time_lock_firing{unit,2}(:,1:twin1+twin3);
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    cross_fixation_sequence_fixation{1} = [cross_fixation_sequence_fixation{1}; firing_rate];
                    
                    %--For Sequence Fixation--%
                    firing_rate =seq_fix{unit};
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    cross_fixation_sequence_fixation{2} = [cross_fixation_sequence_fixation{2}; firing_rate];
                    
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %---For Fixation on cross and List Fixation---%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ((epoch_data.rate_prctile(unit,2) > 95) && (epoch_data.temporalstability_prctile(unit,2) > 95)) ...
                        && (((fixation_info.rate(unit) >  fixation_info_95.rate(unit)) &&  (saccade_info.rate(unit) >  saccade_info_95.rate(unit))) || ...
                        ((fixation_info2.rate(unit) >  fixation_info2_95.rate(unit)) &&  (saccade_info2.rate(unit) >  saccade_info2_95.rate(unit))))
                    
                    %--for fixation on cross hair--%
                    firing_rate = time_lock_firing{unit,2}(:,1:twin1+twin3);
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    cross_fixation_list_fixation{1} = [cross_fixation_list_fixation{1}; firing_rate];
                    
                    %--For Sequence Fixation--%
                    firing_rate =fixation_locked_firing{unit}(fixation_information{unit}(:,4) > 2*twin,:);
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    cross_fixation_list_fixation{2} = [cross_fixation_list_fixation{2}; firing_rate];
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %---For Image On and Sequence Fixation---%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (epoch_data.rate_prctile(unit,3) > 95) && (epoch_data.temporalstability_prctile(unit,3) > 95) ...
                        && ( seq_info.rate_prctile(unit) > 95 ||  seq_info2.rate_prctile(unit) > 95)
                    
                    %--for fixation on cross hair--%
                    firing_rate = time_lock_firing{unit,2}(:,1:twin1+twin3);
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    image_on_sequence_fixation{1} = [image_on_sequence_fixation{1}; firing_rate];
                    
                    %--For Sequence Fixation--%
                    firing_rate =seq_fix{unit};
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    image_on_sequence_fixation{2} = [image_on_sequence_fixation{2}; firing_rate];
                end
                
                
                %---For Image On and Image Off---%
                if (epoch_data.rate_prctile(unit,3) > 95) && (epoch_data.temporalstability_prctile(unit,3) > 95)...
                        && (epoch_data.rate_prctile(unit,4) > 95) && (epoch_data.temporalstability_prctile(unit,4) > 95)
                    
                    %---for image on---%
                    firing_rate = time_lock_firing{unit,3};
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    image_on_image_off{1} = [image_on_image_off{1} ; firing_rate];
                    
                    %---for image off---%
                    firing_rate = time_lock_firing{unit,4};
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    image_on_image_off{2} = [image_on_image_off{2} ; firing_rate];
                end
            end
        end
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Plots Cross on and Image on---%

figure

subplot(2,2,1)
[m,cross_i] = max(cross_on_image_on{1},[],2);
[mm,cross_ii] = sort(cross_i);
imagesc([-twin1:twin3-1],[1:size(cross_on_image_on{1},1)],cross_on_image_on{1}(cross_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(cross_on_image_on{1},1)],'w--');
hold off
xlabel('Time from Crosshair Onset (ms)')
ylabel('Neuron #')
title('Aligned to Crosshair Onset')

subplot(2,2,2)
[m,image_i] = max(cross_on_image_on{2},[],2);
[mm,image_ii] = sort(image_i);
imagesc([-twin1:twin2-1],[1:size(cross_on_image_on{2},1)],cross_on_image_on{2}(image_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(cross_on_image_on{1},1)],'w--');
hold off
xlabel('Time from Image Onset (ms)')
ylabel('Neuron #')
title('Aligned to Image Onset')

subplot(2,2,3)
imagesc([-twin1:twin3-1],[1:size(cross_on_image_on{1},1)],cross_on_image_on{1}(image_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(cross_on_image_on{1},1)],'w--');
hold off
xlabel('Time from Crosshair Onset (ms)')
ylabel('Neuron #')
title('Sorted by Image Onset')

subplot(2,2,4)
imagesc([-twin1:twin2-1],[1:size(cross_on_image_on{2},1)],cross_on_image_on{2}(cross_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(cross_on_image_on{1},1)],'w--');
hold off
xlabel('Time from Image Onset (ms)')
ylabel('Neuron #')
title('Sorted by Crosshair Onset')

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Plots Cross Fixation and Sequence---%

figure

subplot(2,2,1)
[m,crossf_i] = max(cross_fixation_sequence_fixation{1},[],2);
[mm,crossf_ii] = sort(crossf_i);
imagesc([-twin1:twin3-1],[1:size(cross_fixation_sequence_fixation{1},1)],cross_fixation_sequence_fixation{1}(crossf_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(cross_fixation_sequence_fixation{1},1)],'w--');
hold off
xlabel('Time from Fixation (ms)')
ylabel('Neuron #')
title('Aligned to Fixation on Crosshair')

subplot(2,2,2)
[m,seq_i] = max(cross_fixation_sequence_fixation{2},[],2);
[mm,seq_ii] = sort(seq_i);
imagesc([-twin3:twin3-1],[1:size(cross_fixation_sequence_fixation{2},1)],cross_fixation_sequence_fixation{2}(seq_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(cross_fixation_sequence_fixation{2},1)],'w--');
hold off
xlabel('Time from Fixation(ms)')
ylabel('Neuron #')
title('Aligned to Fixation in Sequence')

subplot(2,2,3)
imagesc([-twin1:twin3-1],[1:size(cross_fixation_sequence_fixation{1},1)],cross_fixation_sequence_fixation{1}(seq_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(cross_fixation_sequence_fixation{1},1)],'w--');
hold off
xlabel('Time from Fixation (ms)')
ylabel('Neuron #')
title('Aligned to Sequence')

subplot(2,2,4)
imagesc([-twin3:twin3-1],[1:size(cross_fixation_sequence_fixation{2},1)],cross_fixation_sequence_fixation{2}(crossf_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(cross_fixation_sequence_fixation{2},1)],'w--');
hold off
xlabel('Time from Fixation(ms)')
ylabel('Neuron #')
title('Aligned to Crosshiar')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Plots Cross Fixation and List---%

figure

subplot(2,2,1)
[m,crossf2_i] = max(cross_fixation_list_fixation{1},[],2);
[mm,crossf2_ii] = sort(crossf2_i);
imagesc([-twin1:twin3-1],[1:size(cross_fixation_list_fixation{1},1)],cross_fixation_list_fixation{1}(crossf2_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(cross_fixation_list_fixation{1},1)],'w--');
hold off
xlabel('Time from Fixation (ms)')
ylabel('Neuron #')
title('Aligned to Fixation on Crosshair')

subplot(2,2,2)
[m,list_i] = max(cross_fixation_list_fixation{2}(:,450:850),[],2);
[mm,list_ii] = sort(list_i);
imagesc([-twin3:twin3-1],[1:size(cross_fixation_list_fixation{2},1)],cross_fixation_list_fixation{2}(list_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(cross_fixation_sequence_fixation{2},1)],'w--');
hold off
xlabel('Time from Fixation(ms)')
ylabel('Neuron #')
title('Aligned to Fixation in List')


subplot(2,2,3)
imagesc([-twin1:twin3-1],[1:size(cross_fixation_list_fixation{1},1)],cross_fixation_list_fixation{1}(list_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(cross_fixation_list_fixation{1},1)],'w--');
hold off
xlabel('Time from Fixation (ms)')
ylabel('Neuron #')
title('Aligned to Fixation on List')

subplot(2,2,4)
[m,list_i] = max(cross_fixation_list_fixation{2}(:,450:850),[],2);
[mm,list_ii] = sort(list_i);
imagesc([-twin3:twin3-1],[1:size(cross_fixation_list_fixation{2},1)],cross_fixation_list_fixation{2}(crossf2_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(cross_fixation_sequence_fixation{2},1)],'w--');
hold off
xlabel('Time from Fixation(ms)')
ylabel('Neuron #')
title('Aligned to Fixation to Crosshair')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Plots Image on and Image Off---%

figure

subplot(2,2,1)
[m,imgon_i] = max(image_on_image_off{1}(:,200:end),[],2);
[mm,imgon_ii] = sort(imgon_i);
imagesc([-twin1:twin2-1],[1:size(image_on_image_off{1},1)],image_on_image_off{1}(imgon_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(image_on_image_off{1},1)],'w--');
hold off
xlabel('Time from Image On (ms)')
ylabel('Neuron #')
title('Aligned to Image On')

subplot(2,2,2)
[m,imgoff_i] = max(image_on_image_off{2}(:,500:end),[],2);
[mm,imgoff_ii] = sort(imgoff_i);
imagesc([-twin3:twin3-1],[1:size(image_on_image_off{2},1)],image_on_image_off{2}(imgoff_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(image_on_image_off{2},1)],'w--');
hold off
xlabel('Time from Image OFF (ms)')
ylabel('Neuron #')
title('Aligned to Image OFF')


subplot(2,2,3)
imagesc([-twin1:twin2-1],[1:size(image_on_image_off{1},1)],image_on_image_off{1}(imgoff_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(image_on_image_off{1},1)],'w--');
hold off
xlabel('Time from Image On (ms)')
ylabel('Neuron #')
title('Aligned to Image OFF')

subplot(2,2,4)
imagesc([-twin3:twin3-1],[1:size(image_on_image_off{2},1)],image_on_image_off{2}(imgon_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(image_on_image_off{2},1)],'w--');
hold off
xlabel('Time from Image OFF (ms)')
ylabel('Neuron #')
title('Aligned to Image On')

%%
image_on_sequence_fixation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Plots Image on Sequence Fixation---%

figure

subplot(2,2,1)
[m,imgon2_i] = max(image_on_sequence_fixation{1}(:,200:end),[],2);
[mm,imgon2_ii] = sort(imgon2_i);
imagesc([-twin1:twin2-1],[1:size(image_on_sequence_fixation{1},1)],image_on_sequence_fixation{1}(imgon2_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(image_on_sequence_fixation{1},1)],'w--');
hold off
xlabel('Time from Image On (ms)')
ylabel('Neuron #')
title('Aligned to Image On')

subplot(2,2,2)
[m,seq2_i] = max(image_on_sequence_fixation{2}(:,500:end),[],2);
[mm,seq2_ii] = sort(seq2_i);
imagesc([-twin3:twin3-1],[1:size(image_on_sequence_fixation{2},1)],image_on_sequence_fixation{2}(seq2_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(image_on_image_off{2},1)],'w--');
hold off
xlabel('Time from Fixation (ms)')
ylabel('Neuron #')
title('Aligned to Sequence Fixation')


subplot(2,2,3)
imagesc([-twin1:twin2-1],[1:size(image_on_sequence_fixation{1},1)],image_on_sequence_fixation{1}(seq2_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(image_on_sequence_fixation{1},1)],'w--');
hold off
xlabel('Time from Image On (ms)')
ylabel('Neuron #')
title('Aligned to Sequence')

subplot(2,2,4)
[m,seq2_i] = max(image_on_sequence_fixation{2}(:,500:end),[],2);
[mm,seq2_ii] = sort(seq2_i);
imagesc([-twin3:twin3-1],[1:size(image_on_sequence_fixation{2},1)],image_on_sequence_fixation{2}(imgon2_ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(image_on_image_off{2},1)],'w--');
hold off
xlabel('Time from Fixation (ms)')
ylabel('Neuron #')
title('Aligned to Image')