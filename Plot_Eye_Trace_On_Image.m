%plot eye trace on image

clar
task = 'ListSQ';
twin = 500;

base_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\FlickrImages\LSQ';

for monkey = 2
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
    
    for session = 18%:length(session_data) %ignore 1st 6 sessions for simplicity
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import task and unit data---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
            = get_task_data(session_data{session},task);
        if isempty(task_file)
            warning('No file could be found for specificed task. Exiting function...')
            return
        end
        image_set_number = str2double(item_file(7:8));
        img_dir = [base_image_dir num2str(image_set_number) '\S' num2str(image_set_number) 'I'];
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg');

        
        %set/get some general important info
        [eyechans] = find_desired_channels(cfg,'eye');
        num_trials = length(cfg.trl);
        Fs = data(1).fsample; %should be 1000
        ITIstart_code = 15;
        img_on_code= 23;
        img_off_code = 24;
        imageX = 800; %horizontal size of the image
        imageY = 600; %horizontal size of the image
        
        %set the image duration
        if str2double(task_file(3:8)) < 140805
            imgdur = 7000;
        else
            imgdur = 5000;
        end
        
        %get important task specific information
        [itmlist,sequence_items] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        [which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);
        
        for t = 1:60
            if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                
                imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start+twin; %want to avoid visual response
                imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start;
                
                % if monkey isn't paying attention data isn't probably
                % worth much plus have to cut off somewhere
                if imgoff-imgon > 1.5*imgdur-1
                    imgoff = imgon+1.5*imgdur-1;
                end
                
                xn = round(data(eyechans(1)).values{t}(imgon:imgoff));
                yn = round(data(eyechans(2)).values{t}(imgon:imgoff));
                
                img_index = find(img_cnd == cfg.trl(t).cnd);
                
                
                figure
                if which_img(img_index) < 10
                    imshow(imread([img_dir  '0' num2str(which_img(img_index)) '.bmp']))
                else
                    imshow(imread([img_dir  num2str(which_img(img_index)) '.bmp']))
                end
                hold on
                plot(xn,imageY-yn,'y')
                plot([650 770],[570 570],'w')
                
                hold off
                
            end
        end
    end
end