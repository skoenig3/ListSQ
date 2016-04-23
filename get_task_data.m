function [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,waveform_count] = get_task_data(session_data,task)
%function grabs desired task data

if ~isempty(strfind(lower(session_data.task1),lower(task)))
    
    task_file = session_data.task1_file;
    
    %item and condition file are usually the same but not always in case of
    %mistake or when running a different strcture
    itm = strfind(session_data.task1,'itm');
    cnd = strfind(session_data.task1,'cnd');
    item_file = session_data.task1(itm-9:itm+2);
    cnd_file = session_data.task1(cnd-9:cnd+2);
    
    waveform_count = session_data.task1_unit_waveforCount;
    multiunit = session_data.task1_unit_MultiunitSeperability;
    unit_confidence = session_data.task1_unit_confidence;
    sorting_quality = session_data.task1_unit_cutQuality;
    unit_names =  session_data.task1_unit_names;
    
elseif isfield(session_data,'task2')
    if ~isempty(strfind(lower(session_data.task2),lower(task)))
        
        task_file = session_data.task2_file;
        
        %item and condition file are usually the same but not always in case of
        %mistake or when running a different strcture
        itm = strfind(session_data.task2,'itm');
        cnd = strfind(session_data.task2,'cnd');
        item_file = session_data.task2(itm-9:itm+2);
        cnd_file = session_data.task2(cnd-9:cnd+2);
        
        waveform_count = session_data.task2_unit_waveforCount;
        multiunit = session_data.task2_unit_MultiunitSeperability;
        unit_confidence = session_data.task2_unit_confidence;
        sorting_quality = session_data.task2_unit_cutQuality;
        unit_names =  session_data.task2_unit_names;
    end
    
else
    warning('Could not identify the proper tasks, exiting function')
    task_file = [];
    item_file = [];
    cnd_file = [];
    multiunit = [];
    unit_names = [];
    unit_confidence = [];
    sorting_quality = [];
    waveform_count = [];
end

