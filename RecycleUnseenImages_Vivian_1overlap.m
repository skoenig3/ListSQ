% % Edited on 3/4/15 by SDK to generate item files for items ListSQ R-Z with 1
% overalpping location instead of 2 overlapping locaitons. DO NOT edit 
% file.! Previous version had 1st 3 location overlapping. 1 overlapping
% code was taken from Gen_Two_Fwd_Predict.m. 
% % Written December 29, 2014 by Seth Konig
% % code determines the last image Vivian saw during ListSQ sets and then
% % tracks which images she hasn't seen to reuse these in "new sets".
% %
% % Note code does not rename images and keeps original names in new set folders and in
% % the item files.
% %
% % Uses code from Make_ListSQ_Task_Files.m to generate sequneces with 3
% % items overlapping in location/order (i.e. same as section 3.2). Essentially, 
% % exact same code as Make_ListSQ_Task_Files.m but for 4 sets of novel and repeat
% % blocks instead of 6, so 64 images total instead of 96
% %
% % Run in order
% % [1] Determine which images Vivian saw for each set
% % [2] Sort images into Image Sets. 
% % [3] Generate item and condition file with 2 predictable sequences
% 
% %%
% %%---[1] Determine which images Vivian saw for each set---%%
% cortexfile = {...
%     'PW140724.2','PW140725.3','PW140728.3','PW140729.3','PW140730.3',...
%     'PW140801.3','PW140805.3','PW140806.3','PW140825.3','PW140826.3',...
%     'PW140828.3','PW140829.3','PW140908.3','PW140910.3','PW140915.3',...
%     'PW140917.3','PW141001.3','PW141002.3','PW141006.3','PW141007.3',...
%     'PW141008.3','PW141010.3','PW141013.3','PW141014.3','PW141015.3',...
%     'PW141016.3','PW141022.3','PW141023.3','PW141024.3','PW141027.3',...
%     'PW141028.3','PW141029.3','PW141031.2','PW141103.3','PW141104.3',...
%     'PW141105.3','PW141106.3','PW141107.3','PW141110.3','PW141112.3',...
%     'PW141113.3','PW141114.3','PW141117.3','PW141216.3','PW141217.3',...
%     'PW141218.3','PW141219.2','PW141222.2','PW141223.2','PW141226.2',...
%     'PW141229.3','PW141230.3','PW150107.3'};
% setnum =  1:53; %item file number associated with cortex file
% 
% lastimage = NaN(1,50);
% for set = 1:length(cortexfile)
%    if strcmpi(cortexfile{set}(1:2),'PW')
%       datafile = ['R:\Buffalo Lab\Cortex Data\Vivian\' cortexfile{set}]; 
%    end
%    
%    if setnum(set) < 10
%        itmfile = ['ListSQ0' num2str(setnum(set)) '.itm'];
%    else
%        itmfile = ['ListSQ' num2str(setnum(set)) '.itm'];
%    end
%    
%    [~,event_arr,~,~,~,~] = get_ALLdata(datafile); %import event data from cortex file
%    [itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(itmfile);
%    lastcnd = event_arr((find(event_arr(:,end)>1000,1,'last')),end)-1000;
%    if lastcnd  <= 21 %never got to any of the image trials
%       lastimage(set) = 0;
%       continue
%    end
% 
%    lastitm = itmlist(lastcnd);
%    while lastitm <= sequence_items(end) %means sequence trials
%       lastcnd = lastcnd-1;
%       lastitm =  itmlist(lastcnd);
%    end
%    %the image may be within a repeat block so have to check if there were
%    %any items (i.e. images) that came before it with a larger number
%    largestitem = max(itmlist(1:lastcnd));
%    lastimage(setnum(set)) = largestitem-sequence_items(end);
% end
% %%
% %%---[2]  Copy unseen pictures to a recycled directory---%%
% image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\FlickrImages\';
% recycle_dir = [image_dir 'Recycled_Images\'];
% mkdir(recycle_dir)
% 
% for set = 1:length(setnum)
%    if set < 10 %set names were named with 2 digits, a zero preceded # if set was less than 10
%        set_zero = '0'; 
%    else
%        set_zero = '';
%    end
%    set_dir = [image_dir 'LSQ' set_zero num2str(set) '\'];
%    for img = lastimage(set)+1:96 %won't do anything if seen all images in set
%         if img < 10 %same as with set names 2 digits for image #
%             img_zero = '0';
%         else
%             img_zero = '';
%         end
%         %copy file to recycle dir. DO NOT DELETE original FILE
%         copyfile([set_dir 'S' set_zero num2str(set) 'I' img_zero num2str(img) '.bmp'],...
%             [recycle_dir 'S' set_zero num2str(set) 'I' img_zero num2str(img) '.bmp']);
%    end
% end
% %%
% %%---[2] Sort images into Image Sets---%%
% % derived from section 4 of Make_ListSQ_Task_Files.m
% % secition sorts images from unused library into image sets and keeps track
% % of the original name of files. Original files are put into used libary.
root_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\FlickrImages\';
% unused_dir = 'Recycled_Images\'; %updated by SDK 1/7/15
% %unused_dir = 'Unused\'; %updated by SDK 1/7/15
% used_dir ='Used\';
% 
% cd([root_dir unused_dir]);
% 
% %updated by SDK 1/7/15
% setnums = ['A':'Z'  num2str([1:9]')']; %DO NOT edit line always start add more if needed
% %setnums = [1:51]; %do not write over original sets!!! Start at 1 will skip ones already created if in root_dir
% 
% num_images_per_set = 64; %was 96 updated by SDK 1/7/15
% for set = setnums;
%     d=dir([root_dir unused_dir '*.bmp']);
%     
%     %updated by SDK 1/7/15
%     set_dir = [root_dir 'LSQV' set '\']; %updated by SDK 1/7/15
%     %     if set < 10
%     %          set_dir = [root_dir 'LSQ0' num2str(set) '\']; %shortened from ListSQ or directory name too long for cortex
%     %     else
%     %         set_dir = [root_dir 'LSQ' num2str(set) '\']; %shortened from ListSQ or directory name too long for cortex
%     %     end
%     
%     if exist(set_dir,'dir')
%         disp('Image Set already exists') %do not let it rewrite original folders
%         continue;
%     else
%         mkdir(set_dir)
%     end
%     rr = randperm(length(d));
%     for img = 1:num_images_per_set;
%         original = d(rr(img)).name;
%         new = original;     %updated by SDK 1/7/15
%         
%         %updated by SDK 1/7/15
%         %         if set < 10;
%         %             if img < 10
%         %                 new = ['S0' num2str(set) 'I0' num2str(img) '.bmp'];
%         %             else
%         %                 new = ['S0' num2str(set) 'I' num2str(img) '.bmp'];
%         %             end
%         %         else
%         %             if img < 10
%         %                 new = ['S' num2str(set) 'I0' num2str(img) '.bmp'];
%         %             else
%         %                 new = ['S' num2str(set) 'I' num2str(img) '.bmp'];
%         %             end
%         %         end
%         
%         copyfile([root_dir unused_dir original],[set_dir new]);
%         movefile([root_dir unused_dir original],[root_dir used_dir original],'f');
%     end
% end
%%
%%---[3] Generate item and condition file with 2 predictable sequences---%%
% modified by Seth Konig October 19-20, 2014
% Essentially, exact same code as section 3.2 Make_ListSQ_Task_Files.m but for 4 sets
% of novel and repeat blocks  instead of 6, so 64 images total instead of 96

set_nums = ['R':'Z']; %DO NOT edit line always start add more if needed
number_of_sequences = 2;
number_of_crosshairs = 4; %# of items per sequence
number_sequence_trials_btwn_images = 2;
number_of_images = 64; %was 96 originally updated by SDK 1/7/15
image_spacing = 16; %number of images between novel and repeat presentations

%set to true if you want sequences to be reapeat several times in a row
%during a familrization block
familiarization_block = true; %if want familrization trials before image sets else set to false
familiar_trials = 10; %# of familirization trails if want them
if familiarization_block
    max_sequence_conditions = 256+number_of_sequences*familiar_trials; %was 384 updated by SDK 1/7/15
else
    max_sequence_conditions = 256; %was 384 updated by SDK 1/7/15
end

item_locations = {};
overlap_index = [];
buffer = 5; %minimum distance between 2 items
maxdist = 15; %maximum distance between 2 items
buffer = buffer - 1;% since items are 1 dva grid
sizex = 25;
sizey = 19;
clear x y
[cc,rr] = meshgrid(1:sizex,1:sizey);
rand('seed',300608); %was 150107 now 150304*2
for i = 1:500;
    valid_points = ones(sizey,sizex);
    crosshair_locations{i} = NaN(number_of_crosshairs*2,number_of_sequences);
    seq = 1;
    for cross = 1:number_of_crosshairs
        if cross == 1
            available_points = find(valid_points);
            choosen_ind = available_points(randi(length(available_points)));
            [y,x] = ind2sub([sizey,sizex],choosen_ind);
            C = sqrt((rr-y).^2+(cc-x).^2)<= 5;
            Cind = find(C);
            valid_points(Cind) = 0;
            crosshair_locations{i}(2*cross-1,seq) = x;
            crosshair_locations{i}(2*cross,seq) = y;
        elseif cross == 2
            available_points = find(valid_points);
            dist = sqrt((rr-y).^2+(cc-x).^2);
            C = dist <= 15;
            Cind = find(C);
            Cind = intersect(Cind,available_points);
            choosen_ind = Cind(randi(length(Cind)));
            [y,x] = ind2sub(size(valid_points),choosen_ind);
            C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
            Cind = find(C);
            valid_points(Cind) = 0;
            crosshair_locations{i}(2*cross-1,seq) = x;
            crosshair_locations{i}(2*cross,seq) = y;
        else %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
            %find angle of quickest scan path from previous 2 crosshairs
            dx12 = crosshair_locations{i}(2*cross-3,seq)-crosshair_locations{i}(2*cross-5,seq);
            dy12 = crosshair_locations{i}(2*cross-2,seq)-crosshair_locations{i}(2*cross-4,seq);
            angle12 = atan2d(dy12,dx12);
            angle12(angle12 < 0) = 360+angle12(angle12 < 0);
            
            available_points = find(valid_points);
            dist = sqrt((rr-y).^2+(cc-x).^2);
            C = dist <= 15;
            Cind = find(C);
            Cind = intersect(Cind,available_points);
            [Cy,Cx] = ind2sub(size(valid_points),Cind);
            potential_angles = NaN(1,length(Cy)); %angles between previous crosshair and potential locations
            for pc = 1:length(Cy)
                dx23 = Cx(pc)-crosshair_locations{i}(2*cross-3,seq);
                dy23 = Cy(pc)-crosshair_locations{i}(2*cross-2,seq);
                angle23 = atan2d(dy23,dx23);
                potential_angles(pc) = angle23;
            end
            potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
            angle_difference = potential_angles-angle12;
            angle_difference = abs(angle_difference);
            good_angles = find(angle_difference <=  90);
            if isempty(good_angles);
                break;
            else
                choosen_ind = good_angles(randi(length(good_angles)));
                crosshair_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                crosshair_locations{i}(2*cross,seq) = Cy(choosen_ind);
                C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                x =  Cx(choosen_ind);
                y =  Cy(choosen_ind);
                Cind = find(C);
                valid_points(Cind) = 0;
            end
        end
    end
    
    if all(~isnan(crosshair_locations{i}(:,1)))
        overlapind = rem(i,4)+1;
        overlap_index(i) = overlapind;
        seq = 2;
        switch overlapind
            case 1
                cross = 1;
                x = crosshair_locations{i}(2*cross-1,1);
                y = crosshair_locations{i}(2*cross,1);
                crosshair_locations{i}(2*cross-1,seq) = x;
                crosshair_locations{i}(2*cross,seq) = y;
                for cross = 2:number_of_crosshairs
                    if cross == 2
                        available_points = find(valid_points);
                        dist = sqrt((rr-y).^2+(cc-x).^2);
                        C = dist <= 15;
                        Cind = find(C);
                        Cind = intersect(Cind,available_points);
                        choosen_ind = Cind(randi(length(Cind)));
                        [y,x] = ind2sub(size(valid_points),choosen_ind);
                        C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
                        Cind = find(C);
                        valid_points(Cind) = 0;
                        crosshair_locations{i}(2*cross-1,seq) = x;
                        crosshair_locations{i}(2*cross,seq) = y;
                    else %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                        %find angle of quickest scan path from previous 2 crosshairs
                        dx12 = crosshair_locations{i}(2*cross-3,seq)-crosshair_locations{i}(2*cross-5,seq);
                        dy12 = crosshair_locations{i}(2*cross-2,seq)-crosshair_locations{i}(2*cross-4,seq);
                        angle12 = atan2d(dy12,dx12);
                        angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                        
                        available_points = find(valid_points);
                        dist = sqrt((rr-y).^2+(cc-x).^2);
                        C = dist <= 15;
                        Cind = find(C);
                        Cind = intersect(Cind,available_points);
                        [Cy,Cx] = ind2sub(size(valid_points),Cind);
                        potential_angles = NaN(1,length(Cy)); %angles between previous crosshair and potential locations
                        for pc = 1:length(Cy)
                            dx23 = Cx(pc)-crosshair_locations{i}(2*cross-3,seq);
                            dy23 = Cy(pc)-crosshair_locations{i}(2*cross-2,seq);
                            angle23 = atan2d(dy23,dx23);
                            potential_angles(pc) = angle23;
                        end
                        potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                        angle_difference = potential_angles-angle12;
                        angle_difference = abs(angle_difference);
                        good_angles = find(angle_difference <=  90);
                        if isempty(good_angles);
                            break;
                        else
                            choosen_ind = good_angles(randi(length(good_angles)));
                            crosshair_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                            crosshair_locations{i}(2*cross,seq) = Cy(choosen_ind);
                            x =  Cx(choosen_ind);
                            y =  Cy(choosen_ind);
                            C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                            Cind = find(C);
                            valid_points(Cind) = 0;
                        end
                    end
                end
            case 2
                cross = 2;
                x = crosshair_locations{i}(2*cross-1,1);
                y = crosshair_locations{i}(2*cross,1);
                crosshair_locations{i}(2*cross-1,seq) = x;
                crosshair_locations{i}(2*cross,seq) = y;
                cross = 1;
                available_points = find(valid_points);
                dist = sqrt((rr-y).^2+(cc-x).^2);
                C = dist <= 15;
                Cind = find(C);
                Cind = intersect(Cind,available_points);
                choosen_ind = Cind(randi(length(Cind)));
                [y,x] = ind2sub(size(valid_points),choosen_ind);
                C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
                Cind = find(C);
                valid_points(Cind) = 0;
                crosshair_locations{i}(2*cross-1,seq) = x;
                crosshair_locations{i}(2*cross,seq) = y;
                for cross = 3:4
                    %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                    %find angle of quickest scan path from previous 2 crosshairs
                    dx12 = crosshair_locations{i}(2*cross-3,seq)-crosshair_locations{i}(2*cross-5,seq);
                    dy12 = crosshair_locations{i}(2*cross-2,seq)-crosshair_locations{i}(2*cross-4,seq);
                    angle12 = atan2d(dy12,dx12);
                    angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                    
                    available_points = find(valid_points);
                    dist = sqrt((rr-y).^2+(cc-x).^2);
                    C = dist <= 15;
                    Cind = find(C);
                    Cind = intersect(Cind,available_points);
                    [Cy,Cx] = ind2sub(size(valid_points),Cind);
                    potential_angles = NaN(1,length(Cy)); %angles between previous crosshair and potential locations
                    for pc = 1:length(Cy)
                        dx23 = Cx(pc)-crosshair_locations{i}(2*cross-3,seq);
                        dy23 = Cy(pc)-crosshair_locations{i}(2*cross-2,seq);
                        angle23 = atan2d(dy23,dx23);
                        potential_angles(pc) = angle23;
                    end
                    potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                    angle_difference = potential_angles-angle12;
                    angle_difference = abs(angle_difference);
                    good_angles = find(angle_difference <=  90);
                    if isempty(good_angles);
                        break;
                    else
                        choosen_ind = good_angles(randi(length(good_angles)));
                        crosshair_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                        crosshair_locations{i}(2*cross,seq) = Cy(choosen_ind);
                        x =  Cx(choosen_ind);
                        y =  Cy(choosen_ind);
                        C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                        Cind = find(C);
                        valid_points(Cind) = 0;
                    end
                end
            case 3
                %going to do this in the reverse order as if this were cross #
                %2, then flip this sequence around so 2nd is now the 3rd cross.
                %It might be easier to do something else but I already have the
                %code.
                cross = 3;
                x = crosshair_locations{i}(2*cross-1,1);
                y = crosshair_locations{i}(2*cross,1);
                cross = 2;
                crosshair_locations{i}(2*cross-1,seq) = x;
                crosshair_locations{i}(2*cross,seq) = y;
                available_points = find(valid_points);
                dist = sqrt((rr-y).^2+(cc-x).^2);
                C = dist <= 15;
                Cind = find(C);
                Cind = intersect(Cind,available_points);
                choosen_ind = Cind(randi(length(Cind)));
                [y,x] = ind2sub(size(valid_points),choosen_ind);
                C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
                Cind = find(C);
                valid_points(Cind) = 0;
                cross  = 1;%rewrtie to 2 so that when we flip it will work out
                crosshair_locations{i}(2*cross-1,seq) = x;
                crosshair_locations{i}(2*cross,seq) = y;
                for cross = 3:4
                    %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                    %find angle of quickest scan path from previous 2 crosshairs
                    dx12 = crosshair_locations{i}(2*cross-3,seq)-crosshair_locations{i}(2*cross-5,seq);
                    dy12 = crosshair_locations{i}(2*cross-2,seq)-crosshair_locations{i}(2*cross-4,seq);
                    angle12 = atan2d(dy12,dx12);
                    angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                    
                    available_points = find(valid_points);
                    dist = sqrt((rr-y).^2+(cc-x).^2);
                    C = dist <= 15;
                    Cind = find(C);
                    Cind = intersect(Cind,available_points);
                    [Cy,Cx] = ind2sub(size(valid_points),Cind);
                    potential_angles = NaN(1,length(Cy)); %angles between previous crosshair and potential locations
                    for pc = 1:length(Cy)
                        dx23 = Cx(pc)-crosshair_locations{i}(2*cross-3,seq);
                        dy23 = Cy(pc)-crosshair_locations{i}(2*cross-2,seq);
                        angle23 = atan2d(dy23,dx23);
                        potential_angles(pc) = angle23;
                    end
                    potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                    angle_difference = potential_angles-angle12;
                    angle_difference = abs(angle_difference);
                    good_angles = find(angle_difference <=  90);
                    if isempty(good_angles);
                        break;
                    else
                        choosen_ind = good_angles(randi(length(good_angles)));
                        crosshair_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                        crosshair_locations{i}(2*cross,seq) = Cy(choosen_ind);
                        x =  Cx(choosen_ind);
                        y =  Cy(choosen_ind);
                        C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                        Cind = find(C);
                        valid_points(Cind) = 0;
                    end
                end
                xs = crosshair_locations{i}(1:2:end,seq);
                ys = crosshair_locations{i}(2:2:end,seq);
                xs = xs(end:-1:1);
                ys = ys(end:-1:1);
                crosshair_locations{i}(1:2:end,seq) = xs;
                crosshair_locations{i}(2:2:end,seq) = ys;
            case 4
                %going to do this in the reverse order as if this were cross #
                %2, then flip this sequence around so 1st is now the 4th cross.
                %It might be easier to do something else but I already have the
                %code.
                cross = 4;
                x = crosshair_locations{i}(2*cross-1,1);
                y = crosshair_locations{i}(2*cross,1);
                cross = 1;
                crosshair_locations{i}(2*cross-1,seq) = x;
                crosshair_locations{i}(2*cross,seq) = y;
                for cross = 2:number_of_crosshairs
                    if cross == 2
                        available_points = find(valid_points);
                        dist = sqrt((rr-y).^2+(cc-x).^2);
                        C = dist <= 15;
                        Cind = find(C);
                        Cind = intersect(Cind,available_points);
                        choosen_ind = Cind(randi(length(Cind)));
                        [y,x] = ind2sub(size(valid_points),choosen_ind);
                        C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
                        Cind = find(C);
                        valid_points(Cind) = 0;
                        crosshair_locations{i}(2*cross-1,seq) = x;
                        crosshair_locations{i}(2*cross,seq) = y;
                    else %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                        %find angle of quickest scan path from previous 2 crosshairs
                        dx12 = crosshair_locations{i}(2*cross-3,seq)-crosshair_locations{i}(2*cross-5,seq);
                        dy12 = crosshair_locations{i}(2*cross-2,seq)-crosshair_locations{i}(2*cross-4,seq);
                        angle12 = atan2d(dy12,dx12);
                        angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                        
                        available_points = find(valid_points);
                        dist = sqrt((rr-y).^2+(cc-x).^2);
                        C = dist <= 15;
                        Cind = find(C);
                        Cind = intersect(Cind,available_points);
                        [Cy,Cx] = ind2sub(size(valid_points),Cind);
                        potential_angles = NaN(1,length(Cy)); %angles between previous crosshair and potential locations
                        for pc = 1:length(Cy)
                            dx23 = Cx(pc)-crosshair_locations{i}(2*cross-3,seq);
                            dy23 = Cy(pc)-crosshair_locations{i}(2*cross-2,seq);
                            angle23 = atan2d(dy23,dx23);
                            potential_angles(pc) = angle23;
                        end
                        potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                        angle_difference = potential_angles-angle12;
                        angle_difference = abs(angle_difference);
                        good_angles = find(angle_difference <= 90);
                        if isempty(good_angles);
                            break;
                        else
                            choosen_ind = good_angles(randi(length(good_angles)));
                            crosshair_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                            crosshair_locations{i}(2*cross,seq) = Cy(choosen_ind);
                            x =  Cx(choosen_ind);
                            y =  Cy(choosen_ind);
                            C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                            Cind = find(C);
                            valid_points(Cind) = 0;
                        end
                    end
                end
                xs = crosshair_locations{i}(1:2:end,seq);
                ys = crosshair_locations{i}(2:2:end,seq);
                xs = xs(end:-1:1);
                ys = ys(end:-1:1);
                crosshair_locations{i}(1:2:end,seq) = xs;
                crosshair_locations{i}(2:2:end,seq) = ys;
        end
    end
    
    if all(~isnan(crosshair_locations{i}))
        crosshair_locations{i}(1:2:end,:) = crosshair_locations{i}(1:2:end,:)-13;
        crosshair_locations{i}(2:2:end,:) = crosshair_locations{i}(2:2:end,:)-10;
        for seq = 1:size(crosshair_locations{i},2);
            xs =crosshair_locations{i}(1:2:end,seq);
            ys =crosshair_locations{i}(2:2:end,seq);
            
            %double check to make sure distances are good
            d = pdist([xs,ys]);
            if any(d < 5)
                disp('error item locations too close')
                [min(d),seq]
            end
            dx = diff(xs);
            dy = diff(ys);
            if any(sqrt(dx.^2+dy.^2) > 15)
                disp('error locations too far appart')
                crossnum = find(sqrt(dx.^2+dy.^2) > 15)+1;
                if length(crossnum) > 1
                    [crossnum;seq]
                else
                    [crossnum seq]
                end
            end
            if any(abs(xs) > 12) || any(abs(ys) > 9)
                disp('error locations out of bounds')
            end
            
            dy12 = ys(2)-ys(1);
            dx12 = xs(2)-xs(1);
            angle12 = atan2d(dy12,dx12);
            
            %angle from crosshair 2 to 3
            dy23 = ys(3)-ys(2);
            dx23 = xs(3)-xs(2);
            angle23 = atan2d(dy23,dx23);
            
            %angle from crosshair 3 to 4
            dy34 = ys(4)-ys(3);
            dx34 = xs(4)-xs(3);
            angle34 = atan2d(dy34,dx34);
            
            angle12(angle12 < 0) = 360+angle12(angle12 < 0);
            angle23(angle23 < 0) = 360+angle23(angle23 < 0);
            angle34(angle34 < 0) = 360+angle34(angle34 < 0);
            
            dang123 = angle23-angle12;
            dang234 = angle34-angle23;
            dang123(dang123 < 0) =  dang123(dang123 < 0) + 360;
            dang234(dang234 < 0) =  dang123(dang234 < 0) + 360;
            dang123(dang123 > 180) = 360-dang123(dang123 > 180);
            dang234(dang234 > 180) = 360-dang234(dang234 > 180);
            if dang123 > 90 || dang234 > 90
                disp('locations not in a forward direction')
                seq
            end
        end
    end
end

%in rare case length(ovelap_index) < length(crosshair_locations)
%because the last few locations are aborted because they don't work
crosshair_locations(length(overlap_index)+1:end) = [];

all_nans = [];
for s = 1:length(crosshair_locations);
    if any(any(isnan(crosshair_locations{s}))) 
        all_nans = [all_nans s];
    end
end
crosshair_locations(all_nans) = [];
overlap_index(all_nans) = [];
overlap1 = find(overlap_index == 1);
overlap2 = find(overlap_index == 2);
overlap3 = find(overlap_index == 3);
overlap4 = find(overlap_index == 4);
new_crosshair_locations = [];
new_overlap_index = [];
%ensure that the items that overlaps in a sequences rotates evenly through
%1 to 4
for o = 1:min([length(overlap1),length(overlap2),length(overlap3),length(overlap4)]);
    new_crosshair_locations = [new_crosshair_locations crosshair_locations(overlap1(o))];
    new_crosshair_locations = [new_crosshair_locations crosshair_locations(overlap2(o))];
    new_crosshair_locations = [new_crosshair_locations crosshair_locations(overlap3(o))];
    new_crosshair_locations = [new_crosshair_locations crosshair_locations(overlap4(o))];
    new_overlap_index = [new_overlap_index overlap_index(overlap1(o))];
    new_overlap_index = [new_overlap_index overlap_index(overlap2(o))];
    new_overlap_index = [new_overlap_index overlap_index(overlap3(o))];
    new_overlap_index = [new_overlap_index overlap_index(overlap4(o))];
end
item_locations = new_crosshair_locations; %reassign cuz this is the variable the subsequent code uses
%%
clr=['255 255   0 x';
    '  0 255 255 x';
    '255   0 255 x';
    '  0 255   0 x';
    '  0   0   0 x'];

itemtype = [1 2 9 12 12 14]; %rectangle circle elipse triangle hexagon rosshair
int1 = [NaN NaN NaN 3 6 NaN];
rotation = {[0 45 90],0,[0 45 90],[0 30 180 210],0,45};

for itm = set_nums
    take = randperm(length(itemtype));
    itemtypeorder = itemtype(take);
    
    %updated by SDK 1/7/15
     set = ['ListSQV' itm];
%     if itm < 10
%         set = ['ListSQ0' num2str(itm)];
%     else
%         set = ['ListSQ' num2str(itm)];
%     end
    fid = fopen([set '.itm'],'w+');
    
    line1 ='ITEM TYPE FILLED CENTERX CENTERY BITPAN HEIGHT WIDTH ANGLE OUTER INNER -R- -G- -B- C A ------FILENAME------ int1';
    line2 =' -4   14      1    0.00    0.00      0   0.50  0.50  0.00              255 255 255 x';
    line3 =' -3   14      1    0.00    0.00      0   0.50  0.50  0.00               75  75  75 x';
    line4 =' -2    1      1    0.00    0.00      0   0.20  0.20  0.00  0.50         10  10  10 x';
    line5 =' -1    1      0    0.00    0.00      0   1.00  1.00  0.00              200 200 200 x';
    line6 ='  0    1      1    0.00    0.00      0   0.15  0.15  0.00              150 150 150 x';
    
    color_order = [clr(randperm(4),:); clr(end,:)];
    for line = 1:6
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    
    % updated 1/7/15
    temp_matrix = item_locations{find(itm == set_nums)}; %first 2 are controlled sequences, 2nd 2 random
 
    %one center clrchng trial for offset correction
    str='  1    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        150 150 150 x\r\n';
    fprintf(fid,str);
    str='  2    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        150 150 150 x\r\n';
    fprintf(fid,str);
    str='  3    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        175 175 130 x\r\n';
    fprintf(fid,str);
    
    itmnum = 4;
    for seq = 1:size(temp_matrix,2)
        rot = rotation{take(seq)}(randi(length(rotation{take(seq)})));
        for point = 1:size(temp_matrix,1)/2;
            x_pos(point) = temp_matrix(2*point-1,seq);
            y_pos(point) = temp_matrix(2*point,seq);
        end
        
        for point = 1:length(x_pos);
            if itmnum < 10
                itmspace = '  '; %2 spaces
            else
                itmspace = ' ';%1 space
            end
            if itemtypeorder(seq) >= 10
                type_space = '   ';%3 spaces
            else
                type_space = '    ';%3 spaces
            end
            filled_space = '      ';%5 spaces
            
            if x_pos(point) <= -10;
                x_space = '  ';%2 spaces
            elseif x_pos(point) < 0
                x_space = '   ';%3 spaces
            elseif x_pos(point) >= 10;
                x_space = '   ';%3 spaces
            else
                x_space = '    ';%4 spaces
            end
            if y_pos(point) <= -10;
                y_space = '  ';%2 spaces
            elseif y_pos(point) < 0
                y_space = '   ';%3 spaces
            elseif y_pos(point) >= 10;
                y_space = '   ';%3 spaces
            else
                y_space = '    ';%4 spaces
            end
            
            if rot < 100
                rotationspace = '  ';
                if rot < 10
                    rotation_add = '.00';
                else
                    rotation_add = '.0';
                end
            else
                rotationspace = '   ';
                rotation_add = '';
            end
            
            if any(itemtypeorder(seq) == [1,9]);
                height_width_space = '0.75  0.38';
            elseif any(itemtypeorder(seq) == [2,12]);
                height_width_space = '0.38  0.38';
            else
                height_width_space = '0.75  0.75';
            end
            
            if rem(abs(x_pos(point)),ceil(abs(x_pos(point)))) > 0.01 %is a decimal number, can't use zero its a float
                decimalx = '0';
            else
                decimalx = '.00';
            end
            
            if  rem(abs(y_pos(point)),ceil(abs(y_pos(point)))) > 0.01 %is a decimal number, can't  zero its a float
                decimaly = '0';
            else
                decimaly = '.00';
            end
            
            integer_space = '                        ';
            if any(itemtypeorder(seq) == [1,9,14]);
                inneroutter_space = '              ';
            elseif any(itemtypeorder(seq) == [2,12]);
                inneroutter_space = '  0.38        ';
            else
                inneroutter_space = '  0.75        ';
            end
            
            if isnan(int1(take(seq)))
                str = [itmspace num2str(itmnum) type_space  num2str(itemtypeorder(seq)) filled_space '1' x_space ...
                    num2str(x_pos(point)) decimalx  y_space num2str(y_pos(point)) decimaly '      0   ' height_width_space rotationspace...
                    num2str(rot) rotation_add inneroutter_space color_order(seq,:) '\r\n'];
            else
                str = [itmspace num2str(itmnum) type_space  num2str(itemtypeorder(seq)) filled_space '1' x_space ...
                    num2str(x_pos(point)) decimalx y_space num2str(y_pos(point)) decimaly '      0   ' height_width_space rotationspace...
                    num2str(rot) rotation_add inneroutter_space color_order(seq,:) integer_space num2str(int1(take(seq))) '\r\n'];
            end
            
            fprintf(fid,str);
            itmnum = itmnum+1;
        end
        
        for point = 1:length(x_pos);
            if itmnum < 10
                itmspace = '  '; %2 spaces
            else
                itmspace = ' ';%1 space
            end
            type_space = '   ';%3 spaces
            filled_space = '      ';%6 spaces
            
            if x_pos(point) <= -10;
                x_space = '  ';%2 spaces
            elseif x_pos(point) < 0
                x_space = '   ';%3 spaces
            elseif x_pos(point) >= 10;
                x_space = '   ';%3 spaces
            else
                x_space = '    ';%4 spaces
            end
            if y_pos(point) <= -10;
                y_space = '  ';%2 spaces
            elseif y_pos(point) < 0
                y_space = '   ';%3 spaces
            elseif y_pos(point) >= 10;
                y_space = '   ';%3 spaces
            else
                y_space = '    ';%4 spaces
            end
            
            if rem(abs(x_pos(point)),ceil(abs(x_pos(point)))) > 0.01 %is a decimal number, can't use zero its a float
                decimalx = '0';
            else
                decimalx = '.00';
            end
            
            if  rem(abs(y_pos(point)),ceil(abs(y_pos(point)))) > 0.01 %is a decimal number, can't  zero its a float
                decimaly = '0';
            else
                decimaly = '.00';
            end
            
            str = [itmspace num2str(itmnum) type_space  ' 1' filled_space '1' x_space ...
                num2str(x_pos(point)) decimalx y_space num2str(y_pos(point)) decimaly '      0   5.00  5.00  0.00'...
                '              ' color_order(end,:) '\r\n'];
            
            fprintf(fid,str);
            itmnum = itmnum+1;
        end
    end
    
    base_image_itmnum = itmnum;  % get the starting item number for image 1
    % to use later in condition file
    
    % updated by SDK 1/8/15
    d=dir([root_dir 'LSQV' itm  '\*.bmp']); %get the images in the set folder
     rr = randperm(length(d)); %to randomize the order of images in case something changed over time
    for img = 1:number_of_images;
        if itmnum < 100
            itmspace = ' '; %1 spaces
        else
            itmspace = '';%0 space
        end
        if itm < 10
            setzero = '0';
        else
            setzero = '';
        end
        if img < 10
            imgzero = '0';
        else
            imgzero = '';
        end
        
        % updated 1/8/15
        % was
        %         str = [itmspace num2str(itmnum)...
        %             '    8           0.00    0.00      0                                  '...
        %             '75  75  75 x   C:\\LSQ' setzero num2str(itm) ...
        %             '\\S' setzero num2str(itm) 'I' imgzero num2str(img) '.bmp' '\r\n'];
        % now
        str = [itmspace num2str(itmnum)...
            '    8           0.00    0.00      0                                  '...
            '75  75  75 x   C:\\LSQV' itm ...
            '\\' d(rr(img)).name '\r\n'];
        
        fprintf(fid,str,'%s');
        itmnum = itmnum+1;
    end
    fclose(fid);
    
    %----------------------------------------------------------------%
    %write unique random conditions files for each day to psueorandomize
    %sequence trials
    line1 = 'COND# BCKGND TIMING FIX_ID ---COLOR-PALETTE--- TEST0 TEST1 TEST2 TEST3 TEST4 TEST5 TEST6 TEST7'; %header
    line2 = '  1     -3    1       1                          2     3'; %color change line
    
    %updated by SDK 1/7/15
    set = ['ListSQV' itm];
    %     if itm < 10
%         set = ['ListSQ0' num2str(itm)];
%     else
%         set = ['ListSQ' num2str(itm)];
%     end
    fid = fopen([set '.cnd'],'w+');
    
    for line = 1:2
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    
    all_conditions = cell(1,number_of_sequences);
    trials_per_condition = max_sequence_conditions/number_of_sequences;
    
    %for predictable sequnces
    for seq = 1:size(item_locations{find(itm == set_nums)},2); %updated 1/8/15
        for t = 1:trials_per_condition;
            if seq == 1;
                all_conditions{1,seq} = [all_conditions{1,seq}; 4:11];
            elseif seq == 2;
                all_conditions{1,seq} = [all_conditions{1,seq}; 12:19];
            elseif seq == 3;
                all_conditions{1,seq} = [all_conditions{1,seq}; 20:27];
            elseif seq == 4;
                all_conditions{1,seq} = [all_conditions{1,seq}; 28:35];
            end
        end
    end
    
    %organize conditions as desired
    conditions = [];
    %put into familirization block if desired
    if familiarization_block
        order = randperm(numel(all_conditions));
        for i=1:numel(all_conditions);
            conditions = [conditions;all_conditions{order(i)}(1:familiar_trials,:)];
            all_conditions{order(i)}(1:familiar_trials,:) = [];
        end
    end
    for i=1:numel(all_conditions);
        conditions = [conditions;all_conditions{i}];
    end
    if familiarization_block
        order = 1:size(all_conditions,1)*familiar_trials*number_of_sequences;
        order = [order randperm(size(conditions,1)-length(order))+length(order)];
        conditions = conditions(order,:);
    else
        conditions = conditions(randperm(size(conditions,1)),:);
    end
    
    %write to file
    cndline = 2; %1st one is devoted to clrchng
    seq_cnd = 1;
    %if there is familirization block write it first to the file
    if familiarization_block %if want familrization trials before image sets else set to false
        for seq = 1:number_of_sequences
            for ft = 1:familiar_trials
                if cndline < 10
                    cndspace = '  ';
                elseif cndline < 100
                    cndspace = ' ';
                else
                    cndspace = '';
                end
                btfc = '     -3    2                             '; %background timing, fixid, color palate
                teststr = [];
                for t = 1:8;
                    if conditions(seq_cnd,t) < 10;
                        testspace = '     ';
                    else
                        testspace = '    ';
                    end
                    teststr = [teststr testspace num2str(conditions(seq_cnd,t))];
                end
                fprintf(fid,[cndspace num2str(cndline) btfc teststr '\r\n']);
                cndline = cndline+1;
                seq_cnd = seq_cnd + 1;
            end
        end
    end
    
    %write images with number_sequence_trials_btwn_images trials in between
    for block = 1:number_of_images/image_spacing
        for nov_rep = 1:2; %write 2x for novel then for repeat presentation
            for imgpair = 1:image_spacing/2;
                imgnums = [imgpair*2-1 imgpair*2]+(block-1)*image_spacing;
                if cndline < 10
                    cndspace = '  ';
                elseif cndspace < 100
                    cndspace = ' ';
                else
                    cndspace = '';
                end
                str = [cndspace num2str(cndline) '     -2    3      -4                         ' ...
                    num2str(imgnums(1)+base_image_itmnum-1) '\r\n'];
                fprintf(fid,str);
                cndline = cndline+1;
                for seq = 1:number_sequence_trials_btwn_images
                    if cndline < 10
                        cndspace = '  ';
                    elseif cndline < 100
                        cndspace = ' ';
                    else
                        cndspace = '';
                    end
                    btfc = '     -3    2                             '; %background timing, fixid, color palate
                    teststr = [];
                    for t = 1:8;
                        if conditions(seq_cnd,t) < 10;
                            testspace = '     ';
                        else
                            testspace = '    ';
                        end
                        teststr = [teststr testspace num2str(conditions(seq_cnd,t))];
                    end
                    fprintf(fid,[cndspace num2str(cndline) btfc teststr '\r\n']);
                    cndline = cndline+1;
                    seq_cnd = seq_cnd+1;
                end
                str = [cndspace num2str(cndline) '     -2    3      -4                         ' ...
                    num2str(imgnums(2)+base_image_itmnum-1) '\r\n'];
                fprintf(fid,str);
                cndline = cndline+1;
                for seq = 1:number_sequence_trials_btwn_images
                    if cndline < 10
                        cndspace = '  ';
                    elseif cndline < 100
                        cndspace = ' ';
                    else
                        cndspace = '';
                    end
                    btfc = '     -3    2                             '; %background timing, fixid, color palate
                    teststr = [];
                    for t = 1:8;
                        if conditions(seq_cnd,t) < 10;
                            testspace = '     ';
                        else
                            testspace = '    ';
                        end
                        teststr = [teststr testspace num2str(conditions(seq_cnd,t))];
                    end
                    fprintf(fid,[cndspace num2str(cndline) btfc teststr '\r\n']);
                    cndline = cndline+1;
                    seq_cnd = seq_cnd+1;
                end
            end
        end
    end
    fclose(fid);
end
