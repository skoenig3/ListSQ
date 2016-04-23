% Just briefly looking at KL divergence

dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\';
KLshuff = NaN(2000,5);
KLnorm = NaN(2000,5);
imageX = 800;
imageY = 600;

cd(dir);
a = what;
m = a.mat;
image_number = 0;
for file = 1:size(m,1);
    if ~isempty(strfind(m{file},'preprocessed'))
        fixationstats = [];
        load(m{file},'fixationstats','cfg','item_set');
        if ~isempty(fixationstats) %so is ListsQ
            [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
            [which_img,~] = get_image_numbers(cfg,itmlist,sequence_items,23);
            for img = 1:max(which_img);
                image_number = image_number+1;
                img_ind = find(which_img == img);
                if length(img_ind) == 2 %so novel and repeat
                    %these will contain PDFs (probability distribution functions) of fixation locations
                    novelfixations = cell(1,5);
                    repeatfixations = cell(1,5);
                    shuffled_novelfixations = cell(1,5);
                    shuffled_repeatfixations = cell(1,5);
                    for i = 1:5; %fill with a matrix of zeros
                        novelfixations{i} = zeros(imageY,imageX);
                        repeatfixations{i}=zeros(imageY,imageX);
                        shuffled_novelfixations{i} = zeros(imageY,imageX);
                        shuffled_repeatfixations{i} = zeros(imageY,imageX);
                    end
                    
                    nov_fixations = fixationstats{img_ind(1)}.fixations; %fixation locations from novel presentation
                    rep_fixations = fixationstats{img_ind(2)}.fixations; %fixation locations from repeat presentation
                    
                    %remove 1st fixation if this is a fixation on the cross-hair
                    if nov_fixations(1,1) > imageX/2-100 && nov_fixations(1,1) < imageX/2+100 &&...
                            nov_fixations(2,1) < imageY/2+100 && nov_fixations(2,1) > imageY/2-100
                        nov_fixations(:,1) = [];
                    end
                    nov_fixations = round(nov_fixations);
                    %remove 1st fixation if this is a fixation on the cross-hair
                    if rep_fixations(1,1) > imageX/2-100 && rep_fixations(1,1) < imageX/2+100 &&...
                            rep_fixations(2,1) < imageY/2+100 && rep_fixations(2,1) > imageY/2-100
                        rep_fixations(:,1) = [];
                    end
                    rep_fixations=round(rep_fixations);
                    
                    %we want to take the same number of fixations from the novel and repeat trials
                    maxfixations = min(size(nov_fixations,2),size(rep_fixations,2));
                    maxfixations(maxfixations > 25) = 25; %don't care if there are more than 25 fixations
                    
                    %since puting into groups of five remove the remainder of #/5
                    maxfixations = maxfixations-rem(maxfixations,5);
                    
                    if maxfixations >=5
                        for fixation = 1:maxfixations;
                            nov_fix_x = nov_fixations(1,fixation);%horizontal fixation position for novel presentation
                            nov_fix_y = nov_fixations(2,fixation);%vertical fixation position for novel presentation
                            rep_fix_x = rep_fixations(1,fixation);%horizonal fixation position for repeat presentation
                            rep_fix_y = rep_fixations(2,fixation);%vertical fixation position for repeat presentation
                            shuff_nov_fix_x = randi(imageX);
                            shuff_nov_fix_y = randi(imageY);
                            shuff_rep_fix_x = randi(imageX);
                            shuff_rep_fix_y = randi(imageY);
                            
                            %make sure fixations are within image borders
                            nov_fix_x(nov_fix_x < 1) = 1;
                            nov_fix_x(nov_fix_x > imageX) = imageX;
                            nov_fix_y(nov_fix_y < 1) = 1;
                            nov_fix_y(nov_fix_y > imageY) = imageY;
                            rep_fix_x(rep_fix_x < 1) = 1;
                            rep_fix_x(rep_fix_x > imageX) = imageX;
                            rep_fix_y(rep_fix_y < 1) = 1;
                            rep_fix_y(rep_fix_y > imageY) = imageY;
                            %shuffled may be 0 but not greater than imageX or imageY
                            shuff_nov_fix_x(shuff_nov_fix_x < 1) = 1;
                            shuff_nov_fix_y(shuff_nov_fix_y < 1) = 1;
                            shuff_rep_fix_x(shuff_rep_fix_x < 1) = 1;
                            shuff_rep_fix_y(shuff_rep_fix_y < 1) = 1;
                            
                            %put fixations in their appropriate PDF. Mark matrix with a
                            %1 where there was a fixation
                            if fixation <= 5
                                novelfixations{1}(nov_fix_y,nov_fix_x) = novelfixations{1}(nov_fix_y,nov_fix_x)+1;
                                repeatfixations{1}(rep_fix_y,nov_fix_x) = repeatfixations{1}(rep_fix_y,nov_fix_x)+1;
                                shuffled_novelfixations{1}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{1}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                                shuffled_repeatfixations{1}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{1}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
                            elseif fixation <= 10
                                novelfixations{2}(nov_fix_y,nov_fix_x) = novelfixations{2}(nov_fix_y,nov_fix_x)+1;
                                repeatfixations{2}(rep_fix_y,nov_fix_x) = repeatfixations{2}(rep_fix_y,nov_fix_x)+1;
                                shuffled_novelfixations{2}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{2}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                                shuffled_repeatfixations{2}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{2}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
                            elseif fixation <=15
                                novelfixations{3}(nov_fix_y,nov_fix_x) = novelfixations{3}(nov_fix_y,nov_fix_x)+1;
                                repeatfixations{3}(rep_fix_y,nov_fix_x) = repeatfixations{3}(rep_fix_y,nov_fix_x)+1;
                                shuffled_novelfixations{3}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{3}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                                shuffled_repeatfixations{3}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{3}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
                            elseif fixation <=20
                                novelfixations{4}(nov_fix_y,nov_fix_x) = novelfixations{4}(nov_fix_y,nov_fix_x)+1;
                                repeatfixations{4}(rep_fix_y,nov_fix_x) = repeatfixations{4}(rep_fix_y,nov_fix_x)+1;
                                shuffled_novelfixations{4}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{4}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                                shuffled_repeatfixations{4}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{4}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
                            elseif fixation <=25
                                novelfixations{5}(nov_fix_y,nov_fix_x) = novelfixations{5}(nov_fix_y,nov_fix_x)+1;
                                repeatfixations{5}(rep_fix_y,nov_fix_x) = repeatfixations{5}(rep_fix_y,nov_fix_x)+1;
                                shuffled_novelfixations{5}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{5}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                                shuffled_repeatfixations{5}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{5}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
                            end
                        end
                        
                        for i = 1:5
                             Distance = KL_Divergence(novelfixations{i},repeatfixations{i});
                             KLnorm(image_number,i) = Distance;                        
                             Distance = KL_Divergence(shuffled_novelfixations{i},shuffled_repeatfixations{i});
                             KLshuff(image_number,i) = Distance;
                        end
                    end
                end
            end
        end
    end
end
