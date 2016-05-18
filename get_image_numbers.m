function [which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code)
% define which imagse are novel and which images are repeat
which_img = NaN(1,96*2);
img_cnd = NaN(1,96*2);
img_count = 1;
for t = 1:length(cfg.trl);
    if itmlist(cfg.trl(t).cnd-1000) > sequence_items(end)
        if any(cfg.trl(t).allval == img_on_code);
            if length (cfg.trl(t).cnd) > 1
                which_img(img_count) = itmlist(cfg.trl(t).cnd(1)-1000);
                img_cnd(img_count) = cfg.trl(t).cnd(1);
                warning('More than 1 condition for this trial')
            else
                which_img(img_count) = itmlist(cfg.trl(t).cnd-1000);
                img_cnd(img_count) = cfg.trl(t).cnd;
            end
            img_count = img_count+1;
        end
    end
end
which_img = which_img-which_img(1)+1;

novel_vs_repeat = NaN(1,96*2);
for img = 1:max(which_img)
    imgind = find(which_img == img);
    if ~isempty(imgind);
        novel_vs_repeat(imgind(1)) = 1;
    end
    if length(imgind) > 1
        novel_vs_repeat(imgind(2)) = 2;
    end
end
end