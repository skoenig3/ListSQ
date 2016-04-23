load(['C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\'...
    'PW150127_3-preprocessed.mat'])

[pupil_channel] = find_desired_channels(hdr,data,'pupil');

[itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
[which_img,novel_vs_repeat] = get_image_numbers(cfg,itmlist,sequence_items,23);

pupil_nov = NaN(96,6000);
pupil_rep = NaN(96,6000);

img_num = 1; %more like index
nov_img = 1;
rep_img = 1;

num_trials = length(cfg.trl);
for t = 1:num_trials
    if any(cfg.trl(t).allval == 23); %in which image was displayed or 1st item in sequence was displayed
        if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(end) 
            continue % go to the next trial since this is a sequece trial
        end
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == 15);
        fixation_start = cfg.trl(t).alltim(cfg.trl(t).allval == 23);
        fixation_time = fixation_start(1)-trial_start;
        % multiple fixations may occur within image periods if monkeys looks away from image
        
       
        
        if novel_vs_repeat(img_num) == 1
            pupil_nov(nov_img,:) =data(pupil_channel).values{t}(fixation_time-999:fixation_time+5000);
            nov_img = nov_img+1;
        else
            pupil_rep(rep_img,:) =data(pupil_channel).values{t}(fixation_time-999:fixation_time+5000);
            rep_img = rep_img+1;
        end
        img_num = img_num + 1; %more like index
    end
end
%%
pupil_nov = laundry(pupil_nov);
pupil_rep = laundry(pupil_rep);
%%
figure
hold on
dofill(1:6000,pupil_nov/1000/abs(max(mean(pupil_nov))),'blue',1,120);
dofill(1:6000,pupil_rep(1:end-1,:)/1000/abs(max(mean(pupil_rep))),'red',1,120);
hold off
xlabel('Time from Image Onset')
ylabel('Pupil Diameter (a.u.)')
set(gca,'Xtick',[0:500:6000])
set(gca,'XtickLabel',num2cell([0:500:6000]-1000))
xlim([700 5500])
legend('Novel','Repeat')
%%
n = sampsizepwr('t',[mean(mean(pupil_nov(:,1500:4500))) ...
    std(mean(pupil_nov(:,1500:4500)'))], mean(mean(pupil_rep(1:end-1,1500:4500))),0.9)

%% Peform a power analysis on the pupil data to determine what frequency we need 
% collect the data at

Fs = 1000;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = size(pupil_nov,2);                     % Length of signal
t = (0:L-1)*T;                % Time vector
f = Fs/2*linspace(0,1,NFFT/2+1);

NFFT = 2^nextpow2(L); % Next power of 2 from length of y

Ynov = [];
Yrep = [];

for trial = 1:size(pupil_nov,1);
    Y = fft(pupil_nov(trial ,:),NFFT)/L;
    Ynov = [Ynov; Y];
end


for trial  = 1:size(pupil_rep,1);
    Y = fft(pupil_rep(trial ,:),NFFT)/L;
    Yrep = [Yrep; Y];
end

Ynov = mean(Ynov);
Yrep = mean(Yrep);

figure
plot(f,2*abs(Ynov(1:NFFT/2+1))) 

figure
plot(f,2*abs(Yrep(1:NFFT/2+1))) 

%%
Fs = 1000;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = size(pupil_nov,2);                     % Length of signal
t = (0:L-1)*T;                % Time vector
f = Fs/2*linspace(0,1,NFFT/2+1);

NFFT = 2^nextpow2(L); % Next power of 2 from length of y

Ynov = [];
Yrep = [];

Ynov = fft(nanmean(pupil_nov),NFFT)/L;

Yrep = fft(nanmean(pupil_rep),NFFT)/L;

figure
plot(f,sqrt(2*abs(Ynov(1:NFFT/2+1))))

figure
plot(f,sqrt(2*abs(Yrep(1:NFFT/2+1))))
