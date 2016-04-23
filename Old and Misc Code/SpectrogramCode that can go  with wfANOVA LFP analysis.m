%% Average Spectrogram for difference between sequence 1 and sequence 2 ignoring 1st item
S1 = zeros(40,63);
S2 = zeros(40,63);
trial_count = zeros(1,2);
for t = 1:size(parsed_SQLFP_data,1);
    if item_num(t) ~= 1
         [S F T] = spectrogram(parsed_SQLFP_data(t,1:end),16,8,2:2:80,1000);
        if seq(t) == 1;
            S1 = S1+real(S);
            trial_count(1) = trial_count(1)+1;
        else
            S2 = S2+real(S);
            trial_count(2) = trial_count(2)+1;
        end
    end
end
S1 = S1/trial_count(1);
S2 = S2/trial_count(2);

S1 = S1(end:-1:1,:);
S2 = S2(end:-1:1,:);

figure
subplot(3,1,1)
imagesc(S1);
set(gca,'Xtick',1:10:length(T))
set(gca,'XtickLabel',num2cell(1000*T(1:10:end)))
xlabel('Time (ms)')
set(gca,'Ytick',1:5:length(F))
set(gca,'YtickLabel',num2cell(F(end:-5:1)));
ylabel('Frequence (Hz)')
title('Relative Power in LFP for Sequence 1')

subplot(3,1,2)
imagesc(S2);
set(gca,'Xtick',1:10:length(T))
set(gca,'XtickLabel',num2cell(1000*T(1:10:end)))
xlabel('Time (ms)')
set(gca,'Ytick',1:5:length(F))
set(gca,'YtickLabel',num2cell(F(end:-5:1)));
ylabel('Frequence (Hz)')
title('Relative Power in LFP for Sequence 2 ')

subplot(3,1,3)
imagesc(S2-S1);
set(gca,'Xtick',1:10:length(T))
set(gca,'XtickLabel',num2cell(1000*T(1:10:end)))
xlabel('Time (ms)')
set(gca,'Ytick',1:5:length(F))
set(gca,'YtickLabel',num2cell(F(end:-5:1)));
ylabel('Frequence (Hz)')
title('Difference in Relative Power between Sequence 1 and Sequence 2 ignoring 1st fixation')

for sb = 1:2
    subplot(3,1,sb)
    caxis([min(clims(1,:)) max(clims(2,:))])
end

%% Average Spectrogram for difference between novel and repeat images
S1 = zeros(40,127);
S2 = zeros(40,127);
trial_count = zeros(1,2);
for t = 1:size(parsed_imgLFP_data,1);
    [S F T] = spectrogram(parsed_imgLFP_data(t,1:end),32,16,2:2:80,1000);
    if novel_vs_repeat(t) == 1;
        S1 = S1+real(S);
        trial_count(1) = trial_count(1)+1;
    else
        S2 = S2+real(S);
        trial_count(2) = trial_count(2)+1;
    end
end
S1 = S1/trial_count(1);
S2 = S2/trial_count(2);

S1 = S1(end:-1:1,:);
S2 = S2(end:-1:1,:);

clims = NaN(2,3);
figure
subplot(3,1,1)
imagesc(S1);
set(gca,'Xtick',1:10:length(T))
set(gca,'XtickLabel',num2cell(1000*T(1:10:end)))
xlabel('Time (ms)')
set(gca,'Ytick',1:5:length(F))
set(gca,'YtickLabel',num2cell(F(end:-5:1)));
ylabel('Frequence (Hz)')
title('Relative Power in LFP for Novel Images')
clims(:,1) = caxis;

subplot(3,1,2)
imagesc(S2);
set(gca,'Xtick',1:10:length(T))
set(gca,'XtickLabel',num2cell(1000*T(1:10:end)))
xlabel('Time (ms)')
set(gca,'Ytick',1:5:length(F))
set(gca,'YtickLabel',num2cell(F(end:-5:1)));
ylabel('Frequence (Hz)')
title('Relative Power in LFP for Repeat Images')
clims(:,3) = caxis;

subplot(3,1,3)
imagesc(S2-S1);
set(gca,'Xtick',1:10:length(T))
set(gca,'XtickLabel',num2cell(1000*T(1:10:end)))
xlabel('Time (ms)')
set(gca,'Ytick',1:5:length(F))
set(gca,'YtickLabel',num2cell(F(end:-5:1)));
ylabel('Frequence (Hz)')
title('Difference in Relative Power between Novel and Repeat images')
clims(:,3) = caxis;

for sb = 1:2
    subplot(3,1,sb)
    caxis([min(clims(1,:)) max(clims(2,:))])
end