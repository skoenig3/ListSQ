data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';

file = 'TO160105_3-preprocessed.mat';

load([data_dir file]);
%%
LFPchans = find_desired_channels(cfg,'LFP');

sample = cell(1,4);

for chan = 1:4
    for trl = 70:85;
        sample{chan} =  [sample{chan} data(LFPchans(chan)).values{trl}];
    end
end

thresh = NaN(1,4);

for chan = 1:4
    thresh(chan) =  mean(sample{chan})-2*std(sample{chan});
end
%%

events = cell(1,4);

for chan = 1:4
    events{chan} = find([0 sample{chan}] > thresh(chan) & [sample{chan} 0 ] < thresh(chan));
end
%%

forms = cell(1,4);
for chan = 1:4
    for evt = 1:length(events{chan})
        if events{chan}(evt) > 250 && events{chan}(evt) < length(sample{chan}) - 250
            forms{chan} = [forms{chan};sample{chan}(events{chan}(evt)-250:events{chan}(evt)+250)];
        end
    end
end
%%


Fs = 1000;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = 501;             % Length of signal
t = (0:L-1)*T;        % Time vector
f = Fs*(0:(L/2))/L;

figure(1)
figure(2)
for chan = 1:4
    [coeff,score,latent] = pca(forms{chan});
    
    sil = zeros(1,5); %determines the number of clusters by comparing the ratio
    for numclusts = 2:5
        T = kmeans(score(:,1:3),numclusts,'replicate',5);
        [silh] = InterVSIntraDist(score,T);
        sil(numclusts) = mean(silh);
    end
    sil(sil > 0.9*max(sil)) = 1;
    numclusters = find(sil == max(sil));
    numclusters = numclusters(1);
    
    T = kmeans(score(:,1:3),numclusters,'replicate',5);
    
    figure(1)
    subplot(2,2,chan)
    hold all
    for n = 1:numclusters
        plot(mean(forms{chan}(T == n,:)));
    end
    
    
    figure(2)
    subplot(2,2,chan)
    hold all
    for n = 1:numclusters
        Y = fft(mean(forms{chan}(T == n,:)));
        
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        
        plot(f,P1);
        xlim([2 102])
    end
end

%%
Fs = 1000;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = 501;             % Length of signal
t = (0:L-1)*T;        % Time vector
f = Fs*(0:(L/2))/L;

figure(1)
figure(2)
for chan = 1:4
    [coeff,score,latent] = pca(forms{chan});
    
    sil = zeros(1,5); %determines the number of clusters by comparing the ratio
    for numclusts = 2:5
        T = kmeans(score(:,1:3),numclusts,'replicate',5);
        [silh] = InterVSIntraDist(score,T);
        sil(numclusts) = mean(silh);
    end
    sil(sil > 0.9*max(sil)) = 1;
    numclusters = find(sil == max(sil));
    numclusters = numclusters(1);
    
    T = kmeans(score(:,1:3),numclusters,'replicate',5);
    
    figure(1)
    subplot(2,2,chan)
    hold all
    for n = 1:numclusters
        plot(mean(forms{chan}(T == n,:)));
    end
    
    
    figure(2)
    subplot(2,2,chan)
    hold all
    for n = 1:numclusters
        fo = forms{chan}(T == n,:);
      
        avg = zeros(1,size(fo,2)-250);
        for ff = 1:size(fo,1);
            Y = fft(fo(ff,:));
            P2 = abs(Y/L);
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            avg = avg+P1;
        end
        avg = avg/ff; 
            
        plot(f,avg);
        xlim([2 102])
    end
end



