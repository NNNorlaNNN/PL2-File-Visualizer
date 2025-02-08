% function PLVZ(pl2fullfile)
% input:
% pl2fullfile: .pl2 file directory
% selective output?
%
addpath(genpath("resource"))
addpath("output")


% enter pl2 file name
[pl2.file,pl2.path] = uigetfile("*.pl2");
% pl2fullfile = "data\LFP20240725\LFP\OC393002.pl2";
% pl2fullfile = "data\LFP20240723\FC210001.pl2";
% pl2fullfile = "data\LFP20240628\trkb\shrna-ol-a112001.pl2";

[pl2.path,pl2.name] = fileparts(pl2fullfile)

addpath(pl2.path)

pl2.ind = PL2GetFileIndex(pl2fullfile);

PL2Print(pl2.ind.AnalogChannels) % print out channels name

% Reads A/D data for A/D channel.
%cnl = [4 5 12 13] % enter the number of channel 

% for c = cnl(1:end)
%     pl2.ADdata(c) = PL2Ad(pl2.file, c);
% end

% adfreq:sample rate    n:data points    ts:
chnl = {'FP04','FP05','FP12','FP13'};

% get parameter & create empty vector for mean
[fs, ~, ~, n, ad] = plx_ad_v(pl2fullfile, chnl{1}); % get ad length
format short
t = 1/fs % Sampling period (1 / sampling frequency)
timelen = n/fs; % timelength is 600 seconds, n is the length of data
timevec = (0:n-1)*t % time vector

% =====Downsampling====
dsr = 5; % downsample rate
n = n/dsr;
fs = fs/dsr; % sampling frequncy(fs) change
t = 1/fs % Sampling period (1 / sampling frequency)
timelen = n/fs; % timelength is 600 seconds, n is the length of data
timevec = (0:ceil(n-1))*t % time vector
fprintf('downsample: %i', dsr)
% =====================

matrix_ad = zeros(4,ceil(n)); % specify the sum of ad value for future mean calculation
% Get ad_data for every channel
for x = 1:4
    [~, ~, ~, ~, ad] = plx_ad_v(pl2fullfile, chnl{x});
    if exist('dsr') % if data has been downsampled
        ad = ad(1:dsr:end);
    end
    matrix_ad(x,:) = ad; % store each ad in row
end

% ==mean ad matrix or specific channel==
% OPTION A: average channels
spe_chnl = 'mean'
mad = mean(matrix_ad)
% OPTION B: read specific channel
% spe_chnl = 4; 
% mad = matrix_ad(spe_chnl,:)

% Signal Preprocess

% bandpass
mad = bandpass(mad,[.1 100],fs);
% filter design
% d = designfilt('bandstopiir','FilterOrder',2, ...
%                'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
%                'DesignMethod','butter','SampleRate',fs);
% mad = filtfilt(d, mad);

% Plot ad data
plot(timevec,mad);

xlabel('Time(s)')
ylabel('Voltage(mV)')
hold off;


% short-time fast transfrom

opts = {'Window', hamming(2000), 'OverlapLength',1000, 'FFTLength', 2000, 'FrequencyRange', 'onesided'};

[S,F,T] = stft(mad,fs,opts{:});
% upper frequency range map==
upper_range = [50 100];  
% lower frequency range map==
lower_range = [0.1 50];  

% fre_range = [30 50]; 
upper_start_F = floor(length(F)*upper_range(1)/100);
upper_end_F = floor(length(F)*upper_range(2)/100);
lower_start_F = floor(length(F)*lower_range(1)/100);
lower_end_F = floor(length(F)*lower_range(2)/100);

% ==Prune upper scale heatmap
upper_F = F(upper_start_F:upper_end_F);
upper_S = S((upper_start_F:upper_end_F),:);
imagesc(T,upper_F,abs(upper_S(:,:,1)));
set(gca, 'YDir', 'normal');
colorbar;
title('upper range')
% caxis([8 10])

% Save Heatmap Data
gen_savepath = "output\heatmap\";
% prune_savepath = strcat(gen_savepath,pl2.name,'_chnl',num2str(spe_chnl),'_prune.mat')
upper_savepath = strcat(gen_savepath,pl2.name,'_chnl',num2str(spe_chnl),'_upper.mat');

save(upper_savepath,"T", "upper_F", "upper_S") % save heamap data

upperfig_savepath = erase(upper_savepath,'.mat'); % save heatmap fig

saveas(gca,upperfig_savepath,'fig') % save heatmap figure


% ==Prune lower scale heatmap
lower_F = F(lower_start_F:lower_end_F);
lower_S = S((lower_start_F:lower_end_F),:);
imagesc(T,lower_F,abs(lower_S(:,:,1)));
set(gca, 'YDir', 'normal');
colorbar;
title('lower range')
% caxis([8 10])

% Save Heatmap Data
gen_savepath = "output\heatmap\";
% prune_savepath = strcat(gen_savepath,pl2.name,'_chnl',num2str(spe_chnl),'_prune.mat')

lower_savepath = strcat(gen_savepath,pl2.name,'_chnl',num2str(spe_chnl),'_lower.mat');

save(lower_savepath,"T", "lower_F", "lower_S") % save heamap data

lowerfig_savepath = erase(lower_savepath,'.mat'); % save heatmap fig

saveas(gca,lowerfig_savepath,'fig') % save heatmap figure

% ==Full scale heatmap
imagesc(T,F,abs(S(:,:,1)));
set(gca, 'YDir', 'normal');
colorbar;
% caxis([0 30])

% Save Heatmap Data
% full_savepath = strcat(gen_savepath,pl2.name,'_chnl',num2str(spe_chnl),'_full.mat')
% save(full_savepath,"T", "F", "S") % save heamap data
% fullfig_savepath = erase(full_savepath,'.mat') % save heatmap fig
% saveas(gca,fullfig_savepath,'fig') % save heatmap figure

% discrete fourier transform
dft = fft(mad);
frevec = fs/n*(0:ceil(n-1)); % frequency vector
nyq = fs/2 %  Nyquist point
plot(frevec,abs(dft))
xlim([0,nyq])
xlabel('Frequency(Hz)')
ylabel('magnitude of Fourier coefficients')
title('DFT')

% beautiful waterfall

% frequency range limit==
% fre_range = [30 50]; 
% start_F = floor(length(F)*fre_range(1)/100);
% end_F = floor(length(F)*fre_range(2)/100);
% prune_F = F(start_F:end_F)
% prune_S = S((start_F:end_F),:)
% waterfall(prune_F,T,abs(prune_S(:,:,1))')
