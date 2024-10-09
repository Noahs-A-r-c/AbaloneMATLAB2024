% ApplyButterB.m Apply a butterworth filter to remove high frequecy noise.

function B_ReadFilt = ApplyButterB(B_Read,interpolFact,cutoff,filterOrd)
fs = 5 * interpolFact; % sampling frequency times interpolation

[b,a] = butter(filterOrd, cutoff / (fs/2));
B_ReadFilt = filter(b,a,B_Read);
% freqz(b,a,[],fs)

% % Perform FFT on BxRead 
% Y_fft = fft(By_filt);
% 
% % Plot the magnitude of the FFT
% hold on;
% figure(1);
% scatterlin = linspace(0,n,n);
% stem(scatterlin(2:end),Y_mag(2:end));
% title('FFT of BxRead');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');

% Get the magnitude of the FFT result
% Y_mag = abs(Y_fft);


% % Take the inverse FFT to produce the filtered result
% B_ReadFilt = ifft(Y_fft);
end