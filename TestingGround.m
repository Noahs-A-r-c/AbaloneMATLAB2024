% Adjust the format for displaying numbers
format short
n = 1024;
interpolFact = 70;
amplify = 1;



% Load the data from CSV
T = readtable("Magnetometer Magnitude to Distance - 1 Dimension - 8_23_24 9_54 PM - 2D CSV 1D Movement for MATLAB.csv");

% Amplify the Bx and By readings
BxRead = interp(T.Bx * amplify, interpolFact);
ByRead = interp(T.By * amplify, interpolFact);

figure(2)
hold on;


%%

[b,a] = butter(7, 1 / (fs/2));
By_filt = filter(b,a,ByRead);
% By_filt = ByRead
% freqz(b,a,[],fs)

% Perform FFT on BxRead (using 512-point FFT)
Y_fft = fft(By_filt);

% Get the magnitude of the FFT result
Y_mag = abs(Y_fft);

% Plot the magnitude of the FFT
hold on;
figure(1);
scatterlin = linspace(0,n,n);
stem(scatterlin(2:end),Y_mag(2:end));
title('FFT of BxRead');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

figure(2);
newBy = ifft(Y_fft);
plot(linspace(0,length(newBy),length(newBy)),newBy);
plot(linspace(0,length(ByRead),length(ByRead)),ByRead);
grid on;
