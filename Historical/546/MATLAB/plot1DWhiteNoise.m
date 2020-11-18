%  plot1DWhiteNoise.m
%
%  make noise with mean 0 and standard deviation 1
%  plots the first 5 octave bandpass, k = 1, 2:3, 4:7, 8:15, 16:31

clear
close all
N = 128;

typeOfNoise = 1;   %  1 gives Gaussian,  2 gives binary

switch typeOfNoise
    case 1
        % Gaussian noise
        noise = randn(1,N);
    case 2
        % binary noise with mean set to 0
        noise = (randn(1,N) > 0);
        noise = noise - mean(noise);
end
xrange = 0:N-1;
kRange = 0:N/2-1;

figure(1)
% Make two panels.  The upper shows one white noise signal.   The lower
% shows the power spectrum.

subplot(2,1,1)
plot(xrange,noise, '*-')
title(['white noise signal with mean 0, \sigma^2_n = ' num2str(var(noise))]);
axis( [ 0, N-1, min(noise)*1.05, max(noise)*1.05 ] )
subplot(2,1,2)
data = power(abs(fft(noise)),2);
plot(kRange, data(1:N/2),'*-') ;
hold on;


powerSpectrum = zeros(1,N);
M = 250;
for j = 1:M
    switch typeOfNoise
    case 1
        noise = randn(1,N);
    case 2
        % binary noise with mean set to 0
        noise = (randn(1,N) > 0);
        noise = noise - mean(noise);
    end
    powerSpectrum = powerSpectrum + power(abs(fft(noise)),2);
end
powerSpectrum = powerSpectrum/M;
subplot(2,1,2)
plot(kRange, powerSpectrum(1:N/2),'r*-') ;
title('power spectrum (up to Nyquist frequency)')
legend('one trial (example above)', ['average of ' num2str(M) ' trials']);
print -djpeg whiteNoise.jpg

