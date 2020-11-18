%  contrastSensitiveIdealObserver.m
%
%  Here is the code used to generate the plot in lecture 13 where I show an
%  ideal observer who is trying to decide which of two noisy images
%  contains a cosine function.    It doesn't fit the data.  Its just an illustration.
%
%  Notice if the ideal observer runs a much smaller number of samples, then
%  the curve is much less smooth.

clear
close all
N = 512;

I0 = 100;
sigmaN = 3;
k0  =  6;
numTrials = 2000;
numDeltaI = 15;
correct = zeros(numDeltaI,1);
cosine = cos( 2*pi * k0 / N * (0:N-1));  
scaleDeltaI = .1;

for deltaI = 1:numDeltaI
    for trial = 1:numTrials
        I1 = I0 + sigmaN * randn(1,N);
        I2 = I0 + sigmaN * randn(1,N) + scaleDeltaI * deltaI * cosine;

        I1hat = abs(fft(I1));
        I2hat = abs(fft(I2));
        if  ( max( I1hat(2:N/2)) / mean( I1hat(2:N/2)) < max( I2hat(2:N/2)) / mean( I2hat(2:N/2)))
            correct(deltaI) = correct(deltaI) + 1;
        end
    end
end
figure(1)
plot(scaleDeltaI* (1:numDeltaI), correct/numTrials *100, '*')
xlabel('\Delta I');
ylabel('percent correct');
title('psychometric curve for detecting \Delta I * cos( ) versus \Delta I * cos( ) + noise')
%figure(2)
%plot(I1hat);  hold on;  plot(I2hat)
