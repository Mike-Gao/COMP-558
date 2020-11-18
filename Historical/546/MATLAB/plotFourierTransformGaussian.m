%  plotFourierTransformGaussian.m
%  author:  Michael Langer
%
%  This code shows that the Fourier transform of a Gaussian with
%  standard deviation sigx is a Gaussian (in the frequency domain)
%  with standard deviation sigk such that sigx * sigk =  N/(2*pi).

N    = 128; 
sigx =  3;

close all
%  define the Gaussian on x = -N/2-1, ...,  N/2
%  Note that Matlab indexes this with indices 1,.. N
g = 1/sqrt(2*pi)/sigx * exp( -power(-N/2-1:N/2,2)/ (2*sigx*sigx));

%  fftshift is used to make the function defined on mod N 
plot(real( fft( fftshift(g) )));
%  There is a tiny imaginary component because of roundoff error.
hold on;

% fit it to 
sigk = N/ (2*pi*sigx);
plot( fftshift(exp( -power(-N/2-1:N/2,2)/ (2*sigk*sigk))),'r');
legend('Fourier transform of Gaussian, g(x)', 'Gaussian defined over k (using formula)')
xlabel('frequency k');
title('Fourier transform of a Gaussian');