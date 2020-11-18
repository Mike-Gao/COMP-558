% Here is some very basic code that illustrates a likelihood function.
%
% Suppose we have N samples of an image I.  The samples are  
% an intensity I_0 + noise (Gaussian with mean 0 and variance sig*sig).
%
% You would like to estimate the likelihood of various values I_0.   
% For any hypothetical I_0,  you look at the difference 
% between the I(x) and the hypothetical I_0, and you
% consider the probability of obtaining those image value(s).

clear; close all

N = 400;            % number of pixels
I_0 = 2;             % actual I_0
sigNoise = 3;       % standard deviation of the noise

%  Compute the likelihoods for various hypothetical I_0 in a 
%  neighborhood (plus or minus 2) around the true value of I_0.
%  We want to choose the I_0 that maximizes the likelihood.

rangeOfI_0 = I_0 - 2 : 0.05 : I_0 + 2;
numberI_0s = length(rangeOfI_0);
logLikelihood = zeros(1, numberI_0s);

%  generate the random image which is I_0 plus Gaussian noise

I =  I_0 + sigNoise * randn(N,1);

%  Compute and plot the likelihoods.  
%  Rather than computing the product of the likelihoods for the different
%  noise values, I use the fact that e^a * e^b = e^(a+b) and I compute 
%  the sum of log likelihoods.   Here I am relying on the assumption
%  that the noise is Gaussian.
%  I also ignore constant factors  1/(sqrt(2 pi) sigNoise) since they
%  don't depend on the actual noise values.

for j = 1:numberI_0s
   err = I - rangeOfI_0(j);   % these are the noise values, for each model I_0
   
%  product of likelihoods can be computed by taking sum of log likelihoods i.e.  e.g.  sum of  
   logLikelihood(j) =  sum( - (err .* err)/(2*sigNoise*sigNoise)  ) ;
end

%plot(rangeOfI_0, logLikelihood,'-*');
plot(rangeOfI_0, exp(logLikelihood),'-*');
xlabel('I_0','fontsize',16);
ylabel('p(  I(x)   |   I_0 )','fontsize',16)
title([ 'N = ' num2str(N) ...
       ', I_0 = ' num2str(I_0) ...
    ', \sigma_n= ' num2str(sigNoise), ',  mean is ',  num2str(mean(I),3)  ],'fontsize',16)
 