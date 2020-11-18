%  Plot the 1D blur function B*B*..*B  and its Fourier transform.
%  First B is shown, then B*B, then B*B*B, etc up to a convolution
%  with m=20  B(x)'s.    Each plot is a slightly different grey value from
%  the previous one.

N = 32;

B = zeros(1,N);
B(1) = .5;
B(2) = .25;
B(N) = .25;

figure

h = zeros(1,N);
h(1) = 1;          % start with a delta function

close all
figure
hold on
xvalues = (-N/2:N/2-1);
for m = 1:20
    
  %  circular convolution
    
  h = cconv(h, B, N);
  
  %  choose the color for convolving m times
  
  c = .5-m/40;
  subplot(2,1,1);
  
  plot(xvalues,fftshift(h), '-*', 'color', [c c c]  );
  axis([min(xvalues) max(xvalues) 0 .5]);
  xlabel('x')
  title('B(x) * ...  * B(x),  i.e.   m times');
  hold on

  %   now compute and plot the Fourier transform
  
  subplot(2,1,2);
  
  %   the "abs" (amplitude) below isn't necessary strictly speaking 
  %   since the Fourier transform of B() is always real.  But
  %   Matlab sometimes produces floating point approximations 
  %   that lead to tiny leakage into the imaginary part and then
  %   when you try to plot the value, it gives an error.
  
  plot( 0:N-1, abs(fft(h)), '-*', 'color', [c c c]  );
  axis([0, N-1, 0, 1]);
  xlabel('frequency k');
  title('B(k)^m = (0.5 + 0.5*cos( 2 pi / N k))^m');
  hold on
  pause(1)
end