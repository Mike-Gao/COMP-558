%  sumOfSinusoids.m
%
%  This script plots the sum of cos( 2 pi / N * k (x/M) ) 
%  and sin( 2 pi / N * k (x/M) ) for k=1:N and x/M = 0,1,.. N-1
%  where x = 0,1, ... N*M-1.
%  The sums are plotted in black.
%  The individual cos and sin functions are plotted, in red.
%

N = 128;

sumOfCos = zeros(1,N);
sumOfSin = zeros(1,N);
figure
for k = 0:N-1
    subplot(2,1,1)
    sumOfCos = sumOfCos + cos( 2*pi/N * k* (0:N-1));
    plot( 0:N-1, sumOfCos,  '-k.' );
    hold on
    plot( 0:N-1,0,'k-');
    hold off
    axis([0, N-1, min(-k,-1), max(1,k)])
    title(['sum of cos(2 pi/N * k x) ) from k = 0 to ' int2str(k) ]);
    subplot(2,1,2);
    sumOfSin = sumOfSin + sin( 2*pi/N * k * (0:N-1));
    plot( 0:N-1, sumOfSin,  '-k.' );
    hold on
    plot( 0:N-1,0,'k-');
    hold off
    axis([0, N-1, min(-k,-1), max(1,k)])
    title(['sum of sin(2 pi / N * k x ) from k = 0 to ' int2str(k)  ]);
    pause(.1)
end  