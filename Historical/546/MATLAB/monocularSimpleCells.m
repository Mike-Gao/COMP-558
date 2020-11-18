%  monocularSimpleCells.m
%
%  Michael Langer  COMP 546 Fall 2015

clear all
close all
N = 64;
k = 2;
[cosGabor, sinGabor] = make2DGabor(N,k,0);
Gabor = cosGabor + 1i * sinGabor;

IcosGabor = zeros(N,N,3);
IsinGabor = zeros(N,N,3);
%  The following is a bit annoying, but if I don't make these images RGB
%  then I end up with colored plots when I show them because the colormap
%  has to be for a whole figure, not just for subplots
for c = 1:3
   IcosGabor(:,:,c) = remapImage(real(Gabor));    
   IsinGabor(:,:,c) = remapImage(imag(Gabor));
end

screen_size = get(0, 'ScreenSize');

%This will return a 4 element array: (left, bottom, width, height);
%The ones we are interested are the "width" and "height".
%Now we have the size of the screen, we can make our figure full screen:

f1 = figure(1);
set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );


%  show cosine and sine Gabors

subplot(2,3,1);
subimage(IcosGabor/255);
title('cosine Gabor')
hold on

subplot(2,3,4);
subimage(IsinGabor/255);
title('sine Gabor')
hold on

I  = zeros(N,N,3);
response = zeros(1,N);

degPerRadian = 180/pi;
noiseLevel = 0.2;

for x = 1: N
    
    %  make noisy left image of a line
    
    I  = noiseLevel*rand(N,N);
    I( :, x, : ) = 1;
        
    if (x == N/2)
        Isaveforfinalplot = I;
    end
    
    subplot(2,3,2);
    imagesc(I);
    axis square
    colormap gray
    title('image of a vertical line + noise');
    
    response(x)  = sum(sum(I  .* Gabor));
    theta = 10/degPerRadian;
    
    %  polor plot of (cosGabor, sinGabor) responses
    
    subplot(2,3,3);
    plot(response(x), '*k');
    axis square
    title('polar plot of cos and sine Gabor responses')
    axis([-.05 .05 -.05 .05]);
    hold on
    plot(response(1:x), 'k');
    line([0 0],  [-.05 .05],  'Color','k')
    line([-.05 .05], [0 0],  'Color','k')   
    hold off
    pause(.05);
end

subplot(2,3,5);
%plot(abs(response),'k');
title('responses')
maxresponse = max(abs(response(:)));
axis([0 N -maxresponse 2*maxresponse]);
hold on
%plot(abs(response), 'k','LineWidth',2);
plot(real(response), '-k'); axis square
plot(imag(response), '--k'); axis square
%legend('amplitude', 'cos Gabor', 'sin Gabor','Location','Best');
legend('cos Gabor', 'sin Gabor','Location','Best');
xlabel('x position of vertical line')

subplot(2,3,6);
plot( degPerRadian*angle(response),'k'); axis square
title('phase (degrees)')
axis([0 N -180 180 ]);

subplot(2,3,2);
imagesc(Isaveforfinalplot); axis square
title('image of a vertical line + noise');    
axis square
