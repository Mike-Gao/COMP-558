%  plot2DBands.m
%
%  plots the first 5 octave bandpass, k = 1, 2:3, 4:7, 8:15, 16:31

clear

imageNumber = 3;

switch imageNumber
    case 1
        NX = 128;
        NY = 128;
        I = randn(NY,NX);
    case 2
        NX = 128;
        NY = 128;
        I  = make2Dcosine(NX,6,8);
    case 3
%        imageName = '../../images/forest.jpg';  
%        imageName = '../../images/stones.jpg';
%        imageName = '../../images/landscape1.jpg';    
%        imageName = '../../images/Ansel-Adams-Mountains.jpg'
         imageName = '../../images/McGillCampus1.jpg'
        
        I = imread(imageName);
        %  avoid dimensions which are odd because it will give an indexing
        %  error when I take the Fourier transform
        NY = 2*floor(size(I,1)/2);
        NX = 2*floor(size(I,2)/2);
        I = (double(squeeze(I(1:NY,1:NX,1))) + ...
            double(squeeze(I(1:NY,1:NX,2))) + ...
            double(squeeze(I(1:NY,1:NX,3))))/3;
end

subplot(2,3,1)
imagesc(I); 
colormap gray; 
title(['image, \sigma^2 = ' num2str( var(I(:)) , '%6.0f')] )

%  now filter it with one octave bandpass

IHat = fft2(I);

distFromCenter = sqrt( power( (0:NY-1) - NY/2, 2 )'*ones(1,NX) + ...
                       ones(NY,1) * power( (0:NX-1) - NX/2, 2 ));
for j = 1:5
  %  each band in 2D is defined by an annulus
    
  %  set to 1 for frequencies within band j and 0 otherwise
  %  bands are for |k| in [1,2), [2:4), [4:8), [8:16), [16:32)
  
  bandpass = fftshift( (distFromCenter >= power(2,j-1)) & ...
                       (distFromCenter < power(2,j))   );
  
  %figure(2); image(255*bandpass(1:16,1:16)); colormap gray; %pause(5); figure(1);
  
  %  Use the following instead to look at lowpass instead of bandpass.
  %
  %  bandpass = fftshift(  distFromCenter < power(2,j)   );
  
  %  only plot a given band if the image is large enough
  if (power(2,j) <= min(NX/2, NY/2))
      subplot(2,3,j+1);
      Ifiltered = real( ifft2(bandpass .* IHat) );
      imagesc(Ifiltered);  colormap gray
      title( ['\sigma^2 = ' num2str( var(Ifiltered(:)) , '%6.0f') ]);
  end
end
% print -djpeg plot2DBands-McGillCampus.jpg