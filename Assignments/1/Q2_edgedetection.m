rng(1);
rect = ones(500,500) * 256;
for i=1:200
    height = randi(300);
    width = randi(300);
    height_max = height + randi(180);
    width_max = height + randi(180);
    intensity=randi([0,255]);
    rect(height:height_max,width:width_max) = intensity;
end
rect = uint8(rect(20:269,20:269));
figure(1);
imshow(rect);

Gaussian = fspecial('gaussian',20,3);
x=conv2(Gaussian, [1 -2 1], 'same');
y=conv2(Gaussian, [1 -2 1]', 'same');
Laplace = x + y;
figure(2);
imshow(Laplace, 'DisplayRange', [-0.003 0.001], 'InitialMagnification', 1000);
colorbar

ZeroCrossingImg = CreateZeroCrossingImage(rect,20,3);
figure(3);
imshow(ZeroCrossingImg);

NoisyRect = rect+uint8(0.1*std(double(rect(:)))*randn(size(rect)));
figure(4);
imshow(NoisyRect);

figure(5);
NoisyZeroCrossingImg = CreateZeroCrossingImage(NoisyRect,20,3);
imshow(NoisyZeroCrossingImg);

figure(6);
NewNoisyZeroCrossingImg = CreateZeroCrossingImage(NoisyRect,40,6);
imshow(NewNoisyZeroCrossingImg);


function img = CreateZeroCrossingImage(Img,Width,Stddev)
    Gaussian = fspecial('gaussian', Width, Stddev);
    x=conv2(Gaussian, [1 -2 1], 'same');
    y=conv2(Gaussian, [1 -2 1]', 'same');
    Laplace = x+y;
    LaplaceImg = conv2(Img, Laplace,'same');
    ShiftLeft=(circshift(LaplaceImg,[0 -1]) .* LaplaceImg)<0;
    ShiftRight=(circshift(LaplaceImg,[0 1]) .* LaplaceImg)<0;
    ShiftUp=(circshift(LaplaceImg,[1 0]) .* LaplaceImg)<0;
    ShiftDown=(circshift(LaplaceImg,[-1 0]) .* LaplaceImg)<0;
    ShiftTopLeft=(circshift(LaplaceImg,[1 -1]) .* LaplaceImg)<0;
    ShiftTopRight=(circshift(LaplaceImg,[1 1]) .* LaplaceImg)<0;
    ShiftBottomLeft=(circshift(LaplaceImg,[-1 -1]) .* LaplaceImg)<0;
    ShiftBottomRight=(circshift(LaplaceImg,[-1 1]) .* LaplaceImg)<0;
    img = ShiftLeft+ShiftRight+ShiftUp+ShiftDown+ShiftTopLeft+ShiftTopRight+ShiftBottomLeft+ShiftBottomRight;
end