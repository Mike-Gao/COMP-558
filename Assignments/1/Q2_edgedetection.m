rect = ones(500,500) * 256;
for i=1:150
    p1=[randi([1,500],1),randi([1,500],1)];
    height = randi(150,1);
    p2=[p1+height,p1(2)];
    width = randi(150,1);
    p3=[p1(1),p1(2)+width];
    intensity=randi([0,250],1);
    smallrec=intensity*ones(height,width);
    rect(p1(1):p2(1)-1,p1(2):p3(2)-1) = smallrec;
end
rect = uint8(rect(125:375,125:375));
figure(1);
imshow(rect);

Gaussian = fspecial('gaussian',20,3);
x=conv2(Gaussian, [1 -2 1], 'same');
y=conv2(Gaussian, [1 -2 1]', 'same');
Laplace = x + y;
figure(2);
imshow(Laplace, 'DisplayRange', [-0.003 0.001], 'InitialMagnification', 1000);
colorbar


figure(3);
ZeroCrossingImg = CreateZeroCrossingImage(rect,20,3);
imshow(ZeroCrossingImg);

figure(4);
NoisyRect = rect+uint8(0.1*std(double(rect(:)))*randn(size(rect)));
imshow(NoisyRect);

figure(5)
NoisyZeroCrossingImg = CreateZeroCrossingImage(NoisyRect,20,3);
imshow(NoisyZeroCrossingImg);

figure(6)
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