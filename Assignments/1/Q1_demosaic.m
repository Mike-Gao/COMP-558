im = imread("demosaicing.png");

row = size(im,1);
col = size(im,2);

%initialiating R, G and B channels
R = zeros(row, col);
G = zeros(row, col);
B = zeros(row, col);

% fill red
i = 2:2:row;
j = 2:2:col;
R(i,j) = im(i,j,1);
red = uint8(cat(3, R, zeros(row, col), zeros(row, col)));
subplot(1,5,1);
image(red);
title("red");


% fill green
i = 2:2:row;
j = 1:2:col;
G(i,j) = im(i,j,2);
i = 1:2:row;
j = 2:2:col;
G(i,j) = im(i,j,2);
green = uint8(cat(3, zeros(row, col), G, zeros(row, col)));
subplot(1,5,2);
image(green);
title("green");


% fill blue
i = 1:2:row;
j = 1:2:col;
B(i,j) = im(i,j,3);
blue = uint8(cat(3, zeros(row, col), zeros(row, col), B));
subplot(1,5,3);
image(blue);
title("blue");


R = conv2(R,[1 2 1;2 4 2;1 2 1]/4, 'same');
G = conv2(G,[0 1 0;1 4 1;0 1 0]/4, 'same');
B = conv2(B,[1 2 1;2 4 2;1 2 1]/4, 'same');


output = cat(3, R, G, B);
output=uint8(output);
subplot(1,5,4);
image(output);
title("output");
subplot(1,5,5);
image(im);
title("origional");



