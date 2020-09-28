im = imread("demo.jpg")

Rcfa = repmat([1 0;0 0], size(im)/2);
Gcfa = repmat([0 1;1 0], size(im)/2);
Bcfa = repmat([0 0;0 1], size(im)/2);

Rh = im.*Rcfa
Gh = im.*Gcfa
Bh = im.*Bcfa

R = conv2(Rh,[1 2 1;2 4 2;1 2 1]/4, 'same');
G = conv2(Gh,[0 1 0;1 4 1;0 1 0]/4, 'same');
B = conv2(Bh,[1 2 1;2 4 2;1 2 1]/4, 'same');

output(:,:,1)=R; output(:,:,2)=G; output(:,:,3)=B;