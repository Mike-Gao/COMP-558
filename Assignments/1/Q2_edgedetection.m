%% 
% Read in the image, convert to grayscale, and detect edges.
% Creates an array edges where each row is    (x, y, cos theta, sin theta)   

%im = imread("...");
%im = imresize(rgb2gray(im), 0.5);

Iedges = edge(im,'canny');
%  imgradient computes the gradient magnitude and gradient direction
%  using the Sobel filter.  
[~,grad_dir]=imgradient(im);
%  imgradient defines the gradient direction as an angle in degrees
%  with a positive angle going (CCW) from positive x axis toward
%  negative y axis.   However, the (cos theta, sin theta) formulas from the lectures define theta
%  positive to mean from positive x axis to positive y axis.  For this
%  reason,  I will flip the grad_dir variable:
grad_dir = - grad_dir;

imshow(Iedges)

%Now find all the edge locations, and add their orientations (cos theta,sin theta). 
%  row, col is  y,x
[row, col] = find(Iedges);
% Each edge is a 4-tuple:   (x, y, cos theta, sin theta)   
edges = [col, row, zeros(length(row),1), zeros(length(row),1) ];
for k = 1:length(row)
     edges(k,3) = cos(grad_dir(row(k),col(k))/180.0*pi);
     edges(k,4) = sin(grad_dir(row(k),col(k))/180.0*pi);
end