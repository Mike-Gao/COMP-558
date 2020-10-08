%% 
% Read in the image, convert to grayscale, and detect edges.
% Creates an array edges where each row is    (x, y, cos theta, sin theta)   

im = imread("vanishing_point_demo.png");
im = imresize(rgb2gray(im), 0.5);

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

im_new = zeros(size(im));
horizontal_edges = edges(abs(edges(:,3))<(sqrt(2)/2),:);
vertical_edges = edges(abs(edges(:,3))>=(sqrt(2)/2),:);


for i=1:length(horizontal_edges)
    r = horizontal_edges(i,1) * horizontal_edges(i,3) + horizontal_edges(i,2) * horizontal_edges(i,4);
    for j=1:size(im_new,2)
        y = round((r-(horizontal_edges(i,3)*j)) / horizontal_edges(i,4));
        if y >= 1 && y <= size(im_new,1)
            im_new(y,j) = im_new(y,j) + 1;
        end
    end
end

for i=1:length(vertical_edges)
     r = vertical_edges(i,1) * vertical_edges(i,3) + vertical_edges(i,2) * vertical_edges(i,4);
     for j=1:size(im_new,1)
         x = round((r-(vertical_edges(i,4)*j)) / vertical_edges(i,3));
         if x >= 1 && x <= size(im_new,2)
             im_new(j,x) = im_new(j,x) + 1;
         end
     end
end


im_new = uint8(im_new);
imagesc(im_new);
colorbar;