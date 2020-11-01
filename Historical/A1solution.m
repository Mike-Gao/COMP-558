%  COMP 558  Fall 2010
%  Assignment 1 Solution
%
%  Prepared by Prof. Michael Langer
%  with help from T.A.s Florian and Chatavut (Patrick)

clear

%  Input dialogs use "cell arrays".  Index a cell with [] and index
%  the contents of a cell with {}.

prompt = 'Enter image name ';
userInput = inputdlg(prompt);
imageName = ['./' userInput{1} '.jpg'];
clear userInput

%  Read image and exif file.

I = imread(imageName);
imageInfo = imfinfo(imageName);
NY = imageInfo.DigitalCamera.CPixelYDimension;
NX = imageInfo.DigitalCamera.CPixelXDimension;
f = imageInfo.DigitalCamera.FocalLength;

pixelsPerMM = 63.5;   %  For the Nikon D70, this is the density of pixels per mm.

close all
figure;
imshow(I);
hold on;

%  User must select two pairs of points (four points) which are the 
%   corners of a scene rectangle.  These four points will be used to 
%   compute vanishing points.  The selected rectangle will define an 
%   orthogonal world coordinate system.  The points must be selected in 
%   either clockwise or counterclockwise order.  The first two points 
%   should be a known distance apart in the scene.
%   Read the first pair of points, draw the line, and then ask how far
%   apart they are in 3D.

displacement = 25;  %  label the selected points
[x1, y1, x2, y2] = readPairXY();
text(x1,y1+displacement,'1');
text(x2,y2+displacement,'2');

hold on;
plot([x1,x2], [y1, y2], 'bs');
[a1, b1, c1] = computeLine(x1, y1, x2, y2);
[xc1, xc2, yc1, yc2] = clipLineToImage(a1,b1,c1,NX,NY);
plot([xc1,xc2], [yc1, yc2], 'gs-')

prompt = 'How far apart were those points in 3D (in millimeters) ?';
userInput = inputdlg(prompt);
distP1P2MM = str2double(userInput{1});
clear userInput

%  Read the second pair of points and draw the line.

[x3, y3, x4, y4] = readPairXY();
text(x3,y3+displacement,'3');
text(x4,y4+displacement,'4');
plot([x3,x4], [y3, y4], 'bs');
[a2, b2, c2] = computeLine(x3, y3, x4, y4);
[xc1, xc2, yc1, yc2] = clipLineToImage(a2,b2,c2,NX,NY);
plot([xc1,xc2], [yc1, yc2], 'gs-')

%display(['(x1,y1) = (' num2str(round(x1)) ',' num2str(round(y1)) ') '  ]);
%display(['(x2,y2) = (' num2str(round(x2)) ',' num2str(round(y2)) ')'  ]);
%display(['(x3,y3) = (' num2str(round(x3)) ',' num2str(round(y3)) ')'  ]);
%display(['(x4,y4) = (' num2str(round(x4)) ',' num2str(round(y4)) ')' ]);

%  Compute first vanishing point in pixel coordinates

vp1 = computeVP(a1, b1, c1, a2, b2, c2);
pix = round( vp1');
%display(['First vanishing point is at pixel (x,y) = ( ', num2str(pix(1)), ' , ', num2str(pix(2)), ' ).']);
plot(vp1(1),vp1(2),'gs')

% Next, pair (2nd,3rd) and (1st,4th) points, and compute the second
% vanishing point in pixel coordinates.

[a1, b1, c1] = computeLine(x1, y1, x4, y4);
[xc1, xc2, yc1, yc2] = clipLineToImage(a1,b1,c1,NX,NY);
plot([xc1,xc2], [yc1, yc2], 'gs-');
[a2, b2, c2] = computeLine(x2, y2, x3, y3);
[xc1, xc2, yc1, yc2] = clipLineToImage(a2,b2,c2,NX,NY);
plot([xc1,xc2], [yc1, yc2], 'gs-');

vp2 = computeVP(a1, b1, c1, a2, b2, c2);
pix = round( vp2');
%display(['Second vanishing point is at pixel (x,y) = ( ', num2str(pix(1)), ', ', num2str(pix(2)), ' ).']);
plot(vp2(1),vp2(2),'rs')

%    Draw vanishing line.

[a, b, c] = computeLine(vp1(1), vp1(2), vp2(1), vp2(2) );
  display(['Vanishing line (in pixel space) is ', num2str(a), ' x + ', num2str(b), ' y + ', ...
         num2str(c), ' =  0']);
[xc1, xc2, yc1, yc2] = clipLineToImage(a,b,c,NX,NY);
plot([xc1,xc2], [yc1, yc2], 'b-')  %  changed to blue on Sept 29

%%   Question 1:  Estimate equation of vanishing line on Z=f plane.   
%
%  (vp(1), vp(2), 1) is  homogeneous coordinates
%  vp1 and vp2 are measured in pixel units,  so convert to mm

%  BEGIN SOLUTION
vp1ProjPlane =  [pixelsPerMM,  0,   NX/2; 0,  pixelsPerMM, NY/2; 0, 0, 1 ] \ [vp1(1) vp1(2) 1]' ;
vp2ProjPlane =  [pixelsPerMM,  0,   NX/2; 0,  pixelsPerMM, NY/2; 0, 0, 1 ] \ [vp2(1) vp2(2) 1]' ;
%  END SOLUTION

display(['1.)  First  vanishing point is at ( ', num2str(vp1ProjPlane(1)), ', ', num2str(vp1ProjPlane(2)), ' ).']);
display(['     Second vanishing point is at ( ', num2str(vp2ProjPlane(1)), ', ', num2str(vp2ProjPlane(2)), ' ).']);
[a, b, c] = computeLine(vp1ProjPlane(1), vp1ProjPlane(2), vp2ProjPlane(1), vp2ProjPlane(2) );
%display(['  Vanishing line (on proj plane Z=f) is ', num2str(a), ' x + ', num2str(b), ' y + ', ...
%         num2str(c), ' =  0.']);
     
%%   Question 2:  Estimate the slant and tilt of the scene plane, AX + BY + CZ = D.        
%
%  Hint:  the A,B,C values are the same as those of the plane that passes
%  through the origin and through the vanishing line in the Z=f plane.

%  ---------  BEGIN SOLUTION   -----------------------
%  The vanishing line on Z=f plane is ax + by + c = 0
%
%  Consider a plane, AX + BY + CZ = 0, that passes through 
%  the origin and has the same normal as the scene plane.   
%  Since we know the vanishing line as ax + by + c = 0, we
%  can associate a=A, b=B, c=Cf for the vanishing line.
%  Thus, the plane that passes through the origin and the 
%  vanishing line is  aX + bY + c/f * Z = 0.
%  We can compute the slant and tilt using definition from class.
%

A=a; B=b; C= c/f;
tilt  = angle( - A/C - 1i * B/C );
slant = atan2(  sqrt(A*A+B*B), abs(C) );

%  An alternative solution would be to recall from basic linear 
%  algebra that the normal of a plane AX + BY + CZ = 0 is (A,B,C).
%  One could compute the normal by taking the cross product of the 
%  two vanishing points.


%  ---------   END SOLUTION ---------------------------

display(['2.)  Slant is ' num2str(round(slant/pi*180)) '  degrees.']);
display(['     Tilt is ' num2str(round(tilt/pi*180)) ' degrees' ...
    ' (measured clockwise from horizontal - see comment)']);

%  Tilt is measured clockwise in this program because Matlab indexes image
%  pixels with origin at topleft

%%  Question 3:   Estimate the distance to plane, namely where the optical
%                 axis intersects the plane.  
%
%  Take the first two selected image points (x1,y1), (x2,y2).  
%  The user inputs the Euclidean 3D distance between the scene points 
%  that project to these image points.  
%  The plane is AX + BY + CZ = D.  Thus, D/C is the Z value where the  
%  plane intersects the optical axis (X=0,Y=0).  
%  We solve for D.
%
%  Hint:  
%  - Convert the selected image points from pixels to mm, and create new 
%    2-vectors p1 and p2 that represent these points on a projection plane.
%    (You can use Z=1, for example.) 
%    The rays from the origin through these vectors will intersect the
%    plane at two points that are the given distance apart.

%  -----------  BEGIN SOLUTION -----------------------------
%   Estimate the equation of the plane,  AX+BY+CZ=D.  Project onto the Z=1
%   plane by dividing by Z,  e.g. A p1(1) + B p1(2) + C  = D/Z,  or 
%
%   Z = D / ( A p1(1) + B p1(2) + C )                (*)
%
%   We want to find D and Z(x,y).
%  Take the first two selected image points (x1,y1), (x2,y2).  They are in
%  pixel coords and Z=f plane. First convert to mm.

p1 =  [pixelsPerMM,  0,   NX/2; 0,  pixelsPerMM, NY/2; 0, 0, 1 ] \ [x1 y1 1]'; 
p2 =  [pixelsPerMM,  0,   NX/2; 0,  pixelsPerMM, NY/2; 0, 0, 1 ] \ [x2 y2 1]'; 

%  Then rescale to find the xy coords on the Z=1 plane.
p1 =  [1/f,  0,  0; 0,  1/f, 0; 0, 0, 1 ] * p1; 
p2 =  [1/f,  0,  0; 0,  1/f, 0; 0, 0, 1 ] * p2; 

%  The user has input the Euclidean distance between the backprojected 
%  points in the scene as distP1P2MM.   Thus,
%
%   | p1 * Z1  - p2 * Z2 | = distP1P2MM              (**)
%
%   From (*), we get Z1 = D/(A p1(1) + B p1(2) + C), Z2 = D/(A p2(1) + B p2(2) + C)  
%   so we can plug those into (**) and solve for | D |.
%
absD = distP1P2MM   / norm( p1 / (A*p1(1) + B*p1(2) + C)  - p2 / (A *p2(1) + B *p2(2) + C) );

%   One subtlety is that we have only solved for the absolute value of D.
%   The plane could intersect the optical axis either in front or behind 
%   the observer.  e.g.  if we have a ground plane and the observer is
%   looking upwards rather than downwards, then the intersection with the
%   ground will occur at Z < 0.
%
%   To check if Z0 is positive or negative, we just need to check whether
%   the origin and one of the selected points are on the same side of
%   the vanishing line or not.  (If yes, then Z>0.  If no, then Z<0.)

Zpos = (sign(a* f*p1(1) + b* f*p1(2) + c) == sign(c));

if (Zpos)        %  then D has the same sign as C
    D = absD * sign(C) ;
else             %   D has the opposite sign of C.
    D = - absD * sign(C); 
end
Z0 = D/C;

% -----------  END SOLUTION ---------------------------------

display(['3.)  Distance to plane is ' num2str(round(Z0)) ' millimetres.']); 

%%   Question 4:   Estimate and draw lines where blurwidth is 0, 1.
%
%  The thin lens equation implies: zi = f * Z0 / (Z0 - f)
%  But  zsensor = zi  at the center of the image,  so zsensor is easy to
%  get if you know the depth Z0.

zsensor = f * Z0 / (Z0 - f);
aperture = f/imageInfo.DigitalCamera.FNumber;

%  Let blurwidth  be Delta X_i  in the lecture notes.

%  -----------BEGIN SOLUTION ------------------------------------
%
%  Plane in scene is:  AX + BY + CZ = D
%  Project onto the Z=zsensor plane  by multiplying by zsensor/Z.
%  (Ax + By + C zsensor)/D  = zsensor/Z, ...   
%
%  blurwidth = aperture * (zsensor / zi - 1);  
%            = aperture * ( zsensor * (1/f - 1/Z) - 1 );
%            = aperture * ( zsensor * (1/f - (Ax + By + C zsensor)/(zsensor D)) - 1 );
%
%  so   grad blurwidth = - (aperture /D) * (A, B)  
%  But the change in x with a unit change in blurwidth is |grad blurwidth|^{-1}
%
%  Notice that this is a ratio, so the units of "change in x" just need
%  to be the same as the units of "change in blurwidth".  The units can be
%  mm or pixels or whatever.   Since we want to draw it in pixels and
%  measure blur in pixels, "pixels" is the unit we use.   So for a change
%  in blur width of 1 pixel, we want a blurwidthStep as follows:

blurwidthStep = D/( aperture * sqrt(A*A + B*B) );

%  Draw a line through the origin and parallel to the vanishing line.

[aTmp, bTmp, cTmp] = computeLine( NX/2, NY/2, NX/2 - b, NY/2 + a);
[xc1, xc2, yc1, yc2] = clipLineToImage(aTmp,bTmp,cTmp,NX,NY);
plot([xc1,xc2], [yc1, yc2], 'r-')

% draw another line that is stepped by blurwidthStep in a direction
% perpendicular to the line through the origin.
%
% step in direction (a,b)/sqrt(a*a+b*b
%

stepX = blurwidthStep*a/sqrt(a*a+b*b);
stepY = blurwidthStep*b/sqrt(a*a+b*b);
display(['4.)  Distance from image center to blurwidth=1 line is ' ...
    num2str(round(sqrt(stepX*stepX + stepY*stepY))) ' pixels.']); 

[aTmp, bTmp, cTmp] = computeLine( NX/2 + stepX, NY/2+ stepY, NX/2 - b + stepX, NY/2 + a + stepY);
[xc1, xc2, yc1, yc2] = clipLineToImage(aTmp,bTmp,cTmp,NX,NY);
plot([xc1,xc2], [yc1, yc2], 'r-')
[aTmp, bTmp, cTmp] = computeLine( NX/2 - stepX, NY/2- stepY, NX/2 - b - stepX, NY/2 + a - stepY);
[xc1, xc2, yc1, yc2] = clipLineToImage(aTmp,bTmp,cTmp,NX,NY);
plot([xc1,xc2], [yc1, yc2], 'r-')

% ------------  END SOLUTION -------------------------------------