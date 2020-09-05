%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Introducing Matlab (adapted from http://www.cns.nyu.edu/~eero and
% http://www.cs.dartmouth.edu/~farid/teaching/cs88/matlab.intro.html)
% via http://www-cse.ucsd.edu/%7Esjb/classes/matlab/matlab.intro.html
% Sourced from Morteza Rezanejad, McGill University (2017)
% Modified by Tabish Syed, McGill University (2019)
% Modifed  by Mike Langer, Sept. 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Help and basics

% The symbol "%" is used in front of a comment.

% To get quick help inside command window type "help" (will give list of help topics) or "help topic"
% To get detailed help/ documentation on any command type doc command_name

% If you don't know the exact name of the topic or command you are looking for,
% type "lookfor keyword" (e.g., "lookfor regression")

% When writing a long matlab statement that exceeds a single row use ...
% to continue statement to next row.

% When using the command line, a ";" at the end means matlab will not
% display the result. If ";" is omitted then matlab will display result.

% Use the up-arrow to recall commands without retyping them (and down
% arrow to go forward in commands).  

% Using up(down)-arrow with a prefix recalls commands that start with 
% that prefix.

% Pressing tab key after partial command completes the commands or lists 
% available options if multiple completions are possible.

% Other commands borrowed from emacs and/or tcsh:
% C-a moves to beginning of line (C-e for end), C-f moves forward a
% character (C-b moves back), C-d deletes a character, C-k deletes 
% the line to the right of the cursor, C-p goes back through the
% command history and C-n goes forward (equivalent to up and down arrows),
% tab command completion.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Objects in matlab -- the basic objects in matlab are scalars,
% vectors, and matrices...

echo off;    %  If you want to echo the commands in the command line, then turn it on.
            %  This is arguably useful for demos.

N	= 5				% a scalar
v 	= [1 0 0]			% a row vector
v 	= [1;2;3]			% a column vector
v 	= v'				% transpose a vector 
						%(row to column or column to row)
v	= [1:.5:3]			% a vector in a specified range: 
v	= pi*[-4:4]/4			% 	[start:end] or [start:stepsize:end]
v	= []				% empty vector

m 	= [1 2 3; 4 5 6]		% a matrix: 1ST parameter is ROWS
					% 	    2ND parameter is COLS 
m	= zeros(2,3)   			% a matrix of zeros
v	= ones(1,3)  			% a matrix of ones
m	= eye(3)			% identity matrix
% v	= rand(3,1)			% random matrix with values in [0,1] (see also randn)

pause    %  click in command window and hit any key to continue
          %  ML:  I don't use it much.   I find breaking into sections to
          %  be usually more convenient.

%%   Two ampersands defines a section.   You can "run section" to step through chunks of code.

%Save data to an ASCII file, and view the contents of the file with the type function:

p = rand(1, 10);
q = ones(10);
save('pqfile.txt', 'p', 'q', '-ASCII')
%display contents of file
type pqfile.txt

%Alternatively, use command syntax for the save operation:
save pqfile.txt p q -ASCII
%You can also store variables in mat files
% create a file 'matrix_data' containing:
matrix_data =  [2,     3,     4;
			    5,     6,     7;
				1,     2,     3]

save('matrix_data.mat','matrix_data')   %  does same thing as "save matrix_datamat matrix_data"

clear matrix_data           % clears variables
load matrix_data 			% read data from a file  matrix_data.mat

disp(matrix_data)           % displays values  (prettier than just typing "matrix_data")

%%%%  Indexing vectors and matrices
%%
v	= [1 2 3];			% access a vector element
v(3)					%	vector(number) 
   				  	    % Index starts from 1, not 0.

m 	= [1 2 3; 4 5 6]
m(1,3)					% access a matrix element
					% 	matrix(rownumber, columnnumber)
m(2,:)    				% access a matrix row (2nd row)
m(:,1)    				% access a matrix column (1st row)

size(m)					% size of a matrix
size(m,1)  				% number rows
size(m,2)  				% number of columns

m1	= zeros(size(m))		% create a new matrix with size of m

who					    % list of variables
whos					% list/size/type of variables

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Simple operations on vectors and matrices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (A) Pointwise (element by element) Operations:

% addition of vectors/matrices and multiplication by a scalar
% are done "element by element"
a	= [1 2 3 4];			% vector
2 * a 					% scalar multiplication
a / 4					% scalar multiplication
b	= [5 6 7 8];			% vector
a + b					% pointwise vector addition
a - b					% pointwise vector addition
a .^ 2					% pointise vector squaring (note .)
a .* b					% pointwise vector multiply (note .)
a ./ b					% pointwise vector divide (note .)

log( [1 2 3 4] )			% pointwise arithmetic operation
round( [1.5 2; 2.2 3.1] )		% pointwise arithmetic operation

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (B) Vector Operations (no for loops needed -- this is huge in Matlab)
% Built-in matlab functions operate on vectors.    If a (2D) matrix is given,
% then the function operates on each column of the matrix.

a	= [1 4 6 3]			% vector
sum(a)					% sum of vector elements
mean(a)					% mean of vector elements
var(a)					% variance
std(a)					% standard deviation
max(a)					% maximum

a 	= [1 2 3; 4 5 6]		% matrix
a(:)                 		% vectorized version of the matrix
mean(a)                     % mean of each column
max(a)                      % max of each column

max(max(a))		     		% to obtain max of matrix 
max(a(:))		     		% 	or...
max(a,'all')                %   or
%%
%%%%%%%%%%%%%%%%%%%%%%%%
% (C) Matrix Operations:

[1 2 3] * [4 5 6]'  			% row vector 1x3 times column vector 3x1 
                    			% results in single number, also
                   			% known as dot product or inner product

[1 2 3]' * [4 5 6]  			% column vector 3x1 times row vector 1x3
                    			% results in 3x3 matrix, also
                    			% known as outer product

a	= rand(3,2)			% 3x2 matrix
b	= rand(2,4)			% 2x4 matrix
c	= a * b				% 3x4 matrix

a	= [1 2; 3 4; 5 6]		% 3 x 2 matrix
b	= [5 6 7];			% 1 x 3 vector
b * a					% matrix multiply
a' * b'					% matrix multiply

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(4) Saving your work

save mysession      			% creates mysession.mat with all variables
save mysession a b  			% save only variables a and b

clear all				% clear all variables
clear a b           			% clear variables a and b

load mysession				% load session
a
b

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(5) Relations and control statements

% Example:  Given a vector v, create a new vector with values equal to
% v if they are greater than 0, and equal to 0 if they less than or
% equal to 0.

%  You would naturally do this with a "for loop" as follows.   (But there
%  is another way to do it Matlab which we'll see after.)

v	= [3 5 -2 5 -1 0]		% 1: FOR LOOPS
u 	= zeros( size(v) );		% initialize
for i = 1:size(v,2)			% size(v,2) is the number of columns
	if( v(i) ~= 0 )
		u(i) = v(i); 
	end
end
u

%  We can do this without for loops as follows:

u2	= zeros( size(v) );		% initialize
ind	= find( v>0 )			% index into >0 elements 
u2(ind)	= v( ind )			
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(6) Creating functions using m-files:
% Functions in matlab are written in m-files. Create a file called 
% 'thres.m'.     In this file put the following five lines (without comments):
% (Alternatively, you can define a function at the end of this file, as
% I've done for you.)

v	= [3 5 -2 5 -1 0]
res = thres( v )				% call from command line

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(7) Plotting 

x	= [0 1 2 3 4];			% basic plotting
plot( x );                  %  note that a plot window has been created
plot( x, 2*x );
axis( [0 8 0 8] );

x = pi*[-24:24]/24;
plot( x, sin(x) );
xlabel( 'radians' );
ylabel( 'sin value' );
title( 'sine curve' );
gtext( 'put cursor where you want text and press mouse' );

figure;					% multiple functions in separate graphs within one plotting window 
subplot( 1,2,1 );
plot( x, sin(x) );
axis square;
subplot( 1,2,2 );
plot( x, 2.*cos(x) );
axis square;

figure;					% multiple functions in single graph
plot( x,sin(x) );
hold on;              			% hold on tells matlab to write on top 
plot (x, 2.*cos(x), '--' );		% of the current plot
legend( 'sin', 'cos' );
hold off;

%%
figure;					% matrices as images
m = randn(64,64);       % 64 x 64 matrix of random values form normal distribution
imagesc(m)
colormap autumn;       % check out "doc colormap" to see the colormap that "autumn" corresponds to
axis image
axis off;

%%

function u = thres( v )
   u	= zeros( size(v) );		% initialize
   ind	= find( v>0 )			% index into v>0 elements 
   u(ind)	= v( ind )			
end

