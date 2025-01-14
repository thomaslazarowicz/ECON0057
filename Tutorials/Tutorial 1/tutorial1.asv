%**************************************************************************
%                       Tutorial 1. Intro to Matlab                       %
%                 MSc Advanced Economic Theory (ECON0057)                 %
%                     Thomas Lazarowicz (UCL), Jan 2025                    %
%**************************************************************************


%--------------------------------------------------------------------------
% 1. What do we see when we open Matlab?  
%--------------------------------------------------------------------------
% - Current Folder (on the left): where your files are stored.
% - Editor (where you are reading this): to edit and run Matlab program files 
%   (.m files) such as scripts (like matlab_intro.m; equivalent to do-files 
%   in Stata) and functions.
% - Command Window (below the editor): alternative way to enter commands,
%   also displays output.
% - Workspace (on the right): displays your objects.

% To run code:
% - Either click "Run"  or 
% - Use keyboard shortcuts  (better!)

%--------------------------------------------------------------------------
% 2. Preparation commands 
%--------------------------------------------------------------------------
close all; clear; clc;  % close open images, clear workspace, clear command window

cd("C:\Users\uctpttl\ECON0057\Tutorials") % set working directory
addpath("Tutorial 1") % Additional folder 
addpath("Functions") %Directory for storing functions


%--------------------------------------------------------------------------
% 3. Entering data 
%--------------------------------------------------------------------------

% a. A scalar: 
k_ss = 4.54;
fprintf('Steady state capital: %2.2f\n',k_ss) % display in the command window

% b. A matrix:
C = eye(4); % 4x4 identity matrix using inbuilt function "eye"

% c. Long sum: 
manyCs = C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+... 
         C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C+C;  % ... divides the command in two lines

clear manyCs; % delete object from workspace 

% d. Matrices: 
A = [1 2 3; 4 5 6]; % 2x3 array  
A  % show in command window. Notice we don't add ";", why? 

B = zeros(5,2); % 5x2 matrix of zeros 
D = ones(12,1); % 12x1 vector of ones

% e. Vector of equidistant elements: 
v1 = (0:2:10);           % (start: increments: end)
v2 = linspace(0,10,6);   % start, end, number of values 

%--------------------------------------------------------------------------
% 4. Matrix operations 
%--------------------------------------------------------------------------
clear

A = randn(5);   % 5x5 matrix of pseudorandom standard normal values 
B = rand(5);    % 5x5 matrix of pseudorandom Uniform(0,1) value

C1 = A + B;     % Add A and B
C2 = A - B;     % Subtract B from A
C3 = A*B;       % Right-multiply A by B 
C4 = A/B;       % A*B^(-1), this is preferred to A*inv(B)
C5 = A\B;       % A^(-1)*B, this is preferred to inv(A)*B
C6 = A\eye(5);  % A^(-1), this is preferred to inv(A)

C7 = C1';       % Transpose of C1

C8  = A.*B;     % (i,j)th element of J is (a_ij*b_ij)
C9  = A./B;     % (i,j)th element of K is (a_ij/b_ij)
C10 = A.^3;     % (i,j)th element of L is (a_ij^3)

%--------------------------------------------------------------------------
% 5. Indexing 
%--------------------------------------------------------------------------
clear 

% a. 2x2 array
D          = zeros(5)    % 5x5 matrix of zeros
D(1,1)     = 1           % Replace element in first row and column with 1
D(5,:)     = 2           % Replace elements in fifth row with 2
D(2:3,4:5) = 3           % Replace elements in the 4th and fifth columns 
                         % of the second and third rows with 3
% b. Array of higher dimension 
E = zeros(5,5,3)         % 5x5x3 array of zeros (three 5x5 matrices)
E(:,:,1)   = D           % Replace the first matrix with D'
E(:,:,2)   = D           % Replace the second matrix with D

% c. Generate matrix of standard normal random variables and replace
% negative values with zero.
G          = randn(5)
G(G < 0)   = 0           % Turn every negative value in G into 0

%--------------------------------------------------------------------------
% 6. Inbuilt functions 
%--------------------------------------------------------------------------
[n,m] = size(D);   % Compute number of rows and columns of the matrix D
clear m
[n,~] = size(D)    % ~ tells Matlab not to do anything with the second output

F1    = sum(D,1);  % Sum along rows
F2    = sum(D,2);  % Sum along columns

%--------------------------------------------------------------------------
% 7. User-written functions 
%--------------------------------------------------------------------------

% a. Linear combination of two values with weight lambda
v1     = 5;
v2     = 10;
lambda = 0.5; % must be between 0 and 1 
x      = myfunction(v1, v2, lambda);

y = myfunction(v1, v2, 1.2);

% b. Anonymous function 
a =  1;
b = -3;
c =  2;
fx = @(x) exp(x) + a*x.^2 + b*x + c; % Define anonymous function fx(x)

fx1 = fx(1) % evaluate at x = 1 

% c. Symbolic math 
clear 

syms x;                 % symbolic variable 
fun = x^2;              % define a function
fprime = diff(fun);     % take derivative 
x0 = 2;                 % choose a point
y0 = subs(fprime,x,x0); % evaluate derivative at that point
y0 = double(y0);        % save as numeric value 

%--------------------------------------------------------------------------
% 8. Structures 
%--------------------------------------------------------------------------
% Sometimes we want our workspace to be better organized

clear 
parameters = struct();    % initialize structure 

parameters.alpha = 0.3; 
parameters.beta  = 0.9;
parameters.rho   = 0.85;

A = parameters.beta + parameters.alpha; % we use them as any other object 

%--------------------------------------------------------------------------
% 9. Loops 
%--------------------------------------------------------------------------
clear 

% a. Simple example
% I can repeat an operation may times manually: 
disp(1)
disp(0.8)
disp(0.6)

% Or I can automatize the process with a loop: 
for v = 1.0:-0.2:0.0 % Recall, start:gap:end
   disp(v)
end

% b. Putting loops to work

% Consider the neoclassical growth model. We wish to evaluate period 
% consumption of the representative household for each value of capital 
% today and capital choice tomorrow over some grid 
% for capital.

clear 

% Parameters 
pars.beta  = 0.99;
pars.eta   = 1.3;
pars.alpha = 0.35;
pars.delta = 0.1;

% Steady state values of consumption and capital (take as given for now)
css = 1.27;
kss = 5.93;

% Construct grids 
kgrid = linspace(0.1,2*kss,100)'; % Construct 100-point grid for capital
cgrid = zeros(100);               % Create storage matrix for consumption on grid


tic
for ii = 1:100 % For each value of capital today
    
    for jj = 1:100 % For each choice of capital tomorrow
        
        % Compute implied level of consumption using budget constraint
        % (just a formula for consumption!) 
        cgrid(ii,jj) = kgrid(ii)^pars.alpha + (1-pars.delta)*kgrid(ii) ...
                       - kgrid(jj);
        
    end
      
end
mytime = toc; 

%--------------------------------------------------------------------------
% 10. Plots 
%--------------------------------------------------------------------------
% Example: plot the curve y = x^2 - 3*x - 2 over the range [-20,20].
% plot and line take optional arguments such as 'color' and 'blue'.

clear 

x = -20:20;
y = x.^2 - 3*x - 2;

figure;
plot(x,y,'color','blue','LineWidth',2);
hold on;
line([-20 20],[0 0],'color','black','LineStyle','--');
title('x^2 - 3x - 2');
ylabel('y');
xlabel('x');

