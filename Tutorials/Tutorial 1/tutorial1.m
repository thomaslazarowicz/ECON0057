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

%--------------------------------------------------------------------------
% 11. Importing and Exporting Data 
%--------------------------------------------------------------------------
clear

% a. Import CSV data (e.g., GDP data)
gdp_data = readtable('gdp_quarterly.csv');  % assuming you provide this file
fprintf('First few rows of GDP data:\n')
head(gdp_data)


%--------------------------------------------------------------------------
% 12. Time Series Operations
%--------------------------------------------------------------------------
clear

% a. Create sample quarterly GDP data
dates = (2020:0.25:2024.75)'; % Quarterly dates from 2020Q1 to 2024Q4
gdp = 100 + 0.5*(1:length(dates))' + 2*randn(length(dates),1);

% b. Calculate growth rates
gdp_growth = 100 * (gdp(2:end) - gdp(1:end-1))./gdp(1:end-1);

% c. Moving average (4-quarter)
gdp_ma4 = movmean(gdp, 4);

% d. Plot time series
figure;
plot(dates, gdp, 'b-', 'LineWidth', 1.5)
hold on
plot(dates, gdp_ma4, 'r--', 'LineWidth', 1.5)
title('Quarterly GDP with 4-Quarter Moving Average')
ylabel('GDP Level')
xlabel('Year')
legend('GDP', '4-Quarter MA', 'Location', 'best')

% e. Calculate year-over-year growth rates
gdp_yoy = 100 * (gdp(5:end) - gdp(1:end-4))./gdp(1:end-4);
dates_yoy = dates(5:end);

% Additional plot for growth rates
figure;
plot(dates_yoy, gdp_yoy, 'b-', 'LineWidth', 1.5)
title('Year-over-Year GDP Growth')
ylabel('Percent Change')
xlabel('Year')
grid on

%--------------------------------------------------------------------------
% 13. Basic Simulation
%--------------------------------------------------------------------------
clear

% a. Simulate AR(1) process: y(t) = ρy(t-1) + ε(t)
T = 100;           % number of periods
rho = 0.8;         % persistence parameter
sigma = 0.1;       % standard deviation of shock
eps = sigma*randn(T,1);  % random shocks

% Initialize
y = zeros(T,1);
y(1) = eps(1);     % initial condition

% Generate series
for t = 2:T
    y(t) = rho*y(t-1) + eps(t);
end

% Plot with confidence bands
figure;
plot(y, 'b-', 'LineWidth', 1.5)
hold on
plot(1:T, 2*sigma*ones(T,1)/(1-rho^2)^0.5, 'r--')  % Upper confidence band
plot(1:T, -2*sigma*ones(T,1)/(1-rho^2)^0.5, 'r--') % Lower confidence band
title('Simulated AR(1) Process')
xlabel('Time')
ylabel('Value')
legend('Series', '±2SD Band', 'Location', 'best')

%--------------------------------------------------------------------------
% 14. Optimisation Example
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 14. Solow Model: Analytical and Numerical Solutions
%--------------------------------------------------------------------------
clear

% Parameters for Solow model
alpha = 0.3;    % capital share
delta = 0.1;    % depreciation rate
s = 0.2;        % savings rate
n = 0.02;       % population growth rate
g = 0.02;       % technology growth rate

% Analytical steady state solution
% In steady state: s*k^alpha = (delta + n + g)k
% Therefore: k* = (s/(delta + n + g))^(1/(1-alpha))
k_star_analytical = (s/(delta + n + g))^(1/(1-alpha));

% Numerical solution using fsolve (as a check)
ss_fn = @(k) s*k^alpha - (delta + n + g)*k;
k0 = 1;  % initial guess
options = optimset('Display', 'off');
[k_star_numerical, fval] = fsolve(ss_fn, k0, options);

% Display results
fprintf('Analytical steady state k*: %2.4f\n', k_star_analytical)
fprintf('Numerical steady state k*: %2.4f\n', k_star_numerical)
fprintf('Difference between methods: %2.4e\n', abs(k_star_analytical - k_star_numerical))

% Calculate steady state output and consumption
y_star = k_star_analytical^alpha;
c_star = (1-s)*y_star;
fprintf('Steady state output y*: %2.4f\n', y_star)
fprintf('Steady state consumption c*: %2.4f\n', c_star)

% Plot phase diagram
k_grid = linspace(0, 3.5, 100);
k_next = s*k_grid.^alpha./(delta + n + g);  % Evolution of capital

figure;
plot(k_grid, k_next, 'b-', 'LineWidth', 1.5)
hold on
plot(k_grid, k_grid, 'k--')  % 45-degree line
plot(k_star_analytical, k_star_analytical, 'ro', 'MarkerSize', 10)
title('Phase Diagram - Solow Model')
xlabel('k_t')
ylabel('k_{t+1}')
legend('Capital accumulation', '45° line', 'Steady state')
axis([0 3.5 0 3.5])  % Set consistent axis limits
grid on

% Additional plot: Production, Saving, and Depreciation
figure;
y = k_grid.^alpha;           % production function
sy = s*y;                    % saving
dk = (delta + n + g)*k_grid; % depreciation (including pop growth and tech)

plot(k_grid, sy, 'b-', 'LineWidth', 1.5)
hold on
plot(k_grid, dk, 'r--', 'LineWidth', 1.5)
plot(k_star_analytical, s*k_star_analytical^alpha, 'ko', 'MarkerSize', 10)
title('Solow Diagram - Saving and Depreciation')
xlabel('Capital (k)')
ylabel('Output (y)')
legend('sy', '(δ+n+g)k', 'Steady state')
grid on

%--------------------------------------------------------------------------
% 15. More Plotting
%--------------------------------------------------------------------------
clear

% Generate sample data
t = (1:40)';
gdp = 100 * exp(0.02*t + 0.1*randn(40,1));
unemp = 5 + 0.5*cos(t/4) + 0.2*randn(40,1);

% Create figure with two y-axes
figure;
yyaxis left
plot(t, gdp, 'b-', 'LineWidth', 1.5)
ylabel('GDP')

yyaxis right
plot(t, unemp, 'r-', 'LineWidth', 1.5)
ylabel('Unemployment Rate')

title('GDP and Unemployment Over Time')
xlabel('Quarters')
legend('GDP', 'Unemployment', 'Location', 'best')

%--------------------------------------------------------------------------
% 16. Error Handling and Data Checks
%--------------------------------------------------------------------------
clear

% Function with error handling for growth rate calculation
try
    % Sample data with a negative value
    values = [100; 105; -103; 108; 112];
    
    growth_rates = zeros(length(values)-1, 1);
    for i = 2:length(values)
        if values(i) <= 0 || values(i-1) <= 0
            warning('Negative or zero values found at period %d', i)
            growth_rates(i-1) = NaN;
        else
            growth_rates(i-1) = 100 * (values(i) - values(i-1))/values(i-1);
        end
    end
    
    % Display results
    fprintf('Growth rates calculated successfully\n')
    disp(growth_rates)
    
catch ME
    fprintf('Error in calculation: %s\n', ME.message)
end

%--------------------------------------------------------------------------
% 17. Saving and Loading Results
%--------------------------------------------------------------------------
clear

% Create some results
results = struct();
results.params = struct('alpha', 0.3, 'beta', 0.99, 'delta', 0.025);
results.simulation = randn(100,3);
results.dates = datetime(2020,1,1):calmonths(1):datetime(2024,4,1);

% Save to MAT file
save('macro_results.mat', 'results')

% Clear and reload
clear
load('macro_results.mat')
disp('Loaded parameters:')
disp(results.params)

