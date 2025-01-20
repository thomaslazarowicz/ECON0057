%==========================================================================
%
%  Numerical Solution to McCall Search Model
%

%
%
%==========================================================================

clc;
clear;

cd("C:\Users\uctpttl\ECON0057\Tutorials") % set working directory
addpath("Tutorial 2") % Additional folder 
addpath("Functions") %Directory for storing functions

%% Parameters
beta = 0.9;                     % Discount factor: reflects the agent's patience (closer to 1 = more patient)
c = 2;                          % Unemployment benefits: fallback utility while searching for a job
w = (1:100)';                   % Wage offers (grid of possible wages)

% Define the probability distribution for wage offers
% Uncomment one of the following lines to choose the wage distribution:
% p = 0.01 * ones(100, 1);        % Uniform distribution (all wages equally likely)
p = w.^(-1) / sum(w.^(-1));      % Inverse distribution (higher wages less likely)

% Tolerance for convergence and maximum number of iterations
tol = 1.0e-5;                   % Convergence tolerance: how close successive iterations must be to stop
max_iter = 1000;                % Maximum number of iterations allowed

%% Plot the Wage Distribution
% Visualize the probability density function (PDF) of wage offers
figure(1);
plot(w, p, 'LineWidth', 1.5);
title('Probability Density of Wage Offers');
xlabel('Wage');
ylabel('Probability');
grid on;

%% Initialize Variables
% V0: Initial guess for the value function, starting at zero for all wages
% V1: Placeholder for the updated value function
V0 = zeros(size(w));            
V1 = V0;                        

tic;                            % Start timing the computation
disp('Starting Value Function Iteration...');

%% Value Function Iteration
% This loop iteratively computes the value function until convergence.
% The agent's goal is to determine the maximum value (expected utility)
% for each wage level, considering two choices:
% 1. Accept the current wage and get immediate utility.
% 2. Reject the current wage, receive unemployment benefits, and wait for
%    another wage offer (with future value discounted by beta).

for iter = 1:max_iter
    % The value function is updated using the Bellman equation:
    %   V1(w) = max{ Immediate payoff from working at wage w,
    %                Future payoff from waiting for a better offer }
    % Immediate payoff: w / (1 - beta) (value of accepting the wage offer)
    % Future payoff: c + beta * E[V(w')] (value of rejecting and waiting)
    V1 = max(w / (1 - beta), c + beta * (p' * V0));

    % Convergence check: stop if the maximum difference between iterations
    % is below the specified tolerance (tol)
    if norm(V1 - V0, inf) <= tol
        fprintf('Converged after %d iterations.\n', iter);
        break;
    end

    % Update the guess for the next iteration
    V0 = V1;

    % Display progress every 50 iterations to monitor convergence
    if mod(iter, 50) == 0
        fprintf('Iteration %d: Max Difference = %.6f\n', iter, norm(V1 - V0, inf));
    end
end

% If the loop completes without convergence, display a warning
if iter == max_iter
    disp('WARNING: No convergence after maximum iterations.');
end

% Stop timing and display elapsed time
time_elapsed = toc;
fprintf('Evaluation completed in %.2f seconds.\n', time_elapsed);

%% Analytical Reservation Wage
% The reservation wage is the lowest wage the agent is willing to accept.
% This is calculated analytically by solving the Bellman equation for
% equality between accepting the wage and rejecting it.
res_wage_analytical = (c + beta * (p' * V1)) * (1 - beta);
fprintf('Analytical Reservation Wage: %.2f\n', res_wage_analytical);

%% Plot Value Function
% Visualize the final value function, which shows the maximum expected
% utility for each possible wage level.
figure(2);
plot(w, V1, 'LineWidth', 1.5);
title('Value Function vs. Wage');
xlabel('Wage');
ylabel('Value Function');
grid on;

%% Graphical Reservation Wage
% Alternatively, the reservation wage can be identified graphically by
% finding the point where the value function "jumps" from the unemployment
% utility to the value of working.
for k = 2:(length(w) - 1)
    if V1(k + 1) > V1(k) && abs(V1(k) - V1(k - 1)) < tol
        fprintf('Graphical Reservation Wage: %.2f\n', w(k));
        break;
    end
end

%% Sensitivity Analysis Instructions
disp('-----------------------------------------------------');
disp('Sensitivity Analysis:');
disp('1. To test different probability distributions, modify "p" and rerun.');
disp('2. Observe how changes in "c" and "beta" affect the reservation wage:');
disp('   - As beta approaches 0, reservation wage approaches c.');
disp('   - As beta approaches 1, reservation wage approaches max(w).');
disp('-----------------------------------------------------');
