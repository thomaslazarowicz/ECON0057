%**************************************************************************
%    Tutorial 4. Example of the contraction mapping theorem at work       %
%                 MSc Advanced Economic Theory (ECON0057)                 %
%**************************************************************************

%% 1. Initialization
clear all; close all; clc;

%% 2. Parameters
% Economic Parameters
params.alpha = [2; 1];     % State contingent utilities 
params.beta  = 0.9;        % Discount factor
params.p1    = 0.9;        % Probability of staying in state 1 
params.p2    = 0.5;        % Probability of staying in state 2

% Numerical Parameters
params.maxit = 10000;      % Maximum number of iterations 
params.tol   = 1.0e-5;     % Tolerance
params.print = true;       % Plot option

%% 3. Setup Markov Transition Matrix
% Initialize transition matrix
P = zeros(2, 2);

% Fill in transition probabilities
P(1,1) = params.p1;                
P(2,2) = params.p2;                    
P(1,2) = 1 - params.p1;     % Rows sum up to 1  
P(2,1) = 1 - params.p2;           

%% 4. Direct Solution
% Solve system directly using matrix algebra
x_direct = (eye(2) - params.beta * P) \ params.alpha;

% Display results
fprintf('\nSolution by direct calculation:\n')
disp(x_direct)

%% 5. Iterative Solution
% Initialize
x0 = zeros(2, 1);
x_history = x0;  % Store sequence for visualization

% Start timer
tic;

% Iterate until convergence
for it = 1:params.maxit
    % Calculate next iteration
    x1 = params.alpha + params.beta * P * x0;
    x_history = [x_history, x1];
    
    % Check convergence
    dist = norm(x1 - x0);
    if dist < params.tol
        fprintf('Convergence achieved after %d iterations\n', it)
        break
    elseif it == params.maxit
        warning('Maximum iterations reached without convergence')
    end
    
    % Update guess
    x0 = x1;
end

% End timer
computation_time = toc;

%% 6. Results Display
% Print summary statistics
fprintf('\nSolution by successive approximation:\n')
disp(x1)
fprintf('Computation time: %.2f seconds\n', computation_time)

% Display sequence evolution
fprintf('\nFirst 8 iterations:\n')
disp(x_history(:, 2:min(9, it)))

fprintf('\nLast 8 iterations:\n')
disp(x_history(:, max(1, it-7):it))

%% 7. Visualization
if params.print
    % Create figure with improved styling
    figure('Name', 'Contraction Mapping Convergence', 'Position', [100 100 800 500])
    
    % Plot convergence path
    it_seq = 0:it;
    plot(it_seq, x_history, '.', 'MarkerSize', 10)
    
    % Add labels and title
    xlabel('Number of Iterations', 'FontSize', 12)
    ylabel('Value', 'FontSize', 12)
    title('Convergence of Sequence {x_n} by Successive Approximation', 'FontSize', 14)
    
    % Add legend and grid
    legend({'x_n(1)', 'x_n(2)'}, 'Location', 'southeast', 'FontSize', 12)
    grid on
    
    % Adjust axes for better visualization
    ylim([min(x_history(:))*0.9, max(x_history(:))*1.1])
end