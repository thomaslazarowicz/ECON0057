%**************************************************************************
%               Tutorial 4. McCall Search Model with Layoffs              %
%                 MSc Advanced Economic Theory (ECON0057)                 %
%                       Problem set 3. Question 1.1.                      %
%**************************************************************************
%  
% Numerical Solution to McCall Search Model with Layoffs
%
% The dynamic problem:
%
% V(w) = max{ w + beta [(1-alpha)V(w) + alpha*(c+beta*E[V(w')]) ], 
%             c + beta*E[V(w')] }

%% 1. Initialization
clear all; close all; clc;

%% 2. Parameters
% Economic Parameters
params.beta  = 0.9;      % Discount factor
params.c     = 2;        % Unemployment benefits
params.alpha = 0.1;      % Probability of being fired

% Grid Parameters
params.n     = 100;      % Length of the distribution
params.wmin  = 1;        % Minimum wage
params.wmax  = 100;      % Maximum wage

% Numerical Parameters
params.maxit = 1000;     % Maximum iterations
params.tol   = 1.0e-5;   % Convergence tolerance
params.print = true;     % Plot option

%% 3. Grid Construction
% Create wage grid and probability distribution
w_grid = linspace(params.wmin, params.wmax, params.n)';

% Default uniform PDF for wage offers
p = ones(params.n, 1) / params.n;

% Alternative PDF (commented out)
% p = w_grid.^(-1);
% p = p / sum(p);

%% 4. Value Function Iteration
% Initialize
V0 = zeros(params.n, 1);
tic;

% Main iteration loop
for it = 1:params.maxit
    % Expected value of future wages
    EV = p' * V0;
    
    % Value of accepting vs rejecting offer
    V_accept = w_grid + params.beta * ((1-params.alpha) * V0 + ...
               params.alpha * (params.c + params.beta * EV));
    V_reject = params.c + params.beta * EV;
    
    % Compute value function
    VF = max(V_accept, V_reject);
    
    % Check convergence
    dist = norm(VF - V0);
    if dist < params.tol
        fprintf('Convergence achieved after %d iterations\n', it)
        break
    elseif it == params.maxit
        warning('Maximum iterations reached without convergence')
    end
    
    V0 = VF;
end

computation_time = toc;

%% 5. Results Analysis
% Display summary statistics
fprintf('\nValue Function Iteration Results:\n')
fprintf('--------------------------------\n')
fprintf('Total iterations: %d\n', it)
fprintf('Computation time: %.4f seconds\n', computation_time)
fprintf('Value function range: [%.2f, %.2f]\n', VF(1), VF(end))

%% 6. Reservation Wage Calculation
% Analytical reservation wage
res_wage_analytical = (params.c + params.beta * (p' * VF)) * (1 - params.beta);
fprintf('\nAnalytical Results:\n')
fprintf('--------------------------------\n')
fprintf('Reservation wage: %.2f\n', res_wage_analytical)

% Numerical reservation wage (first wage where value function increases)
VF_diff = diff(VF);
res_wage_idx = find(VF_diff > 0, 1);
res_wage_numerical = w_grid(res_wage_idx);

fprintf('\nNumerical Results:\n')
fprintf('--------------------------------\n')
fprintf('Reservation wage: %.2f\n', res_wage_numerical)

%% 7. Visualization
if params.print
    figure('Name', 'McCall Search Model Results', 'Position', [100 100 800 500])
    
    % Plot value function
    plot(w_grid, VF, 'LineWidth', 2)
    hold on
    
    % Add reservation wage vertical line
    plot([res_wage_numerical res_wage_numerical], [min(VF) max(VF)], ...
         '--r', 'LineWidth', 1.5)
    
    % Formatting
    xlabel('Wage (w)', 'FontSize', 12)
    ylabel('Value Function V(w)', 'FontSize', 12)
    title('Value Function and Reservation Wage', 'FontSize', 14)
    legend({'Value Function', 'Reservation Wage'}, 'Location', 'southeast', ...
           'FontSize', 12)
    grid on
    
    % Add text annotation for reservation wage
    text(res_wage_numerical + 2, min(VF) + 0.1 * (max(VF) - min(VF)), ...
         sprintf('w_r = %.2f', res_wage_numerical), 'FontSize', 12)
end