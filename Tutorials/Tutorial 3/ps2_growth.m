%**************************************************************************
%                  Tutorial 3. Neoclassical Growth Model                  %
%                 MSc Advanced Economic Theory (ECON0057)                 %
%                             Problem set 2                               %
%**************************************************************************
%  
% Numerical Solution to the Neoclassical Growth Model 
%
% The dynamic problem:
% 
%                V(k) = max U[f(k) - k'] + beta*V(k')
%
% where: U(c) = log(c)
%        f(k) = A*k^alpha
%


%% 1. Initialisation
clear all; close all; clc;


%% 2. Model Parameters
% Economic Parameters
params.beta  = 0.96;                  % Discount factor
params.alpha = 0.25;                  % Capital output elasticity 
params.A     = (params.alpha * params.beta)^(-1);  % TFP (implies kss = 1)

% Numerical Parameters
params.kmin  = 0.8;                   % Minimum capital stock
params.kmax  = 1.2;                   % Maximum capital stock
params.n     = 41;                    % Grid size
params.maxit = 1000;                  % Maximum iterations
params.tol   = 1.0e-5;               % Convergence tolerance
params.print = true;                  % Plot option

%% 3. Discretisation of the state space (aka capital grid) 
%--------------------------------------------------------------------------
% Create evenly-spaced capital grid
k_grid = linspace(params.kmin, params.kmax, params.n)';

%% 4. Utility matrix
%--------------------------------------------------------------------------
% a. Compute utility matrix
k      = repmat(k_grid, 1, params.n);        % capital today matrix (auxiliary)
kprime = repmat(k_grid', params.n, 1);       % capital tomorrow matrix (auxiliary)

% Check for negative consumption
consumption = params.A * k.^params.alpha - kprime;
if any(consumption(:) <= 0)
    error('Negative consumption detected. Adjust grid or parameters.');
end
U = log(consumption);

clear k kprime consumption

% b. Illustration of utility matrix structure
%    |----------------------------------|
%    |u(k(1),k'(1)), ... ,u(k(1),k'(n)) |
% U =|  :                       :       |
%    |u(k(n),k'(1)), ... ,u(k(n),k'(n)) |                                                |
%    |----------------------------------|

%% 5. Value Function Iteration 
%--------------------------------------------------------------------------

% a. Initialise
V0 = zeros(params.n, 1);        % VF initial guess for RHS
V_history = zeros(params.n, params.maxit);  % Store VF evolution

% b. Iterate on VF 
for it = 1:params.maxit
    % Compute value function matrix
    VF = U + params.beta * repmat(V0', params.n, 1);
    
    % Maximise along rows to get optimal k' and value
    [Vmax, I] = max(VF, [], 2);
    
    % Check convergence
    dist = norm(Vmax - V0);
    V_history(:,it) = Vmax;  % Store current iteration
    
    if dist < params.tol
        fprintf('Convergence achieved after %d iterations\n', it);
        break
    elseif it == params.maxit
        warning('Maximum iterations reached without convergence');
    end
    
    V0 = Vmax;  % Update value function for next iteration
end

% c. Illustration of the VF matrix structure
%          |---------------------------------|         |-----------------------|
%          |u(k(1),k'(1)), ..., u(k(1),k'(n))|         |V(k'(1)), ..., V(k'(n))|
% VF = max |  :                       :      | + beta *|    :             :    |
%          |u(k(n),k'(1)), ..., u(k(n),k'(n))|         |V(k'(1)), ..., V(k'(n))|      
%          |---------------------------------|         |-----------------------| 

clear V0 VF dist U

%% 6. Policy Functions 
%--------------------------------------------------------------------------

% a. Capital policy function 
k_pol = k_grid(I);                    % numerical solution (elements of k_grid in position I)
k_pol_analytical = params.alpha * params.beta * params.A * k_grid.^params.alpha;  % analytical solution

% b. Consumption policy function
c_pol = params.A * k_grid.^params.alpha - k_pol;              % numerical solution
c_pol_analytical = params.A * k_grid.^params.alpha - k_pol_analytical;  % analytical solution

% c. Steady state values
k_ss = 1;  % By construction due to choice of A
c_ss = params.A * k_ss^params.alpha - k_ss;

%% 7. Visualisation
%--------------------------------------------------------------------------
if params.print
    % Create figure with subplots for better presentation
    figure('Name', 'Neoclassical Growth Model Results', 'Position', [100 100 1200 500])
    
    % a. Capital policy function
    subplot(1,2,1)
    plot(k_grid, k_pol, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Approximation at n=%d', params.n))
    hold on
    plot(k_grid, k_pol_analytical, 'r--', 'LineWidth', 2, 'DisplayName', 'Analytical policy function')
    plot(k_grid, k_grid, 'k:', 'LineWidth', 1.5, 'DisplayName', '45 degree line')
    plot(k_ss, k_ss, 'k*', 'MarkerSize', 10, 'DisplayName', 'Steady State')
    hold off
    
    xlabel('k_{t}', 'FontSize', 12)
    ylabel('Policy k_{t+1}', 'FontSize', 12)
    title('Capital Policy Function', 'FontSize', 14)
    legend('Location', 'northwest', 'FontSize', 10)
    grid on
    xlim([k_grid(1) k_grid(end)])
    ylim([k_grid(1) k_grid(end)])
    
    % b. Consumption policy function
    subplot(1,2,2)
    plot(k_grid, c_pol, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Approximation at n=%d', params.n))
    hold on
    plot(k_grid, c_pol_analytical, 'r--', 'LineWidth', 2, 'DisplayName', 'True policy function')
    plot([k_ss k_ss], [0 max(c_pol)*1.1], 'k:', 'LineWidth', 1.5, 'DisplayName', 'Steady State k')
    plot([k_grid(1) k_grid(end)], [c_ss c_ss], 'k:', 'LineWidth', 1.5, 'DisplayName', 'Steady State c')
    hold off
    
    xlabel('k_{t}', 'FontSize', 12)
    ylabel('Policy c_{t}', 'FontSize', 12)
    title('Consumption Policy Function', 'FontSize', 14)
    legend('Location', 'northeast', 'FontSize', 10)
    grid on
    
   
end
