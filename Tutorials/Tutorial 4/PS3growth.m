%**************************************************************************
%                   Tutorial 4. One-sector Growth Model                   %
%                 MSc Advanced Economic Theory (ECON0057)                 %
%                             Problem set 3                               %
%**************************************************************************
%  
% Numerical Solution to the Neoclassical Growth Model
%
% The dynamic problem:
% 
%          V(k,z) = max U[exp(z)*f(k) - k'] + beta*E[V(k',z')]
%
% where: U(c) = log(c)
%        f(k) = A*k^alpha

%% 1. Initialization
clear all; close all; clc;

%% 2. Parameters
% Economic Parameters
params.beta  = 0.96;                  % Discount factor
params.alpha = 0.25;                  % Capital output elasticity 
params.sigma = 0.5;                   % Standard deviation of shocks
params.A     = (params.alpha * params.beta)^(-1);  % TFP (implies kss = 1)

% Grid Parameters
params.kmin  = 0.0875;               % Lowest level of capital
params.kmax  = 3.5875;               % Highest level of capital
params.n     = 50;                   % Capital grid length
params.zmin  = -3 * params.sigma;    % Lowest TFP shock
params.zmax  = 3 * params.sigma;     % Highest TFP shock 
params.m     = 7;                    % TFP shock grid length

% Numerical Parameters
params.maxit = 1000;                 % Maximum iterations
params.tol   = 1.0e-5/(1-params.beta); % Convergence tolerance
params.print = true;                 % Plot option

%% 3. Grid Construction
% Create state space grids
k_grid = linspace(params.kmin, params.kmax, params.n)';
z_grid = linspace(params.zmin, params.zmax, params.m)';

% Compute shock density
q_dens = normpdf(z_grid/params.sigma, 0);
q_dens = q_dens / sum(q_dens);

% Visualize grids if requested
if params.print
    figure('Name', 'Shock Distributions', 'Position', [100 100 1200 500])
    
    % Normal distribution of shocks
    subplot(1,2,1)
    plot(z_grid, q_dens, 'LineWidth', 2)
    xlabel('z', 'FontSize', 12)
    ylabel('Pr(z)', 'FontSize', 12)
    title('Density of Shocks', 'FontSize', 14)
    grid on
    
    % Log-normal distribution
    subplot(1,2,2)
    plot(exp(z_grid), q_dens, 'LineWidth', 2)
    xlabel('e^z', 'FontSize', 12)
    ylabel('Pr(e^z)', 'FontSize', 12)
    title('Density of Log-normal Shocks', 'FontSize', 14)
    grid on
end

%% 4. Utility Computation
% Create 3D matrices for computation
[K, KPRIME, Z] = ndgrid(k_grid, k_grid, z_grid);

% Compute production and utility
F = exp(Z) .* params.A .* K.^params.alpha;
U = log(max(F - KPRIME, 0));

clear K KPRIME Z F

%% 5. Value Function Iteration
% Initialize
V0 = zeros(params.n, params.m);
tic;

% Main iteration loop
for it = 1:params.maxit
    % Compute expected value
    EV = V0 * q_dens;
    EV1 = repmat(EV', params.n, 1, params.m);
    V1 = U + params.beta * EV1;
    
    % Find optimal policy
    [VF, I] = max(V1, [], 2);
    VF = reshape(VF, params.n, params.m);
    
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

%% 6. Policy Function Construction
% Reshape policy function
kindex = reshape(I, params.n, params.m);
k_star = k_grid(kindex);
k_star = reshape(k_star, params.n, params.m);

%% 7. Visualization
if params.print
    % Figure 3: Policy functions for all z_t
    figure('Name', 'Policy Functions', 'Position', [100 100 1200 400])
    
    subplot(1,2,1)
    plot(k_grid, k_star, 'LineWidth', 2)
    xlabel('k_t', 'FontSize', 12)
    ylabel('k_{t+1}', 'FontSize', 12)
    title('Capital Policy Functions for Different z_t', 'FontSize', 14)
    grid on
    
    % Figure 4: Policy function for specific z_t
    zval = floor(4*params.m/5);
    subplot(1,2,2)
    plot(k_grid, k_star(:,zval), 'LineWidth', 2)
    xlabel('k_t', 'FontSize', 12)
    ylabel('k_{t+1}', 'FontSize', 12)
    title(sprintf('Capital Policy Function (z_t = %.2f)', z_grid(zval)), 'FontSize', 14)
    grid on
    
    % Figure 5: Value function
    figure('Name', 'Value Function', 'Position', [100 100 600 400])
    plot(k_grid, VF(:,zval), 'LineWidth', 2)
    xlabel('k_t', 'FontSize', 12)
    ylabel('V(k,z)', 'FontSize', 12)
    title(sprintf('Value Function (z_t = %.2f)', z_grid(zval)), 'FontSize', 14)
    grid on
end

% Display computation time
fprintf('Computation time: %.2f seconds\n', computation_time)