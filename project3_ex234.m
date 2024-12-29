% Define the function and gradient
f = @(x1, x2) (1/3)*x1.^2 + 3*x2.^2; % Objective function
grad_f = @(x1, x2) [(2/3)*x1; 6*x2]; % Gradient of the function

% Problem constraints
x1_min = -10; x1_max = 5;
x2_min = -8; x2_max = 12;

% Projection function
projection = @(x) [min(max(x(1), x1_min), x1_max); min(max(x(2), x2_min), x2_max)];

% Parameters
epsilon = 0.01; % Tolerance for convergence
max_iter = 1250; % Maximum iterations

% Define test cases
cases = {
    struct('x0', [5; -5], 'sk', 5, 'gamma', 0.5, 'title', 'Case 2'),
    struct('x0', [-5; 10], 'sk', 15, 'gamma', 0.1, 'title', 'Case 3'),
    struct('x0', [8; -10], 'sk', 0.1, 'gamma', 0.2, 'title', 'Case 4')
};

% Iterate through cases
for c = 1:length(cases)
    % Extract parameters for the current case
    x0 = cases{c}.x0;
    sk = cases{c}.sk;
    gamma = cases{c}.gamma;
    title_str = cases{c}.title;
    
    % Initialize variables
    xk = x0;
    iter = 0;
    results = [xk', f(xk(1), xk(2))]; % Store initial results
    
    % Projected Gradient Descent Loop
    while iter < max_iter
        % Compute gradient
        grad = grad_f(xk(1), xk(2));
        
        % Update x̄ (unconstrained step)
        x_bar = xk - sk * grad;
        
        % Project onto feasible region
        xk_new = xk + gamma * (projection(x_bar) - xk);
        
        % Check convergence
        if norm(xk_new - xk) < epsilon
            break;
        end
        
        % Update variables
        xk = xk_new;
        iter = iter + 1;
        results = [results; xk', f(xk(1), xk(2))];
    end
    
    % Display results for the current case
    fprintf('%s:\n', title_str);
    fprintf('  Final Point: x1 = %.6f, x2 = %.6f\n', xk(1), xk(2));
    fprintf('  Function Value: f(x) = %.6f\n', f(xk(1), xk(2)));
    fprintf('  Iterations: %d\n\n', iter);
    
    % Plot convergence
    iterations = 1:size(results, 1);
    figure;
    plot(iterations, results(:, 3), 'o-', 'LineWidth', 2);
    xlabel('Iteration');
    ylabel('f(x)');
    title(sprintf('Convergence for %s', title_str));
    grid on;
    
    % Plot trajectory on contour plot
    x1_range = linspace(x1_min, x1_max, 100);
    x2_range = linspace(x2_min, x2_max, 100);
    [X1, X2] = meshgrid(x1_range, x2_range);
    Z = f(X1, X2);
    figure;
    contour(X1, X2, Z, 50);
    hold on;
    plot(results(:, 1), results(:, 2), 'o-', 'LineWidth', 2, ...
         'DisplayName', sprintf('%s Trajectory', title_str));
    xlabel('x1');
    ylabel('x2');
    title(sprintf('Optimization Path for %s', title_str));
    legend;
    grid on;
    hold off;
end
