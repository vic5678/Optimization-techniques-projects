% Define the function, gradient
f = @(x1, x2) (1/3)*x1.^2 + 3*x2.^2;
grad_f = @(x1, x2) [(2/3)*x1; 6*x2];

% Parameters
tol = 0.001; % Convergence tolerance
max_iter = 150; % Maximum iterations
initial_point = [8; -10]; % Example initial point
gammas = [0.1, 0.3, 3, 5]; % Different gamma values for testing

% Results storage
results_all = cell(length(gammas), 1);

% Define the grid for contour plots
x_range = linspace(-10, 10, 100);
y_range = linspace(-10, 10, 100);
[X, Y] = meshgrid(x_range, y_range);
Z = f(X, Y);

% Iterate for each gamma
for g = 1:length(gammas)
    gamma = gammas(g);
    xk = initial_point; % Starting point
    iter = 0; % Iteration counter
    results = [xk', f(xk(1), xk(2))]; % Store initial point and function value
    
    % Steepest Descent Loop
    while iter < max_iter
        % Compute gradient
        grad = grad_f(xk(1), xk(2));
        
        % Check stopping criterion
        if norm(grad) < tol
            break;
        end
        
        % Update xk
        xk = xk - gamma * grad;
        
        % Store results
        iter = iter + 1;
        results = [results; xk', f(xk(1), xk(2))];
    end
    
    % Save results
    results_all{g} = results;
end

% Plot convergence and isocontours for each gamma in a subplot
figure;
for g = 1:length(gammas)
    results = results_all{g};
    
    % Create a subplot for each gamma
    subplot(ceil(length(gammas)/2), 2, g); % Arrange subplots in 2 columns
    
    % Plot the contour of the function
    contour(X, Y, Z, 50); % 50 contour levels
    hold on;
    
    % Overlay the optimization trajectory
    plot(results(:, 1), results(:, 2), 'o-', 'LineWidth', 2, 'DisplayName', sprintf('Gamma = %.1f', gammas(g)));
    xlabel('x_1');
    ylabel('x_2');
    title(sprintf('Convergence for Gamma = %.1f', gammas(g)));
    legend('show');
    grid on;
end
sgtitle('Convergence of Steepest Descent with Isocontours'); % Overall title
% Display results
for g = 1:length(gammas)
    fprintf('Gamma: %.1f\n', gammas(g));
    final_result = results_all{g}(end, :);
    fprintf('  Final Point: x1 = %.6f, x2 = %.6f\n', final_result(1), final_result(2));
    fprintf('  Function Value: f(x) = %.6f\n', final_result(3));
    fprintf('  Iterations: %d\n', size(results_all{g}, 1) - 1);
end

% Plot convergence for each gamma in a subplot
figure;
for g = 1:length(gammas)
    results = results_all{g};
    iterations = 1:size(results, 1);
    
    % Create a subplot for each gamma
    subplot(ceil(length(gammas)/2), 2, g); % Arrange subplots in 2 columns
    plot(iterations, results(:, 3), 'o-', 'LineWidth', 2);
    xlabel('Iteration');
    ylabel('f(x)');
    title(sprintf('Convergence for Gamma = %.1f', gammas(g)));
    grid on;
end
sgtitle('Convergence of Steepest Descent for Different Gamma'); 