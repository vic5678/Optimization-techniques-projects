
% Define the function and gradient
f = @(x, y) x.^5 .* exp(-x.^2 - y.^2);
grad_f = @(x, y) [5*x.^4.*exp(-x.^2 - y.^2) - 2*x.^6.*exp(-x.^2 - y.^2), ...
                  -2*x.^5.*y.*exp(-x.^2 - y.^2)];

% Parameters for the methods
tol = 1e-6; % Convergence tolerance
max_iter = 100; % Maximum iterations
initial_points = [0, 0; -1, 1; 1, -1]; % Initial points (i), (ii), (iii)
methods = {'constant', 'exact', 'armijo'}; % Step size strategies
alpha_constant = 0.1; % Fixed step size for constant method
c1 = 1e-4; % Armijo parameter
beta = 0.5; % Reduction factor for Armijo rule

% Golden Section parameters
golden_tol = 1e-6; % Tolerance for Golden Section Search

% Golden Section Search function
function gamma = golden_section(f_line, a, b, tol)
    % Golden ratio
    phi = (1 + sqrt(5)) / 2;
    rho = phi - 1;

    % Initial points
    c = b - rho * (b - a);
    d = a + rho * (b - a);

    % Function evaluations
    fc = f_line(c);
    fd = f_line(d);

    % Iterative loop
    while (b - a) > tol
        if fc < fd
            b = d;
            d = c;
            c = b - rho * (b - a);
            fd = fc;
            fc = f_line(c);
        else
            a = c;
            c = d;
            d = a + rho * (b - a);
            fc = fd;
            fd = f_line(d);
        end
    end

    % Return the midpoint as the result
    gamma = (a + b) / 2;
end

% Iterate over all initial points and methods
results_all = cell(size(initial_points, 1), length(methods));
for p = 1:size(initial_points, 1)
    x0 = initial_points(p, :)'; % Current initial point
    for m = 1:length(methods)
        method = methods{m};
        
        % Initialize variables
        xk = x0;
        iter = 0;
        results = [xk', f(xk(1), xk(2))]; % Store results for analysis
        
        % Steepest Descent Loop
        while iter < max_iter
            % Compute gradient
            grad = grad_f(xk(1), xk(2))';
            
            % Check stopping criterion
            if norm(grad) < tol
                break;
            end
            
            % Compute direction
            d = -grad; % Negative gradient direction
            
            % Step size selection
            switch method
                case 'constant'
                    % Fixed step size
                    alpha = alpha_constant;
                
                case 'exact'
                    % Exact line search (using Golden Section Search)
                    f_line = @(gamma) f(xk(1) + gamma * d(1), xk(2) + gamma * d(2));
                    alpha = golden_section(f_line, 0, 1, golden_tol);
                
                case 'armijo'
                    % Armijo rule
                    alpha = 1; % Start with initial step size
                    while f(xk(1) + alpha * d(1), xk(2) + alpha * d(2)) > ...
                          f(xk(1), xk(2)) + c1 * alpha * d'* (grad)
                        alpha = beta * alpha; % Reduce step size iteratively
                    end
            end
            
            % Update xk
            xk = xk + alpha * d;
            
            % Store results
            iter = iter + 1;
            results = [results; xk', f(xk(1), xk(2))];
        end
        
        % Save results
        results_all{p, m} = results;
    end
end

% Display results for each initial point and method
for p = 1:size(initial_points, 1)
    fprintf('Initial Point: (%.1f, %.1f)\n', initial_points(p, 1), initial_points(p, 2));
    for m = 1:length(methods)
        fprintf('  Method: %s\n', methods{m});
        final_result = results_all{p, m}(end, :);
        fprintf('    Minimum at: x = %.6f, y = %.6f\n', final_result(1), final_result(2));
        fprintf('    Function value: f(x, y) = %.6f\n', final_result(3));
        fprintf('    Iterations: %d\n', size(results_all{p, m}, 1) - 1);
    end
end

% Plot convergence for each initial point and method
figure;
for p = 1:size(initial_points, 1)
    for m = 1:length(methods)
        results = results_all{p, m};
        iterations = 1:size(results, 1);
        subplot(size(initial_points, 1), length(methods), (p-1)*length(methods)+m);
        plot(iterations, results(:, 3), 'o-');
        xlabel('Iteration');
        ylabel('f(x, y)');
        title(sprintf('Point: (%.1f, %.1f), Method: %s', ...
                      initial_points(p, 1), initial_points(p, 2), methods{m}));
        grid on;
    end
end



% Δημιουργία πλέγματος για τις ισοϋψείς καμπύλες
x_range = linspace(-2, 2, 100); % Εύρος για τον άξονα x
y_range = linspace(-2, 2, 100); % Εύρος για τον άξονα y
[X, Y] = meshgrid(x_range, y_range); % Δημιουργία πλέγματος
Z = f(X, Y); % Υπολογισμός της συνάρτησης στο πλέγμα

% Plot ισοϋψών καμπυλών και διαδρομή βελτιστοποίησης
figure;
for p = 1:size(initial_points, 1)
    for m = 1:length(methods)
        results = results_all{p, m};
        subplot(size(initial_points, 1), length(methods), (p-1)*length(methods)+m);
        
        % Σχεδίαση ισοϋψών καμπυλών
        contour(X, Y, Z, 50); % 50 ισοϋψείς καμπύλες
        hold on;
        
        % Σχεδίαση διαδρομής βελτιστοποίησης
        plot(results(:, 1), results(:, 2), 'o-', 'LineWidth', 2, ...
             'DisplayName', 'Optimization Path');
        
        % Εμφάνιση λεπτομερειών
        xlabel('x');
        ylabel('y');
        title(sprintf('Point: (%.1f, %.1f), Method: %s', ...
                      initial_points(p, 1), initial_points(p, 2), methods{m}));
        grid on;
        legend('show');
        hold off;
    end
end
