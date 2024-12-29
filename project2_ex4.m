clear;

% Define the function, gradient, and Hessian
f = @(x, y) x.^5 .* exp(-x.^2 - y.^2);
grad_f = @(x, y) [5*x.^4.*exp(-x.^2 - y.^2) - 2*x.^6.*exp(-x.^2 - y.^2), ...
                  -2*x.^5.*y.*exp(-x.^2 - y.^2)];
hessian_f = @(x, y) [exp(-x.^2 - y.^2) .* (20*x.^3 - 14*x.^5 + 4*x.^7), ...
                     -2*x.^5.*y.*exp(-x.^2 - y.^2) .* (2*x.^2 - 1); ...
                     -2*x.^5.*y.*exp(-x.^2 - y.^2) .* (2*x.^2 - 1), ...
                     -2*x.^5.*exp(-x.^2 - y.^2) .* (1 - 2*x.^2 - 2*y.^2)];

% Parameters for the algorithm
tol = 1e-6; % Convergence tolerance
max_iter = 100; % Maximum iterations
initial_points = [0, 0; -1, 1; 1, -1]; % Initial points (i), (ii), (iii)
methods = {'constant', 'exact', 'armijo'}; % Step size strategies
mu = 1; % Initial regularization parameter
beta_mu = 2; % Scaling factor for mu adjustment
alpha_constant = 0.1; % Fixed step size for constant method
c1 = 1e-4; % Armijo parameter
tau = 0.5; % Reduction factor for Armijo rule
golden_tol = 1e-6; % Tolerance for Golden Section Search

% Golden Section Search function
function gamma = golden_section(f_line, a, b, tol)
    phi = (1 + sqrt(5)) / 2;
    rho = phi - 1;

    c = b - rho * (b - a);
    d = a + rho * (b - a);

    fc = f_line(c);
    fd = f_line(d);

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

    gamma = (a + b) / 2;
end

% Criteria Testing Function
function [criteria_flag] = criteria_testing(gamma, d, xk, grad_f, f)
    % Criteria testing for Levenberg-Marquardt
    % Returns true if both criteria 3 and 4 are satisfied
    criteria_flag = false;
    
    % Compute new point x_{k+1}
    xk_new = xk + gamma * d;
    
    % Criterion 3: Directional derivative check
    z_1 = d' * grad_f(xk_new(1), xk_new(2))';
    z_2 = d' * grad_f(xk(1), xk(2))';
    for beta = linspace(0.1, 1, 50) % Adjust beta range
        if z_1 > beta * z_2
            % Criterion 4: Sufficient decrease condition
            w_1 = f(xk_new(1), xk_new(2));
            w_2 = f(xk(1), xk(2));
            w_3 = gamma * d' * grad_f(xk(1), xk(2))';
            
            for alpha = linspace(0.001, beta, 50) % Adjust alpha range
                if w_1 <= w_2 + alpha * w_3
                    criteria_flag = true;
                    return;
                end
            end
        end
    end
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
        
        % Levenberg-Marquardt method loop
        while iter < max_iter
            % Compute gradient and Hessian
            grad = grad_f(xk(1), xk(2))';
            hess = hessian_f(xk(1), xk(2));
             % Compute eigenvalues of Hessian to initialize mu
            eigenvalues = eig(hess);
            lambda_max = max(abs(eigenvalues));
            if lambda_max <= 0
                lambda_max = 1e-3; % Safeguard for non-positive definite Hessians
            end
            mu = lambda_max + 0.1; % Ensure mu is slightly larger than max eigenvalue
            % Check stopping criterion
            if norm(grad) < tol
                break;
            end

            % Modify the Hessian for Levenberg-Marquardt
            hess_mod = hess + mu * eye(2);

            % Solve for direction d_k
            d = -hess_mod \ grad;

            % Step size selection
            switch method
                case 'constant'
                    gamma = alpha_constant;
                case 'exact'
                    f_line = @(gamma) f(xk(1) + gamma * d(1), xk(2) + gamma * d(2));
                    gamma = golden_section(f_line, 0, 1, golden_tol);
                case 'armijo'
                    gamma = 1;
                    while f(xk(1) + gamma * d(1), xk(2) + gamma * d(2)) > ...
                          f(xk(1), xk(2)) + c1 * gamma * grad' * d
                        gamma = tau * gamma;
                    end
            end

            % Ensure criteria are satisfied
            while ~criteria_testing(gamma, d, xk, grad_f, f)
                gamma = gamma / 2; % Reduce gamma iteratively
                if gamma < 1e-8
                    warning('Gamma reached minimum threshold. Stopping criteria adjustment.');
                    break;
                end
            end

            % Update xk
            xk_new = xk + gamma * d;

            % Update mu based on progress
            if f(xk_new(1), xk_new(2)) < f(xk(1), xk(2)) % Progress made
                mu = mu / beta_mu; % Decrease regularization
                xk = xk_new; % Accept the step
            else
                mu = mu * beta_mu^2; % Increase regularization more aggressively
            end

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
