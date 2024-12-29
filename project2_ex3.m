
% Define the function, gradient, and Hessian
f = @(x, y) x.^5 .* exp(-x.^2 - y.^2);
grad_f = @(x, y) [5*x.^4.*exp(-x.^2 - y.^2) - 2*x.^6.*exp(-x.^2 - y.^2), ...
                  -2*x.^5.*y.*exp(-x.^2 - y.^2)];
hessian_f = @(x, y) [exp(-x.^2 - y.^2) .* (20*x.^3 - 14*x.^5 + 4*x.^7), ...
                     -2*x.^5.*y.*exp(-x.^2 - y.^2) .* (2*x.^2 - 1); ...
                     -2*x.^5.*y.*exp(-x.^2 - y.^2) .* (2*x.^2 - 1), ...
                     -2*x.^5.*exp(-x.^2 - y.^2) .* (1 - 2*x.^2 - 2*y.^2)];

% Parameters
tol = 1e-6; % Convergence tolerance
max_iter = 100; % Maximum iterations
initial_points = [0, 0; -1, 1; 1, -1]; % Initial points (i), (ii), (iii)
methods = {'constant', 'exact', 'armijo'}; % Step size strategies
alpha_constant = 0.1; % Fixed step size for constant method
c1 = 1e-4; % Armijo parameter
tau = 0.5; % Armijo scaling factor
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

% Iterate over all initial points and methods
results_all = cell(size(initial_points, 1), length(methods));
for p = 1:size(initial_points, 1)
    x0 = initial_points(p, :)'; % Initial point
    for m = 1:length(methods)
        method = methods{m};
        
        xk = x0;
        iter = 0;
        results = [xk', f(xk(1), xk(2))]; % Store results for analysis

        % Newton's method loop
        while iter < max_iter
            grad = grad_f(xk(1), xk(2))';
            hess = hessian_f(xk(1), xk(2));

            if norm(grad) < tol
                break;
            end

            dx = -hess \ grad;

            switch method
                case 'constant'
                    alpha = alpha_constant;
                case 'exact'
                    f_line = @(gamma) f(xk(1) + gamma * dx(1), xk(2) + gamma * dx(2));
                    alpha = golden_section(f_line, 0, 1, golden_tol);
                case 'armijo'
                    alpha = 1;
                    while f(xk(1) + alpha * dx(1), xk(2) + alpha * dx(2)) > ...
                          f(xk(1), xk(2)) + c1 * alpha * grad' * dx
                        alpha = tau * alpha;
                    end
            end

            xk = xk + alpha * dx;
            iter = iter + 1;
            results = [results; xk', f(xk(1), xk(2))];
        end

        results_all{p, m} = results;
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
