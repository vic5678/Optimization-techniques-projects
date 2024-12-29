a = -1;
b = 3;
% Define the three functions
f1 = @(x) (x - 2)^2 + x * log(x + 3);
f2 = @(x) exp(-2*x) + (x - 2)^2;
f3 = @(x) exp(x) * (x^3 - 1) + (x - 1) * sin(x);
functions = {f1, f2, f3};
titles = {'f1(x) = (x - 2)^2 + x * log(x + 3)', ...
          'f2(x) = exp(-2x) + (x - 2)^2', ...
          'f3(x) = exp(x) * (x^3 - 1) + (x - 1) * sin(x)'};

 % Experiment: Vary epsilon for Dichotomous Search with constant l
% Define a range of epsilon values to test
epsilon_values = [0.002,0.001, 0.0005, 0.0001];
l=0.01;
figure;
for j = 1:3
    f = functions{j};
    iterations_epsilon = zeros(size(epsilon_values));
    for i = 1:length(epsilon_values)
        epsilon = epsilon_values(i);
        [~, ~, iterations_epsilon(i), ~, ~] = dichotomous_search(f, a, b, epsilon, l);
    end
    subplot(3, 1, j);%sublot for each function
    plot(epsilon_values, iterations_epsilon, '-o');
    title(['Dichotomous Search: Iterations vs Epsilon for ', titles{j}]);
    xlabel('Epsilon');
    ylabel('Number of Iterations');
    grid on;
end





 % Experiment: Vary l for Dichotomous Search with constant epsilon
epsilon=0.001;
% Define a range of l values to test
l_values = [0.1,0.08, 0.05,0.02, 0.01,0.005];
% Create a figure for the iterations vs l plot
figure;
for j = 1:3
    f = functions{j};
    iterations_l = zeros(size(l_values));
    for i = 1:length(l_values)
        l = l_values(i);
        [~, ~, iterations_l(i), ~, ~] = dichotomous_search(f, a, b, epsilon, l);
    end 
    subplot(3, 1, j);
    plot(l_values, iterations_l, '-o');
    title(['Dichotomous Search: Iterations vs Final Interval Width (l) for ', titles{j}]);
    xlabel('Final Interval Width (l)');
    ylabel('Number of Iterations');
    grid on;
end





% Create a figure for the interval boundaries plot
% Define a range of l values to test
l_values = [0.1, 0.05, 0.01];
figure;
for j = 1:3
    f = functions{j};
    % Subplot setup for each function with different l values
    subplot(3, 1, j);
    hold on;
    for i = 1:length(l_values)
        l = l_values(i);
        [~, ~, ~, a_vals, b_vals] = dichotomous_search(f, a, b, epsilon, l);
        iterations_index = 1:length(a_vals);
        plot(iterations_index, a_vals, '-o', 'DisplayName', ['a, l = ' num2str(l)]);
        plot(iterations_index, b_vals, '-x', 'DisplayName', ['b, l = ' num2str(l)]);
    end
    hold off;
    title(['[a, b] Boundaries over Iterations for ', titles{j}]);
    xlabel('Iteration');
    ylabel('Boundary values');
    legend show;
    grid on;
end


function [x_min, f_min, iterations, a_vals, b_vals] = dichotomous_search(f, a, b, epsilon, l)
    % Initialize iteration counter and boundary storage
    iterations = 0;
    a_vals = [a];
    b_vals = [b];

    % Dichotomous Search Loop
    while (b - a) >= l
        iterations = iterations + 2;% Increase iteration counter
        x1 = (a + b) / 2 - epsilon;
        x2 = (a + b) / 2 + epsilon;
        f1 = f(x1);
        f2 = f(x2);
        
        % Update the interval based on f1 and f2
        if f1 < f2
            b = x2;    
        else 
            a = x1;     
        end
        % Store updated boundaries
        a_vals = [a_vals, a];
        b_vals = [b_vals, b];
    end
     x_min = (a + b) / 2;
    f_min = f(x_min);
end







