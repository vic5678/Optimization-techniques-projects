% Set interval [a, b]
a = -1;
b = 3;
l_values = [0.1, 0.05, 0.01, 0.005, 0.001];  % Different l values to test

% Define the three functions
f1 = @(x) (x - 2)^2 + x * log(x + 3);
f2 = @(x) exp(-2*x) + (x - 2)^2;
f3 = @(x) exp(x) * (x^3 - 1) + (x - 1) * sin(x);
functions = {f1, f2, f3};
titles = {'f1(x) = (x - 2)^2 + x * log(x + 3)', ...
          'f2(x) = exp(-2*x) + (x - 2)^2', ...
          'f3(x) = exp(x) * (x^3 - 1) + (x - 1) * sin(x)'};

figure;
for j = 1:3
    f = functions{j};
    iterations_l = zeros(size(l_values));
    for i = 1:length(l_values)
        l = l_values(i);
        [iterations_l(i), ~, ~] = golden_section_search(f, a, b, l);
    end
    subplot(3, 1, j);
    plot(l_values, iterations_l, '-o');
    title([' ', titles{j}]);
    xlabel('Final Interval Width (l)');
    ylabel('Number of Iterations');
    grid on;
end



% Create a figure for interval boundaries with different l values
figure;
for j = 1:3
    f = functions{j};
    % Subplot for each function
    subplot(3, 1, j);
    hold on;
    for i = 1:length(l_values)
        l = l_values(i);
        [ ~, a_vals, b_vals] = golden_section_search(f, a, b, l);
        % Plot [a, b] boundaries for current l with iteration index
        iterations_index = 1:length(a_vals);
        plot(iterations_index, a_vals, '-o', 'DisplayName', ['a, l = ' num2str(l)]);
        plot(iterations_index, b_vals, '-x', 'DisplayName', ['b, l = ' num2str(l)]);
    end
    title(['[a, b] Boundaries over Iterations for ', titles{j}]);
    xlabel('Iteration');
    ylabel('Boundary values');
    legend show;
    grid on;
    hold off;
end



function [counter, a_vals, b_vals] = golden_section_search(f, a, b, l)
    gamma = 0.618;
    % Initialize iteration counter and boundaries
    counter= 0;
    a_vals = [a];
    b_vals = [b];

    % Initialize x1,x2
    x1 = a + (1 - gamma) * (b - a);
    x2 = a + gamma * (b - a);
    f1 = f(x1);
    f2 = f(x2);
    counter=counter+2;

    while (b - a) >= l
        counter = counter + 1;
        if f1 < f2
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + (1 - gamma) * (b - a);
            f1 = f(x1);
        else
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + gamma * (b - a);
            f2 = f(x2);
        end
        % Store updated boundaries
        a_vals = [a_vals, a];
        b_vals = [b_vals, b];
    end
 
end

