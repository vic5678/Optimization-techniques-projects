% Set interval [a, b]
a = -1;
b = 3;
epsilon = 0.001;  
l_values = [0.1, 0.05, 0.01, 0.005, 0.001];  % Different l values to test

% Define the three functions and their derivatives
f1 = @(x) (x - 2)^2 + x * log(x + 3);
df1 = @(x) 2*(x - 2) + log(x + 3) + x * (1 / (x + 3));

f2 = @(x) exp(-2*x) + (x - 2)^2;
df2 = @(x) -2*exp(-2*x) + 2 * (x - 2);

f3 = @(x) exp(x) * (x^3 - 1) + (x - 1) * sin(x);
df3 = @(x) exp(x) .* (3 * x.^2 + x.^3 - 1) + (x - 1) .* cos(x) + sin(x);

functions = {f1, f2, f3};
derivatives = {df1, df2, df3};
titles = {'f1(x) = (x - 2)^2 + x * log(x + 3)', ...
          'f2(x) = exp(-2*x) + (x - 2)^2', ...
          'f3(x) = exp(x) * (x^3 - 1) + (x - 1) * sin(x)'};

% Create a figure for the iterations plot vs l
figure;
for j = 1:3
    f = functions{j};
    df = derivatives{j};
    iterations_l = zeros(size(l_values));
    for i = 1:length(l_values)
        l = l_values(i);
        [ iterations_l(i), ~, ~] = dichotomous_search_derivative(f, df, a, b, epsilon, l);
    end
    % Subplot for Number of Iterations vs. l
    subplot(3, 2, 2*j - 1);
    plot(l_values, iterations_l, '-o');
    title([' ', titles{j}]);
    xlabel('Final Interval Width (l)');
    ylabel('Number of Iterations');
    grid on;
end



derivative={df1,df2,df3};
% Create a figure for interval boundaries with different l values
figure;
for j = 1:3
    f = functions{j};
    df = derivative{j};
    % Subplot for each function
    subplot(3, 1, j);
    hold on;
    for i = 1:length(l_values)
        l = l_values(i);
        [ ~, a_vals, b_vals] =dichotomous_search_derivative(f, df, a, b, epsilon, l);
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


function [ iterations, a_vals, b_vals] = dichotomous_search_derivative(f, df, a, b, epsilon, l)
   
    % Initialize iteration counter and boundaries
    iterations = 0;
    a_vals = [a];
    b_vals = [b];
    n=log2((b-a)/l); %initial condition for n iterations
    k=1;

   while k<=n && b-a> l
        iterations = iterations + 1;
        % Evaluate derivative at midpoint
        midpoint = (a + b) / 2;
        df_mid = df(midpoint);
        % Update interval based on derivative
        if df_mid==0
            x_min=midpoint;
            break;
        elseif df_mid < 0
            a = midpoint;
        else
            b = midpoint;
        end
        % Store updated boundaries
        a_vals = [a_vals, a];
        b_vals = [b_vals, b];
    end
   
end

