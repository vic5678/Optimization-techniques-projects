
function [counter, a_vals, b_vals] = fibonacci_search(f, a, b, l)
    % Generate Fibonacci numbers up to the desired accuracy
    fib_numbers = [1, 1];
    while fib_numbers(end) < (b - a) / l
        fib_numbers = [fib_numbers, fib_numbers(end) + fib_numbers(end - 1)];
    end
    epsilon=0.001;
    n = length(fib_numbers);
    counter = 0;
    a_vals = [a];
    b_vals = [b];
    
    % Initial interior points
    x1 = a + fib_numbers(n - 2) / fib_numbers(n) * (b - a);
    x2 = a + fib_numbers(n - 1) / fib_numbers(n) * (b - a);
    f1 = f(x1);
    f2 = f(x2);
    counter=2;
    % Fibonacci Search Loop
    for k = 1:(n - 2)  % Adjust loop range to avoid out-of-bounds indexing
        if k > n - 2
            break;  % Prevent any invalid indexing
        end
        if k==n-2
            x2=x1+epsilon;
            f2=f(x2);
            counter=counter+1;
             if f1 < f2
              b = x2;
              break;
             end
            if f1>f2
            a = x1;
            break;
            end
        else
          if f1 < f2
              b = x2;
              x2 = x1;
              f2 = f1;
              x1 = a + fib_numbers(n - k - 2) / fib_numbers(n - k) * (b - a);
              f1 = f(x1);
              counter = counter + 1;
          else
              a = x1;
              x1 = x2;
              f1 = f2;
              x2 = a + fib_numbers(n - k - 1) / fib_numbers(n - k) * (b - a);
              f2 = f(x2);
              counter = counter + 1;        
          end
        end
        % Store updated boundaries
        a_vals = [a_vals, a];
        b_vals = [b_vals, b];  
    end
end


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

% Create a figure for the iterations plot vs l
figure;
for j = 1:3
    f = functions{j};
    iterations_l_fib = zeros(size(l_values));
    for i = 1:length(l_values)
        l = l_values(i);
        [iterations_l_fib(i), ~, ~] = fibonacci_search(f, a, b, l);
    end    
    subplot(3, 1, j);
    plot(l_values, iterations_l_fib, '-o');
    title(['Fibonacci Search- ', titles{j}]);
    xlabel('Final Interval Width (l)');
    ylabel('Number of Iterations');
    grid on;
end

% Create a figure for the interval boundaries plot
%figure;
%for j = 1:3
   % f = functions{j};
    %l = 0.01;  % Use a fixed value of l for this plot
    
    % Run Fibonacci Search to capture [a, b] values
   % [ ~, a_vals, b_vals] = fibonacci_search(f, a, b, l);
    % Subplot for [a, b] boundaries over iterations
    %subplot(3, 1, j);
    %iterations_index = 1:length(a_vals);
   % plot(iterations_index, a_vals, '-o', iterations_index, b_vals, '-x');
   % title(['[a, b] Boundaries over Iterations - ', titles{j}]);
    %xlabel('Iteration');
   % ylabel('Boundary values');
    %legend('a_k', 'b_k');
   % grid on;
%end


% Create a figure for interval boundaries with different l values
figure;
for j = 1:3
    f = functions{j};
    % Subplot for each function
    subplot(3, 1, j);
    hold on;
    for i = 1:length(l_values)
        l = l_values(i);
        [ ~, a_vals, b_vals] = fibonacci_search(f, a, b, l);
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