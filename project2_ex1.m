% Define the function
f = @(x, y) x.^5 .* exp(-x.^2 - y.^2);

% Define the range for x and y
x = linspace(-5, 5, 100); % Range for x
y = linspace(-5, 5, 100); % Range for y

% Create a meshgrid for plotting
[X, Y] = meshgrid(x, y);

% Evaluate the function
Z = f(X, Y);

% Plot the function
figure;
surf(X, Y, Z); % 3D surface plot
shading interp; % Smooth shading
colormap jet;
colorbar;
xlabel('x');
ylabel('y');
zlabel('f(x, y)');
title('Plot of f(x, y) = x^5 * exp(-x^2-y^2)');
    
grid on;