% Define the function
f = @(x1, x2) (1/3)*x1.^2 + 3*x2.^2;

% Define the grid for x1 and x2
x1 = linspace(-10, 10, 100); % Range for x1
x2 = linspace(-10, 10, 100); % Range for x2
[X1, X2] = meshgrid(x1, x2); % Create the meshgrid

% Evaluate the function
Z = f(X1, X2);

% Plot the surface
figure;
surf(X1, X2, Z); % 3D surface plot
shading interp; % Smooth shading
colormap jet; % Color map
colorbar; % Color bar
xlabel('x_1');
ylabel('x_2');
zlabel('f(x_1, x_2)');
title('3D Plot of f(x_1, x_2) = (1/3)x_1^2 + 3x_2^2');
grid on;

% Plot the contour (optional, if needed)
figure;
contour(X1, X2, Z, 50); % Contour plot with 50 levels
xlabel('x_1');
ylabel('x_2');
title('Contour Plot of f(x_1, x_2) = (1/3)x_1^2 + 3x_2^2');
colorbar;
grid on;
