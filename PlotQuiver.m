% PlotQuiver.m
% Function to plot normalized colored quiver

function quiv = PlotQuiver(X,Y,Bx,By)
clf; close all;
% Normalize the vectors for better visualization
% Calculate the magnitude
magnitude = sqrt(Bx.^2 + By.^2);
% Normalize only if magnitude is not zero
Bx_normalized = Bx ./ (magnitude + eps);  % 'eps' avoids division by zero
By_normalized = By ./ (magnitude + eps);

% Normalize the magnitude for coloring (0 to 1 range)
max_magnitude = max(magnitude(:));  % Get the maximum magnitude
norm_mag = magnitude / max_magnitude;  % Normalize magnitudes

% Scale down the normalized vectors for better visualization
scale_factor = 0.1;  % Adjust this value as needed
Bx_scaled = Bx_normalized * scale_factor;
By_scaled = By_normalized * scale_factor;

% Create a downsampled grid for quiver plot (optional)
% For example, you can take every nth point (e.g., every 5th point)
downsample_factor = 1;  % Adjust this value as needed
quiver_X = X(1:downsample_factor:end, 1:downsample_factor:end);
quiver_Y = Y(1:downsample_factor:end, 1:downsample_factor:end);
quiver_Bx = Bx_scaled(1:downsample_factor:end, 1:downsample_factor:end);
quiver_By = By_scaled(1:downsample_factor:end, 1:downsample_factor:end);

% Prepare the colormap
colormap_name = 'jet';  % Change 'jet' to any colormap you prefer
color_map = colormap(colormap_name);  % Get colormap data

% Plot the vector field with colored arrows
figure; hold on;

for i = 1:numel(quiver_X)
    % Get the color based on the normalized magnitude at the current position
    color_index = round(norm_mag(quiver_Y(i) == Y(:,1) & quiver_X(i) == X(1,:)) * (size(color_map, 1) - 1)) + 1;  % Adjust index
    color_index = max(min(color_index, size(color_map, 1)), 1);  % Clamp index to valid range
    % Plot each arrow separately with its color
    quiver(quiver_X(i), quiver_Y(i), quiver_Bx(i), quiver_By(i), 0, 'Color', color_map(color_index, :), 'LineWidth', 1.5, "MaxHeadSize", 3);
end

hold off;
xlabel('x (m)');
ylabel('y (m)');
title('Magnetic Field of a Dipole with Color-coded Magnitude');
axis equal;  % Keep the aspect ratio equal for better visualization
colorbar;  % Optional: add a colorbar to represent magnitude scale
end