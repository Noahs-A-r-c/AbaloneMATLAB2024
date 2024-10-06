T = readtable("Magnetometer Magnitude to Distance 1 Dimension - 8_23_24 9_54 PM Data Raw.csv");

% Save the original scatter data
xOrig = T.Var2;
yOrig = T.Var1;

% Cut the flat part out for best fit; max magnetism is undeterminable
xCut = T.Var2(26:end);
yCut = T.Var1(26:end);

xInv3rdRoot = 1./(xCut.^3);
% Make the best fit line for a inverse cube law fitting line to Magnitudes
% and Distances
coeffs = polyfit(xInv3rdRoot,yCut,1);

% Plot the original data points
figure;
s = scatter(xOrig, yOrig, 'b', 'DisplayName', 'Original Data');
s.SizeData = 10;
hold on;

% Calculate the fitted y values using the coefficients from the fit
yFit = coeffs(1) * (1./(xCut.^3)) + coeffs(2);

% Plot the best fit line
plot(xCut, yFit, 'r-', 'DisplayName', 'Best Fit Line');

% Add labels and a legend
xlabel('x (Distance)');
ylabel('y (Magnetometer Magnitude)');
title('Inverse Cube Law Fit: y vs x');
legend('show');

% Hold off the figure
hold off;

% Print out the equation in the format: y = a * 1/x^3 + b
fprintf('Fitted equation: y = %.4f * 1/x^3 + %.4f\n', coeffs(1), coeffs(2));


% Find and Print the minimum vertex of the magnitude and it's distance
% This is the value at which the Earth's magnetic field takes over
% and our measurement is no longer useful
maxUsefulDist = xOrig(find(yOrig==min(yOrig)))
