function yaxisPlot = ExtractYaxisMiddleAnalysis(X,Y,Bx,By)
close all;

middleIndex = ceil(length(X)/2);

Xmid = X(middleIndex:end,middleIndex);
Ymid = Y(middleIndex:end,middleIndex);
Bxmid = Bx(middleIndex:end,middleIndex);
Bymid = By(middleIndex:end,middleIndex);

midMagnitude = sqrt(Bxmid.^2 + Bymid.^2);

hold on;
scatterPlot = scatter(Ymid(2:end,:), midMagnitude(2:end,:), "DisplayName", "Dipole Magnitude");

% Set up the type of polynomial for the best fit
yInv3rdRoot = 1./(Ymid(2:end,:).^3);
% Make the best fit line for a inverse cube law fitting line to Magnitudes
% and Distances 1/ymid^3
coeffs = polyfit(yInv3rdRoot,midMagnitude(2:end,:),1);

% Calculate the fitted magnitude values using the coefficients from the fit
magFit = coeffs(1) * (1./(Ymid(2:end,:).^3)) + coeffs(2);

% Plot the best fit line
yaxisPlot = plot(Ymid(2:end), magFit, 'r-', 'DisplayName', 'Best Fit Line');

% Add labels and a legend
xlabel('y axis (Samples)');
ylabel('Simulated Magnetic Field Magnitude');
title('Inverse Cube Law Fit: simulated dipole magnitude vs y samples within 1 m');
legend('show');

% Print out the equation in the format: y = a * 1/x^3 + b
fprintf('Fitted equation: y = %.4f * 1/x^3 + %.4f\n', coeffs(1), coeffs(2));
end