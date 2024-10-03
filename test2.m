% Declare constants
u0 = 4 * pi * 10^-7;

% Assuming N42 and a 3" magnet of diameter 0.5mm and thickness of 75mm
dpMoment = 9.98; % A*m^2

x = linspace(-1, 1, 20); % -1 to 1 m in 100 intervals
y = linspace(-1, 1, 20); % 1 to -1 m in 100 intervals
[X, Y] = meshgrid(x, y);

% Initialize magnetic field components
Bx = zeros(size(X));
By = zeros(size(Y));

% Loop through grid and compute the magnetic field using the dipole formula
for i = 1:numel(X)
    % Position vector from dipole center (assuming dipole at origin)
    r = [X(i); Y(i)];  
    r_mag = norm(r);  % Magnitude of position vector
    theta = 90 - atan2(Y(i), X(i));
    
    % Store the field components
    const = u0 * dpMoment / (4 * pi * r_mag^3);
    Bx(i) = real(const * 3 * sin(theta) * cos(theta));  % Use sind and cosd for degrees
    By(i) = real(const * (3 * cos(theta)^2 - 1));
end

% Replace all Inf and NaN with zero
Bx(isnan(Bx) | ~isfinite(Bx)) = 0;
By(isnan(By) | ~isfinite(By)) = 0;

%PlotQuiver(X,Y,Bx,By)
