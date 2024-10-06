% DipoleMake.m
% Create Dipole Model given a scaling factor (Dipole Moment)
% ouptuts the Bx and By vectors given the parameters for coordinate planes
% X an Y in meters range and a number of samples, n

function [X,Y,Bx,By] = DipoleMake(Xstart,Xend,Ystart,Yend,n,dpMoment)

% A script to create a 2D Matlab model of a magnetic dipole and display if
% needed.

% Declare constants
u0 = 4 * pi * 10^-7;

x = linspace(Xstart, Xend, n); % 0.4 to -0.4 m in 21 intervals
y = linspace(Ystart, Yend, n); % 0.4 to 0 m in 21 intervals
[X, Y] = meshgrid(x, y);

% Initialize magnetic field components
Bx = zeros(size(X));
By = zeros(size(Y));

% Loop through grid and compute the magnetic field using the dipole formula
for i = 1:numel(X)
    % Position vector from dipole center (assuming dipole at origin)
    r = [X(i); Y(i)];  
    r_mag = norm(r);  % Magnitude of position vector
    theta = atan2(X(i), Y(i));
    
    % Store the field components
    const = u0 * dpMoment / (4 * pi * r_mag^3);
    Bx(i) = real(const * 3 * sin(theta) * cos(theta));  % Use sind and cosd for degrees
    By(i) = real(const * (3 * cos(theta)^2 - 1));
end

% Replace all Inf and NaN with zero
Bx(isnan(Bx) | ~isfinite(Bx)) = 0;
By(isnan(By) | ~isfinite(By)) = 0;

% Replace the center value for Bx and By with the next value as the origin
% point is not accounted for in the model
Bx(1,ceil(length(Bx)/2)) = Bx(2,ceil(length(Bx)/2));
By(1,ceil(length(Bx)/2)) = By(2,ceil(length(Bx)/2));

% Find the highest magnitude for reference
mag = sqrt(By.^2 + Bx.^2);
maxMag = max(max(mag));
[maxMagIndexY,maxMagIndexX] = find(mag == maxMag,1);
fprintf("The maximum grid magnitude for a vector is: %d at [%d,%d]\n",maxMag,maxMagIndexX,maxMagIndexY)

end