%% DipoleModel.m computation section
% Last modified 10/1/24 1:34 AM - Noah Nguyen UCD ECE Jlab
clc;

% A script to create a 2D Matlab model of a magnetic dipole and display if
% needed.

% Declare constants
u0 = 4 * pi * 10^-7;

% Assuming N42 and a 3" magnet of diameter 0.5mm and thickness of 75mm
%dpMoment = 9.98; % A*m^2 
dpMoment =  9.98 * 9.2050494347e15; % A*m^2  * conversion to IMU data number to match OneDimMagnitude.m

x = linspace(-0.4, 0.4, 21); % -1 to 1 m in 100 intervals
y = linspace(0, 0.4, 21); % 1 to -1 m in 100 intervals
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

%% Plot the Quiver
%PlotQuiver(X,Y,Bx,By)

%% Extract the values on the positive Y axis only
%yaxisExtracted = ExtractYaxisMiddleAnalysis(X,Y,Bx,By);

%% Determine Inital Position based on input

% Combine Bx and By for analysis, X on top of Y
BxByCat = cat(3,Bx,By);
% Manually set the middle value to the one above it so it doesn't stay at 0
BxByCat(1,ceil(length(Bx)/2),:) = BxByCat(2,ceil(length(Bx)/2),:);

% Define some arbitrary values within the DipoleModel grid
format shortE;
BxInit = -2.3794e+11;
ByInit = -33.1775e+10;

% Determine the index of the initial Bx and By values from a Y(Forward)
% facing device
[indexXInit, indexYInit] = MagIndInit(BxByCat,BxInit, ByInit);

% WORKS!

% initialize Next function
% we will rename these just for funs (but really to keep it modular for
% later)
indexXFirst = indexXInit;
indexYFirst = indexYInit;

[indexXNext,indexYNext] = MagIndNext(BxByCat,indexXFirst,indexYFirst);
