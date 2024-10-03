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
format shortE;
initialBx = -2.3794e+11;
initialBy = -33.1775e+10;
initialBxByCat = cat(3,initialBx,initialBy);

BxByCat = cat(3,Bx,By);
% Manually set the middle value to the one above it so it doesn't stay at 0
BxByCat(1,ceil(length(Bx)/2),:) = BxByCat(2,ceil(length(Bx)/2),:);

% Subtract the initial position from all other positions
initialDifferences = BxByCat - repmat(initialBxByCat, size(BxByCat,1), size(BxByCat,2));

% Find the magnitude of the differences so that the minimum may be found
differencesMagnitudes = sqrt(initialDifferences(:,:,1).^2 + initialDifferences(:,:,2).^2);

% Find the minimum
minDifference = min(min(differencesMagnitudes));
[indexY, indexX] = find(differencesMagnitudes == minDifference);
fprintf("the initial coordinate is: [%d,%d]\n",indexX,indexY) 
%BxByCat(minRowIndexY,minColIndexX,:);

% WORKS!
% Now we need to index the 8 surrounding cells to find which is closest to
% a second coordinate from a first coordinate
firstBx = -2.3794e+11;
firstBy = -33.1775e+10;
IndexX1 = indexX;
IndexY1 = indexY;
firstBxByCat = cat(3,firstBx,firstBy);

secondBx = -2.3794e+11;
secondBy = -3.1775e+10;
secondBxByCat = cat(3,secondBx,secondBy);

% create an array of each of the surrounding indeces for comparison
surroundingArray = zeros(3,3,2);
for i = 1:9
    switch i
        case 1
            try 
                indexY2 = IndexY1 + 1;
                indexX2 = indexX1 + -1;
                surroundingArray(1,1,:) = BxByCat(indexY2,indexX2,:);
            catch
                fprintf("the coordinate: [%d,%d] does not exist\n",indexX2,indexxY2) 
            end

        case 2
            try 
                indexY2 = IndexY1 + 1;
                indexX2 = indexX1 + 0;
                surroundingArray(1,2,:) = BxByCat(indexY2,indexX2,:);
            catch
                fprintf("the coordinate: [%d,%d] does not exist\n",indexX2,indexxY2) 
            end

        case 3
            try 
                indexY2 = IndexY1 + 1;
                indexX2 = indexX1 + 1;
                surroundingArray(1,3,:) = BxByCat(indexY2,indexX2,:);
            catch
                fprintf("the coordinate: [%d,%d] does not exist\n",indexX2,indexxY2) 
            end

        case 4
            try 
                indexY2 = indexY1 + 0;
                indexX2 = indexX1 + -1;
                surroundingArray(2,1,:) = BxByCat(indexY2,indexX2,:);
            catch
                fprintf("the coordinate: [%d,%d] does not exist\n",indexX2,indexxY2) 
            end

        case 5
            surroundingArray(2,2,:) = firstBxByCat;

        case 6
            try 
                indexY2 = indexY1 + 0;
                indexX2 = indexX1 + 1;
                surroundingArray(2,3,:) = BxByCat(indexY2,indexX2,:);
            catch
                fprintf("the coordinate: [%d,%d] does not exist\n",indexX2,indexxY2) 
            end

        case 7 
            try 
                indexY2 = indexY1 + -1;
                indexX2 = indexX1 + -1;
                surroundingArray(3,1,:) = BxByCat(indexY2,indexX2,:);
            catch
                fprintf("the coordinate: [%d,%d] does not exist\n",indexX2,indexxY2) 
            end

        case 8
            try 
                indexY2 = indexY1 + -1;
                indexX2 = indexX1 + 0;
                surroundingArray(3,2,:) = BxByCat(indexY2,indexX2,:);
            catch
                fprintf("the coordinate: [%d,%d] does not exist\n",indexX2,indexxY2) 
            end

        case 9
            try 
                indexY2 = indexY1 + -1;
                indexX2 = indexX1 + 1;
                surroundingArray(3,3,:) = BxByCat(indexY2,indexX2,:);
            catch
                fprintf("the coordinate: [%d,%d] does not exist\n",indexX2,indexxY2) 
            end
    end
end

% find the differences between the Bx and By elements from the first and
% second 
surroudingDiff = surroundingArray - repmat(secondBxByCat, size(surroundingArray,1), size(surroundingArray,2));

surroundingDiffMag = sqrt(surroudingDiff(:,:,1).^2 + surroudingDiff(:,:,2).^2);

surroundDiffMagMin = min(min(surroundingDiffMag));
[tempInd] = find(surroundingDiffMag == surroundDiffMagMin);

% Set index2 based on which surrounding array element (or the original) was
% closest to the new secondBy and secondBx
for i = 1:9
    switch tempInd
        case 1
                indexY2 = indexY1 + 1;
                indexX2 = indexX1 + -1;

        case 2
                indexY2 = indexY1 + 1;
                indexX2 = indexX1 + 0;

        case 3
                indexY2 = indexY1 + 1;
                indexX2 = indexX1 + 1;

        case 4
                indexY2 = indexY1 + 0;
                indexX2 = indexX1 + -1;

        case 5
                indexY2 = indexY1 + 0;
                indexX2 = indexX1 + 0;

        case 6
                indexY2 = indexY1 + 0;
                indexX2 = indexX1 + 1;

        case 7 
                indexY2 = indexY1 + -1;
                indexX2 = indexX1 + -1;

        case 8
                indexY2 = indexY1 + -1;
                indexX2 = indexX1 + 0;

        case 9
                indexY2 = indexY1 + -1;
                indexX2 = indexX1 + 1;
    end
end

% print the new coord
fprintf("the next coordinate is: [%d,%d]\n",indexX2,indexY2) 