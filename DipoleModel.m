%% DipoleModel.m computation section
% Last modified 10/1/24 1:34 AM - Noah Nguyen UCD ECE Jlab
    clc; format short;

% Set up prameters
    % Assuming N42 and a 3" magnet of diameter 0.5mm and thickness of 75mm
    %dpMoment = 9.98; % A*m^2 
    amplify = 1; % amplification to possibly drown out low amplitude noise. Did not work.
    dpMoment =  9.98 * 9.2050494347e15 * -1.426593504e-11 * 1.00796677702 * amplify; % A*m^2  * conversion to IMU data number to match OneDimMagnitude.m

    % Set up physical boundary size of grid
    Xstart = -0.4;
    Xend = 0.4;
    Ystart = 0;
    Yend = 0.4;

    % number of samples per axis
    n = 401;

% Call the function to create the magnetic vector field
    [X,Y,Bx,By] = DipoleMake(Xstart,Xend,Ystart,Yend,n,dpMoment);

%% Plot the Quiver
    %PlotQuiver(X,Y,Bx,By)

%% Extract the values on the positive Y axis only
    %yaxisExtracted = ExtractYaxisMiddleAnalysis(X,Y,Bx,By);

%% Determine Inital Position based on input
% 
% Combine Bx and By for analysis, X on top of Y
    BxByCat = cat(3,Bx,By);
% Manually set the middle value to the one above it so it doesn't stay at 0
    BxByCat(1,ceil(length(Bx)/2),:) = BxByCat(2,ceil(length(Bx)/2),:);

%% Try the location tracking with moving the device in a straight line

% Read the table from the 1 Dim data
T = readtable("Magnetometer Magnitude to Distance - 1 Dimension - 8_23_24 9_54 PM - 2D CSV 1D Movement for MATLAB.csv");
%T = readtable("test straight line.csv");
% Save the original scatter data for the X and Y coordinates.
interpolFact = 3;
BxRead = interp(T.Bx * amplify,interpolFact);
ByRead = interp(T.By * amplify,interpolFact);

[indexXInit, indexYInit] = MagIndInit(BxByCat,BxRead(1),ByRead(1));

% Initialize a zeros block to grab the coordinate for each table entry
indexLog = zeros(length(BxRead),2);
% The first coordinate is 
indexLog(1,:) = [indexXInit, indexYInit];

% kick of the itterations by assigning the first indeces with the init
% indices
indexX1 = indexXInit;
indexY1 = indexYInit;
for i = 1:length(BxRead)-1
    disp(i)
    % Find and print the next coordinate decided by a surrounding magnitude
    % match
    [indexX2,indexY2] = MagIndNext(BxByCat,indexX1,indexY1,BxRead(i+1),ByRead(i+1));

    % save it to the next row of the index log
    indexLog(i+1,:) = [indexX2,indexY2];

    % the current X2 and Y2 are now the X1 and Y1 for the next itteration
    indexX1 = indexX2;
    indexY1 = indexY2;
end

%figure 1;
plot(indexLog(:,1),indexLog(:,2));
hold on;
scatter(indexLog(:,1),indexLog(:,2));
axis equal;
xlim([0 n]);
ylim([0 n]);
%hold off;