%% DipoleModel.m computation section, this is the main file
% Last modified 10/1/24 1:34 AM - Noah Nguyen UCD ECE Jlab
    clc; format short;

%% Set up prameters

    % amplification to possibly drown out low amplitude noise. Did not work.
        amplify = 1; 
    % interpolation factor set
        interpolFact = 70;
    % the absolute value boundary from one side for the surrounding values check.
        absBound = 2; 
        surroundRange = -absBound:absBound;
    % cuttoff at cutoff/interpolationFactor/2
        cutoff = 1; 
    % the order of the applied butterorth filter
        filterOrd = 1; 

    % Assuming N42 and a 3" magnet of diameter 0.5mm and thickness of 75mm
        %dpMoment = 9.98; % A*m^2 
        dpMoment =  9.98 * 9.2050494347e15 * -1.426593504e-11 * 1.00796677702 * amplify; % A*m^2  * conversion to IMU data number to match OneDimMagnitude.m

    % Set up physical boundary size of grid
        Xstart = -0.4;
        Xend = 0.4;
        Ystart = 0;
        Yend = 0.4;

    % number of samples per axis
        n = 2001;

% Call the function to create the magnetic vector field
        [X,Y,BxByCat] = DipoleMake(Xstart,Xend,Ystart,Yend,n,dpMoment);


%% Plot the Quiver
    %PlotQuiver(X,Y,Bx,By)

%% Extract the values on the positive Y axis only
    %yaxisExtracted = ExtractYaxisMiddleAnalysis(X,Y,Bx,By);

%% Try the location tracking with moving the device in a straight line

% % Read the table from the 1 Dim data
T = readtable("Magnetometer Magnitude to Distance - 1 Dimension - 8_23_24 9_54 PM - 2D CSV 1D Movement for MATLAB.csv");
%T = readtable("test straight line.csv");

% Save the original scatter data for the X and Y coordinates.
BxRead = interp(T.Bx * amplify,interpolFact);
ByRead = interp(T.By * amplify,interpolFact);


% Apply the filter to reduce high frequency noise.
BxRead = ApplyButterB(BxRead,interpolFact,cutoff,filterOrd);
ByRead = ApplyButterB(ByRead,interpolFact,cutoff,filterOrd);

% Determine the first point based on vector subtractions considering the
% device starts in y direction
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
    % display which iteration we are in
    %disp(i)
    % Find and print the next coordinate decided by a surrounding magnitude
    % match
    [indexX2,indexY2] = MagIndNext3(BxByCat,indexX1,indexY1,BxRead(i+1),ByRead(i+1),surroundRange);

    % save it to the next row of the index log
    indexLog(i+1,:) = [indexX2,indexY2];

    % the current X2 and Y2 are now the X1 and Y1 for the next itteration
    indexX1 = indexX2;
    indexY1 = indexY2;
end

%% Plot the path and add a legend for each run

% Create a plot of the path taken by the algorithm
hold on;
% Construct the DisplayName string with multiple lines
% Construct the DisplayName string with multiple lines using sprintf
legendString = sprintf('Interpolation = %d\nSearch Grid Size = %dx%d\nButterworth: Order = %d, Cutoff (Hz) = %d', ...
                        interpolFact, 1+absBound*2, 1+absBound*2, filterOrd, cutoff);

% Assuming indexLog(:,1) and indexLog(:,2) are your data for the scatter plot
plot(indexLog(:,1), indexLog(:,2), 'DisplayName', legendString);
axis equal;

% Set x and y limits as specified
xlim([0 n]);
ylim([0 n]);

% Add labels and title
xlabel('X-axis (Index)');
ylabel('Y-axis (Index)');
title('Path and Position Scatter Plot With Similar Min Randomization');

% Add the legend to explain the colors/markers
legend('show');
