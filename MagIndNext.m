% MagIndNext.m
% A function to take in the next By and Bx values along with a starting set
% of coordinates. The function outputs the next set of coordinates that
% matches closest to the next Bx and By magnitude, even if that is the same
% coordinate that it started at

function [indexXNext, indexYNext] = MagIndNext(BxByCat,indexXFirst,indexYFirst,secondBx,secondBy)
format shortE;
% Now we need to index the 8 surrounding cells to find which is closest to
% a second coordinate from a first coordinate
firstBx = BxByCat(indexYFirst,indexXFirst,1);
firstBy = BxByCat(indexYFirst,indexXFirst,2);
firstBxByCat = cat(3,firstBx,firstBy);

secondBxByCat = cat(3,secondBx,secondBy);

% Create an array of each of the surrounding indices for comparison
surroundingArray = zeros(5,5,2);
surroundingArray(surroundingArray == 0) = Inf;

% Loop through the 5x5 surrounding cells
for i = 1:25
    switch i
        case 1
            indexYNext = indexYFirst + 2;
            indexXNext = indexXFirst - 2;
            try
                surroundingArray(1,1,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 2
            indexYNext = indexYFirst + 2;
            indexXNext = indexXFirst - 1;
            try
                surroundingArray(1,2,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 3
            indexYNext = indexYFirst + 2;
            indexXNext = indexXFirst + 0;
            try
                surroundingArray(1,3,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 4
            indexYNext = indexYFirst + 2;
            indexXNext = indexXFirst + 1;
            try
                surroundingArray(1,4,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 5
            indexYNext = indexYFirst + 2;
            indexXNext = indexXFirst + 2;
            try
                surroundingArray(1,5,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        % Similarly for row 2
        case 6
            indexYNext = indexYFirst + 1;
            indexXNext = indexXFirst - 2;
            try
                surroundingArray(2,1,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 7
            indexYNext = indexYFirst + 1;
            indexXNext = indexXFirst - 1;
            try
                surroundingArray(2,2,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 8
            indexYNext = indexYFirst + 1;
            indexXNext = indexXFirst + 0;
            try
                surroundingArray(2,3,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 9
            indexYNext = indexYFirst + 1;
            indexXNext = indexXFirst + 1;
            try
                surroundingArray(2,4,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 10
            indexYNext = indexYFirst + 1;
            indexXNext = indexXFirst + 2;
            try
                surroundingArray(2,5,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        % Row 3, where the center (3,3) is the current position
        case 11
            indexYNext = indexYFirst + 0;
            indexXNext = indexXFirst - 2;
            try
                surroundingArray(3,1,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 12
            indexYNext = indexYFirst + 0;
            indexXNext = indexXFirst - 1;
            try
                surroundingArray(3,2,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 13
            % Center point is the current position (no change in indices)
            surroundingArray(3,3,:) = firstBxByCat;

        case 14
            indexYNext = indexYFirst + 0;
            indexXNext = indexXFirst + 1;
            try
                surroundingArray(3,4,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 15
            indexYNext = indexYFirst + 0;
            indexXNext = indexXFirst + 2;
            try
                surroundingArray(3,5,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        % Rows 4 and 5 similarly expanded
        case 16
            indexYNext = indexYFirst - 1;
            indexXNext = indexXFirst - 2;
            try
                surroundingArray(4,1,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 17
            indexYNext = indexYFirst - 1;
            indexXNext = indexXFirst - 1;
            try
                surroundingArray(4,2,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 18
            indexYNext = indexYFirst - 1;
            indexXNext = indexXFirst + 0;
            try
                surroundingArray(4,3,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 19
            indexYNext = indexYFirst - 1;
            indexXNext = indexXFirst + 1;
            try
                surroundingArray(4,4,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 20
            indexYNext = indexYFirst - 1;
            indexXNext = indexXFirst + 2;
            try
                surroundingArray(4,5,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        % Last row (5th row)
        case 21
            indexYNext = indexYFirst - 2;
            indexXNext = indexXFirst - 2;
            try
                surroundingArray(5,1,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 22
            indexYNext = indexYFirst - 2;
            indexXNext = indexXFirst - 1;
            try
                surroundingArray(5,2,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 23
            indexYNext = indexYFirst - 2;
            indexXNext = indexXFirst + 0;
            try
                surroundingArray(5,3,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 24
            indexYNext = indexYFirst - 2;
            indexXNext = indexXFirst + 1;
            try
                surroundingArray(5,4,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 25
            indexYNext = indexYFirst - 2;
            indexXNext = indexXFirst + 2;
            try
                surroundingArray(5,5,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                %fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end
    end
end


% find the differences between the first Bx By elements' magnitude and the
% magnitude of the surrounding vectors. We use magnitude here instead of
% vector subtraction because the device may rotate in the 2D plane. This
% still needs to be tested for accuracy.
surroundingArrayMag = sqrt(surroundingArray(:,:,1).^2 + surroundingArray(:,:,2).^2);

secondBxByCatMag = sqrt(secondBxByCat(:,:,1).^2 + secondBxByCat(:,:,2).^2);
secondBxByCatMagArray = repmat(secondBxByCatMag,size(surroundingArray,1), size(surroundingArray,2));
magDiffArray = abs(surroundingArrayMag - secondBxByCatMagArray);

% Determine the minimum measureable difference
magDiffArrayMin = min(min(magDiffArray));
[tempInd] = find(magDiffArray == magDiffArrayMin,1);

% % Assume that the the device stays straight instead
%     % Subtract the initial position from all other positions
%     nextDiffs = surroundingArray - repmat(secondBxByCat, size(surroundingArray,1), size(surroundingArray,2));
% 
%     % Find the magnitude of the differences so that the minimum may be found
%     nextDiffsMags = sqrt(nextDiffs(:,:,1).^2 + nextDiffs(:,:,2).^2);
% 
%     % Find the minimum
%     minDifference = min(min(nextDiffsMags));
%     [tempInd] = find(nextDiffsMags == minDifference,1);

format short g
% % Display the computed arrays
% disp('Surrounding Array Magnitudes:');
% disp(surroundingArrayMag);
% 
% disp('Second BxByCat Magnitude Array:');
% disp(secondBxByCatMagArray);
% 
% disp('Magnitude Difference Array:');
% disp(magDiffArray);



% Set index2 based on which surrounding array element (or the original) was
% closest to the new secondBy and secondBx
for i = 1:9
    switch tempInd
        case 1
                indexYNext = indexYFirst + 1;
                indexXNext = indexXFirst + -1;

        case 4
                indexYNext = indexYFirst + 1;
                indexXNext = indexXFirst + 0;

        case 7
                indexYNext = indexYFirst + 1;
                indexXNext = indexXFirst + 1;

        case 2
                indexYNext = indexYFirst + 0;
                indexXNext = indexXFirst + -1;

        case 5
                indexYNext = indexYFirst + 0;
                indexXNext = indexXFirst + 0;

        case 8
                indexYNext = indexYFirst + 0;
                indexXNext = indexXFirst + 1;

        case 3 
                indexYNext = indexYFirst + -1;
                indexXNext = indexXFirst + -1;

        case 6
                indexYNext = indexYFirst + -1;
                indexXNext = indexXFirst + 0;

        case 9
                indexYNext = indexYFirst + -1;
                indexXNext = indexXFirst + 1;
    end
end

if indexXNext ~= indexXFirst || indexYNext ~= indexYFirst
    % print the new coord and corresponding distance if it has changed
    fprintf("the next coordinate is: [%d,%d]\n",indexXNext,indexYNext) 
    
    % distX = 200 - (indexXNext/length(BxByCat)) * 400;
    % distY = (indexYNext/length(BxByCat)) * 400;
    % % display the current coordinates and computed distances
    % fprintf("next coordinate: [%d,%d]   next distance: [%.2f mm,%.2f mm]\n",indexXNext,indexYNext,distX,distY)
end

end

