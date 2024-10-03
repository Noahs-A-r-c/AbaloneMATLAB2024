% MagIndNext.m
% A function to take in the next By and Bx values along with a starting set
% of coordinates. The function outputs the next set of coordinates that
% matches closest to the next Bx and By magnitude, even if that is the same
% coordinate that it started at

function [indexXNext, indexYNext] = MagIndNext(BxByCat,indexXFirst,indexYFirst)

% Now we need to index the 8 surrounding cells to find which is closest to
% a second coordinate from a first coordinate
firstBx = BxByCat(indexYFirst,indexXFirst,1);
firstBy = BxByCat(indexYFirst,indexXFirst,2);
firstBxByCat = cat(3,firstBx,firstBy);

secondBx = -2.3794e+11;
secondBy = -3.1775e+10;
secondBxByCat = cat(3,secondBx,secondBy);

% create an array of each of the surrounding indeces for comparison
surroundingArray = zeros(3,3,2);
for i = 1:9
    switch i
        case 1 
                indexYNext = indexYFirst + 1;
                indexXNext = indexXFirst + -1;
            try
                surroundingArray(1,1,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 2 
                indexYNext = indexYFirst + 1;
                indexXNext = indexXFirst + 0;
            try
                surroundingArray(1,2,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 3
                indexYNext = indexYFirst + 1;
                indexXNext = indexXFirst + 1;
            try
                surroundingArray(1,3,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 4
                indexYNext = indexYFirst + 0;
                indexXNext = indexXFirst + -1;
            try 
                surroundingArray(2,1,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 5
            surroundingArray(2,2,:) = firstBxByCat;

        case 6
                indexYNext = indexYFirst + 0;
                indexXNext = indexXFirst + 1;
            try
                surroundingArray(2,3,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 7 
                indexYNext = indexYFirst + -1;
                indexXNext = indexXFirst + -1;
            try
                surroundingArray(3,1,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 8 
                indexYNext = indexYFirst + -1;
                indexXNext = indexXFirst + 0;
            try
                surroundingArray(3,2,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
            end

        case 9
                indexYNext = indexYFirst + -1;
                indexXNext = indexXFirst + 1;
            try
                surroundingArray(3,3,:) = BxByCat(indexYNext,indexXNext,:);
            catch
                fprintf("the coordinate: [%d,%d] does not exist\n",indexXNext,indexYNext) 
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
magDiffArray = surroundingArrayMag - secondBxByCatMagArray;

magDiffArrayMin = min(min(magDiffArray));
[tempInd] = find(magDiffArray == magDiffArrayMin);

% Set index2 based on which surrounding array element (or the original) was
% closest to the new secondBy and secondBx
for i = 1:9
    switch tempInd
        case 1
                indexYNext = indexYFirst + 1;
                indexXNext = indexXFirst + -1;

        case 2
                indexYNext = indexYFirst + 1;
                indexXNext = indexXFirst + 0;

        case 3
                indexYNext = indexYFirst + 1;
                indexXNext = indexXFirst + 1;

        case 4
                indexYNext = indexYFirst + 0;
                indexXNext = indexXFirst + -1;

        case 5
                indexYNext = indexYFirst + 0;
                indexXNext = indexXFirst + 0;

        case 6
                indexYNext = indexYFirst + 0;
                indexXNext = indexXFirst + 1;

        case 7 
                indexYNext = indexYFirst + -1;
                indexXNext = indexXFirst + -1;

        case 8
                indexYNext = indexYFirst + -1;
                indexXNext = indexXFirst + 0;

        case 9
                indexYNext = indexYFirst + -1;
                indexXNext = indexXFirst + 1;
    end
end

% print the new coord
fprintf("the next coordinate is: [%d,%d]\n",indexXNext,indexYNext) 
end

