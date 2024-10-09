% MagIndNext2 is the next itteration of the next index search function
% which allows for a changing grid size with magnitude based search.
function [indexXNext, indexYNext] = MagIndNext2(BxByCat, indexXFirst, indexYFirst, secondBx, secondBy,surroundRange)

    [shiftY, shiftX] = meshgrid(surroundRange, surroundRange);
    shifts = [shiftY(:),shiftX(:)];

    surroundingArray = Inf(size(shifts,1), 1, 2);
    
    [rows, cols, ~] = size(BxByCat);

    for i = 1:size(shifts, 1)
        tempY = indexYFirst + shifts(i, 1);
        tempX = indexXFirst + shifts(i, 2);

        if tempY >= 1 && tempY <= rows && tempX >= 1 && tempX <= cols
        surroundingArray(i,1,:) = BxByCat(tempY, tempX, :);
        end 
    end

    surroundingArrayMag = sqrt(surroundingArray(:,:,1).^2 + surroundingArray(:,:,2).^2);
    secondBxByCatMag = sqrt(secondBx^2 + secondBy^2);
    magDiffArray = abs(surroundingArrayMag - secondBxByCatMag);
    
    % Find all indices where the minimum difference occurs
    minValue = min(magDiffArray(:));
    minIndices = find(magDiffArray(:) == minValue);
    
    % Randomly select one index among the minima
    randomMinIndex = minIndices(randi(length(minIndices)));
    
    % Use the corresponding shift to update the next indices
    indexYNext = indexYFirst + shifts(randomMinIndex, 1);
    indexXNext = indexXFirst + shifts(randomMinIndex, 2);


    if indexXNext ~= indexXFirst || indexYNext ~= indexYFirst
    % print the new coord and corresponding distance if it has changed
    fprintf("the next coordinate is: [%d,%d]\n",indexXNext,indexYNext) 
    
    % distX = 200 - (indexXNext/length(BxByCat)) * 400;
    % distY = (indexYNext/length(BxByCat)) * 400;
    % % display the current coordinates and computed distances
    % fprintf("next coordinate: [%d,%d]   next distance: [%.2f mm,%.2f mm]\n",indexXNext,indexYNext,distX,distY)
    end
end
