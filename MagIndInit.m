% MagIndInit.m 
% A function to take in the beginning two coordinates for a Y (forward)
% facing device and spit out which index this applies to. The next step
% would be translating that to real distances

function [indexXInit,indexYInit] = MagIndInit(BxByCat,initialBx, initialBy)
    initialBxByCat = cat(3,initialBx,initialBy);
    
    % Subtract the initial position from all other positions
    initialDifferences = BxByCat - repmat(initialBxByCat, size(BxByCat,1), size(BxByCat,2));
    
    % Find the magnitude of the differences so that the minimum may be found
    differencesMagnitudes = sqrt(initialDifferences(:,:,1).^2 + initialDifferences(:,:,2).^2);
    
    % Find the minimum
    minDifference = min(min(differencesMagnitudes));
    [indexYInit, indexXInit] = find(differencesMagnitudes == minDifference);
    fprintf("the initial coordinate is: [%d,%d]\n",indexXInit,indexYInit) 
    %BxByCat(minRowIndexY,minColIndexX,:);
end