dir = '/Users/brianross/Documents/GitHub/BXC-ORCA-data/WT3kb/';
inputFile = '/ORCA_WT_3kb_rpt1_reduced.csv';
minSpotFraction = 0.9;

sourceData = csvread([ dir inputFile ], 1, 0);
minNoSpots = ceil(minSpotFraction*max(sourceData(:, 7)));

newCellRows = [ 0 find(diff(sourceData(:, end)) ~= 0)' size(sourceData, 1) ] + 1;
numCells = length(newCellRows) - 1;

numSavedCells = 0;
allspots = zeros(0, 8);
for loopCell = 1:numCells
    cellIdx = newCellRows(loopCell):(newCellRows(loopCell+1)-1);
    numRows = find(diff(sourceData(cellIdx, 7)) < 1, 1, 'first');
    if isempty(numRows)
        numRows = length(cellIdx);
    end
    
    if numRows >= minNoSpots
        numSavedCells = numSavedCells + 1;
        savedIdx = newCellRows(loopCell) + (0:(numRows-1));
        allspots = [ allspots; [ numSavedCells*ones(length(savedIdx), 1), ...
            sourceData(savedIdx, 7), ...
            sourceData(savedIdx, 1:6) / 1e3 ] ];     % convert nm --> um
    end
end

dlmwrite([ dir 'allspots.txt' ], allspots, '\t');