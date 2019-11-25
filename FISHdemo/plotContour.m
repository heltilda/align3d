pathname = '/Users/brianross/Desktop/align3d/experimentalDemo/results/FISH10spot/3color/chr4/';
groupName = '4x1';
fileNo = 'ch12';

axesToPlot = [ 1 2 3 ];
interpMode = 'pchip';   %pchip, linear
numMidpoints = 100;
numEndpoints = 3;
colors = ['b'; 'g'; 'r'];
spotHeight = 150;
viewElevation = 45;

loci = load([pathname 'loci' groupName '.txt']);
spots = load([pathname 'spots' groupName fileNo '.txt']);
% l2s = load([pathname 'l2s' groupName fileNo '.txt']);
l2s = [ 4 5 0 6 1 2 0 3 7 8 ];


thisfig = figure(1);
clf, hold on

sampleAt = (-numEndpoints*2:numMidpoints+numEndpoints*2-1) * (loci(end, 1)-loci(1, 1))/(numMidpoints-1) + loci(1, 1);
ptsToPlot = interp1(loci(l2s ~= 0, 1), spots(l2s(l2s ~= 0), axesToPlot), sampleAt, interpMode);
midpoints = (numEndpoints+1):(numMidpoints+3*numEndpoints);
leftendpoints = (1:numEndpoints+1);
rightendpoints = numMidpoints+3*numEndpoints+(0:numEndpoints);
zBase = min(ptsToPlot(:, 3));

%plot3(
plot3(ptsToPlot(midpoints, 1), ptsToPlot(midpoints, 2), ptsToPlot(midpoints, 3)*0+zBase, 'Color', .8*[1 1 1], 'LineWidth', 13);
plot3(ptsToPlot(leftendpoints, 1), ptsToPlot(leftendpoints, 2), ptsToPlot(leftendpoints, 3), 'Color', .7*[1 1 1], 'LineWidth', 8);
plot3(ptsToPlot(rightendpoints, 1), ptsToPlot(rightendpoints, 2), ptsToPlot(rightendpoints, 3), 'Color', .7*[1 1 1], 'LineWidth', 8);
plot3(ptsToPlot(midpoints, 1), ptsToPlot(midpoints, 2), ptsToPlot(midpoints, 3), 'Color', .3*[1 1 1], 'LineWidth', 10);
for loopSpot = 1:size(spots, 1)
    plot3(spots(loopSpot, 1)+[0 0], spots(loopSpot, 2)+[0 0], [ zBase spots(loopSpot, 3)+spotHeight ], 'k')
    plot3(spots(loopSpot, 1), spots(loopSpot, 2), spots(loopSpot, 3)+spotHeight, 'ok', 'MarkerFaceColor', colors(spots(loopSpot, 7)+1))
end
axis equal, axis off

set(thisfig, 'Color', 'w')
set(thisfig, 'Position', [50, 200, 200, 200])
view(0, viewElevation)