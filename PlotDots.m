pathname = '/Users/brianross/Desktop/align3d/experimentalDemo/results/FISH10spot/3color/chr4/';
%pathname = '/Users/brianross/Documents/GitHub/BXC-ORCA-data/WT2kb/';
lociFileNum = '4x1';
fileNum = [ lociFileNum 'ch12' ];
% pathname = '/Users/brianross/Desktop/align3d/seriesExpansion/results/small100/';
% fileNum = '1';
xax = 1;
yax = 2;
zax = 3;
units = 'nm';

plotContour = true;
plot_NoMSE = true;
plotStraight = false;
plotRealAnswer = false;
interpolateLoci = false;
plot_endpoints = false;
mark_FPs = false;
projectToPlane = false;

thisfig = figure; clf, hold on

loci = load([pathname, 'loci' lociFileNum '.txt']);
numColors = max(loci(:, 2))+1;
if ~plotStraight
    spots = load([pathname, 'spots' fileNum '.txt']);
end
if mark_FPs
    s2l = load([pathname, 's2l' fileNum '.txt']);
end

if interpolateLoci
    l2s = load([pathname, 'l2s' fileNum '.txt']);
    spots = spots(l2s(l2s~=0), :);
    RealChain = spots(:, 1:3);
    noMSE = RealChain;
    fnLoci = zeros(0, 7);
elseif plotRealAnswer && ~plotStraight
    RealChain = load([pathname, 'RealChain' fileNum '.txt']);
    fnLoci = load([pathname, 'FNSpots' fileNum '.txt']);
    if plot_NoMSE == true
        noMSE = load([pathname, 'noMSE' fileNum '.txt']);
    end
end

colors = ['b'; 'g'; 'r'; 'c'; 'm'; 'y'; 'k' ];
%colors = ['g'; 'r'; 'b'; 'c'; 'm'; 'y'; 'k' ];
if numColors > length(colors)
    mapSize = size(colormap, 1);
    cMap = hsv;
    colors = cMap(floor((1:numColors)*mapSize/numColors), :);
end

if plotStraight
    if plotRealAnswer
        L = 0;
        lastXYZ = RealChain(1, :);
        for c1 = 1:size(RealChain,1)
            L = L + sqrt(sum( (RealChain(c1, :) - lastXYZ).^2 ));
            lastXYZ = RealChain(c1, :);
            RealChain(c1, :) = 0;
            RealChain(c1, xax) = L;
        end
    else
        RealChain = zeros(2, 3);
        RealChain(:, xax) = [ min(loci(:, 1)); max(loci(:, 1)) ];
        plotRealAnswer = true;
    end
    
    spots = zeros(size(loci,1), size(spots,2));
    for c1 = 1:size(loci,1)
        spots(c1, xax) = loci(c1, 1);
        spots(c1, 7) = loci(c1, 2);
    end
end

axs = ['x' 'y' 'z'];

if plotRealAnswer
    ChainEnd = [ RealChain(1, :); RealChain(size(RealChain, 1), :); ];
    if plotContour == true
       if plot_endpoints == true
           plot3(ChainEnd(1, xax), ChainEnd(1, yax), ChainEnd(1, zax), ...
              'x', 'MarkerSize', 18, 'LineWidth', 3)
           plot3(ChainEnd(2, xax), ChainEnd(2, yax), ChainEnd(2, zax), ...
              'o', 'MarkerSize', 18, 'LineWidth', 3)
       end
       plot3(RealChain(:,xax), RealChain(:,yax), RealChain(:,zax), ...
           'k', 'LineWidth', 1)
    end
end

zplane = 1.3*min(spots(:, zax)) - 0.3*max(spots(:, zax));
colorOffset = 1-min(spots(:, 7));
for c1 = [1:size(spots, 1)]
    dot_symbol = 'o';
    dot_edge_color = 'k';
    dot_size = 7;
    q = colors(spots(c1,7)+colorOffset, :);
    if plotRealAnswer && mark_FPs
        if s2l(c1) == 0
            dot_size = 6;
            dot_symbol = 'd';
        end
    end
    plot3(spots(c1,xax), spots(c1,yax), spots(c1,zax), dot_symbol,...
                'MarkerFaceColor', q,...
                'MarkerEdgeColor', dot_edge_color,...
                'MarkerSize', dot_size)
    if projectToPlane
        plot3(spots(c1,xax)+[0 0], spots(c1,yax)+[0 0], [spots(c1,zax) zplane], 'k')
    end
    if plotRealAnswer && mark_FPs
        if s2l(c1) == 0
            plot3(spots(c1,xax), spots(c1,yax), spots(c1,zax), '+',...
                'MarkerEdgeColor', colors(spots(c1,7)+colorOffset, :),...
                'MarkerFaceColor', 'k',...
                'LineWidth', 1,...
                'MarkerSize', 12)
        end
    end
    if plotRealAnswer && plot_NoMSE && plotContour && ~plotStraight
        plot3([spots(c1,xax) noMSE(c1,xax)], [spots(c1,yax) noMSE(c1,yax)],...
                    [spots(c1,zax) noMSE(c1,zax)], 'b');
    end
end
if plotRealAnswer && ~plotStraight
    for c1 = [1:size(fnLoci, 1)]
        plot3(fnLoci(c1,xax), fnLoci(c1,yax), fnLoci(c1,zax), 'o',...
                    'MarkerEdgeColor', colors(fnLoci(c1,7)+colorOffset, :),...
                    'MarkerFaceColor', 'w',...
                    'MarkerSize', 7)
    end
end

if plotStraight == false
    xlabel([axs(xax) ' (', units, ')'])
    ylabel([axs(yax) ' (', units, ')'])
    zlabel([axs(zax) ' (', units, ')'])
else
    xlabel('Contour position (Mbp)')
end

if plotStraight == false
    axis equal
end

if projectToPlane
    axs = gca;
    for c1 = 1:length(axs.XTick)
        plot3(axs.XTick(c1)*[1 1], [axs.YTick(1) axs.YTick(end)], zplane*[1 1], 'k:')
    end
    for c1 = 1:length(axs.YTick)
        plot3([axs.XTick(1) axs.XTick(end)], axs.YTick(c1)*[1 1], zplane*[1 1], 'k:')
    end
end

hold off
set(thisfig, 'Color', 'w')

allDists2 = zeros(size(spots, 1));
for c1 = 2:size(spots, 1)
    for c2 = 1:c1
        allDists2(c1, c2) = sum((spots(c1, 1:3)-spots(c2, 1:3)).^2);
    end
end
allDists2 = allDists2(allDists2>0);

disp(['indistinguishable fraction:  ', num2str(mean(allDists2<0.1^2)) ...
    ' to ', num2str(mean(allDists2<0.2^2))])