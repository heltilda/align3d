%pathname = '/Users/brianross/Desktop/align3d/FISH_Demo/superres_072115/';
% pathname = '/Users/brianross/Desktop/align3d/Results/';
%dotsFile = [ pathname 'bross_072115_100x_02_R3D_D3D_TIFFS/xyzPts_ch2_calibrated.txt' ];
%dotsFile = [ pathname 'bross_072115_100x_02_R3D_D3D_TIFFS/allSpotsDbg.txt' ];
%pathname = '/Users/brianross/Desktop/align3d/Results/HundredChains/';
%pathname = '/Users/brianross/Desktop/align3d/Results/BigChain/';
pathname = '/Users/brianross/Desktop/align3d/Results/large4/';
%pathname = '/Users/brianross/Desktop/GEM/conformations/chr18/toPlot/';

plot_contour = true;
plot_NoMSE = true;
plot_straight = false;
plotRealAnswer = true;
plot_endpoints = false;
mark_FPs = false;
projectToPlane = false;
doAutoRainbow = true;

thisfig = figure; clf, hold on

xax = 3;
yax = 2;
zax = 1;
fileNum = '4';
units = 'nm';
if plotRealAnswer
    spots = load([pathname, 'spots' fileNum '.txt']);
    FNSpots = load([pathname, 'FNSpots' fileNum '.txt']);
    if plot_NoMSE == true
        noMSE = load([pathname, 'noMSE' fileNum '.txt']);
    end
    s2l = load([pathname, 's2l' fileNum '.txt']);
    RealChain = load([pathname, 'RealChain' fileNum '.txt']);
else
    spots = load(dotsFile);
end

if doAutoRainbow
    numColors = max(spots(:, 7))+1;
    mapSize = size(colormap, 1);
    cMap = hsv;
    colors = cMap(floor((1:numColors)*mapSize/numColors), :);
else
    colors = ['g'; 'r'; 'b'; 'c'; 'm'; 'y'; 'k' ];
end

if plot_straight
    loci = load([pathname, 'loci' fileNum '.txt']);
    
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
    if plot_contour == true
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
    if plotRealAnswer && plot_NoMSE && plot_contour && ~plot_straight
        plot3([spots(c1,xax) noMSE(c1,xax)], [spots(c1,yax) noMSE(c1,yax)],...
                    [spots(c1,zax) noMSE(c1,zax)], 'b');
    end
end
if plotRealAnswer && ~plot_straight
    for c1 = [1:size(FNSpots, 1)]
        plot3(FNSpots(c1,xax), FNSpots(c1,yax), FNSpots(c1,zax), 'o',...
                    'MarkerEdgeColor', colors(FNSpots(c1,7)+colorOffset, :),...
                    'MarkerFaceColor', 'w',...
                    'MarkerSize', 7)
    end
end

if plot_straight == false
    xlabel([axs(xax) ' (', units, ')'])
    ylabel([axs(yax) ' (', units, ')'])
    zlabel([axs(zax) ' (', units, ')'])
else
    xlabel('Contour position (Mbp)')
end

if plot_straight == false
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