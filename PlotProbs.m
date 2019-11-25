%pathname = '/Users/brianross/Desktop/align3d/experimentalDemo/results/oligopaintsRecolored/chr21/';
pathname = '/Users/brianross/Desktop/align3d/experimentalDemo/results/FISH10spot/3color/chr4/analysis/';
pathname = '/Users/brianross/Desktop/probsControl0.txt';

plotRealAnswer = false;
plotSecondRealAnswer = false;
plotTestAnswer = false;
useS2L = true;
chromosomeName = '12';
plotindex = 3;       % 2-4 = x-z
colorful = false;

%AllProbs = load([pathname, 'chr4x1cell' chromosomeName 'Probs.txt']);
AllProbs = load('/Users/brianross/Desktop/probsControl1.txt');
if plotRealAnswer
    if useS2L
        loci = load([pathname, 'loci2kbx1.txt']);
        spots = load([pathname, 'spots2kbx1ch2.txt']);
        l2s = load([pathname, 'l2s2kbx1ch2.txt']);
        TrueProbs = [ loci(l2s~=0, 1), spots(l2s(l2s~=0), 1:3) ];
    else
        TrueProbs = load([pathname, 'analysis/series1chain1TrueProbs.txt']);
    end
    if plotSecondRealAnswer
        secondTrueProbs = load([pathname, '21TrueProbs.txt']);
    end
elseif plotTestAnswer
    testAnswer = [ 3 5 4 1 1 2 2 3 2 1 ];   % cell 2 ch 2, alt 1
    endLength = 300;
    circleRadius = 150;
    circleNumPts = 30;
    contourWidth = 30;
    pinOffset = 300;
    pinSize = 13;
    
    startIdx = [ 1 find(diff(AllProbs(:, 1)) ~= 0)'+1 ];
    testIdx = startIdx+testAnswer-1;
    testIdx = testIdx(~isnan(testIdx));
end

thisfig = figure; clf, hold on

c1_min = min(AllProbs(:, plotindex));
c1_max = max(AllProbs(:, plotindex));
fn_level = c1_max + 0.2*(c1_max-c1_min);

if colorful == true
    colors = ['r'; 'g'; 'b'; 'c'; 'm'; 'y'; 'k'; 'r'; 'g' ];
else
    colors = repmat('b', 1, max(AllProbs(:, 6)));
end
yax_labels = [ 'q', 'x', 'y', 'z' ];
p_sum = 0;

for i = 1:size(AllProbs, 1)
    dsize = AllProbs(i, 5).^.5*20;
    if dsize > .49
        plot(AllProbs(i, 1), AllProbs(i, plotindex), ...
            ['o' colors(AllProbs(i, 6))], 'MarkerSize', dsize)
    end
    p_sum = p_sum + AllProbs(i, 5);
    
    fn_dot = 0;
    if i == size(AllProbs, 1), fn_dot = 1;
    else if AllProbs(i+1, 1) ~= AllProbs(i, 1), fn_dot = 1; end
    end
    
    if fn_dot == 1
        dsize = (1-p_sum).^.5*20;
        if dsize > .49
            plot(AllProbs(i, 1), fn_level, 'om', 'MarkerSize', dsize)
            if plotRealAnswer
            if sum(TrueProbs(:, 1) == AllProbs(i, 1)) == 0
                plot(AllProbs(i, 1), fn_level, 'xm', 'MarkerSize', dsize)
            end, end
        end
        p_sum = 0;
    end
end
if plotRealAnswer
    plot(TrueProbs(:, 1), TrueProbs(:, plotindex), '-')
    if plotSecondRealAnswer
        plot(secondTrueProbs(:, 1), secondTrueProbs(:, plotindex), 'r-')
    end
    plot(TrueProbs(:, 1), 0*TrueProbs(:, 1)+fn_level, ':')
elseif plotTestAnswer
    plot(AllProbs(testIdx, 1), AllProbs(testIdx, plotindex), '-.')
end

xlabel('Contour position (Mbp)')
ylabel(['R_' yax_labels(plotindex) '(nm)'])
set(thisfig, 'Color', 'w')

if plotTestAnswer
    zplane = 1.3*min(AllProbs(:, 4)) - 0.3*max(AllProbs(:, 4));
    firstL = AllProbs(testIdx(1), 1);
    lastL = AllProbs(testIdx(end), 1);
    
    lstep = (lastL-firstL)/400;
    ls = firstL:lstep:lastL;
    
    xs = [1;1;1]*ls;
    for c3 = 1:3
        xs(c3, :) = spline(AllProbs(testIdx, 1), ...
            AllProbs(testIdx, 1+c3), ls);
    end
    
    clrList = 'grb';
    figure, hold on
    plot3(xs(1, :), xs(2, :), xs(3, :), ...
        'Color', .3*[1 1 1], 'LineWidth', contourWidth)
    
    n1 = xs(:, 1)-xs(:, 2);
    n2 = xs(:, end)-xs(:, end-1);
    n1 = n1/sqrt(n1'*n1);
    n2 = n2/sqrt(n2'*n2);
    xend1 = [ xs(:, 1)+endLength*n1, xs(:, 1) ];
    xend2 = [ xs(:, end), xs(:, end)+endLength*n2 ];
    plot3(xend1(1, :), xend1(2, :), xend1(3, :), ...
        'Color', .5*[1 1 1], 'LineWidth', contourWidth/2)
    plot3(xend2(1, :), xend2(2, :), xend2(3, :), ...
        'Color', .5*[1 1 1], 'LineWidth', contourWidth/2)
    plot3(xs(1, :), xs(2, :), 0*xs(3, :)+zplane, ...
        'Color', .8*[1 1 1], 'LineWidth', 1.6*contourWidth)
    for cd = 1:length(testIdx)
        plot3(AllProbs(testIdx(cd), 2), AllProbs(testIdx(cd), 3), ...
             AllProbs(testIdx(cd), 4)+pinOffset, 'o', ...
             'MarkerSize', pinSize, 'MarkerEdgeColor', 'k', ...
             'MarkerFaceColor', clrList(AllProbs(testIdx(cd), 6)-1))
        plot3(AllProbs(testIdx(cd), 2)+[0 0], ...
              AllProbs(testIdx(cd), 3)+[0 0], ...
             [AllProbs(testIdx(cd), 4)+pinOffset zplane], ...
             'k', 'LineWidth', contourWidth/10)
    end
    axis equal
    axis off
    view(0, 50)
end

try
    controls = load([pathname, 'controls', chromosomeName, '.txt']);
    figure, hold on
    plot([1 size(controls, 1)-1], controls(1, 1)*[1 1], 'r')
    plot(sort(controls(2:end, 1)), 'b')
    xlabel('sorted control #')
    ylabel('entropy (bits)')
    legend('true mapping', 'control mappings')
end