pathname = '/Users/brianross/Desktop/align3d/Results/large4/';
fileNum = '4';

AllProbs = load([pathname, 'analysis/series2chain' fileNum 'Probs.txt']);
TrueProbs = load([pathname, 'analysis/series2chain' fileNum 'TrueProbs.txt']);
RealChain = load([pathname, 'RealChain' fileNum '.txt']);

figure, hold on
plotindex = 4;       % 2-4 = x-z

p_sum = 0;
bestP = 0;
GuessContour = zeros(0, 3);
for i = 1:size(AllProbs, 1)
    
    if AllProbs(i, 5) > bestP
        bestP = AllProbs(i, 5);
        bestPIdx = i;
    end
    p_sum = p_sum + AllProbs(i, 5);
    
    fn_dot = false;
    if i == size(AllProbs, 1)
        fn_dot = true;
    else
        if AllProbs(i+1, 1) ~= AllProbs(i, 1)
            fn_dot = true;
        end
    end
    
    if fn_dot
        if bestP > 1-p_sum
            GuessContour = [ GuessContour; AllProbs(bestPIdx, 2:4) ];
        end
        bestP = 0;
        p_sum = 0;
    end
end

plot(RealChain(:, 3), RealChain(:, 2), 'Color', [ .8 .8 .8 ], 'LineWidth', 1)

plot(GuessContour(:, 3), GuessContour(:, 2), 'r')
plot(TrueProbs(:, 4), TrueProbs(:, 3), 'b')

axis equal
title([ 'conformation ', fileNum ])
legend('DNA contour', 'reconstruction errors', 'ideal reconstruction')
% plot3(GuessContour(:, 1), GuessContour(:, 2), GuessContour(:, 3))
% plot3(TrueProbs(:, 2), TrueProbs(:, 3), TrueProbs(:, 4))
set(gcf, 'Position', [0, 3000, 300, 200])