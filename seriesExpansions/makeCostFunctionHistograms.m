costDir = '/Users/brianross/Desktop/align3d/Results/large100/analysis/';
sampleNos = 1:100;
costThreshold = 1;
preOptStrings = { '_noPO', '' };
WFstrings = { 'w', 'f' };
numSamples = length(sampleNos);

numIterations = { nan(numSamples, 2), nan(numSamples, 2) };


    % load each cost function curve and extract the number of steps to
    % convergence

for preOpt = 1:2
    for sample = sampleNos
        doSkipNextWF = false;
        for wf = 1:2
            if ~doSkipNextWF
                costTable = dlmread([costDir 'chain' num2str(sample) '_' WFstrings{wf} 'Cs' preOptStrings{preOpt} '.txt']);
                if costTable(end, 1) <= costThreshold
                    numIterations{preOpt}(sample, wf) = size(costTable, 1)-1;
                else
                    doSkipNextWF = true;
                end
            end
        end
    end
end
convergedTotals = [ [ numSamples sum(~isnan(numIterations{1}), 1) ]; [ numSamples sum(~isnan(numIterations{2}), 1) ] ];


    % generate histograms of convergence time

figure
for wf = 1:2
    subplot(1, 2, wf), hold on
    hist([ numIterations{1}(:, wf), numIterations{2}(:, wf) ])
    legend([ 'no preOpt:  ' num2str(convergedTotals(1, wf+1)) '/' num2str(convergedTotals(1, wf)) ], ...
           [ 'preOpt:  ' num2str(convergedTotals(2, wf+1)) '/' num2str(convergedTotals(2, wf)) ])
    xlabel('# iterations'), ylabel('count')
    ylim([0 numSamples])
    title(WFstrings(wf))
end
set(gcf, 'Position', [0, 3000, 500, 200])