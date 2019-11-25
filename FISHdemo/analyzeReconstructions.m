exptNo = 1;
filterExpt1 = true;
filterCriteria = [ 200 8 12 ];
numColors = [ 3 5 10 20 ];
numColoringsPerExpt = [ 1 3 2 4 4 4 ];
numLoci = [ nan, 28, 9, 52, 54, 70 ];
plotLimits = [ nan 4 3 6 6 6 ];
controlPcutoff = 0.05;

exptNames = { 'chr4', 'chr21', 'chr22', 'WT2kb', 'WT3kb', 'WT10kb' };
if exptNo == 1
    exptFolder = '~/Desktop/align3d/experimentalDemo/results/FISH10spot/';
    numColors = 3;
    widths = load([ exptFolder, num2str(numColors) 'color/', exptNames{exptNo}, '/widths.txt' ]);
    numsOfSpots = load([ exptFolder, num2str(numColors) 'color/', exptNames{exptNo}, '/numsOfSpots.txt' ]);
elseif exptNo <= 3
    exptFolder = '~/Desktop/align3d/experimentalDemo/results/Oligopaints/';
else
    exptFolder = '~/Documents/GitHub/BXC-ORCA-data/';
end

S0idx = 5;
I0idx = 6;
itsIdx = 9;
logZidx = 11;
logZAdjustedIdx = 12;
Sidx = 13;
Iidx = 14;

figure(1), clf
figure(2), clf
numColorings = numColoringsPerExpt(exptNo);
maxColorings = max(numColoringsPerExpt);
for loopColoring = 1:numColorings
    
    summaryFile = load([ exptFolder, num2str(numColors(loopColoring)) 'color/', exptNames{exptNo}, '/analysis/' exptNames{exptNo} 'summary.txt' ]);
    if exptNo ~= 1
        mismatchesFile = load([ exptFolder, num2str(numColors(loopColoring)) 'color/', exptNames{exptNo}, ...
            '/analysis/' exptNames{exptNo} 'mismatches.txt' ]);
    end
    numLabelings = max(summaryFile(:, 1));
    numChromosomes = max(summaryFile(:, 2));
    numControlsPerExperiment = max(summaryFile(:, 3));
    
    if exptNo == 1 && filterExpt1
        criterion1 = (widths > filterCriteria(1));
        criterion2 = (numsOfSpots >= filterCriteria(2)) & (numsOfSpots <= filterCriteria(3));
        usedChromosomes = (criterion1 & criterion2);
    else
        usedChromosomes = true(numChromosomes, 1);
    end
    numUsedChromosomes = sum(usedChromosomes);
    
    
        % figure 1 -- null hypothesis test:  control ranks
    
    exptIdx = find(summaryFile(:, 3) == 0);
    usedExptIdx = exptIdx(usedChromosomes);
    statIdx = [ itsIdx, Sidx, logZidx, logZAdjustedIdx ];
    statNames = { '# iterations', 'S', 'log Z', 'adjusted log Z' };
    targetLowRank = [ true, true, false, false ];
    
    allRanks = (1:(numControlsPerExperiment+1))/(numControlsPerExperiment+1);
    logP = log(allRanks);
    meanLogP = mean(logP);
    stdevLogP = sqrt(mean(log(allRanks).^2) - meanLogP^2);
    
    figure(1)
    rankIdx = zeros(numChromosomes, 4);
    for c4 = 1:4
        ranks = zeros(1, numControlsPerExperiment+1);
        for loopChromosome = 1:numChromosomes
            rankIdx(loopChromosome, c4) = sum(summaryFile(exptIdx(loopChromosome)+(1:numControlsPerExperiment), statIdx(c4)) ...
                    < summaryFile(exptIdx(loopChromosome), statIdx(c4))) + 1;
            if usedChromosomes(loopChromosome)
                ranks(rankIdx(loopChromosome, c4)) = ranks(rankIdx(loopChromosome, c4)) + 1;
            end
        end
        
        if targetLowRank(c4)
            ourRank = ranks.*logP;
        else
            ourRank = fliplr(ranks).*logP;
        end
        stdevsBelowMean = (numUsedChromosomes*meanLogP - sum(ourRank)) / ((2*numUsedChromosomes)^.5 * stdevLogP);
        pValue = 0.5 * erfc(stdevsBelowMean);
        
        subplot(maxColorings, 4, 4*(loopColoring-1) + c4), hold on
        bar(1:(numControlsPerExperiment+1), ranks)
        xlim([ 0.5 numControlsPerExperiment + 1.5 ])
        title([ statNames{c4} ' rank; p = ' num2str(pValue) ])
        xlabel([ 'rank (1 + ' num2str(numControlsPerExperiment) ' controls)' ])
        ylabel('counts')
    end
    set(gcf, 'Position', [0, 2000, 1200, 700])
    
    
    figure(2)
    if exptNo == 1
        hold on
        plot(rankIdx(~criterion2, 3), widths(~criterion2), 'o', 'MarkerSize', 3, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r')
        plot(rankIdx(criterion2, 3), widths(criterion2), 'o', 'MarkerEdgeColor', 'b')
        plot([ 1 numControlsPerExperiment+1 ], filterCriteria(1) + [0 0], 'k-.')
        title('Z rank vs. chromosome thickness'), xlabel('Z rank'), ylabel('thickness (nm)')
        legend('wrong # spots', 'correct # spots', 'thickness cutoff')
    else
        controlS = sort(summaryFile(summaryFile(:, 3) ~= 0, Sidx));
        controlI = sort(summaryFile(summaryFile(:, 3) ~= 0, Iidx));
        numControls = length(controlS);
        Sbounds = controlS(round([ controlPcutoff 1-controlPcutoff ]*numControls));
        Ibounds = controlI(round([ controlPcutoff 1-controlPcutoff ]*numControls));
        
        colorFrac = (loopColoring-1) / (maxColorings-1);
        colorSub = [ 1-colorFrac, 1, colorFrac ];
        
        subplot(1, 3, 1), hold on
        plot(summaryFile(usedExptIdx, Sidx), summaryFile(usedExptIdx, Iidx), 'o', ...
            'MarkerEdgeColor', 1-colorSub/3, 'MarkerSize', (3e3/numUsedChromosomes)^.5)
        plot(summaryFile(usedExptIdx, S0idx), summaryFile(usedExptIdx, I0idx), 'd', 'MarkerFaceColor', 1-colorSub, 'MarkerEdgeColor', 'none')
        if Ibounds(2) > plotLimits(exptNo)
            plot(Sbounds([1 1 2 2]), min(plotLimits(exptNo), Ibounds([2 1 1 2])), '-.', 'Color', 1-colorSub, 'LineWidth', 1)
        else
            plot(Sbounds([1 2 2 1 1]), Ibounds([1 1 2 2 1]), '-.', 'Color', 1-colorSub, 'LineWidth', 1)
        end
        
        subplot(1, 3, 2), hold on
        plot(summaryFile(usedExptIdx, Sidx), mismatchesFile(1:numChromosomes) / numLoci(exptNo), 'o', ...
            'MarkerEdgeColor', 1-colorSub/3, 'MarkerSize', (3e3/numUsedChromosomes)^.5)
        plot(Sbounds(1)+[0 0], [0 1], '-.', 'Color', 1-colorSub, 'LineWidth', 1)
        
        subplot(1, 3, 3), hold on
        plot(summaryFile(usedExptIdx, Iidx), mismatchesFile(1:numChromosomes) / numLoci(exptNo), 'o', ...
            'MarkerEdgeColor', 1-colorSub/3, 'MarkerSize', (3e3/numUsedChromosomes)^.5)
    end
    set(gcf, 'Position', [0, 0, 600, 150])
end

if exptNo == 1
    set(gcf, 'Position', [0, 0, 300, 200])
else
    set(gcf, 'Position', [0, 0, 600, 150])
    
    subplot(1, 3, 1), title('S vs. I')
    plot([ 0 plotLimits(exptNo) ], [ 0 plotLimits(exptNo) ], ':')
    xlim([0 plotLimits(exptNo)]), ylim([0 plotLimits(exptNo)])
    xlabel('S (bits)'), ylabel('I (bits)')
    
    subplot(1, 3, 2)
    title('S vs. mismatch rate'), xlabel('S (bits)'), ylabel('mismatch rate')
    
    subplot(1, 3, 3)
    xlim([0 plotLimits(exptNo)])
    title('Mismatch rate vs. I'), xlabel('I (bits)'), ylabel('mismatch rate')
end
