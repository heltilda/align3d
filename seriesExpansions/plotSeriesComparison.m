global seriesSamples logNMax clrs


    % adjust these parameters

seriesDir = '/Users/brianross/Desktop/align3d/seriesExpansion/Results/large100/analysis/';
sampleNos = 1:100;
binSpacing = .25;
logNMax = 3;
infoMax = 0.4;
clrs = [ .2 .2 .5; .5 .2 .2; .5 .2 .2 ];


    % variable definitions

idxN = 1;
idxNterms = 2;
idxZ = 3;
idxS = 4;
idxI = 5;
idxT = 6;
idxSit = 7;
idxIit = 8;
idxTit = 9;
idxSI = 10;
idxSIit = 11;
idxIdiff = 12;
idxIitdiff = 13;
numCharacteristics = 13;

numSamples = length(sampleNos);
numBins = length(0:binSpacing:logNMax);
seriesSamples = cell(3, numBins);
sampleZ0times = nan(2, numSamples);
sampleTimeBins = nan(2, numSamples);
doIndividualPlots = (numSamples <= 10);
numNans = zeros(2, 2);
numDataPoints = [ 0 0 ];


    % make the plots

figure
for sampleCounter = 1:numSamples
    sampleStr = num2str(sampleNos(sampleCounter));
    series1 = dlmread([seriesDir 'series1chain' sampleStr '.txt']);
    series2 = dlmread([seriesDir 'series2chain' sampleStr '.txt']);
    
    bothSeries = { series1(series1(:, idxN) >= 0, :), series2(series2(:, idxN) >= 0, :) };
    for loopSeries = 1:2
        bothSeries{loopSeries}(:, [ idxSI idxSIit ]) = ...
            bothSeries{loopSeries}(:, [ idxS idxSit ]) - bothSeries{loopSeries}(:, [ idxI idxIit ]);
        whichZ0 = find(bothSeries{loopSeries}(:, idxNterms) == 1, 1);
        bothSeries{loopSeries}(:, idxIdiff) = bothSeries{loopSeries}(:, idxI) - bothSeries{loopSeries}(whichZ0, idxIit);
        bothSeries{loopSeries}(:, idxIitdiff) = bothSeries{loopSeries}(:, idxIit) - bothSeries{loopSeries}(whichZ0, idxIit);
        
        upperTimeBin = find(bothSeries{loopSeries}(:, idxT) > bothSeries{loopSeries}(whichZ0, idxTit), 1, 'first');
        if ~isempty(upperTimeBin)
            sampleZ0times(loopSeries, sampleCounter) = bothSeries{loopSeries}(whichZ0, idxTit);
            if upperTimeBin > 1
                binFrac = (log10(bothSeries{loopSeries}(whichZ0, idxTit)) - log10(bothSeries{loopSeries}(upperTimeBin-1, idxT))) ...
                        / diff(log10(bothSeries{loopSeries}(upperTimeBin-1:upperTimeBin, idxT)));
                sampleTimeBins(loopSeries, sampleCounter) = (1-binFrac)*log10(bothSeries{loopSeries}(upperTimeBin-1, idxNterms)) ...
                        + binFrac*log10(bothSeries{loopSeries}(upperTimeBin, idxNterms));
            end
        end
        
        numNans(:, loopSeries) = numNans(:, loopSeries) + sum(isnan(bothSeries{loopSeries}(:, [ idxZ idxTit ])), 1)';
        numDataPoints(loopSeries) = numDataPoints(loopSeries) + size(bothSeries{loopSeries}, 1);
    end
    
    
        % first tally up (S, I, log Z) for each sample
    
    if ~doIndividualPlots
        lastSamples = zeros(numBins, numCharacteristics);
        for loopBin = 1:numBins
            for loopSeries = 1:2
                binVal = (loopBin-1)*binSpacing;
                sample2Idx = find(log10(bothSeries{loopSeries}(:, idxNterms)) > binVal, 1, 'first');
                if ~isempty(sample2Idx)
                    sample2frac = (binVal - log10(bothSeries{loopSeries}(sample2Idx-1, idxNterms))) ...
                        / diff(log10(bothSeries{loopSeries}(sample2Idx-1:sample2Idx, idxNterms)));
                else
                    sample2Idx = size(bothSeries{loopSeries}, idxN);
                    sample2frac = 1.;
                end
                
                if sample2Idx > 1
                    newSample = (1-sample2frac)*bothSeries{loopSeries}(sample2Idx-1, :) + sample2frac*bothSeries{loopSeries}(sample2Idx, :);
                    seriesSamples{loopSeries, loopBin} = [ seriesSamples{loopSeries, loopBin}; newSample ];
                    lastSamples(loopBin, :) = newSample;
                end
            end
            seriesSamples{3, loopBin} = seriesSamples{2, loopBin}-seriesSamples{1, loopBin};
        end
    end
    
    
        % if it's time to draw a plot, do that
    
    if doIndividualPlots
        
        
            % find median/Nth percentile of (S, I, log Z) over all samples
        
        seriesXs = { log10(bothSeries{1}(:, idxNterms))', log10(bothSeries{2}(:, idxNterms))' };
        
        
            % draw the I/S plot
        
        subplot(2, numSamples, sampleCounter), hold on
        %colormap(cMap)
        
        plot(seriesXs{1}, bothSeries{1}(:, idxI), '-', 'LineWidth', .5, 'Color', clrs(1, :))
        plot(seriesXs{2}, bothSeries{2}(:, idxI), '-', 'LineWidth', .5, 'Color', clrs(2, :))
        plot(seriesXs{1}, bothSeries{1}(:, idxIit), '-.', 'LineWidth', .5, 'Color', clrs(1, :))
        plot(seriesXs{2}, bothSeries{2}(:, idxIit), '-.', 'LineWidth', .5, 'Color', clrs(2, :))
        plot(seriesXs{1}, bothSeries{1}(:, idxS), ':', 'Color', clrs(1, :))
        plot(seriesXs{2}, bothSeries{2}(:, idxS), ':', 'Color', clrs(2, :))
        
        plot(seriesXs{1}, bothSeries{1}(:, idxI), 'o', 'MarkerEdgeColor', clrs(1, :), 'MarkerFaceColor', clrs(1, :))
        plot(seriesXs{2}, bothSeries{2}(:, idxI), 'o', 'MarkerEdgeColor', clrs(2, :), 'MarkerFaceColor', clrs(2, :))
        
        axis([0 logNMax 0 infoMax])
        title(['Conformation #' sampleStr])
        xlabel('log_{10} # terms')
        if sampleCounter == 1
            ylabel('missing info')
            plotLabels = { 'series 1', 'series 2', 'S1 iterated', 'S2 iterated' };
            legend(plotLabels)
        end
        %caxis([1 numSamples+1])
        
        
            % draw the log Z plot
        
        subplot(2, numSamples, sampleCounter+numSamples), hold on
%        colormap(cMap)
        
        plot(seriesXs{1}, bothSeries{1}(:, idxZ), '-', 'LineWidth', .5, 'Color', clrs(1, :))
        plot(seriesXs{2}, bothSeries{2}(:, idxZ), '-', 'LineWidth', .5, 'Color', clrs(2, :))
        
        plot(seriesXs{1}, bothSeries{1}(:, idxZ), 'o', 'MarkerEdgeColor', clrs(1, :), 'MarkerFaceColor', clrs(1, :))
        plot(seriesXs{2}, bothSeries{2}(:, idxZ), 'o', 'MarkerEdgeColor', clrs(2, :), 'MarkerFaceColor', clrs(2, :))
        
        xlim([0 logNMax])
        xlabel('log_{10} # terms')
        if sampleCounter == 1
            ylabel('log Z')
        end
%        caxis([1 numSamples+1])
        set(gcf, 'Position', [0, 3000, 250*numSamples, 400])
    end
end


    % draw the entropy error plot

if ~doIndividualPlots
    
    
        % make the 4 main-text plots
    
    subplot(4, 1, 1)
    l1 = makePercentilePlot(1, idxI, 1, false);
    l2 = makePercentilePlot(2, idxI, 1, false);
    l3 = makePercentilePlot(1, idxIit, 1, true); l3{1} = [ l3{1} ' iterated' ];
    l4 = makePercentilePlot(2, idxIit, 1, true); l4{1} = [ l4{1} ' iterated' ];
    ylabel('missing info (I)'), ylim([0 infoMax]), legend([l1 l2 l3 l4])
    
    subplot(4, 1, 2)
    l1 = makePercentilePlot(3, idxI, 3, false);
    hold on, plot([0 logNMax], [0 0], 'Color', [ .2 .2 .5 ], 'LineWidth', 2);
    ylabel('I_2 - I_1'), legend(l1)
    
    subplot(4, 1, 3), thisPlot = gca;
    l1 = makePercentilePlot(2, idxIdiff, 3, false);
    tBinsAv = median(sampleTimeBins(2, ~isnan(sampleTimeBins(2, :))));
    lineBounds = thisPlot.YLim+[.05 -.05]*diff(thisPlot.YLim);
    hold on, plot(tBinsAv*[1 1], lineBounds, ':');
    text(tBinsAv, lineBounds(2), [ 'Z_{0it} comp time = ' num2str(median(sampleZ0times(2, ~isnan(sampleZ0times(2, :))))) ' s' ]);
    plot([0 logNMax], [0 0], 'k-.');
    ylabel('I - I_{it0}'), legend(l1)
    
    subplot(4, 1, 4)
    l1 = makePercentilePlot(2, idxIitdiff, 3, false);
    hold on, plot([0 logNMax], [0 0], 'k:');
    ylabel('I_{it} - I_{it0}'), legend(l1), set(gcf, 'Position', [0, 3000, 250, 860])
    
    
        % make the supplemental plots for the appendix
    
    figure
    subplot(3, 1, 1)
    l1 = makePercentilePlot(1, idxZ, 1, false);
    l2 = makePercentilePlot(2, idxZ, 1, false);
    ylabel('log Z'), legend([l1 l2])
    
    subplot(3, 1, 2)
    l1 = makePercentilePlot(3, idxZ, 3, false);
    hold on, plot([0 logNMax], [0 0], 'Color', [ .2 .2 .5 ], 'LineWidth', 2);
    ylabel('\Delta log Z'), legend(l1)
    
    subplot(3, 1, 3)
    l1 = makePercentilePlot(2, idxSI, 3, false);
    l2 = makePercentilePlot(1, idxSI, 3, false);
    hold on, plot([0 logNMax], [0 0], 'k-.');
    ylabel('S - I'), legend([l1 l2]), set(gcf, 'Position', [0, 3000, 250, 600])
end

disp(['unconverged w fractions (series 1/2):    ', num2str(numNans(1, :)./numDataPoints)])
disp(['unconverged f fractions (series 1/2):    ', num2str(numNans(2, :)./numDataPoints)])
