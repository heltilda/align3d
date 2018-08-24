function legendString = makePercentilePlot(whichSeries, whichIdx, numContours, plotDotted)

global seriesSamples logNMax clrs

percentiles = [ 10 25 50 75 90 ];
fillSaturation = [ .25 .5 1 .5 .25 ];
binSpacing = .25;
seriesTitles = { '1', '2', 'diff' };

numPercentiles = length(percentiles);
midPercentile = (numPercentiles+1)/2;
firstPercentile = midPercentile+1-numContours;
if whichSeries == 1
    lastPercentile = midPercentile+numContours-1;
    lw = 3*ones(1, numPercentiles);
    lw(midPercentile) = 2;
else
    lastPercentile = midPercentile;
    lw = ones(1, numPercentiles);
end

binXs = 0:binSpacing:logNMax;
numBins = length(binXs);

lineType = '-';
if plotDotted
    lineType = '-.';
    lw(:) = 1;
end

seriesPercentiles = nan(numPercentiles, numBins);


    % find the percentile range for each bin

for loopPercentile = 1:numPercentiles
    for loopBin = 1:numBins
        nonNans = ~isnan(seriesSamples{whichSeries, loopBin}(:, whichIdx));
        if sum(nonNans) ~= 0
            seriesPercentiles(loopPercentile, loopBin) = ...
                prctile(seriesSamples{whichSeries, loopBin}(nonNans, whichIdx), percentiles(loopPercentile), 1);
        end
    end
end


    % plot the percentile contours

legendString = {};
hold on
for loopPercentile = firstPercentile:lastPercentile
    if whichSeries == 1 || loopPercentile == midPercentile
        plot(binXs, seriesPercentiles(loopPercentile, :), lineType, 'LineWidth', lw(loopPercentile), ...
            'Color', 1 - fillSaturation(loopPercentile)*(1-clrs(whichSeries, :)))
    else
        nonNans = ~isnan(seriesPercentiles(loopPercentile, :)) & ~isnan(seriesPercentiles(end+1-loopPercentile, :));
        fill( [ binXs(nonNans) fliplr(binXs(nonNans)) ], ...  
            [ seriesPercentiles(loopPercentile, nonNans) fliplr(seriesPercentiles(end+1-loopPercentile, nonNans)) ], ...
            clrs(whichSeries, :)*fillSaturation(loopPercentile) + [1 1 1]*(1-fillSaturation(loopPercentile)), ...
            'EdgeColor', 'none' );
    end
    if loopPercentile < midPercentile
        legendString{end+1} = [ num2str(percentiles(loopPercentile)) '% - ' ...
            num2str(percentiles(numPercentiles+1-loopPercentile)) '%' ];
    elseif loopPercentile == midPercentile
        legendString{end+1} = [ '50% series ' seriesTitles{whichSeries} ];
    end
end

xlim([0 logNMax])
xlabel('log_{10} # terms')

hold off

end