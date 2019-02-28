pathname = '/Users/brianross/Desktop/align3d/seriesExpansion/results/large4superres/';
fileNum = '3';
doNaive = false;

AllProbs = load([pathname, 'analysis/series2chain' fileNum 'Probs.txt']);
TrueProbs = load([pathname, 'analysis/series2chain' fileNum 'TrueProbs.txt']);
RealChain = load([pathname, 'RealChain' fileNum '.txt']);
loci = load([pathname 'loci' fileNum '.txt']);
spots = load([pathname 'spots' fileNum '.txt']);
l2s = load([pathname 'l2s' fileNum '.txt'])';

numColors = max(spots(:, 7))+1;
numLoci = size(loci, 1);
numSpots = size(spots, 1);
numSpotsOfColor = zeros(1, numColors);
for loopColor = 1:numColors
    numSpotsOfColor(loopColor) = sum(spots(:, 7) == loopColor-1);
end
colorIdx = [ 0 cumsum(numSpotsOfColor) ];

APidxBase = zeros(1, numLoci);
guessContourIndices = zeros(1, numLoci);
guessSpotMappings = zeros(1, numLoci);
locusMapped = false(1, numLoci);
spotUse = zeros(1, numSpots);
spotInPlay = true(1, numSpots);

APidx = 0;
firstTime = true;
while sum(~locusMapped) > 0
    for locus = 1:numLoci
        if ~locusMapped(locus)
            locusColor = loci(locus, 2) + 1;
            
            bestP = 0;
            p_sum = 0;
            if firstTime
                APidxBase(locus) = APidx;
            else
                APidx = APidxBase(locus);
            end
            for spot = colorIdx(locusColor)+1:colorIdx(locusColor+1)
                
                APidx = APidx + 1;
                if AllProbs(APidx, 5) > bestP && spotInPlay(spot)
                    bestP = AllProbs(APidx, 5);
                    bestAPIdx = APidx;
                    bestSpot = spot;
                end
                p_sum = p_sum + AllProbs(APidx, 5);
                
            end
            
            if bestP > 1-p_sum
                guessContourIndices(locus) = bestAPIdx;
                guessSpotMappings(locus) = bestSpot;
                spotUse(bestSpot) = spotUse(bestSpot)+1;
            end
        end
    end
    
    locusMapped(:) = true;
    for loopSpot = 1:numSpots
        if spotUse(loopSpot) > 1
            
            overlappingLoci = find(guessSpotMappings == loopSpot);
            locusColor = loci(overlappingLoci(1), 2) + 1;
            APindices = 1:colorIdx(locusColor+1)-colorIdx(locusColor);
            APindices = APindices(spotUse(colorIdx(locusColor)+1:colorIdx(locusColor+1)) > 0);
            
            minPslack = 2;
            for locus = overlappingLoci
                pSlack = 1 - sum(AllProbs(APidxBase(locus) + APindices, 5));
                if pSlack < minPslack
                    minPslack = pSlack;
                    locusToKeep = locus;
                end
            end
            
            overlappingLoci = overlappingLoci(overlappingLoci ~= locusToKeep);
            if ~doNaive
                locusMapped(overlappingLoci) = false;
                guessContourIndices(overlappingLoci) = 0;
                guessSpotMappings(overlappingLoci) = 0;
            end
            spotUse(loopSpot) = 1;
        end
    end
    
    spotInPlay = (spotUse == 0);
    firstTime = false;
end

guessContourIndices = guessContourIndices(guessContourIndices > 0);

figure, hold on
plot(RealChain(:, 3), RealChain(:, 2), 'Color', [ .8 .8 .8 ], 'LineWidth', 1)
plot(AllProbs(guessContourIndices, 4), AllProbs(guessContourIndices, 3), 'r')
plot(TrueProbs(:, 4), TrueProbs(:, 3), 'b')

axis equal
title([ 'conformation ', fileNum ])
legend('DNA contour', 'reconstruction errors', 'ideal reconstruction')
set(gcf, 'Position', [0, 3000, 300, 200])

disp(['error:  ' num2str(lev(guessSpotMappings, l2s)) ' / ' num2str(numLoci) ' loci'])