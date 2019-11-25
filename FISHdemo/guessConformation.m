function [ guess_l2s, guessAllProbsIndices, numMismatches ] = guessConformation(allProbs, locusColors, spotColors, true_l2s)

numColors = max(locusColors)+1;
numLoci = length(locusColors);
numSpotsOfColor = zeros(1, numColors);
for loopColor = 1:numColors
    numSpotsOfColor(loopColor) = sum(spotColors == loopColor-1);
end
colorIdx = [ 0 cumsum(numSpotsOfColor) ];
numSpots = colorIdx(end);

APidxBase = zeros(1, numLoci);
guessAllProbsIndices = zeros(1, numLoci);
guess_l2s = zeros(1, numLoci);
locusMapped = false(1, numLoci);
spotUse = zeros(1, numSpots);
spotInPlay = true(1, numSpots);

APidx = 0;
firstTime = true;
while sum(~locusMapped) > 0
    for locus = 1:numLoci
        if ~locusMapped(locus)
            locusColor = locusColors(locus) + 1;
            
            bestP = 0;
            p_sum = 0;
            if firstTime
                APidxBase(locus) = APidx;
            else
                APidx = APidxBase(locus);
            end
            for spot = colorIdx(locusColor)+1:colorIdx(locusColor+1)
                
                APidx = APidx + 1;
                if allProbs(APidx, 5) > bestP && spotInPlay(spot)
                    bestP = allProbs(APidx, 5);
                    bestAPIdx = APidx;
                    bestSpot = spot;
                end
                p_sum = p_sum + allProbs(APidx, 5);
                
            end
            
            if bestP > max(0, 1-p_sum)
                guessAllProbsIndices(locus) = bestAPIdx;
                guess_l2s(locus) = bestSpot;
                spotUse(bestSpot) = spotUse(bestSpot)+1;
            end
        end
    end
    
    locusMapped(:) = true;
    for loopSpot = 1:numSpots
        if spotUse(loopSpot) > 1
            
            overlappingLoci = find(guess_l2s == loopSpot);
            locusColor = locusColors(overlappingLoci(1)) + 1;
            APindices = 1:colorIdx(locusColor+1)-colorIdx(locusColor);
            APindices = APindices(spotUse(colorIdx(locusColor)+1:colorIdx(locusColor+1)) > 0);
            
            minPslack = 2;
            for locus = overlappingLoci
                pSlack = 1 - sum(allProbs(APidxBase(locus) + APindices, 5));
                if pSlack < minPslack
                    minPslack = pSlack;
                    locusToKeep = locus;
                end
            end
            
            overlappingLoci = overlappingLoci(overlappingLoci ~= locusToKeep);
            locusMapped(overlappingLoci) = false;
            guessAllProbsIndices(overlappingLoci) = 0;
            guess_l2s(overlappingLoci) = 0;
            spotUse(loopSpot) = 1;
        end
    end
    
    spotInPlay = (spotUse == 0);
    firstTime = false;
end

if exist('true_l2s', 'var')
    numMismatches = sum(guess_l2s ~= true_l2s);
end

end
