whichSet = 1;

if whichSet == 1
    pathname = '/Users/brianross/Desktop/align3d/experimentalDemo/results/Oligopaints/';
    folders = { '3color/chr21/', '5color/chr21/', '10color/chr21/', '3color/chr22/', '5color/chr22/' };
    header = 'chr';
    groups = { '21x1', '21x1', '21x1', '22x1', '22x1' };
    numFiles = [ 120, 120, 120, 151, 151 ];
    numPanels = [ 2 3 ];
else
    pathname = '/Users/brianross/Documents/GitHub/BXC-ORCA-data/';
    folders = { '3color/WT2kb/', '5color/WT2kb/', '10color/WT2kb/', '20color/WT2kb/', ...
                '3color/WT3kb/', '5color/WT3kb/', '10color/WT3kb/', '20color/WT3kb/', ...
                '3color/WT10kb/', '5color/WT10kb/', '10color/WT10kb/', '20color/WT10kb/' };
    header = 'WT';
    groups = { '2kbx1', '2kbx1', '2kbx1', '2kbx1', '3kbx1', '3kbx1', '3kbx1', '3kbx1', '10kbx1', '10kbx1', '10kbx1', '10kbx1' };
    numFiles = [ 259*ones(1, 4) 1744*ones(1, 4) 445*ones(1, 4) ];
    numPanels = [ 4 3 ];
end

lineColor = [1 1 1]*0.3;
ebColor = [1 1 1]*0.8;
nSigma = 3;
numBins = 100;

numFolders = length(folders);

thisfig = figure(1); clf;

for loopFolder = 1:numFolders
    loci = load([ pathname folders{loopFolder} 'loci' groups{loopFolder} '.txt' ]);
    numLoci = size(loci, 1);
    numColors = max(loci(:, 2))+1;
    lastSpot = zeros(1, numColors+1);
    
    binnedProbs = zeros(numBins, 2);
    for loopFile = 1:numFiles(loopFolder)
        l2s = load([ pathname folders{loopFolder} 'l2s' groups{loopFolder} 'ch' num2str(loopFile) '.txt' ]);
        spots = load([ pathname folders{loopFolder} 'spots' groups{loopFolder} 'ch' num2str(loopFile) '.txt' ]);
        lastSpot(1) = 0;
        for loopColor = 1:numColors
            oneLastSpot = find(spots(:, 7) == loopColor-1, 1, 'last');
            if isempty(oneLastSpot)
                oneLastSpot = lastSpot(loopColor);
            end
            lastSpot(loopColor+1) = oneLastSpot;
        end
        numSpots = diff(lastSpot);
        
        probs = load([ pathname folders{loopFolder} 'analysis/' header groups{loopFolder} 'cell' num2str(loopFile) 'Probs.txt' ]);
        whichBin = max(1, min(numBins, ceil(numBins*probs(:, 5))));
        loopP = 0;
        for loopLocus = 1:numLoci
            locusColor = loci(loopLocus, 2)+1;
            if l2s(loopLocus) ~= 0
                correctP = loopP + (l2s(loopLocus) - lastSpot(locusColor));
                binnedProbs(whichBin(correctP), 2) = binnedProbs(whichBin(correctP), 2) + 1;
            end
            for loopSpot = 1:numSpots(locusColor)
                loopP = loopP+1;
                binnedProbs(whichBin(loopP), 1) = binnedProbs(whichBin(loopP), 1) + 1;
            end
        end
    end
    
    binsize = 1/size(binnedProbs, 1);
    x = binsize/2:binsize:1-binsize/2;
    
    subplot(numPanels(1), numPanels(2), loopFolder), hold on
    
    for i = 1:size(x, 2)
        if binnedProbs(i, 1) > 0
            p = min(max(binnedProbs(i, 2), 0.5), binnedProbs(i, 1)-0.5) / binnedProbs(i, 1);
            dy = nSigma * (p*(1-p)/binnedProbs(i, 1))^.5;
            rectangle('Position', [ x(i)-binsize/2, ...
                binnedProbs(i, 2)/binnedProbs(i, 1)-dy, binsize, 2*dy ], ...
                'EdgeColor', 'none', 'FaceColor', ebColor)
        end
    end
    
    plot(x, binnedProbs(:, 2)./binnedProbs(:, 1), '.', 'MarkerSize', 10, 'Color', lineColor)
    plot(0:1, 0:1, ':')
    hold off
    set(thisfig, 'Color', 'w')
    axis([0 1 0 1])
    xlabel('assigned probability'), ylabel('measured probability')
    title(folders{loopFolder}(1:end-1))
end

set(gcf, 'Position', [0, 3000, 200*numPanels(2)+50, 200*numPanels(1)])