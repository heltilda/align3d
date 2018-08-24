pathname = '/Users/brianross/Desktop/align3d/Results/';
folders = { 'small100/analysis/', 'med100/analysis/', 'large100/analysis/' };
fileNames = { 'series1BinnedProbs.txt', 'series2BinnedProbs.txt', ...
              'series1BinnedIteratedProbs.txt', 'series2BinnedIteratedProbs.txt' };
% folders = { 'small100/analysis/' };
% fileNames = { 'exactBinnedProbs.txt' };

lineColor = [1 1 1]*0.3;
ebColor = [1 1 1]*0.8;
nSigma = 3;

numFolders = length(folders);
numFiles = length(fileNames);

thisfig = figure(2);

for loopFolder = 1:numFolders
    for loopFile = 1:numFiles
        
        binnedProbs = load([ pathname folders{loopFolder} fileNames{loopFile} ]);
        
        binsize = 1/size(binnedProbs, 1);
        x = binsize/2:binsize:1-binsize/2;
        
        subplot(numFolders, numFiles, (loopFolder-1)*numFiles+loopFile), hold on
        
        for i = 1:size(x, 2)
            p = min(max(binnedProbs(i, 2), 0.5), binnedProbs(i, 1)-0.5) / binnedProbs(i, 1);
            dy = nSigma * (p*(1-p)/binnedProbs(i, 1))^.5;
            rectangle('Position', [x(i)-binsize/2, ...
                binnedProbs(i, 2)/binnedProbs(i, 1)-dy, binsize, 2*dy], ...
                'EdgeColor', 'none', 'FaceColor', ebColor)
        end
        
        plot(x, binnedProbs(:, 2)./binnedProbs(:, 1), '.', 'MarkerSize', 10, 'Color', lineColor)
        plot(0:1, 0:1, ':')
        hold off
        axis([0 1 0 1])
        
        if loopFolder == numFolders
            xlabel('assigned probability')
        end
        if loopFile == 1
            ylabel('measured probability')
            title([ folders{loopFolder} fileNames{loopFile} ])
        else
            title(fileNames{loopFile})
        end
        set(thisfig, 'Color', 'w')
    end
end

set(gcf, 'Position', [0, 3000, 250*numFiles, 200*numFolders])