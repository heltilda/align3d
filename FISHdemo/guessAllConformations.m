numColorings = [ 3, 2, 4, 4, 4 ];
numLabelings = [ 2, 1, 2, 1, 1 ];
numCells = [ 120, 151, 259, 1744, 445 ];

numColors = { '3', '5', '10', '20' };
exptDirectories = [ 1, 1, 2, 2, 2 ];
exptPaths = { '/Users/brianross/Desktop/align3d/experimentalDemo/results/Oligopaints/', '/Users/brianross/Documents/GitHub/BXC-ORCA-data/' };
exptPrefixes = { 'chr', 'chr', 'WT', 'WT', 'WT' };
exptLabels = { '21', '22', '2kb', '3kb', '10kb' };

for loopExpt = 1:5
    for loopColoring = 1:numColorings(loopExpt)
        pathname = [ exptPaths{exptDirectories(loopExpt)} numColors{loopColoring} 'color/' exptPrefixes{loopExpt} exptLabels{loopExpt} '/' ];
        
        mismatches = zeros(1, numLabelings(loopExpt)*numCells(loopExpt));
        for loopLabeling = 1:numLabelings(loopExpt)
            groupNum = [ exptLabels{loopExpt} 'x' num2str(loopLabeling) ];
            
            loci = load([pathname 'loci' groupNum '.txt']);
            for loopCell = 1:numCells(loopExpt)
                fileNum = num2str(loopCell);
                allProbs = load([pathname, 'analysis/' exptPrefixes{loopExpt} groupNum 'cell' fileNum 'Probs.txt']);
                spots = load([pathname 'spots' groupNum 'ch' fileNum '.txt']);
                l2s = load([pathname 'l2s' groupNum 'ch' fileNum '.txt'])';
                
                [ guess_l2s, ~, numMismatches ] = guessConformation(allProbs, loci(:, 2), spots(:, 7), l2s);
                mismatches((loopLabeling-1)*numCells(loopExpt)+loopCell) = numMismatches;
                
                dlmwrite([pathname, 'analysis/' exptPrefixes{loopExpt} groupNum 'cell' fileNum 'guessl2s.txt'], guess_l2s, '\n');
            end
            disp([ 'finished ' exptPrefixes{loopExpt} exptLabels{loopExpt} ' ' numColors{loopColoring} ' colors x' num2str(loopLabeling) ])
        end
        dlmwrite([pathname, 'analysis/' exptPrefixes{loopExpt} exptLabels{loopExpt} 'mismatches.txt'], mismatches, '\n');
    end
end