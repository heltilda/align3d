exptFolders = { ...
    '~/Desktop/align3d/experimentalDemo/results/Oligopaints/', ...
    '~/Documents/GitHub/BXC-ORCA-data/'  };
exptNames = { 'chr21', 'WT2kb' };
exptColors = { [ 3, 5, 10 ], [ 3, 5, 10, 20 ] };
Iidx = 14;
Icutoff = [ 5 5 5 nan; 6 6 5 3 ];

figure(1), clf
for loopExperiment = 1:length(exptNames)
    for loopColor = 1:length(exptColors{loopExperiment})
        colorString = num2str(exptColors{loopExperiment}(loopColor));
        summaryFile = load([ exptFolders{loopExperiment} colorString 'color/' ...
            exptNames{loopExperiment} '/analysis/' exptNames{loopExperiment} 'summary.txt' ]);
        
        labeling1 = summaryFile(summaryFile(:, 1) == 1 & summaryFile(:, 3) == 0, Iidx);
        labeling2 = summaryFile(summaryFile(:, 1) == 2 & summaryFile(:, 3) == 0, Iidx);
        
        plottable1 = labeling1(labeling1 < Icutoff(loopExperiment, loopColor));
        plottable2 = labeling2(labeling2 < Icutoff(loopExperiment, loopColor));
        plottableBothIdx = (labeling1 < Icutoff(loopExperiment, loopColor) & labeling2 < Icutoff(loopExperiment, loopColor));
        [ ~, histEdges ] = histcounts([ plottable1; plottable2; abs(labeling1(plottableBothIdx)-labeling2(plottableBothIdx)) ]);
        histX = [ histEdges(1) - (histEdges(3)-histEdges(2))/2 ...
            (histEdges(1:end-1)+histEdges(2:end))/2 ...
            histEdges(end) - (histEdges(end-2)-histEdges(end-1))/2 ];
        hist1 = [ 0 histcounts(plottable1, histEdges) 0 ];
        hist2 = [ 0 histcounts(plottable2, histEdges) 0 ];
        
        finiteIdx = (isfinite(labeling1) & isfinite(labeling2));
        finite1 = labeling1(finiteIdx);
        finite2 = labeling2(finiteIdx);
        expectedSigma = (std(finite1)^2 + std(finite2)^2)^.5 / length(finiteIdx)^.5;
        Nsigmas = abs(mean(finite1) - mean(finite2)) / expectedSigma;
        p = erfc(Nsigmas/2^.5);
        diff_hist = [ histcounts(abs(finite1-finite2), histEdges) 0 ];
        RMSdiff = sum((finite1-finite2).^2)^.5;
        
        subplot(length(exptNames), 4, (loopExperiment-1)*4 + loopColor), hold on
        plot(histX, hist1)
        plot(histX, hist2)
        plot(histX(2:end), diff_hist, '-.')
        xlim([ 0 max(histX) ])
        title([ exptNames{loopExperiment} ' ' colorString ' colors; p = ' num2str(p) ])
        xlabel('I'), ylabel('frequency')
        if loopExperiment == 2 && loopColor == 4
            legend('labeling 1', 'labeling 2', '|I_1 - I_2|')
        end
    end
end
set(gcf, 'Position', [0, 0, 800, 150*length(exptNames)])
