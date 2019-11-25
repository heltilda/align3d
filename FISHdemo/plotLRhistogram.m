exptNo = 1;
chrNo = 22;
nkb = 10;
fileRange = 1:100;
msErr = 0.05;
if exptNo == 2 && nkb == 2
    residualBoundaries = [ 0 0.45; 1 inf ];
else
    residualBoundaries = [ 0 inf ];
end

if exptNo == 1
    pathname = [ '/Users/brianross/Desktop/align3d/experimentalDemo/results/Oligopaints/5color/chr' num2str(chrNo) '/' ];
    lociFileNum = [ num2str(chrNo) 'x1' ];
    if chrNo == 21
        Lslices = (.35:.1:1.5) * 1e6;
    else
        Lslices = (.35:.2:1.5) * 1e6;
    end
    RbinMax = 3;
    RbinSeparation = .1;
    gaussPeakAverageNum = 3;
else
    pathname = [ '/Users/brianross/Documents/GitHub/BXC-ORCA-data/3color/WT' num2str(nkb) 'kb/' ];
    lociFileNum = [ num2str(nkb) 'kbx1' ];
    Lslices = nkb*(0.5:1:15)*1000;
    RbinMax = 2.5;
    RbinSeparation = .01;
    gaussPeakAverageNum = 7;
end


    % load all adjacent-locus spot distances from the files

dLs = zeros(1, 0);
dRs = dLs;
for loopFile = fileRange
    fileNum = [ lociFileNum 'ch' num2str(loopFile) ];
    
    loci = load([pathname, 'loci' lociFileNum '.txt']);
    spots = load([pathname, 'spots' fileNum '.txt']);
    l2s = load([pathname, 'l2s' fileNum '.txt']);
    
    LS = [ loci(l2s ~= 0, 1), spots(l2s(l2s~=0), 1:3) ];
    numL = size(LS, 1);
    
    dists = cell(1, 4);
    for Lxyz = 1:4
        dists{Lxyz} = tril(repmat(LS(:, Lxyz), 1, numL) - repmat(LS(:, Lxyz)', numL, 1), -1);
    end
    idx = (dists{1} ~= 0);
    dLs = [ dLs, dists{1}(idx)' ];
    dRs = [ dRs, (dists{2}(idx).^2 + dists{3}(idx).^2 + dists{4}(idx).^2)'.^.5 ];
end


    % 0) bin the adjacent spot separations by contour length

numSlices = length(Lslices)-1;
sliceCenters = (Lslices(1:end-1)+Lslices(2:end))/2;

R_binEdges = 0:RbinSeparation:RbinMax;
R_binCenters = (R_binEdges(1:end-1)+R_binEdges(2:end))/2;
numRbins = length(R_binCenters);

mean_dRs = nan(1, numSlices);
idx = cell(1, numSlices);
numInBin = zeros(numSlices, numRbins);

fitM = @(x,y,w) ((sum(w)*sum(w.*x.*y) - sum(w.*x)*sum(w.*y)) / (sum(w)*sum(w.*x.^2) - sum(w.*x)^2));
linFit = @(x,y,w) ([ fitM(x,y,w), (sum(w.*y) - fitM(x,y,w)*sum(w.*x)) / sum(w) ]);

for slice = 1:numSlices
    idx{slice} = find((dLs >= Lslices(slice)) & (dLs < Lslices(slice+1)));
    numInBin(slice, :) = histcounts(dRs(idx{slice}), R_binEdges);
    [ maxNumInBin, maxIdx ] = max(numInBin(slice, :));
    [ ~, highestIdx ] = sort(numInBin(slice, :), 'descend');
    if maxNumInBin > 0
        mean_dRs(slice) = (median(R_binCenters(highestIdx(1:gaussPeakAverageNum))).^2 - 2*msErr^2).^.5;
    end
end


    % 1) find a linear fit to <R>(L), and plot it

usedSliceIdx = ~isnan(mean_dRs);
x = log10(sliceCenters(usedSliceIdx));
y = log10(mean_dRs(usedSliceIdx));
lf = linFit(x, y, mean_dRs(usedSliceIdx) ./ abs(y));
m = lf(1); b = lf(2);
pred_dRs = 10^b*sliceCenters.^m;
expConsts = 1 ./ (pred_dRs.^2 + 2*msErr^2);

figure(1), clf
set(gcf, 'Position', [0, 0, 400, 200])

subplot(1, 2, 1), hold on
plot(sliceCenters, mean_dRs, 'o')
plot(sliceCenters, 10^b*sliceCenters.^m)
xlabel('L (bp)'), ylabel('R_{peak} (\mum)')
title([ 'R_{peak} - MS error = ' num2str(10^b) ' * L^{' num2str(m) '}' ])
disp([ 'R_peak = ' num2str(10^b) ' * L^' num2str(m) ])

subplot(1, 2, 2), hold on
plot(x, y, 'o')
plot(x, m*x+b)
xlabel('log_{10} L'), ylabel('log_{10} R_{peak}')
title([ 'y = ' num2str(m) ' * x + ' num2str(b) ])


    % 2) now go through slice by slice to find the residual part of the
    % distribution

numTails = size(residualBoundaries, 1);
MBA = cell(1, numTails);
for loopTail = 1:numTails
    MBA{loopTail} = nan(numSlices, 3);
end
rescaledBinDensities = cell(1, numSlices);
smoothedBinDensities = cell(1, numSlices);
smoothedBinCenters = cell(1, numSlices);
residuals = cell(numSlices, 1);
residualIdx = cell(numSlices, 1);
tailIdx = cell(numSlices, numTails);
gaussFit = cell(1, numSlices);
w = cell(1, numSlices);

for slice = 1:numSlices
    
    
        % first smooth the p(dR) distribution across bins with no samples
        
    nonemptyBinIdx = find(numInBin(slice, :) ~= 0);
    binsToSmoothIdx = nonemptyBinIdx(1):nonemptyBinIdx(end);
    smoothedBinEdges = R_binEdges([ binsToSmoothIdx binsToSmoothIdx(end)+1 ]);
    smoothedBinMasses = numInBin(slice, nonemptyBinIdx) / sum(numInBin(slice, :));
    emptyBinIdx = find(numInBin(slice, binsToSmoothIdx) == 0);
    while ~isempty(emptyBinIdx)
        lastEmpty = length(emptyBinIdx);
        firstEmpty = lastEmpty;
        while firstEmpty > 1
            if emptyBinIdx(firstEmpty-1) ~= emptyBinIdx(firstEmpty)-1
                break;
            end
            firstEmpty = firstEmpty-1;
        end
        
        removedEdgeIdx = emptyBinIdx(firstEmpty:lastEmpty);
        smoothedBinEdges(emptyBinIdx(lastEmpty)+1) = ...
            smoothedBinEdges(emptyBinIdx(lastEmpty)+1) - RbinSeparation * length(removedEdgeIdx)/2;
        smoothedBinEdges(removedEdgeIdx) = [];
        
        emptyBinIdx(firstEmpty:lastEmpty) = [];
    end
    
    smoothedBinCenters{slice} = (smoothedBinEdges(1:end-1) + smoothedBinEdges(2:end)) / 2;
    smoothedBinDensities{slice} = smoothedBinMasses ./ diff(smoothedBinEdges);
    
    
        % next, construct the predicted Gaussian distribution, subtract it
        % and fit the residuals to a decaying exponential
    
    rescaledBinDensities{slice} = smoothedBinDensities{slice} ./ (4*pi*smoothedBinCenters{slice}.^2);
    gaussFit{slice} = (expConsts(slice)/pi)^1.5 * exp(-expConsts(slice)*smoothedBinCenters{slice}.^2);
    residuals{slice} = rescaledBinDensities{slice} - gaussFit{slice};
    w{slice} = residuals{slice} ./ rescaledBinDensities{slice}.^.5;
    
    tooManyIdx = find(residuals{slice} > 0);
    tooFewIdx = [ 0 find(residuals{slice} < 0) ];
    [ ~, residualPeakIdx ] = max(residuals{slice}(tooFewIdx(end)+1:end));
    residualPeakIdx = residualPeakIdx + tooFewIdx(end);
    residualIdx{slice} = tooManyIdx(tooManyIdx > residualPeakIdx);
    
    for loopTail = 1:numTails
        tailIdx{slice, loopTail} = residualIdx{slice}( smoothedBinCenters{slice}(residualIdx{slice}) > residualBoundaries(loopTail, 1) ...
                & smoothedBinCenters{slice}(residualIdx{slice}) < residualBoundaries(loopTail, 2) );
        if length(tailIdx{slice, loopTail}) > 1
            oneTailRs = smoothedBinCenters{slice}(tailIdx{slice, loopTail});
            oneTailResiduals = residuals{slice}(tailIdx{slice, loopTail});
            for loopTail2 = 1:loopTail-1
                oneTailResiduals = oneTailResiduals - 10.^(MBA{loopTail-1}(slice, 2) - MBA{loopTail-1}(slice, 1)*oneTailRs);
            end
            if min(oneTailResiduals) > 0
                lf = linFit(oneTailRs, log10(oneTailResiduals), w{slice}(tailIdx{slice, loopTail}));
                MBA{loopTail}(slice, :) = [ -lf(1), lf(2), sum(w{slice}(tailIdx{slice, loopTail})) ];
            end
        end
    end
end

figure(2), clf
set(gcf, 'Position', [0, 300, 400, 200])

MBAtitles = { '-m', 'b' };
mb = zeros(numTails, 2);
for loopTail = 1:numTails
    nonNans = find(~isnan(MBA{loopTail}(:, 1)))';
    MBA{loopTail}(nonNans, 3) = MBA{loopTail}(nonNans, 3) / mean(MBA{loopTail}(nonNans, 3));
    mb(loopTail, :) = median(MBA{loopTail}(nonNans, 1:2), 1);
    
    for c2 = 1:2
        subplot(numTails, 2, 2*(loopTail-1) + c2), hold on
        for slice = nonNans
            plot(sliceCenters(slice), MBA{loopTail}(slice, c2), 'ob', 'MarkerSize', 6*MBA{loopTail}(slice, 3)^.5)
        end
        plot(sliceCenters([1 end]), mb(loopTail, c2)*[1 1], 'r')
        title([ MBAtitles{c2} ' = ' num2str(mb(loopTail, c2)) ]), xlabel('L (bp)'), ylabel(MBAtitles{c2})
    end
end

ms = mb(:, 1);
bs = mb(:, 2);
As = 8*pi*10.^bs./(log(10)*ms).^3;
expNorms = (ms*log(10)).^3 / (8*pi);
gaussAmplitude = 1 - sum(As);


    % 3) sum the Gaussian and exp-tail distributions

figure(3), clf
F3dimX = 6; %3*ceil((3*numSlices)^.5/3);
F3dimY = ceil(3*numSlices/F3dimX);
set(gcf, 'Position', [0, 2000, 200*F3dimX, 150*F3dimY+50])

for slice = 1:numSlices
    
    r = smoothedBinCenters{slice};
    Vshell = (4*pi*r.^2);
    gaussWithTail = gaussAmplitude * (expConsts(slice)/pi)^1.5 * exp(-expConsts(slice)*r.^2);
    for loopTail = 1:numTails
        gaussWithTail = gaussWithTail + 10.^(bs(loopTail) - ms(loopTail)*r);
    end
    
    subplot(F3dimY, F3dimX, slice*3-2), hold on
    bar(r, smoothedBinDensities{slice})
    plot(r, gaussAmplitude * Vshell .* gaussFit{slice}, 'r-.');
    for loopTail = 1:numTails
        plot(r, As(loopTail) * expNorms(loopTail) * Vshell .* 10.^(-ms(loopTail)*r), 'g');
    end
    plot(r, Vshell .* gaussWithTail, 'r');
    title(['p(|R|; L = ' num2str(sliceCenters(slice)/1e6) ' Mb)']), xlabel('|R| (\mum)'), ylabel('freq'), xlim([0 RbinMax])
    
    subplot(F3dimY, F3dimX, slice*3-1), hold on
    barPlot = bar(r, log10(rescaledBinDensities{slice})); barPlot.BaseValue = min(log10(rescaledBinDensities{slice}))-0.5; yl = ylim;
    plot(r, log10(gaussAmplitude * gaussFit{slice}), 'r-.');
    for loopTail = 1:numTails
        plot(r, bs(loopTail) - ms(loopTail)*r, 'g');
    end
    plot(r, log10(gaussWithTail), 'r');
    title('log_{10} p(R)'), xlabel('R (\mum)'), ylabel('log_{10} freq'), ylim([ barPlot.BaseValue yl(2) ]), xlim([0 RbinMax]);
    
    subplot(F3dimY, F3dimX, slice*3), hold on
    w{slice} = w{slice} / mean(w{slice}(residualIdx{slice}));
    for cy = 1:length(residualIdx{slice})
        plot(smoothedBinCenters{slice}(residualIdx{slice}(cy)), log10(residuals{slice}(residualIdx{slice}(cy))), ...
            'ob', 'MarkerSize', 4*w{slice}(residualIdx{slice}(cy))^.5);
    end
    for loopTail = 1:numTails
        plot(smoothedBinCenters{slice}(tailIdx{slice, loopTail}), ...
            -ms(loopTail)*smoothedBinCenters{slice}(tailIdx{slice, loopTail})+bs(loopTail), 'g')
    end
    residualR = smoothedBinCenters{slice}(residualIdx{slice});
    plot(residualR, log10(sum(10.^(-ms*residualR+bs*ones(1, length(residualR))), 1)), 'k')
    title('log_{10} residual'), xlabel('R (\mum)'), ylabel('log_{10} residual'), xlim([0 RbinMax])
end

disp([ 'p(Gauss) = ' num2str(gaussAmplitude) ])
for loopTail = 1:numTails
    disp([ 'p(tail ' num2str(loopTail) ') = ' num2str(As(loopTail)) ...
        '; p(R | tail ' num2str(loopTail) ') = ' num2str(expNorms(loopTail)) '*e^(-' num2str(log(10)*ms(loopTail)) '*R)' ])
end
