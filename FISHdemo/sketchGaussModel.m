LsToLabel = [ 1.5 2 4 ] * 1e5;
L1a = [ 0 2e5 ];
L1b = 2e5 + [ 0 .8e5 ];
L2a = 2e5 - 2e4*(0:7);
L2b = 2e5 + 2e4*(0:20);
Rs = 0:.05:1.2;

f1 = @(L) ((1/3e5)*L);
f2 = @(L) (.00149071*L.^0.5);
f = @(L) (min(f1(L), f2(L)));
fGauss = @(R, R_star) ((3/(2*pi*R_star^2)) .* exp(-3*R.^2/(2*R_star^2)));

sampleColor = @(n) ([ 1-(n-1)/2, 0, (n-1)/2 ]);

figure
set(gcf, 'Position', [0, 300, 400, 150])

subplot(1, 2, 1), hold on
ls = cell(1, 3);
for c3 = 1:3
    plot(Rs, fGauss(Rs, f(LsToLabel(c3))), 'Color', sampleColor(c3))
    ls{c3} = ['L = ' num2str(LsToLabel(c3)/1e6) ' Mb'];
%     plot(Rs, 4*pi*Rs.^2 .* fGauss(Rs, f(LsToLabel(c3))), ':')
%     plot(f(LsToLabel(c3)) + [0 0], [0 2], ':')
end
xlim([ 0 Rs(end) ])
xlabel('|{\bfR}| (\mum)'), ylabel('\rho({\bfR} | L)')
ax = gca;
holdXTick = ax.XTick;
holdXTickLabel = ax.XTickLabel;
ax.YTickLabel = {};
legend(ls)

subplot(1, 2, 2), hold on
plot(L1a/1e6, f1(L1a), 'k'), plot(L2b/1e6, f2(L2b), 'k')
plot(L1b/1e6, f1(L1b), 'k:'), plot(L2a/1e6, f2(L2a), 'k:')
for c3 = 1:3
    plot(LsToLabel(c3)/1e6, f(LsToLabel(c3)), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', sampleColor(c3))
end
ax = gca;
ax.YTick = holdXTick;
ax.YTickLabel = holdXTickLabel;
ylim([ 0 f2(L2b(end)) ])
xlabel('L (Mb)'), ylabel('|{\bfR}|_{RMS} (\mum)')
legend('R \propto L', 'R \propto L^{1/2}')