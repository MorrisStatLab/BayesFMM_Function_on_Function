function make_sim_heatmaps(results, outPath)
%% Estimated beta surface
figure
colormap(hot)
imagesc(results.bhat')
colorbar
set(gca, 'Ydir', 'normal')
set(gca,'FontSize',14)
xlabel('Time, t','FontSize',14) %note this is really "v" but for my presentation I change to time "t"
ylabel('CpG site, s', 'FontSize',14) 
title('Estimated Surface','FontSize',16)
saveas(gcf, sprintf('%sBetahat.png',outPath))
%% Bayesian false discovery rate 
figure
colormap(hot)
imagesc(results.psi')
colorbar
caxis([0 1])
set(gca, 'Ydir', 'normal')
set(gca,'FontSize',14)
xlabel('Time, t','FontSize',14) %note this is really "v" but for my presentation I change to time "t"
ylabel('CpG site, s', 'FontSize',14) 
title('Bayesian False Discovery Rate','FontSize',16)
saveas(gcf, sprintf('%sBFDR.png',outPath))
%% Simultaneous Band Scores
figure
colormap(hot)
imagesc(results.MAPs')
colorbar
caxis([0 0.5])
set(gca, 'Ydir', 'normal')
set(gca,'FontSize',14)
xlabel('Time, t','FontSize',14) %note this is really "v" but for my presentation I change to time "t"
ylabel('CpG site, s', 'FontSize',14) 
title('SimBaS','FontSize',16)
saveas(gcf, sprintf('%sSimBaS.png',outPath))
%% Results for scalar covariates        
if results.model.nScalarCov > 0    
    figure
    x = [1 size(results.betahatScalar',2)];
    colormap(hot)
    imagesc(results.betahatScalar')
    colorbar
    set(gca, 'Ydir', 'normal')
    set(gca,'FontSize',14)
    xlabel('Scalar Covariates','FontSize',14) %note this is really "v" but for my presentation I change to time "t"
    ylabel('CpG site, s', 'FontSize',14) 
    title('Estimated Scalar Beta Surface','FontSize',16)
    saveas(gcf, sprintf('%sScalar_beta.png', outPath))
end

fprintf('\n Done saving heatmaps.\n \n');
end