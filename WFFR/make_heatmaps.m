function make_heatmaps(results, Scen, chrNum, startProbe, endProbe, outPath)
    %% save images 
%     %% MZ added: Heatmap of true beta surface, no error
%     figure
%     x = [180 270];
%     y = [startProbe endProbe];
%     colormap(hot)
% %     imagesc(results.b1', [-0.1, 1.5])
%     imagesc(x,y,results.b1')
%     colorbar
%     set(gca, 'Ydir', 'normal')
%     set(gca,'FontSize',14)
%     xlabel('Gestational Age (days)','FontSize',14) %note this is really "v" but for my presentation I change to time "t"
%     ylabel('Genomic Position','FontSize',14) % note this is really "t" but for my pres. I change to site "s"
%     title('True Peak Surface','FontSize',16)
%     saveas(gcf, sprintf('%sS%d_true.jpg',outPath,Scen))
    %% MZ added: Estimated beta surface
    figure
    % x = [180 270];
    x = [90 0];
    y = [startProbe endProbe];
    colormap(hot)
%     imagesc(results.bhat', [-0.1, 1.5])
    imagesc(x,y,results.bhat')
    colorbar
    set(gca, 'Ydir', 'normal')
    set(gca,'xdir','reverse')
    set(gca,'FontSize',14)
    xlabel('Days Before Delivery','FontSize',14) %note this is really "v" but for my presentation I change to time "t"
    ylabel(sprintf('Chr%d: CpG Postion Number',chrNum), 'FontSize',14) 
    title('FFR: Estimated Surface','FontSize',16)
    saveas(gcf, sprintf('%sS%d_ff.jpg',outPath,Scen))
    %print('-dtiff', sprintf('%sS%d_ff.tif',outPath,Scen), '-r350');
        %% MZ: make map of psi matrix 
    figure
    x = [90 0];
    y = [startProbe endProbe];
    colormap(hot)
    imagesc(x,y,results.psi')
    colorbar
    caxis([0 1])
    set(gca, 'Ydir', 'normal')
    set(gca,'xdir','reverse')
    set(gca,'FontSize',14)
    xlabel('Days Before Delivery','FontSize',14) %note this is really "v" but for my presentation I change to time "t"
    ylabel(sprintf('Chr%d: CpG Postion Number',chrNum), 'FontSize',14) 
    title('FFR: BFDR','FontSize',16)
    saveas(gcf, sprintf('%sS%d_BFDR.jpg',outPath,Scen))
    %print('-dtiff', sprintf('%sS%d_BFDR.tif',outPath,Scen), '-r350');
    %% MZ: make map of MAPs matrix
    figure
    x = [90 0];
    y = [startProbe endProbe];
    %colormap(flipud(hot)) % inverting the color bar
    colormap(flipud(hot))    
    imagesc(x,y,results.MAPs')
    colorbar
    caxis([0 0.5])
    set(gca, 'Ydir', 'normal')
    set(gca,'xdir','reverse')
    set(gca,'FontSize',14)
    xlabel('Days Before Delivery','FontSize',14) %note this is really "v" but for my presentation I change to time "t"
    ylabel(sprintf('Chr%d: CpG Postion Number',chrNum), 'FontSize',14) 
    title('FFR: SimBaS','FontSize',16)
    saveas(gcf, sprintf('%sS%d_SimBaS.jpg',outPath,Scen))
    %print('-dtiff', sprintf('%sS%d_SimBaS.tif',outPath,Scen), '-r350');
        %% MZ added: scalar beta surface
        
    if results.model.nScalarCov > 0    
        figure
        x = [1 size(results.betahatScalar',2)];
        y = [startProbe endProbe];
        colormap(hot)
        imagesc(x,y,results.betahatScalar')
        colorbar
        set(gca, 'Ydir', 'normal')
        set(gca,'FontSize',14)
        xlabel('Scalar Covariates','FontSize',14) %note this is really "v" but for my presentation I change to time "t"
        ylabel(sprintf('Chr%d: CpG Postion Number',chrNum), 'FontSize',14) 
        title('Estimated Scalar Beta Surface','FontSize',16)
        saveas(gcf, sprintf('%sS%d_scalarB.jpg',outPath,Scen))
        %print('-dtiff', sprintf('%sS%d_scalarB.tif',outPath,Scen), '-r350');
    end
%     %% MZ added: Combined beta scalar and time-varying covariate surface
%     figure
%     x = [180 270];
%     y = [startProbe endProbe];
%     colormap(hot)
%     imagesc(x,y,[results.betahatScalar' results.bhat'])
%     colorbar
%     set(gca, 'Ydir', 'normal')
%     set(gca, 'FontSize', 14)
%     xlabel('Scalar Covariates + Time-Varying Exposure', 'FontSize', 14)
%     ylabel('Genomic Position', 'FontSize', 16)
%     saveas(gcf, sprintf('%sS%d_allCov.jpg',outPath,Scen))     
        %% MZ: make map of LMAPs matrix
%     figure
%     %colormap(flipud(hot)) % inverting the color bar
%     colormap(hot)    
%     imagesc(results.LMAPs')
%     colorbar
%     set(gca, 'Ydir', 'normal')
%     set(gca,'FontSize',14)
%     xlabel('Gestational Age (days)','FontSize',14) %note this is really "v" but for my presentation I change to time "t"
%     ylabel('Genomic Position','FontSize',14) % note this is really "t" but for my pres. I change to site "s"
%     title('LMAPs','FontSize',16)
%             %% MZ: make map of UMAPs matrix
%     figure
%     %colormap(flipud(hot)) % inverting the color bar
%     colormap(hot)    
%     imagesc(results.UMAPs')
%     colorbar
%     set(gca, 'Ydir', 'normal')
%     set(gca,'FontSize',14)
%     xlabel('Gestational Age (days)','FontSize',14) %note this is really "v" but for my presentation I change to time "t"
%     ylabel('Genomic Position','FontSize',14) % note this is really "t" but for my pres. I change to site "s"
%     title('UMAPs','FontSize',16)
%                 %% MZ: make map of MAPs matrix
%     figure
%     %colormap(flipud(hot)) % inverting the color bar
%     colormap(hot)    
%     imagesc(results.MAPs')
%     colorbar
%     set(gca, 'Ydir', 'normal')
%     set(gca,'FontSize',14)
%     xlabel('Gestational Age (days)','FontSize',14) %note this is really "v" but for my presentation I change to time "t"
%     ylabel('Genomic Position','FontSize',14) % note this is really "t" but for my pres. I change to site "s"
%     title('MAPs','FontSize',16)
fprintf('\n Done saving heatmaps.\n \n');
end