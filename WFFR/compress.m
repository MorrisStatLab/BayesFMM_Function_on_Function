xticks=[500,1000,2000,5000,10000,20000,50000,100000,200000];
alpha=[.9,.95,.977,.989,.994]
[results,settings,D]=Wavelet_compress(X2,wavespecs);


%%% Figure 1 for paper
color_wheel=[1 0 0;0 1 0;1 0 1;0 1 1;1 1 0];
figure(4)
ncoeffs=settings(:,2);
plot(ncoeffs,settings(:,1),'b-','LineWidth',4)
ylabel('Minimum % Energy Preserved','Fontsize',16)
xlabel('Number of Coefficients (K*)','Fontsize',16)
title('Minimum % Energy Preserved vs. Number of Coefficients','Fontsize',16)
set(gca,'XLim',[min(ncoeffs)*0.9,max(ncoeffs)*1.05],'Xscale','log','Xtick',xticks,'Fontsize',16)
hold on
set(0,'DefaultAxesLineStyleOrder','--|--|--|--|--|--|--')
alpha_row=repmat(0,length(alpha),1);
for (i=1:length(alpha))
    alpha_row(i)=sum(settings(:,1)<alpha(i))+1;
    a=[ncoeffs(1:alpha_row(i));wrev(repmat(ncoeffs(alpha_row(i)),alpha_row(i),1))];
    b=[repmat(alpha(i),alpha_row(i),1);wrev(settings(1:alpha_row(i),1))];
    plot(a,b,'Linewidth',2,'Color',color_wheel(mod(i-1,size(color_wheel,1))+1,:))
end;
text(min(ncoeffs),0.98,['T*=10,634'],'FontSize',16)
hold off


%%% Figure 2 for paper
%%%% Plot image #1 with wavelet threshold compressed version.
i=0;
i=i+1;
selected=repmat(0,1,wavespecs.K);
selected(int32(wavespecs.DIndex)+1)=1;
[a,b]=wavedec2(reshape(Y(i,:),wavespecs.t1,wavespecs.t2),6,'db4');
Yreg=waverec2(a.*selected,b,'db4');


Effect_size_indicator=2; %%% 1.25, 1.5, 1.75, 2
choose_effect=3;
zlim=[0,5];
figure(1)
colormap(jet)
subplot(2,4,1)
imagesc(reshape(2.^Y1,t1,t2),zlim);
%title('Uncompressed Gel Image, T=556,206','FontSize',12)
title('Uncompressed Image','FontSize',16)
text(50,50,'T=556,216','FontSize',16,'Color','white')
ylabel('Gel Image 1','FontSize',16)
subplot(2,4,2)
imagesc(2.^Yreg_99,zlim);
%title('Wavelet Compressed (P=0.99), T*=26,520','FontSize',12)
title('P=0.99','FontSize',16)
text(50,50,'T*=26,520','FontSize',16,'Color','white')
subplot(2,4,3)
imagesc(2.^Yreg_975,zlim);
%title('Wavelet Compressed (P=0.975), T*=10,634','FontSize',12)
title('P=0.975','FontSize',16)
text(50,50,'T*=10,634','FontSize',16,'Color','white')
subplot(2,4,4)
imagesc(2.^Yreg_95,zlim);
%title('Wavelet Compressed (P=0.95), T*=4,958','FontSize',12)
title('P=0.95','FontSize',16)
text(50,50,'T*=4,958','FontSize',16,'Color','white')
subplot(2,4,6)
%%% posterior probability of at least 1.5 fold in effect.
local_fdr=reshape(beta_peffects_99(2,((choose_effect-1)*wavespecs.t1*wavespecs.t2+1):(choose_effect*wavespecs.t1*wavespecs.t2)),wavespecs.t1,wavespecs.t2);
imagesc(local_fdr,[0,1])
title('P=0.99','FontSize',16)
subplot(2,4,7)
%%% posterior probability of at least 1.5 fold in effect.
local_fdr=reshape(beta_peffects_975(2,((choose_effect-1)*wavespecs.t1*wavespecs.t2+1):(choose_effect*wavespecs.t1*wavespecs.t2)),wavespecs.t1,wavespecs.t2);
imagesc(local_fdr,[0,1])
title('P=0.975','FontSize',16)
subplot(2,4,8)
%%% posterior probability of at least 1.5 fold in effect.
local_fdr=reshape(beta_peffects_95(2,((choose_effect-1)*wavespecs.t1*wavespecs.t2+1):(choose_effect*wavespecs.t1*wavespecs.t2)),wavespecs.t1,wavespecs.t2);
imagesc(local_fdr,[0,1])
title('P=0.95','FontSize',16)

%%%% Figure 3 for paper
Ymean=reshape(mean(Y),t1,t2);
Ymean_wfmm975=reshape(beta_mean(1,:),t1,t2);
zlim=[-0.7,3];
figure(3)
subplot(1,2,1)
imagesc(Ymean,zlim)
hold on
plot(spotlist_rev(:,1),spotlist_rev(:,2),'wx','MarkerSize',ms)
plot(spotlist_rev(flagged==1,1),spotlist_rev(flagged==1,2),'wo','MarkerSize',ms)
hold off
title('Raw Mean gel','FontSize',16) 
subplot(1,2,2)
imagesc(Ymean_wfmm975,zlim)
title('Model-Based Mean gel: M(t_1, t_2)','FontSize',16)
hold on
plot(spotlist_rev(:,1),spotlist_rev(:,2),'wx','MarkerSize',ms)
plot(spotlist_rev(flagged==1,1),spotlist_rev(flagged==1,2),'wo','MarkerSize',ms)
hold off

%%% Figure 4 for paper

%%% Mean gel and effect gel
figure(4)
subplot(2,2,1)
mean_gel=reshape(beta_mean_975(1,:),wavespecs.t1,wavespecs.t2);
imagesc(mean_gel,[0,5])
%hold on
%plot(spotlist_rev(:,1),spotlist_rev(:,2),'wx','MarkerSize',ms)
%plot(spotlist_rev(flagged==1,1),spotlist_rev(flagged==1,2),'wo','MarkerSiz
%e',ms)
%hold off
title('M(t_1,t_2)','FontSize',16)
subplot(2,2,2)
choose_effect=3;
effect_gel=reshape(beta_mean_975(choose_effect,:),wavespecs.t1,wavespecs.t2);
imagesc(effect_gel,[-1,1])
title('C_{13}(t_1,t_2)','FontSize',16)
subplot(2,2,3)
Effect_size_indicator=2; %%% 1.25, 1.5, 1.75, 2
local_fdr=reshape(beta_peffects_975(Effect_size_indicator,((choose_effect-1)*wavespecs.t1*wavespecs.t2+1):(choose_effect*wavespecs.t1*wavespecs.t2)),wavespecs.t1,wavespecs.t2);
imagesc(local_fdr,[0,1])
title('Probability of 1.5-fold, p^{1.5}(t_1,t_2)','FontSize',16)
subplot(2,2,4)
temp=(sort(1-reshape(local_fdr,1,wavespecs.t1*wavespecs.t2)));
seq=1:wavespecs.t1*wavespecs.t2;
FDR_alpha=.10
FDR_thresh=temp(sum(cumsum(temp)./seq<FDR_alpha))
imagesc(local_fdr>1-FDR_thresh)
title('Flagged Regions, FDR<0.10, 1.5-fold','FontSize',16)

%%% significance matches detected spot
xl=[390,420];
yl=[465,500];
ms=8;
figure(5)
subplot(2,2,1)
mean_gel=reshape(beta_mean_975(1,:),wavespecs.t1,wavespecs.t2);
imagesc(mean_gel,[0,3])
hold on
plot(spotlist_rev(:,1),spotlist_rev(:,2),'wx','MarkerSize',ms)
plot(spotlist_rev(flagged==1,1),spotlist_rev(flagged==1,2),'wo','MarkerSize',ms)
hold off
title('M(t_1,t_2)','FontSize',16)
xlim(xl),ylim(yl)
subplot(2,2,2)
choose_effect=3;
effect_gel=reshape(beta_mean_975(choose_effect,:),wavespecs.t1,wavespecs.t2);
imagesc(effect_gel)
hold on
plot(spotlist_rev(:,1),spotlist_rev(:,2),'wx','MarkerSize',ms)
plot(spotlist_rev(flagged==1,1),spotlist_rev(flagged==1,2),'wo','MarkerSize',ms)
hold off
title('C_{13}(t_1,t_2)','FontSize',16)
xlim(xl),ylim(yl)
subplot(2,2,3)
Effect_size_indicator=2; %%% 1.25, 1.5, 1.75, 2
local_fdr=reshape(beta_peffects_975(Effect_size_indicator,((choose_effect-1)*wavespecs.t1*wavespecs.t2+1):(choose_effect*wavespecs.t1*wavespecs.t2)),wavespecs.t1,wavespecs.t2);
imagesc(local_fdr,[0,1])
hold on
plot(spotlist_rev(:,1),spotlist_rev(:,2),'wx','MarkerSize',ms)
plot(spotlist_rev(flagged==1,1),spotlist_rev(flagged==1,2),'wo','MarkerSize',ms)
hold off
title('Probability of 1.5-fold, p^{1.5}(t_1,t_2)','FontSize',16)
xlim(xl),ylim(yl)
subplot(2,2,4)
temp=(sort(1-reshape(local_fdr,1,wavespecs.t1*wavespecs.t2)));
seq=1:wavespecs.t1*wavespecs.t2;
FDR_alpha=.10
FDR_thresh=temp(sum(cumsum(temp)./seq<FDR_alpha))
imagesc(local_fdr>1-FDR_thresh)
hold on
plot(spotlist_rev(:,1),spotlist_rev(:,2),'wx','MarkerSize',ms)
plot(spotlist_rev(flagged==1,1),spotlist_rev(flagged==1,2),'wo','MarkerSize',ms)
hold off
title('Flagged Regions, FDR<0.10, 1.5-fold','FontSize',16)
xlim(xl),ylim(yl)

%%% significance on part of spot
xl=[120,200];
yl=[390,460];
xlim(xl),ylim(yl)
ms=8;
figure(6)
subplot(2,2,1)
mean_gel=reshape(beta_mean_975(1,:),wavespecs.t1,wavespecs.t2);
imagesc(mean_gel,[0,5])
hold on
plot(spotlist_rev(:,1),spotlist_rev(:,2),'wx','MarkerSize',ms)
plot(spotlist_rev(flagged==1,1),spotlist_rev(flagged==1,2),'wo','MarkerSize',ms)
hold off
title('M(t_1,t_2)','FontSize',16)
xlim(xl),ylim(yl)
subplot(2,2,2)
choose_effect=3;
effect_gel=reshape(beta_mean_975(choose_effect,:),wavespecs.t1,wavespecs.t2);
imagesc(effect_gel)
hold on
plot(spotlist_rev(:,1),spotlist_rev(:,2),'kx','MarkerSize',ms)
plot(spotlist_rev(flagged==1,1),spotlist_rev(flagged==1,2),'ko','MarkerSize',ms)
hold off
title('C_{13}(t_1,t_2)','FontSize',16)
xlim(xl),ylim(yl)
subplot(2,2,3)
Effect_size_indicator=2; %%% 1.25, 1.5, 1.75, 2
local_fdr=reshape(beta_peffects_975(Effect_size_indicator,((choose_effect-1)*wavespecs.t1*wavespecs.t2+1):(choose_effect*wavespecs.t1*wavespecs.t2)),wavespecs.t1,wavespecs.t2);
imagesc(local_fdr,[0,1])
hold on
plot(spotlist_rev(:,1),spotlist_rev(:,2),'wx','MarkerSize',ms)
plot(spotlist_rev(flagged==1,1),spotlist_rev(flagged==1,2),'wo','MarkerSize',ms)
hold off
title('Probability of 1.5-fold, p^{1.5}(t_1,t_2)','FontSize',16)
xlim(xl),ylim(yl)
subplot(2,2,4)
temp=(sort(1-reshape(local_fdr,1,wavespecs.t1*wavespecs.t2)));
seq=1:wavespecs.t1*wavespecs.t2;
FDR_alpha=.10
FDR_thresh=temp(sum(cumsum(temp)./seq<FDR_alpha))
imagesc(local_fdr>1-FDR_thresh)
hold on
plot(spotlist_rev(:,1),spotlist_rev(:,2),'wx','MarkerSize',ms)
plot(spotlist_rev(flagged==1,1),spotlist_rev(flagged==1,2),'wo','MarkerSize',ms)
hold off
title('Flagged Regions, FDR<0.10, 1.5-fold','FontSize',16)
xlim(xl),ylim(yl)


%%% another region where significance is near spot
xl=[680,800];
yl=[260,360];
xlim(xl),ylim(yl)

xl=[90,150];
yl=[550,575];
xlim(xl),ylim(yl)
