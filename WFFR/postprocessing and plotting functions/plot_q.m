function [FDR_res]=plot_q(x,q,delta,alpha,beta0,xlab,ylab,main)
%%%
if nargin<3
    delta=1;
end;
if nargin<4
    alpha=0.05;
end;
if nargin<6
    xlab='x';
end
if nargin<7
    ylab='Prob(Exceed|data)';
end
if nargin<8
    main='Posterior Prob of Exceeding';
end

[phi,regions,res]=FDR_analysis(q,alpha,beta0,x);
summary2.x=x;
summary2.q=q;
summary2.phi=phi;
summary2.regions=regions;
summary2.FDR=res.FDR;
summary2.FNR=res.FNR;
summary2.sens=res.sens;
summary2.spec=res.spec;

%These folowing lines I'm adding to get the indices of the q>phi 
 
indge=[1:length(q)];
indge=indge(q>=phi);
plot(x(indge),logit(1-q(indge)),'b');
%plot(x,logit(1-q),'b')
hold on
for (j=1:size(summary2.regions,1))
     plot(x(summary2.regions(j,1):summary2.regions(j,2)),logit(1-summary2.q(summary2.regions(j,1):summary2.regions(j,2))),'g')
end;
%plot(x(summary2.x_peaks),logit(summary2.postprobs(i-1,summary2.x_peaks)),'r.')
plot(x,repmat(logit(1-summary2.phi),size(x)),'y')
xlabel(xlab,'FontSize',12);
ylabel(ylab,'FontSize',12)
ylim([-7,7])
set(gca,'YTick',logit([0.001,0.01,0.1,0.5,0.9,0.99,0.999]));
set(gca,'YTickLabel',['0.001';' 0.01';' 0.10';' 0.50';' 0.90';' 0.99';'0.999'],'FontSize',12);
title(main,'FontSize',16)
%set(gca,'XTick',xticks2);
%set(gca,'XTickLabel',xlabels2);
text(x(1),logit(1-min(summary2.q))-1.5,['\phi=',num2str(summary2.phi)],'FontSize',12)
hold off

FDR_res=summary2;
