function [phi,regions,res]=FDR_analysis(q,alpha,beta0,x)

%%%% Feed in q as local FDR, so SMALL q <--> significant effect
if nargin<4
    x=1:length(beta0)
end;

temp=sort(q); %%% q-values?
temp2=cumsum(temp)./(1:length(temp)); %%% FDR as function of # sig
nsig=sum(temp2<=alpha);
%%% next 2 lines modified from Jeff's
%%%temp=[0 temp]; %%% in case nsig=0
if nsig > 0
phi=temp(nsig);
%%%% Identify significant regions on curves.
regions=get_regions(q,phi,x,beta0);
res=get_table(q,phi);
else
 error('Unexpected situation: no significant regions')
end
%ROC=get_ROC(1-q);

%plot(ROC.one_minus_spec,ROC.sens)
%xlabel('1-Specificity','FontSize',16)
%ylabel('Sensitivity','FontSize',16)
%title('ROC Curve','FontSize',24)
%text(0.5,0.5,['AUC=',num2str(ROC.AUC)],'FontSize',24)