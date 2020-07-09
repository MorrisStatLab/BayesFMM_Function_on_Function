function res=get_table(p,phi)

%%%%% Get table of FDR, FNR, sensitivity, and specificity given set of
%%%%% posterior probabilities p and threshold phi such that p>phi flagged
%%%%%
%%%%% plus ROC curve, with AUC of ROC curve and AUC_10 and AUC_20;
%%%%%
res.FNR=mean(1-p(p>=phi));
res.FDR=mean(p(p<=phi));
res.sens=sum(p(p>=phi))/sum(p);
res.spec=sum(1-p(p<=phi))/sum(1-p);

phi=(0:1000)/1000;
n=length(phi);
ROC=NaN(n,2);
for (i=1:n)
    ROC(n-i+1,1)=1-(sum(1-p(p<=phi(i))))/sum(1-p);
    ROC(n-i+1,2)=sum(p(p>=phi(i)))/sum(p);
end;
ROC=[0,0;ROC;1,1];
n=size(ROC,1);
res.ROC=ROC;
%%% Triangle rule
res.AUC=0;res.AUC_10=0;res.AUC_20=0;
for (i=1:(n-1))
    dx=ROC(i+1,1)-ROC(i,1);
    dy=ROC(i+1,2)-ROC(i,2);
    delta=ROC(i,2)*dx+.5*dy*dx;
    res.AUC=res.AUC+delta;
    if ROC(i+1,1)<.1
        ctr10=res.ROC(i+1,1);
        res.AUC_10=res.AUC_10+delta;
    end;
    if ROC(i+1,1)<.2
        ctr20=res.ROC(i+1,1);
        res.AUC_20=res.AUC_20+delta;
    end;
end;
res.AUC_10=res.AUC_10/ctr10;
res.AUC_20=res.AUC_20/ctr20;

    

