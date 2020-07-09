function [results,settings,D]=Wavelet_compress(X_raw,wavespecs,alpha,xticks,Image_matrix);

tic
if nargin<2
    wavespecs.wavelet='db4';
    wavespecs.nlevels=floor(log(min(size(X_raw{1})))/log(2));
    wavespecs.ndim=2;
end;
    
if nargin<3
    alpha=[.9,.95,.99,.995,.999];
end;

if nargin<4
    xticks=[200,300,400,500,1000,2000,3000,4000,5000,7500,10000,20000];
end;

x_min=.8;
x_max=.999;
x_step=0.001;
x_length=floor((x_max-x_min)/x_step)+1;
x_tick=(0:x_length)*x_step+x_min;



%%%%% Compute DWT for all images/functions
if wavespecs.ndim==2
    if Image_matrix==1
        X=reshape(X_raw(1,:),wavespecs.t1,wavespecs.t2);
        n=size(X_raw,1);
    else
        n=length(X_raw);
        X=X_raw{1}; 
    end;
    [C,S]=wavedec2(X,wavespecs.nlevels,wavespecs.wavelet);
    D=repmat(0,n,length(C));    
    for i=1:n
        if Image_matrix==1
            X=reshape(X_raw(i,:),wavespecs.t1,wavespecs.t2);
        else
            X=X_raw{i};
        end;
        [D(i,:),S]=wavedec2(X,wavespecs.nlevels,wavespecs.wavelet);
        i;
    end;
elseif wavespecs.ndim==1
        n=size(X_raw,1);
        X=X_raw(1,:);
        [C,S]=wavedec(X,wavespecs.nlevels,wavespecs.wavelet);
        D=repmat(0,n,length(C));    
    for i=1:n
        [D(i,:),S]=wavedec(X_raw(i,:),wavespecs.nlevels,wavespecs.wavelet);
        i;
    end;
end;
'Done with Wavelet Transforms ',toc

clear C S X X_raw;

%%%% Compute En= % energy retained 
%%%%   En_{ij}= % energy retained for function i if all coefficients D_{ij'}
%%%%                with |D_{ij}'|>|D_{ij}| kept and all others set to0
En=single(D);
total_Energy=repmat(0,1,size(D,1));
for i=1:n
    [Csort,ix]=sort(-abs(D(i,:)));
    total_Energy(i)=sum(Csort.^2);
    energy=cumsum(Csort.^2)/total_Energy(i);
    En(i,ix)=energy;
    i;
end;
'Done Computing Relative Energy ', toc

%%%% Compute Ncoeffs_{n*,j}=# coefficients with \sum_{i=1}^n (En_{ij}<p_j)
%%%%    for grid of proportion total energies p_j, j=1,...,J
%%%%
Ncoeffs=repmat(0,size(D,1),x_length+1);
for (j=1:(x_length+1));
    temp=sum(En<x_tick(j));
    for (i=1:n);    
        Ncoeffs(i,j)=sum(temp>(i-1));
    end;
    j;
end;

'Done Computing # coefficients ', toc

%%%% Now for each setting of Ncoeffs, compute the minimum total energy
%%%% preserved for 
En_min=repmat(0,size(Ncoeffs));
for j=1:(x_length+1)
    temp=sum(En<x_tick(j));
    for (i=1:n)
        keep=(temp>(i-1));
        if (sum(keep)>1)
            En_min(i,j)=min(sum(D(:,keep)'.^2)./total_Energy);
        end;
        if (sum(keep)==1)
            En_min(i,j)=min(D(:,keep)'.^2./total_Energy);
        end;
        if (sum(keep)==0)
            En_min(i,j)=0;
        end;
    end;
    j;
end;

'Done Computing Minimum Energies', toc

settings=[x_tick',x_tick',x_tick',x_tick'];
for (k=1:length(x_tick))
    min_energy=x_tick(k);
    rows=1:n;
    i=rows(sum(En_min'>min_energy)>0);
    j=(length(x_tick)+1)-sum(En_min(i,:)'>min_energy);
    compression=(diag(Ncoeffs(i,j)));
    istar=min(rows(compression==min(compression)));
    settings(k,3)=istar;
    settings(k,4)=x_tick(j(istar));
    settings(k,2)=compression(istar);
    k;
end;

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
text(min(ncoeffs),0.98,['K*=',num2str(size(D,2))],'FontSize',16)

hold off

toc

results=settings(alpha_row,:);
