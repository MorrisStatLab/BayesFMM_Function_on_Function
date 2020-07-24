function [results,settings,D,C,keep,wavespecs]=Wavelet_compress(X_raw,wavespecs,Image_matrix,alpha,xticks,graph);

tic
    
if nargin<3
    Image_matrix=1;
end;

if nargin<2
    wavespecs.wavelet='db4';
    wavespecs.wtmode='sym';
    wavespecs.nlevels=floor(log(min(size(X_raw{1})))/log(2));
    wavespecs.ndim=2;
end;

if nargin<4
    alpha=[.9,.95,.99,.995,.999];
end;

if nargin<5
    xticks=[200,300,400,500,1000,2000,3000,4000,5000,7500,10000,20000];
end;

if nargin<6
    graph=0;
end;

x_min=min(alpha);
x_max=.999;
x_step=0.001;
x_length=floor((x_max-x_min)/x_step)+1;
x_tick=(0:x_length)*x_step+x_min;

x_tick=[.95+(0:49)*.001,(.999+(1:9)*.0001)];
yticks=[0.95,0.975,0.99,0.995,0.999,0.9999];
x_length=length(x_tick)-1;

if wavespecs.ndim==3
    if Image_matrix==1
        X=reshape(X_raw(1,:),wavespecs.t1,wavespecs.t2,wavespecs.t3);
        n=size(X_raw,1);
    else
        n=length(X_raw);
        X=X_raw{1}; 
    end;
    [C]=wavedec3(X,wavespecs.nlevels,wavespecs.wavelet);
    S=C.sizes;
    d=vectorize_3d_wavelet(C);
    D=repmat(0,n,length(d));  
    for i=1:n
        if Image_matrix==1
            X=reshape(X_raw(i,:),wavespecs.t1,wavespecs.t2,wavespecs.t3);
        else
            X=X_raw{i};
        end;
        [D(i,:)]=vectorize_3d_wavelet(wavedec3(X,wavespecs.nlevels,wavespecs.wavelet));
        i;
    end;
elseif wavespecs.ndim==2
    if Image_matrix==1
        X=reshape(X_raw(1,:),wavespecs.t1,wavespecs.t2);
        n=size(X_raw,1);
    else
        n=length(X_raw);
        X=X_raw{1}; 
    end;
    if wavespecs.rectangular==1  %%% if rectangular, then wavespecs.nlevels
        %%%% and wavespecs.boundary refer to rows, and then .nlevels2 and
        %%%% .boundary2 refer to columns;
        dwtmode(wavespecs.boundary);
        [C11,S1]=wavedec(X(1,:),wavespecs.nlevels,wavespecs.wavelet);
        C1=repmat(0,size(X,1),size(C11,2));
        C1(1,:)=C11;
        for (j=2:size(X,1));
            [C1(j,:)]=wavedec(X(j,:),wavespecs.nlevels,wavespecs.wavelet);
        end;
        dwtmode(wavespecs.boundary2);
        [C21,S2]=wavedec(C1(:,1),wavespecs.nlevels2,wavespecs.wavelet2);
        C2=repmat(0,size(C21,1),size(C1,2));
        C2(:,1)=C21;
        for (j=2:size(C1,2));
            [C2(:,j)]=wavedec(C1(:,j),wavespecs.nlevels2,wavespecs.wavelet2);
        end;
        %%%% Should be reordered in blocks of K1 x K2 
        [reorder,wavespecs.Kj]=reorder_rectangular(S1(1:(end-1)),S2(1:(end-1)));
        wavespecs.reorder=reorder;
        wavespecs.Kj1=S1(1:(end-1));
        wavespecs.Kj2=S2(1:(end-1));
        wavespecs.J1=length(wavespecs.Kj1);
        wavespecs.J2=length(wavespecs.Kj2);
        wavespecs.T1=S1(end);
        wavespecs.T2=S2(end);
        wavespecs.T=S1(end)*S2(end);
        C=reshape(C2,1,size(C2,1)*size(C2,2));
        D=repmat(0,n,length(C));  
    else
        [C,S]=wavedec2(X,wavespecs.nlevels,wavespecs.wavelet);
        D=repmat(0,n,length(C));    
    end
    for i=1:n
        if Image_matrix==1
            X=reshape(X_raw(i,:),wavespecs.t1,wavespecs.t2);
        else
            X=X_raw{i};
        end;
        if wavespecs.rectangular==1  %%% if rectangular, then wavespecs.nlevels
            dwtmode(wavespecs.boundary);
            C1=repmat(0,size(X,1),size(C2,2));
            for (j=1:size(X,1));
                [C1(j,:)]=wavedec(X(j,:),wavespecs.nlevels,wavespecs.wavelet);
            end;
            C2=repmat(0,size(C2,1),size(C2,2));
            dwtmode(wavespecs.boundary2);
            for (j=1:size(C1,2));
                [C2(:,j)]=wavedec(C1(:,j),wavespecs.nlevels2,wavespecs.wavelet2);
            end;
            %%%% Should be reordered in blocks of K1 x K2 but leave alone for
            %%%% now
            C=reshape(C2,1,size(C2,1)*size(C2,2));
            D(i,:)=C(reorder);
        else
            [D(i,:)]=wavedec2(X,wavespecs.nlevels,wavespecs.wavelet);
        end;
        i
    end;
elseif wavespecs.ndim==1
        dwtmode(wavespecs.wtmode);
        n=size(X_raw,1);
        X=X_raw(1,:);
        [C,S]=wavedec(X,wavespecs.nlevels,wavespecs.wavelet);
        D=repmat(0,n,length(C));
        wavespecs.Kj=S(1:(end-1));
        wavespecs.T=S(end);
    for i=1:n
        [D(i,:),S]=wavedec(X_raw(i,:),wavespecs.nlevels,wavespecs.wavelet);
    end;
end;
'Done with Wavelet Transforms ',toc

clear S X X_raw;

En=single(D);
total_Energy=repmat(0,1,size(D,1));
for i=1:n
    [Csort,ix]=sort(-abs(D(i,:)));
    total_Energy(i)=sum(Csort.^2);
    energy=cumsum(Csort.^2)/total_Energy(i);
    En(i,ix)=energy;
    i;
end;
'Done Computing Relative Energy ',


Ncoeffs=repmat(0,size(D,1),x_length+1);
for (j=1:(x_length+1));
    temp=sum(En<x_tick(j));
    for (i=1:n);    
        Ncoeffs(i,j)=sum(temp>(i-1));
    end;
    j;
end;

'Done Computing # coefficients ', toc

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

if (graph)
    color_wheel=[1 0 0;0 1 0;1 0 1;0 1 1;1 1 0];
    figure(4)
    ncoeffs=settings(:,2);
    plot(ncoeffs,settings(:,1),'b-','LineWidth',4)
    ylabel('Minimum % Energy Preserved','Fontsize',16)
    xlabel('Number of Coefficients (K*)','Fontsize',16)
    title('Minimum % Energy Preserved vs. Number of Coefficients','Fontsize',16)
    set(gca,'XLim',[min(ncoeffs)*0.9,max(ncoeffs)*1.05],'Xscale','log','Xtick',xticks,'Fontsize',16)
    set(gca,'YLim',[min(x_tick)-.001,max(x_tick)],'Yscale','log','Ytick',yticks,'Fontsize',16)
    hold on
    set(0,'DefaultAxesLineStyleOrder','--|--|--|--|--|--|--')
end;

alpha_row=repmat(0,length(alpha),1);
for (i=1:length(alpha))
    alpha_row(i)=sum(settings(:,1)<alpha(i))+1;
    if (graph)
        a=[ncoeffs(1:alpha_row(i));wrev(repmat(ncoeffs(alpha_row(i)),alpha_row(i),1))];
        b=[(repmat(alpha(i),alpha_row(i),1));wrev(settings(1:alpha_row(i),1))];
        plot(a,b,'Linewidth',2,'Color',color_wheel(mod(i-1,size(color_wheel,1))+1,:))
    end;
end;

if (graph)
    text(min(ncoeffs),0.98,['K*=',num2str(size(D,2))],'FontSize',16)
    hold off
end;

toc

results=settings(alpha_row,:);
keep=zeros(size(results,1),size(D,2));
for (i=1:size(results,1))
    keep(i,:)=sum(En<results(i,4))>=results(i,3)==1;
end;
