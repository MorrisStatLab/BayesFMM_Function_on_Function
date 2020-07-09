function [Wstar,Kjstar]=junk(W,Kj,n);

%%% Look at what happens if we compute full W matrix using all
%%% coefficients, then just keep the floor(n/2) middle ones at
%%% each level.
%%%


nn=n;
Kjstar=Kj;
J=length(Kj)-1;
for j=1:J;
    Kjstar(J+2-j)=floor(nn/2);
    nn=ceil(nn/2)
end;
Kjstar(1)=nn;

Wstar=repmat(0,n,n);
newctr=1;
isctr=1;
ictr=1;
for (j=1:J+1);
    Wstar(isctr:(isctr+Kjstar(j)-1),:)=W((ictr+floor((Kj(j)-Kjstar(j))/2)):(ictr+floor((Kj(j)-Kjstar(j))/2)+Kjstar(j)-1),:);
    isctr=isctr+Kjstar(j);
    ictr=ictr+Kj(j);
end;

