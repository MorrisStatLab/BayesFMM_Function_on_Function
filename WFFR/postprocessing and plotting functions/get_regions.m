function [regions,m]=get_regions(q,phi,x,mean_spectrum)
%%%% regions=get_regions(q,phi)
%%%% Get regions for which posterior probabilities exceed threshold.
%%%%
%%%% Input: q: 1 x T vector containing posterior
%%%%            probabilities for some null hypothesis on a functional grid
%%%%            of length T
%%%%        phi: threshold below which posterior probabilities are flagged.
%%%%
%%%% Output: regions: m x 6 matrix containing 
%%%%        1-2) left and right endpoints of sig. regions (clock tick)
%%%%        3-4) left and right endpoints of sig. regions (m/z)
%%%%        5) mean posterior probability in interval
%%%%        6) maximum value of mean spectrum in interval
%%%%        m: number of significant regions.
%%%%

sig=1.0*(q<phi);
diff_sig=diff([0,sig],1);
n=length(q)
seq=1:length(q);
m=sum(diff_sig==1);

regions=repmat(0,m,6);
%%%The following 5 lines are changed from Jeff's 2 lines
regions(:,1)=seq(diff_sig==1);
if sig(n)==1   
regions(:,2)=[seq(diff_sig==-1)-1,n];
else regions(:,2)=seq(diff_sig==-1)-1;
end
%%%jeff's original 2 lines:  
%regions(:,1)=seq(diff_sig==1)+1; 
%regions(:,2)=seq(diff_sig==-1); 
%%%% Put m/z values for endpoints as third/fourth rows.
regions(:,3:4)=x(regions(:,1:2));

%%%% Put mean posterior probability as fifth row
for (i=1:m)
    regions(i,5)=mean(q(regions(i,1):regions(i,2)));
end;

%%%% Put height of log mean curve as sixth row
for (i=1:m)
    regions(i,6)=max(mean_spectrum(regions(i,1):regions(i,2)));
end;

