function [aa,bb]=InvGamma_param(mode1,var,plt,arg)
% This function compute the parameters of InverGamma(a,b) according to its
% mode and variance. 
% Input:
%       mode1: a matrix containing the mode of the distributions.
%       var: a matrix of same size of mode1, or a scalar, containing the variances.
%       plt: 1 or 0, plot the density or not.
%       arg: the argument to evaluate densities.
[m,n]=size(mode1);
ratios=mode1.^2./var; % ratios>=0
a0=-2-ratios; % a0<0
a1=5-2*ratios; 
a2=-4-ratios; % a2<0

q=a1/3-a2.^2/9;
r=(a1.*a2-3*a0)/6-a2.^3/27;

discrim=q.^3+r.^2;

aa=NaN(m,n);
for i=1:m
    for j=1:n        
        if discrim(i,j)>0
            s1=(r(i,j)+sqrt(discrim(i,j))).^(1/3);
            s2=sign((r(i,j)-sqrt(discrim(i,j))))*abs((r(i,j)-sqrt(discrim(i,j))).^(1/3)); % we do this because in Matlab cubroot(-8)!=-2.       
            aa(i,j)=(s1+s2)-a2(i,j)/3;            
        else
            s1=(r(i,j)+sqrt(discrim(i,j))).^(1/3);
            s2=(r(i,j)-sqrt(discrim(i,j))).^(1/3);            
            root1=(s1+s2)-a2(i,j)/3;            
            root2=-(s1+s2)/2-a2(i,j)/3+sqrt(-3)*(s1-s2)/2;
            root3=-(s1+s2)/2-a2(i,j)/3-sqrt(-3)*(s1-s2)/2;
            aa(i,j)=max([root1,root2,root3]);
        end                    
    end
end

bb=(aa+1).*mode1;

if plt==1
    figure()
    for i=1:m
         for j=1:n
             dens=bb(i,j)^aa(i,j)/gamma(aa(i,j))*arg.^(-(aa(i,j)+1)).*exp(-bb(i,j)./arg);             
             clf;
             plot(arg,dens,'-');            
             hold on
             densmode=bb(i,j)^aa(i,j)/gamma(aa(i,j))*mode1(i,j)^(-(aa(i,j)+1))*exp(-bb(i,j)/mode1(i,j));
             plot([mode1(i,j),mode1(i,j)],[0,densmode],'-k');
             meanij=bb(i,j)/(aa(i,j)-1);
             plot([meanij,meanij],[0,bb(i,j)^aa(i,j)/gamma(aa(i,j))*meanij^(-(aa(i,j)+1)).*exp(-bb(i,j)./meanij)],'--k');
             title(['i=',int2str(i),'; j=',int2str(j)]);
             pause;
         end
     end
end



