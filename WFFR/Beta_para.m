function [a1,b1]=Beta_para(mode0,var0,noninform,plt,arg)
% This function compute the parameters for Beta-prior based on the mode
% mode0, and variance var0.
% Input:
%       mode0: a matrix of size m, n, indicates the mode of Beta prior.
%       var0: a matrix of same size of PI, or a scalar. indicates the
%       variance of Beta distribution. Numerically work when var0 in(0,0.09].
%       when uniform, var0=1/12=0.0833. 
%       noninform: 1 indicates using non-informative prior, i.e. uniform
%       prior. a1==1, b1==1; 0 indicates using informative prior.
%       plt: 1 plot the prior density over arg. 0: does not plot.
%       arg: a vector of values over which to compute the densities when
%       plt=1. 


[n1,n2]=size(mode0);

if noninform==1
    a1=ones(n1,n2);
    b1=ones(n1,n2);
    
else
    c1=2*mode0-1;
    c2=-4*mode0+2;

    a0=c2.^2./var0/4;
    a1=c1.*c2./var0/2;
    a2=(c1.^2-1)./var0/4+1;

    q=a1/3-a2.^2/9;
    r=(a1.*a2-3*a0)/6-a2.^3/27;

    discrim=q.^3+r.^2;

    yy=NaN(n1,n2);
    for i=1:n1
        for j=1:n2        
            if discrim(i,j)>0
                s1=sign(r(i,j)+sqrt(discrim(i,j)))*abs((r(i,j)+sqrt(discrim(i,j))).^(1/3));
                s2=sign((r(i,j)-sqrt(discrim(i,j))))*abs((r(i,j)-sqrt(discrim(i,j))).^(1/3)); % we do this because in Matlab cubroot(-8)!=-2.       
                yy(i,j)=(s1+s2)-a2(i,j)/3;            
            else
                s1=(r(i,j)+sqrt(discrim(i,j))).^(1/3);
                s2=(r(i,j)-sqrt(discrim(i,j))).^(1/3);            
                root1=(s1+s2)-a2(i,j)/3;            
                root2=-(s1+s2)/2-a2(i,j)/3+sqrt(-3)*(s1-s2)/2;
                root3=-(s1+s2)/2-a2(i,j)/3-sqrt(-3)*(s1-s2)/2;
                yy(i,j)=max([root1,root2,root3]);
            end                    
        end
    end

    xx=var0.*(yy.^3+yy.^2);

    a1=(yy+sqrt(yy.^2-4.*xx))/2;
    b1=xx./a1;
    
    if plt==1
          figure()
          for i=1:n1
             for j=1:n2
                dens=gamma(a1(i,j)+b1(i,j))/gamma(a1(i,j))/gamma(b1(i,j))*arg.^(a1(i,j)-1).*(1-arg).^(b1(i,j)-1);             
                clf;
                plot(arg,dens,'-');            
                hold on
                densmode=gamma(a1(i,j)+b1(i,j))/gamma(a1(i,j))/gamma(b1(i,j))*mode0(i,j)^(a1(i,j)-1).*(1-mode0(i,j)).^(b1(i,j)-1);
                plot([mode0(i,j),mode0(i,j)],[0,densmode],'-k');
                meanij=a1(i,j)/(a1(i,j)+b1(i,j));
                plot([meanij,meanij],[0,gamma(a1(i,j)+b1(i,j))/gamma(a1(i,j))/gamma(b1(i,j))*meanij^(a1(i,j)-1).*(1-meanij).^(b1(i,j)-1)],'--k');
                title(['i=',int2str(i),'; j=',int2str(j)]);
                pause;
             end
          end
      
   end 
        
end


