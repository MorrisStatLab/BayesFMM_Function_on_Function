function U=UpdateU(beta,theta,model,D,wavespecs)
% This function update the random effect coefficients U conditional on Beta
% and Theta.
Zcat=cat_cell(model.Z);
Dstar=D-model.X*beta;
M=model.M;
U=NaN(M,wavespecs.K);
for j=1:wavespecs.K
   Mjk=diag(model.C*theta(model.H+1:model.H+model.c,j).^(-1));
   Njk=diag(repvec(theta(1:model.H,j).^(-1),model.m));
   
   A11=Zcat'*Mjk*Zcat+Njk;
   A12=Zcat'*Mjk*Dstar(:,j);
   A22=Dstar(:,j)'*Mjk*Dstar(:,j);   
   A=[A11,A12;A12',A22];          
   [junk i_sort]=sort(diag(A11)); %%%% sweep largest diagonals first    
   for ii=1:M   %#zhu#% default sort is to sort in ascending order. 
        A=sweep(A,i_sort(M-ii+1));        
   end         
   vj=A(1:M,1:M);
   muj=A(1:M,M+1);
   %U(:,j)=mvnrnd(muj',vj)';  % both random number generator works, but this one sometimes require "numerical non-singularity". 
   %U(:,j)=randmvn(muj,vj,1); 
  if (max(max(abs(vj)))<1e-8)&&(max(abs(muj))<1e-8)  %% there are cases that both muj and vj are extremely small, we do't sample in those cases.
     U(:,j)=muj;       
  else U(:,j)=randmvn(muj,vj,1); 
  end
end
