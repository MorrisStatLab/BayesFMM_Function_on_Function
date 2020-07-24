function [a,b]=EmpBayes_Theta(theta,model,MCMCspecs)

% EmpBayes: Implement empirical Bayes method of Clyde and George (1998) 
%           to select wavelet shrinkage parameteres.
%
%  Input: 
%           theta = (H+c) x * matrix of variance component estimates (MOM)
%           model = structure containing model information, with elements:
%               X = (n x p) design matrix for fixed effects functions;
%               Z = cell array of H matrices, each containing (n x m_h) design
%                       matrix for set of random effects functions; 
%
%           delta = multiple to use in prior (= # of "datasets" of
%                                               information in prior)
%           wlevels = (J x 1) vector containing the number of wavelet
%                             coefficients per level
%
%  Output: (a,b) = inverse Gamma hyperparameters

delta=MCMCspecs.delta_theta;
H=model.H;
Z=model.Z;


a=repmat(0,size(theta));
b=repmat(0,size(theta));
for h=1:H
   a(h,:)=repmat(delta*rank(Z{h}),1,size(theta,2)); %#zhu#% As explained by Morris 2003 JASA,page 5. delta<<1 means vague prior. 
   b(h,:)=theta(h,:).*(a(h,:)+1); %#zhu#% note that even when a is very close to zero, the mean and variance does not exist, it is still proper.
end

%%%% Begin changes made 1/13/05
cc=0;
for h=(H+1):size(theta,1);
    cc=cc+1;   
    a(h,:)=repmat(delta*sum(model.C(:,cc)==1),1,size(theta,2));  %#zhu#% I am not sure about this here. 
    b(h,:)=theta(h,:).*(a(h,:)+1);
end
%%% end changes made 1/13/05





