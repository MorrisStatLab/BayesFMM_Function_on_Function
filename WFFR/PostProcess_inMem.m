function [betahat,beta_025CI,beta_975CI,alphahat,accept_rate_theta,ghat,ghatns,Q025_ghat,Q975_ghat,thetahat,theta_025CI,theta_975CI,tauhat,pihat,...
    Sigma,uhat,u_025CI,u_975CI,g_Ut,g_Ut025,g_Ut975]=...
    PostProcess_inMem(MCMC_beta,MCMC_alpha,MCMC_flag_theta,MCMC_theta,MCMC_tau,MCMC_pi,betans,theta,model,wavespecs,get_sigma,sampleU,MCMC_U)

p=model.p;
T=wavespecs.T;
K=wavespecs.K;
H=model.H;

betahat=reshape(mean(MCMC_beta)',K,p)';
beta_025CI=reshape(quantile(MCMC_beta,0.025)',K,p)';
beta_975CI=reshape(quantile(MCMC_beta,0.975)',K,p)';
alphahat=reshape(mean(MCMC_alpha)',K,p)';
accept_rate_theta=mean(MCMC_flag_theta);

tauhat=reshape(mean(MCMC_tau),p,wavespecs.J);
pihat=reshape(mean(MCMC_pi)',p,wavespecs.J);

Q025_ghat=repmat(0,p,T);
Q975_ghat=repmat(0,p,T);
ghat=NaN(p,T);
for i=1:p     
    gsamples=idwt_rows(MCMC_beta(:,((i-1)*K+1):(i*K)),wavespecs);  %%% 18 seconds for B=100
    gtemp=quantile(gsamples,[0.025, 0.975],1);
    ghat(i,:)=mean(gsamples);
    Q025_ghat(i,:)=gtemp(1,:);
    Q975_ghat(i,:)=gtemp(2,:);
    fprintf('\n Done with function coefficient B_[i](t), i=%d (of p=%d) \n',[i,p]);
end
ghatns=idwt_rows(betans,wavespecs);%#zhu#% note that betans is not changed during MCMC, it is only an initial values.


if sampleU==1
    M=model.M;
    uhat=reshape(mean(MCMC_U)',K,M);
    u_025CI=reshape(quantile(MCMC_U,0.025)',K,M)';
    u_975CI=reshape(quantile(MCMC_U,0.975)',K,M)';
    g_Ut=NaN(M,T);
    g_Ut025=NaN(M,T);
    g_Ut975=NaN(M,T);
    for j=1:M
        gUtsample=idwt_rows(MCMC_U(:,((j-1)*K+1):(j*K)),wavespecs);
        g_Ut(j,:)=mean(gUtsample);
        g_Ut025(j,:)=quantile(gUtsample,0.025);
        g_Ut975(j,:)=quantile(gUtsample,0.975);
    end
end

thetahat=reshape(mean(MCMC_theta),size(theta,1),size(theta,2));
theta_025CI=reshape(quantile(MCMC_theta,0.05),size(theta,1),size(theta,2));
theta_975CI=reshape(quantile(MCMC_theta,0.95),size(theta,1),size(theta,2));

if (get_sigma==1)
    Sigma.td=NaN(T,T,H+model.c); 
    Sigma.diagtd=NaN(T,H+model.c);
    Sigma.Rhotd=NaN(T,T,H+model.c);    
    for h=1:model.H
        [Sigma.td(:,:,h),Sigma.diagtd(:,h),Sigma.Rhotd(:,:,h)]=GetSigma(thetahat(h,:),wavespecs);
    end
    for cc=1:model.c
    [Sigma.td(:,:,model.H+cc),Sigma.diagtd(:,model.H+cc),Sigma.Rhotd(:,:,model.H+cc)]=GetSigma(thetahat(model.H+cc,:),wavespecs);
    end
else
    Sigma=[];    
end

