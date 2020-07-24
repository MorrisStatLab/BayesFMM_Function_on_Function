function res = opwavfm_eb(Y,model,wavespecsy,MCMCspecs) 
%OPWAVFM Performs ordinal probit wavelet based function-on-scalar 
%           regression. Currently only designed for three level response.
%           Allows for any number of scalar covariates.
%
%   Created: 3/29/2014
%   By: Mark John Meyer

%% get model parameters %%
[ n, p ]    = size(model.X);   %%% n=# functions, p=# covariates
model.c     = size(model.C,2);
model.p     = p;
model.n     = n;
O           = size(Y,2); %% dim of Y
model.O     = O;

%% set initial values and priors for latent variables, cuts, and beta
beta        = zeros(p,O);
cuts        = MCMCspecs.cuts;
a           = MCMCspecs.a;

%% extract MCMC specs
B           = MCMCspecs.B;
burnin      = MCMCspecs.burnin;
thin        = MCMCspecs.thin;

%% matrix to save beta samples
MCMC_beta   = NaN(B,p*O);
MCMC_cuts   = NaN(B,1);

%%
ii          = 0; % counter for MCMC samples to output
tic

%% MCMC Loop: B is desired # samples
for i = 1:(B*thin+burnin)
    
    %% (1) update latent variable
    [ yStar ]   = UpdateLV(beta, Y, model, cuts);
    
    %% (2) update cutpoints
    [ b ]       = UpdateCuts(yStar, Y, a);
    cuts(2)     = b;
    
    %% (3) reshape latent variable, project into wavelet-space
    %       at first iteration, get starting values
    %       at later iterations, update W
    yStar       = reshape(yStar,n,O);
    if i == 1
        % (3.1) use WavInitVals to get initial values on first iteration
        [ initial ]     = WavInitVals(yStar, wavespecsy, model, MCMCspecs);
        bstar           = initial.bstar;
        Vbetans         = initial.Vbetans;
        PiMat           = initial.PiMat;
        TauMat          = initial.TauMat;
        Wv              = initial.Wv;
        wavespecsy      = initial.wavespecsy;
        theta           = initial.theta;
        D               = initial.D;
        W               = initial.W;
        prior_Theta_a   = initial.prior_Theta_a;
        prior_Theta_b   = initial.prior_Theta_b;
        propsd_Theta    = initial.propsd_Theta;
        K               = wavespecsy.K;
        MCMC_alpha      = NaN(B,p*K); 

%         a_tau           = initial.a_tau;
%         b_tau           = initial.b_tau;
%         a_pi            = initial.a_pi;
%         b_pi            = initial.b_pi;
%         meanop          = initial.meanop;
%         expandop        = initial.expandop;
%         MCMC_tau        = NaN(B,p*wavespecsy.J);
%         MCMC_pi         = NaN(B,p*wavespecsy.J);
    else
        % (3.2) update W and D using new Y*
        [D, ~]          = dwt_rows(yStar, wavespecsy);
        W               = GetW(model,D);
    end
    
    %% (4) update model parameters
        % (4.1) Update beta
%         [bstar,gamma,alpha]     = UpdateBetaNoOrthog(bstar,Vbetans,PiMat,TauMat,Wv,model,wavespecsy,MCMCspecs); 
        [bstar,~,alpha]         = UpdateBetaNoOrthog(bstar,Vbetans,PiMat,TauMat,Wv,model,wavespecsy,MCMCspecs); 
        Wv.L2                   = Get_L2(bstar,Wv,wavespecsy); %#zhu#% beta is updated here, so does L2. 

        %% (4.2) Update theta. q_jk, and s_jk
        [theta,~,Vbetans,Wv]    = UpdateTheta(bstar,theta,Vbetans,Wv,D,W,model,prior_Theta_a,prior_Theta_b,propsd_Theta,wavespecsy,MCMCspecs);   

        %% (4.3) Update tau_ijk(as well as TauMat) and PI_ij(as well as PiMat)
%         tau     = Update_tau(bstar,gamma,a_tau,b_tau,meanop);
%         TauMat  = tau*expandop./Vbetans;    
%         PI      = Update_pi(gamma,a_pi,b_pi,wavespecsy);
%         PiMat   = PI*expandop;
    
    %% (5) project beta* into data-space
%     beta        = proBeta(bstar,wavespecsy,wavespecsx,pcaspecsx);
    beta        = idwt_rows(bstar,wavespecsy);
    
    %%  Record MCMC samples %%
    if   (i > burnin) && (mod(i-burnin,thin) == 0)   %% Save MCMC samples of beta in single matrix.
            ii                      = ii+1; % this is the real row number among the B samples.
            is                      = mod(ii-1,B)+1; % this is the row number in the block.

            %% save betas in data-space for easier post-processing
            MCMC_beta(is,:)         = reshape(beta',1,p*O); %#zhu#% each row=[beta(1,j,k),j,k|beta(2,j,k),j,k|...|beta(p,j,k),j,k]
            MCMC_cuts(is)           = cuts(2);

            %% save wavelet-space parameters
            MCMC_alpha(is,:)        = reshape(alpha',1,p*K);
%             MCMC_tau(is,:)          = reshape(tau,1,p*wavespecsy.J);   % each block is one j, contains p values.
%             MCMC_pi(is,:)           = reshape(PI,1,p*wavespecsy.J); % [pi_{11},...pi(1J1);pi_{21},...,pi{2J2};..],J blocks, each block has p values.
    end
    if mod(i, MCMCspecs.time_update) == 0
       fprintf('\n %d \n',i),toc;
    end
    
end

fprintf('\n Done with MCMC \n');

%% perform posterior functional inference
% [psi, pst]                  = BFDR(Beta, model.delt, model.alf);
% [MAPS, upper_CI, lower_CI]  = jointband_maps(MCMC_beta, model.alf);

%% save MCMC samples
res.MCMC_beta   = MCMC_beta;
res.MCMC_alpha  = MCMC_alpha;
res.MCMC_cuts   = MCMC_cuts;
% res.MCMC_tau    = MCMC_tau;
% res.MCMC_pi     = MCMC_pi;

%% put results into "res" %%
res.model       = model;
res.wavespecsy  = wavespecsy;
res.MCMCspecs   = MCMCspecs;

%% reshape and store model estimates %%
res.bhat        = reshape(mean(MCMC_beta),p,O);
res.Q025_bhat   = reshape(quantile(MCMC_beta,0.025),p,O);
res.Q975_bhat   = reshape(quantile(MCMC_beta,0.975),p,O);
res.cuts        = [cuts(1) mean(MCMC_cuts)];
% res.tauhat      = reshape(mean(MCMC_tau),p,wavespecsy.J);
% res.pihat       = reshape(mean(MCMC_pi)',p,wavespecsy.J);

fprintf('\n Done storing model estimates \n');

%% reshape and store PFI results %%
% res.MAPs        = reshape(MAPS,G,O);
% res.UMAPs       = reshape(upper_CI,G,O);
% res.LMAPs       = reshape(lower_CI,G,O);
% res.psi         = reshape(psi,G,O);
% res.pst         = reshape(pst,G,O);
% 
% fprintf('\n Done with PFI \n');

%% latent variable and initial values %%
res.Wv          = Wv;       % last Wv
res.yStar       = yStar;    % last yStar
res.Vbetans     = Vbetans;  % last Vbetans
res.D           = D;        % last dwt(yStar)
res.initial     = initial;
