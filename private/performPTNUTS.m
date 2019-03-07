function res = performPTNUTS( logPostHandleFunction, par, opt )
   
   % performPTNUTS.m uses an parallel tempered No-U-Turn sampler to sample
   % from an objective function
   % 'logPostHandleFunction'. The samples are proposed based on Hamiltonian
   % equations. The algorithm is based on the No-U-Turn sampler with dual 
   % averaging from Hoffman and Gelman (2014) and the parrallel tempering
   % algorithm of the PESTO framework from Stapor et al. (2018).
   %
   % The options 'opt' and 'par' cover:
   % opt.theta0                  : The initial parameter points for each of the
   %                               tempered chains
   % par.min and par.max         : The lower and upper bounds for the
   %                               parameters. Proposed points outside this
   %                               area are getting rejected
   % par.number                  : Number of parameters
   % opt.nIterations             : Number of desired sampling iterations
   % opt.PTNUTS.nTemps           : Number of tempered temperatures
   % opt.PTNUTS.exponentT        : The exponent of the power law for initial
   %                               temperatures. Higher Values lead to more
   %                               separated initial temperatures.
   % opt.PTNUTS.temperatureNu    : Control parameter for adaption decay of the
   %                               temperature adaption. 
   % opt.PTNUTS.temperatureEta   : Scaling factor for temperature adaptation
   % opt.PTNUTS.maxT             : Highest temperature considered in tempering
   % opt.PTNUTS.nIterations_adapt: Number of iterations for which an step size
   %                               adaptation should be performed
   % opt.PTNUTS.delta            : Desired acceptance probability for step size 
   %                               adaptation
   % opt.PTNUTS.g                : Regulization scale for step size  
   %                               size adaptation gamma>0
   % opt.PTNUTS.kappa            : Relaxation exponent for step size adaptation
   % opt.PTNUTS.t0               : Iteration offset for step size adaptation  
   % opt.PTNUTS.epsilon          : for a given step size the step size will be  
   %                               fixed and no adaptation will be performed
   % opt.PTNUTS.maxdepth         : maximal tree depth for each temperature
   % opt.PTNUTS.prior            : A prior distribution can be provided 
   %                               seperately to evaluate the 
   %                               the tempered posterior. If nothing 
   %                               is provided, it will an uniform 
   %                               distributed prior between par.min and 
   %                               par.max will be assumed.
   % It returns a struct 'res' covering:
   % res.par               : The Markov chain of the parameters for each 
   %                         temperature
   % res.logPost           : The objective value corresponding to parameter
   %                         vector
   % res.temperedlogPost   : The tempered objective value correspnding to 
   %                         the parameter vector for each temperature
   % res.temperedlogPostGrad: The gradient of the objective value
   % res.epsilon           : The stepsize considered for to sample
   % res.depth             : The depth of the generated tree for each temperature
   % res.momentum          : The initial momentum given the sample
   % res.temperatures      : The temperatures of all tempered chains
   % res.accSwap           : The acceptance rate of swaps between tempered chains

 
% Initialization
    global logPostHandle R parmax parmin prior
    doDebug                 = opt.debug;
    saveEach                = opt.saveEach;
    saveFileName            = opt.saveFileName;   
    nTemps                  = opt.PTNUTS.nTemps;
    nIter                   = opt.nIterations;
    theta0                  = opt.theta0;
    exponentT               = opt.PTNUTS.exponentT;
    temperatureNu           = opt.PTNUTS.temperatureNu;
    nPar                    = par.number;
    temperatureEta          = opt.PTNUTS.temperatureEta;
    prior                   = opt.PTNUTS.prior;
    nIter_adapt             = opt.PTNUTS.nIterations_adapt;
    delta                   = opt.PTNUTS.delta;
    R                       = 1;
    parmax                  = par.max;
    parmin                  = par.min;
    g                       = opt.PTNUTS.g;
    kappa                   = opt.PTNUTS.kappa;
    t0                      = opt.PTNUTS.t0;
    logPostHandle           = logPostHandleFunction;
    epsilon                 = opt.PTNUTS.epsilon;
	maxdepth                = opt.PTNUTS.maxdepth;
   

% if only one theta0 is provided, for all temperatures the same starting 
% sample will be applied
    switch size(theta0,2)
       case 1
            theta           = repmat(theta0,1,nTemps);
       case nTemps
            theta           = theta0;
       otherwise
            error('Dimension of options.theta0 is incorrect.');
    end

% if only one maximal depth is provided, for all temperatures the same depth 
% will be applied    
    switch size(maxdepth,1)
       case 1
            maxdepth        = repmat(maxdepth,1,nTemps);
       case nTemps

       otherwise
            error('Dimension of options.PTNUTS.maxdepth is incorrect.');
    end
   
    if doDebug     
        res.par                     = nan(nPar, nIter, nTemps);
        res.logPost                 = nan(nIter, nTemps);
        res.temperedlogPost         = nan(nIter, nTemps);
        res.temperedlogPostGrad     = nan(nPar, nIter, nTemps);
        res.accSwap                 = nan(nIter, nTemps-1);
        res.temperatures            = nan(nIter, nTemps);
        res.epsilon 	            = nan(nIter,nTemps);
        res.depth                   = nan(nIter,nTemps);
        res.momentum                = nan(nPar,nIter,nTemps);
    else
        res.par                     = nan(nPar, nIter);
        res.logPost                 = nan(nIter, 1);      
    end
   

% generating temperature ladder
   maxT         = opt.PTNUTS.maxT;
   T            = linspace(1,maxT^(1/exponentT),nTemps).^exponentT;
   beta         = 1./T;
   % Special case of NUTS: necessary due to linspace behavior
   if nTemps == 1
      T         = 1;
      beta      = 1;
   end

% Initialization   
   accSwap              = zeros(1,nTemps-1);
   propSwap             = zeros(1,nTemps-1);
   logPost              = nan(1,nTemps);
   temperedlogPost      = nan(1,nTemps);
   gradtemperedlogPost  = nan(nPar,nTemps);
   r0=nan(nPar,nTemps);
   for l = 1:nTemps
      [temperedlogPost(l),gradtemperedlogPost(:,l)] = ...
                PostCalculation(theta(:,l),beta(l));
      logPost(l)        = temperedlogPost(l)/beta(l);
   end


% Determine initial step size if none is provided
   if isempty(epsilon)
       initeps          = false;
       epsilon          = nan(1,nTemps);
       for l=1:nTemps
       epsilon(l)       = FindReasonableEpsilon(theta(:,l), ...
                temperedlogPost(l),gradtemperedlogPost(:,l),beta(l));
       end
       epsilon_dual     = ones(1,nTemps);
   else
       initeps          = true;
       epsilon          = repmat(epsilon,1,nTemps);
       epsilon_dual     = epsilon;
   end
   mu                   = repmat(log(10*epsilon(1)),1,nTemps);
   H                    = zeros(1,nTemps);
   msg                  = '';
   timer                = tic; 
   dspTime              = toc(timer);
   i                    = 1;
   
   % Reset Progress
   if (saveEach ~= 0) && ~isempty(saveFileName) && ...
         exist([saveFileName '.mat'],'file')==2
      switch opt.mode
         case {'visual','text'}
            disp('Restoring aborted run...')
      end
      try
         load(saveFileName);
      catch
         disp('File corrupt.');
      end
   end   

   
   % Perform MCMC
   while i <= nIter
      
      % Save progress - nice for grid computing
      if (saveEach ~= 0) && ~isempty(saveFileName) && mod(i,saveEach)==0
         save(saveFileName);
      end      
      
      % Reporting Progress
      switch opt.mode
         case {'visual','text'}
            if toc(timer)-dspTime > 0.5 
               fprintf(1, repmat('\b',1,numel(msg)-2)) ;
               msg = ['Progress: ' num2str(i/(nIter)*100,'%2.2f') ' %%\n'];
               fprintf(1,msg);
               dspTime = toc(timer); 
            end
         case 'silent'
      end

      % Do MCMC step for each temperature
      for l = nTemps:-1:1

% Tree building initialization, set left-most and right-most leaf is current one
            r0(:,l)         = R*mvnrnd(zeros(1,nPar),eye(nPar),1)';
            u(l)            = unifrnd(0,exp(logPost(l)-dot(r0(:,l),r0(:,l))/2));
            theta_minus     = theta(:,l);
            temperedlogPost_minus       = temperedlogPost(l);
            gradtemperedlogPost_minus   = gradtemperedlogPost(:,l);
            theta_plus      = theta(:,l);
            temperedlogPost_plus        = temperedlogPost(l);
            gradtemperedlogPost_plus    = gradtemperedlogPost(:,l);
            theta_0         = theta(:,l);
            temperedlogPost_0           = temperedlogPost(l);
            gradtemperedlogPost_0       = gradtemperedlogPost(:,l);
            r_minus         = r0(:,l);
            r_plus          = r0(:,l);
            j(l)            = 0;
            n               = 1;
            s               = 1;


% increase the size of the tree until maxdepth, U-Turn or too large integration
% error is reached
            while s==1 && j(l)<=maxdepth(l)
                v(j(l)+1)   = (unidrnd(2)-1)*2-1;
                if v(j(l)+1)==-1
% build tree to the left side
                [theta_minus, temperedlogPost_minus, gradtemperedlogPost_minus,...
                r_minus,~,~,~,~,theta_temp, temperedlogPost_temp, gradtemperedlogPost_temp,...
                n_temp,s_temp,alphaN,n_alpha] = ...
                BuildTree(theta_minus,temperedlogPost_minus,...
                gradtemperedlogPost_minus,r_minus,u(l),v(j(l)+1),j(l),...
                epsilon(l),theta_0,temperedlogPost_0,gradtemperedlogPost_0,...
                r0(:,l),beta(l));
                else
% build tree to the right side
                [~,~,~,~,theta_plus, temperedlogPost_plus, gradtemperedlogPost_plus,...
                r_plus,theta_temp, temperedlogPost_temp, gradtemperedlogPost_temp,...
                n_temp,s_temp,alphaN,n_alpha] = ...
                BuildTree(theta_plus,temperedlogPost_plus,...
                gradtemperedlogPost_plus,r_plus,u(l),v(j(l)+1),j(l),...
                epsilon(l),theta_0,temperedlogPost_0,gradtemperedlogPost_0,...
                r0(:,l),beta(l));
                end

 % accept chosen state of subtree, update number of accepted states, check for 
% U-Turn and increase tree depth     
                if s_temp==1 && (rand)<=min(1,n_temp/n)
                    theta(:,l)                  = theta_temp;
                    temperedlogPost(l)          = temperedlogPost_temp;
                    logPost(l)                  = temperedlogPost(l)/beta(l);
                    gradtemperedlogPost(:,l)    = gradtemperedlogPost_temp;
                end
                n       = n+n_temp;
                s       = s_temp*(dot((theta_plus-theta_minus),r_minus)>=0)* ...
                        (dot((theta_plus-theta_minus),r_plus)>=0)* ...
                        (temperedlogPost_minus>-inf)*...
                        (temperedlogPost_plus>-inf);
                j(l)    = j(l)+1;
            end

% if none step size is provided, update it     
            if ~initeps
                if i <= nIter_adapt
                    H(l) = (1-1./(i+t0))*H(l)+1./(i+t0)*(delta-alphaN/n_alpha);
                    epsilon(l) = exp(mu(l)-sqrt(i)/g*H(l));
                    epsilon_dual(l) = exp (i ^ (-kappa) * log(epsilon(l))+ ...
                            (1-i ^ (-kappa))* log(epsilon_dual(l)));
                else
                    epsilon(l)=epsilon_dual(l);
                end
            end
      end
      % Swaps between all adjacent chains as in Vousden16
      if nTemps > 1
            dBeta = beta(1:end-1) - beta(2:end);
            for l = nTemps:-1:2
                pAccSwap(l-1) = dBeta(l-1) .* (logPost(l)-logPost(l-1))';
                A(l-1) = log(rand) < pAccSwap(l-1);
                propSwap(l-1) = propSwap(l-1) + 1;
                accSwap(l-1) = accSwap(l-1) + A(l-1);
                if A(l-1)
                    theta(:,[l,l-1])    = theta(:,[l-1,l]);
                    logPost([l,l-1])    = logPost([l-1,l]);
                    temporary_post      = temperedlogPost(l-1)/beta(l-1)*beta(l);
                    temperedlogPost(l-1)= temperedlogPost(l)/beta(l)*beta(l-1);
                    temperedlogPost(l)  =temporary_post;
                    temporary_grad      = gradtemperedlogPost(l-1)/beta(l-1)*beta(l);
                    gradtemperedlogPost(l-1)= gradtemperedlogPost(l)/beta(l)*beta(l-1);
                    gradtemperedlogPost(l)=temporary_grad;
                end
            end
       end
      
      % Adaptation of the temperature values (Vousden 2016)
      if nTemps > 1
         kappaT = temperatureNu / ( i + 1 + temperatureNu ) / temperatureEta;
         dS = kappaT*(A(1:end-1)-A(2:end)); 
         dT = diff(1./beta(1:end-1));
         dT = dT .* exp(dS);
         beta(1:end-1) = 1./cumsum([1,dT]); 
      end
      
      % Store iteration
      if doDebug     
            res.par(:,i,:)                  = theta;
            res.logPost(i,:)                = logPost;
            res.temperedlogPost(i,:)        = temperedlogPost;
            res.temperedlogPostGrad(:,i,:)  = gradtemperedlogPost;
            res.accSwap(i,:)                = accSwap;
            res.temperatures(i,:)           = 1./beta;
            res.epsilon(i,:)                = epsilon;
            res.depth(i,:)                  = j-1;
            res.momentum(:,i,:)             = r0;
      else
            res.par(:,i)                    = theta(:,1);
            res.logPost(i)                  = logPost(1);         
      end
      i = i + 1; 
   end
   
   switch opt.mode
      case {'visual','text'}
         fprintf(1, repmat('\b',1,numel(msg)-2)) ;
      case 'silent'
   end
end





function epsilon=FindReasonableEpsilon(theta,temperedlogPost,gradtemperedlogPost,beta)
% Algorithm to find a reasonable step size, provided by Hoffman and Gelman (2014).
% For a given initial theta a step size epsilon will be returned, such that 
% the tempered acceptance probability crosses 0.5
    [nPar,nTemps]                     =size(theta);
    epsilon                           = ones(nTemps,1);
    r                                 = mvnrnd(zeros(1,nPar),eye(nPar))';
    [~,temperedlogPost_temp,~,r_temp] = Leapfrog(theta,gradtemperedlogPost,r,epsilon,beta);
    a = 2*(exp(temperedlogPost_temp-dot(r_temp,r_temp)/2)/exp(temperedlogPost-dot(r,r)/2)>0.5)-1;
    while (exp(temperedlogPost_temp-dot(r_temp,r_temp)/2)/exp(temperedlogPost-dot(r,r)/2))^a>2^(-a)
        epsilon                       = 2^a*epsilon;
        [~,temperedlogPost_temp,~,r_temp] = Leapfrog(theta,gradtemperedlogPost,r,epsilon,beta);
    end
end




function [theta_minus, temperedlogPost_minus, gradtemperedlogPost_minus,...
        r_minus,theta_plus, temperedlogPost_plus, gradtemperedlogPost_plus, ...
        r_plus,theta_temp, temperedlogPost_temp, gradtemperedlogPost_temp,...
        n_temp,s_temp,alpha,n_alpha] = BuildTree(theta,temperedlogPost, ...
        gradtemperedlogPost,r,u,v,j,epsilon,theta_0,temperedlogPost_0, ...
        gradtemperedlogPost_0,r_0,beta)

% Algorithm to build tree in NUTS. For a given theta with momentum r a tree of 
% depth j will be build in the direction v with step size epsilon. 
% For the step size adaptation the starting parameter of the tree with its 
% momentum is taken into account.
% The function returns the sample of the most left theta_ minus and right sample 
% theta_plus with its momentums of the built subtree, as well es the parameter 
% chosen theta_temp in this range and the number of accepted states n_temp, 
% the number of all visited states n_alpha, the cummulative acceptance 
% probability alpha, and state of stopping criterion s_temp.
% For a better performance are the known tempered posterior and its gradient 
% part of the input and output.

    if j==0
% Base Case
        [theta_temp,temperedlogPost_temp,gradtemperedlogPost_temp,r_temp] = ...
            Leapfrog(theta,gradtemperedlogPost,r,v*epsilon,beta);
        n_temp=(u<=exp(temperedlogPost_temp-dot(r_temp,r_temp)/2));
        LambdaM=1000;
        s_temp=(u<exp(temperedlogPost_temp-dot(r_temp,r_temp)/2+LambdaM));
        theta_minus=theta_temp;
        temperedlogPost_minus=temperedlogPost_temp;
        gradtemperedlogPost_minus=gradtemperedlogPost_temp;
        r_minus=r_temp;
        theta_plus=theta_temp;
        temperedlogPost_plus=temperedlogPost_temp;
        gradtemperedlogPost_plus=gradtemperedlogPost_temp;
        r_plus=r_temp;
        alpha=min(1,exp(temperedlogPost_temp-dot(r_temp,r_temp)/2 - ...
            temperedlogPost_0+dot(r_0,r_0)/2));
        n_alpha=1;
    else
% Recursion to build left and right subtrees
        [theta_minus, temperedlogPost_minus, gradtemperedlogPost_minus,...
        r_minus,theta_plus, temperedlogPost_plus, gradtemperedlogPost_plus, ...
        r_plus,theta_temp, temperedlogPost_temp, gradtemperedlogPost_temp,...
        n_temp,s_temp,alpha_temp,n_alpha_temp] = ...
            BuildTree(theta,temperedlogPost,gradtemperedlogPost,r,u,v,j-1, ...
            epsilon,theta_0,temperedlogPost_0,gradtemperedlogPost_0,r_0,beta);
        if s_temp==1
            if v==-1
    [theta_minus, temperedlogPost_minus, gradtemperedlogPost_minus,...
    r_minus,~,~,~,~,theta_temp2, temperedlogPost_temp2, gradtemperedlogPost_temp2,...
    n_temp2,s_temp2,alpha_temp2,n_alpha_temp2] = ...
        BuildTree(theta_minus,temperedlogPost_minus,gradtemperedlogPost_minus, ...
        r_minus,u,v,j-1,epsilon,theta_0,temperedlogPost_0,gradtemperedlogPost_0,r_0,beta);
            else
    [~,~,~,~,theta_plus, temperedlogPost_plus, gradtemperedlogPost_plus,...
    r_plus,theta_temp2, temperedlogPost_temp2, gradtemperedlogPost_temp2,...
    n_temp2,s_temp2,alpha_temp2,n_alpha_temp2] = ...
        BuildTree(theta_plus,temperedlogPost_plus,gradtemperedlogPost_plus, ...
        r_plus,u,v,j-1,epsilon,theta_0,temperedlogPost_0,gradtemperedlogPost_0,r_0,beta);
            end

            if (rand)<=n_temp2/(n_temp+n_temp2)
                theta_temp  = theta_temp2;
                temperedlogPost_temp = temperedlogPost_temp2;
                gradtemperedlogPost_temp = gradtemperedlogPost_temp2;
            end
            alpha_temp      = alpha_temp+alpha_temp2;
            n_alpha_temp    = n_alpha_temp+n_alpha_temp2;
            s_temp = s_temp2*(dot((theta_plus-theta_minus),r_minus)>=0)* ...
                        (dot((theta_plus-theta_minus),r_plus)>=0)* ...
                        (temperedlogPost_minus>-inf)* ...
                        (temperedlogPost_plus>-inf);
            n_temp          = n_temp+n_temp2;
        end
        alpha               = alpha_temp;
        n_alpha             = n_alpha_temp;
    end
end







function [theta_temp,logPost_temp,gradientnew,r_temp] = ...
    Leapfrog(theta,gradientold,r,epsilon,beta)
% the Leapfrog integrator, which is based on Hockney 1970 returns for a given
% state theta and momentum r the new state theta_temp and momentum r_temp 
% after traveling the step size epsilon. For a better performance the known
% tempered posteriors and its gradients are handled in the in- and output
    r_temp                      = r+(epsilon/2)*gradientold;
    theta_temp                  = theta+epsilon*r_temp;
    [logPost_temp,gradientnew]  = PostCalculation(theta_temp,beta);
    r_temp                      = r_temp+epsilon/2*gradientnew;
end




function [varargout]=PostCalculation(theta,beta)
% This function is introduced to make sure that the borders are handled
% correctly and the possibility to provide a prior, which is not part of the
% tempering process can be provided.
    global logPostHandle parmax parmin prior
    if sum((theta<parmin)+(theta>parmax))>0
        if(nargout > 1)
            varargout{1} = -inf;
            varargout{2} =zeros(length(theta),1);
        else
            varargout{1} = -inf;
        end
    else
        if isempty(prior)
            if(nargout > 1)
                [llh,sllh]=logPostHandle(theta);
                varargout{1} = llh*beta;
                varargout{2} = sllh*beta;
            else
                llh=logPostHandle(theta);
                varargout{1} = llh*beta;
            end
        else
            if(nargout > 1)
                [llh,sllh]=logPostHandle(theta);
                [lprior,slprior]=prior(theta);
                varargout{1} = llh*beta+lprior;
                varargout{2} = sllh*beta+slprior;
            else
                llh=logPostHandle(theta);
                lprior=prior(theta);
                varargout{1} = llh*beta+lprior;
            end
        end
    end
end

