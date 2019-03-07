function res = performNUTS(logPostHandleFunction, par, opt)

   % performNUTS.m uses an No-U-Turn sampler to sample
   % from an objective function
   % 'logPostHandleFunction'. The samples are proposed based on Hamiltonian
   % equations. The algorithm is based on the No-U-Turn sampler with dual 
   % averaging from Hoffman and Gelman (2014).
   %
   %
   % The options 'opt' and 'par' cover:
   % opt.theta0                  : The initial parameter points
   % par.min and par.max         : The lower and upper bounds for the
   %                               parameters. Proposed points outside this
   %                               area are getting rejected
   % par.number                  : Number of parameters
   % opt.nIterations             : Number of desired sampling iterations
   % opt.NUTS.nIterations_adapt  : Number of iterations for which an step size
   %                               adaptation should be performed
   % opt.NUTS.delta              : Desired acceptance probability for step size 
   %                               adaptation
   % opt.NUTS.g                  : Regulization scale for step size  
   %                               size adaptation gamma>0
   % opt.NUTS.kappa              : Relaxation exponent for step size adaptation
   % opt.NUTS.t0                 : Iteration offset for step size adaptation  
   % opt.NUTS.epsilon            : for a given step size the step size will be  
   %                               fixed and no adaptation will be performed
   % opt.NUTS.maxdepth           : maximal tree depth.  
   % It returns a struct 'res' covering:
   % res.par               : The Markov chain of the parameters
   % res.logPost           : The objective value corresponding to parameter
   %                         vector
   % res.epsilon           : The stepsize considered for to sample
   % res.depth             : The depth of the generated tree
   % res.momentum          : The initial momentum given the sample



% Initialization
    doDebug                 = opt.debug;
    nPar                    = par.number;
    nIter                   = opt.nIterations;
    nIter_adapt             = opt.NUTS.nIterations_adapt;
    delta                   = opt.NUTS.delta;
    theta                   = nan(nPar,nIter+1);
    theta(:,1)              = opt.theta0;
    nPar                    = par.number;
    global logPostHandle parmax parmin
    R                       = 1;
    parmax                  = par.max;
    parmin                  = par.min;
    g                       = opt.NUTS.g;
    kappa                   = opt.NUTS.kappa;
    t0                      = opt.NUTS.t0;
    logPostHandle           = logPostHandleFunction;
    epsilon                 = opt.NUTS.epsilon;
    maxdepth                = opt.NUTS.maxdepth;
   
    if doDebug     
        res.par             = nan(nPar, nIter);
        res.logPost         = nan(nIter,1);
        res.momentum        = nan(nPar, nIter);
        res.epsilon 	    = nan(nIter,1);
        res.depth           = nan(nIter,1);
    else
        res.par             = nan(nPar, nIter);
        res.logPost         = nan(nIter, 1);      
    end


% Determine initial step size if none is provided
    if isempty(epsilon)
        initeps             = false;
        epsilon             = nan(nIter+1,1);
        epsilon_dual        = nan(nIter+1,1);
        epsilon(1)          = FindReasonableEpsilon(theta(:,1));
        epsilon_dual(1)     = 1;
    else
        initeps             = true;
        epsilon             = repmat(epsilon,nIter+1,1);
        epsilon_dual        = epsilon;
    end
% Initialization    
    H                       = nan(nIter+1,1);
    mu                      = log(10*epsilon(1));
    H(1)                    = 0;
    r0                      = R*mvnrnd(zeros(1,nPar),eye(nPar),nIter)';
    logPost                 = NaN(1,nIter+1);
    logPost(1)              = PostCalculation(theta(:,1));
    msg                     = '';
    timer                   = tic; dspTime      = toc(timer);

for m=1:nIter

% Reporting Progress
    switch opt.mode
         case {'visual','text'}
            if toc(timer)-dspTime > 0.5 
               fprintf(1, repmat('\b',1,numel(msg)-2)) ;
               msg = ['Progress: ' num2str(m/(nIter)*100,'%2.2f') ' %%\n'];
               fprintf(1,msg);
               dspTime = toc(timer); 
            end
         case 'silent'
    end

% Tree building initialization, set left-most and right-most leaf is current one
    u               = unifrnd(0,exp(logPost(m)-dot(r0(:,m),r0(:,m))/2));
    theta_minus     = theta(:,m);
    theta_plus      = theta(:,m);
    r_minus         = r0(:,m);
    r_plus          = r0(:,m);
    theta(:,m+1)    = theta(:,m);
    logPost(m+1)    = logPost(m);
    j               = 0;
    n               = 1;
    s               = 1;

% increase the size of the tree until maxdepth, U-Turn or too large integration
% error is reached
    while s==1 && j<=maxdepth
        v(j+1)=(unidrnd(2)-1)*2-1;
        if v(j+1)==-1
% build tree to the left side
            [theta_minus,r_minus,~,~,theta_temp,logPost_temp,n_temp,s_temp,alpha, ...
                n_alpha]=BuildTree(theta_minus,r_minus,u,v(j+1),j, ...
                epsilon(m),theta(:,m),logPost(m),r0(:,m));
        else
% build tree to the right side
            [~,~,theta_plus,r_plus,theta_temp,logPost_temp,n_temp,s_temp,alpha, ...
                n_alpha]=BuildTree(theta_plus,r_plus,u,v(j+1),j, ...
                epsilon(m),theta(:,m),logPost(m),r0(:,m));
        end

% accept chosen state of subtree, update number of accepted states, check for 
% U-Turn and increase tree depth
        if s_temp==1 && (rand)<=min(1,n_temp/n)
            theta(:,m+1)    = theta_temp;
	        logPost(m+1)    = logPost_temp;
        end
        n                   = n+n_temp;
        s                   = s_temp*(dot((theta_plus-theta_minus),r_minus)>=0)* ...
                                (dot((theta_plus-theta_minus),r_plus)>=0);
        j                   = j+1;
    end

% if none step size is provided, update it 
    if ~initeps
        if m <= nIter_adapt
            H(m+1)          = (1-1./(m+t0))*H(m)+1./(m+t0)*(delta-alpha/n_alpha);
            epsilon(m+1)    = exp(mu-sqrt(m)/g*H(m+1));
            epsilon_dual(m+1) = exp (m ^ (-kappa) * log(epsilon(m+1))+ ...
                (1-m ^ (-kappa))* log(epsilon_dual(m)));
        else
            epsilon(m+1)    = epsilon_dual(nIter_adapt+1);
        end
    end

% Store iteration    
    if doDebug
        res.par(:,m)        = theta(:,m+1);
		res.logPost(m)      = logPost(m+1);
        res.epsilon(m)      = epsilon(m+1);
        res.depth(m)        = j-1;
        res.momentum(:,m)   = r0(:,m);
    else
        res.par(:,m)        = theta(:,m+1);
        res.logPost(m)      = logPost(m+1);
    end
end

   switch opt.mode
      case {'visual','text'}
         fprintf(1, repmat('\b',1,numel(msg)-2)) ;
      case 'silent'
   end
   
end


function epsilon=FindReasonableEpsilon(theta)
% Algorithm to find a reasonable step size, provided by Hoffman and Gelman (2014).
% For a given initial theta a step size epsilon will be returned, such that 
% the acceptance probability crosses 0.5
    epsilon = 1;
    nPar    = length(theta);
    logPost = PostCalculation(theta);
    r       = mvnrnd(zeros(1,nPar),eye(nPar))';
    [ ~,r_temp,logPost_temp] = Leapfrog(theta,r,epsilon);
    a       = 2*(exp(logPost_temp-dot(r_temp,r_temp)/2)/exp(logPost-dot(r,r)/2)>0.5)-1;
    while (exp(logPost_temp-dot(r_temp,r_temp)/2)/exp(logPost-dot(r,r)/2))^a>2^(-a)
        epsilon                 = 2^a*epsilon;
        [~,r_temp,logPost_temp] = Leapfrog(theta,r,epsilon);
    end
end


function [theta_minus,r_minus,theta_plus,r_plus,theta_temp,logPost_temp, ...
    n_temp,s_temp,alpha,n_alpha] = BuildTree(theta,r,u,v,j,epsilon,theta_0,logPost,r_0)
% Algorithm to build tree in NUTS. For a given theta with momentum r a tree of 
% depth j will be build in the direction v with step size epsilon. 
% For the step size adaptation the starting parameter of the tree with its 
% momentum is taken into account.
% The function returns the sample of the most left theta_ minus and right sample 
% theta_plus with its momentums of the built subtree, as well es the parameter 
% chosen theta_temp in this range and the number of accepted states n_temp, 
% the number of all visited states n_alpha, the cummulative acceptance 
% probability alpha, and state of stopping criterion s_temp

    if j==0
% Base Case
        [theta_temp,r_temp,logPost_temp] = Leapfrog(theta,r,v*epsilon);
        n_temp      = (u<=exp(logPost_temp-dot(r_temp,r_temp)/2));
        LambdaM     = 1000;
        s_temp      = (u<exp(logPost_temp-dot(r_temp,r_temp)/2+LambdaM));
        theta_minus = theta_temp;
        r_minus     = r_temp;
        theta_plus  = theta_temp;
        r_plus      = r_temp;
        alpha       = min(1,exp(logPost_temp-dot(r_temp,r_temp)/2 - ...
            logPost+dot(r_0,r_0)/2));
        n_alpha     = 1;
    else
% Recursion to build left and right subtrees
        [theta_minus,r_minus,theta_plus,r_plus,theta_temp,logPost_temp,n_temp,...
        s_temp, alpha_temp,n_alpha_temp] = ...
        BuildTree(theta,r,u,v,j-1,epsilon,theta_0,logPost,r_0);
        if s_temp==1
            if v==-1
                [theta_minus,r_minus,~,~,theta_temp2,logPost_temp2,n_temp2, ...
                s_temp2,alpha_temp2,n_alpha_temp2] = ...
                BuildTree(theta_minus,r_minus,u,v,j-1,epsilon,theta_0,logPost,r_0);
            else
                [~,~,theta_plus,r_plus,theta_temp2,logPost_temp2,n_temp2, ...
                s_temp2,alpha_temp2,n_alpha_temp2] = ...
                BuildTree(theta_plus,r_plus,u,v,j-1,epsilon,theta_0,logPost,r_0);
            end
        
            if (rand)<=n_temp2/(n_temp+n_temp2)
                theta_temp      = theta_temp2;
                logPost_temp    = logPost_temp2;
            end

            alpha_temp          = alpha_temp+alpha_temp2;
            n_alpha_temp        = n_alpha_temp+n_alpha_temp2;
            s_temp=s_temp2*(dot((theta_plus-theta_minus),r_minus)>=0)* ...
                (dot((theta_plus-theta_minus),r_plus)>=0);
            n_temp              = n_temp+n_temp2;
        end
        alpha                   = alpha_temp;
        n_alpha                 = n_alpha_temp;
    end
end

function [theta_temp,r_temp,logPost_temp] = Leapfrog(theta,r,epsilon)
% the Leapfrog integrator, which is based on Hockney 1970 returns for a given
% state theta and momentum r the new state theta_temp and momentum r_temp 
% after traveling the step size epsilon
    [~,gradientold]             = PostCalculation(theta);
    r_temp                      = r+(epsilon/2)*gradientold;
    theta_temp                  = theta+epsilon*r_temp;
    [logPost_temp,gradientnew]  = PostCalculation(theta_temp);
    r_temp                      = r_temp+epsilon/2*gradientnew;
end

function [varargout]=PostCalculation(theta)
% This function is introduced to make sure that the borders are handled correctly
    global logPostHandle parmax parmin
    if sum((theta<parmin)+(theta>parmax))>0
        if(nargout > 1)
            varargout{1} = -inf;
            varargout{2} = zeros(length(theta),1);
        else
            varargout{1} = -inf;
        end
    else
        if(nargout > 1)
            [llh,sllh]=logPostHandle(theta);
            varargout{1} = llh;
            varargout{2} = sllh;
        else
            llh=logPostHandle(theta);
            varargout{1} = llh;
        end
    end
end

