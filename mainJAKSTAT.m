function parameters=mainJAKSTAT(k_str,sampler,sensitivityCalculation,niteration)
    dbstop if error
    k=str2num(k_str);
    rng(k);
    [resultFolderName,~,~]=fileparts(which('mainJAKSTAT.m'));
    options                          = PestoOptions();
    options.mode                     = 'text';
    options.MCMC                     = PestoSamplingOptions();
    options.MCMC.objOutNumber        = 1;
    options.MCMC.nIterations         = str2num(niteration);
    options.MCMC.mode                = 'text';
    options.MCMC.debug               = true;
    datatable         = xlsread('pnas_data_original.xls');
    amiData.t         = datatable(:,1);       % time points
    amiData.Y         = datatable(:,[2,4,6]); % measurement
    amiData.condition = [1.4,0.45];           % initial conditions 
    amiData.Sigma_Y   = NaN(size(amiData.Y)); % preallocation of variances

    % objective function
    logP = @(theta) logLikelihoodJakstat(theta, amiData,sensitivityCalculation);

    % Set required sampling options.MCMC for Parallel Tempering
    par.min     = -5 * ones(17,1);
    par.max     =  3 * ones(17,1);
    par.max(4)  =  6;
    par.max(2)  =  6;
    par.min(10) = -6;
    par.min(4)  = -3;
    par.min(2)  = -3;
    par.number  = length(par.min);
    par.name    = {'log_{10}(p1)','log_{10}(p2)','log_{10}(p3)','log_{10}(p4)','log_{10}(init_{STAT})',...
     'log_{10}(sp1)','log_{10}(sp2)','log_{10}(sp3)','log_{10}(sp4)','log_{10}(sp5)',...
     'log_{10}(offset_{tSTAT})','log_{10}(offset_{pSTAT})','log_{10}(scale_{tSTAT})','log_{10}(scale_{pSTAT})',...
     'log_{10}(\sigma_{pSTAT})','log_{10}(\sigma_{tSTAT})','log_{10}(\sigma_{pEpoR})'};
    optimum = [0.602656039696963;5.99158941975455;-0.954928696723120;... 
        -0.0111612709630796;2.99159087292026;-2.80956590680809;... 
        -0.255716320541754;-0.0765445346531297;-0.407313978699970;...
        -5.46184329322403;-0.731536104114366;-0.654123977718441;...
        -0.108667272925215;0.0100555269616438;-1.42650133555338;...
        -1.34879659859495;-1.16004385000543];

    switch sampler
        case 'PT'
            options.MCMC.samplingAlgorithm   = 'PT';
            options.MCMC.PT.nTemps           = 10;
            options.MCMC.PT.exponentT        = 1000;
            options.MCMC.PT.maxT             = 4000;
            options.MCMC.PT.alpha            = 0.51;
            options.MCMC.PT.temperatureNu    = 1e4;
            options.MCMC.PT.memoryLength     = 1;
            options.MCMC.PT.regFactor        = 1e-8;
            options.MCMC.PT.temperatureEta   = 10;
            options.MCMC.theta0              = repmat(optimum,1,options.MCMC.PT.nTemps);
            options.MCMC.sigma0              = 1e6*diag(ones(1,par.number));

        case 'NUTS'
            options.MCMC.samplingAlgorithm   = 'NUTS';
            options.MCMC.NUTS.maxdepth       = 6;                
            options.MCMC.theta0              = optimum;

        case 'PTNUTS'
            options.MCMC.samplingAlgorithm           = 'PTNUTS';
            options.MCMC.PTNUTS.nTemps               = 10;
            options.MCMC.PTNUTS.exponentT            = 1000;
            options.MCMC.PTNUTS.maxT                 = 4000;
            options.MCMC.PTNUTS.temperatureNu        = 1e4;
            options.MCMC.PTNUTS.temperatureEta       = 10;
            options.MCMC.PTNUTS.nIterations_adapt    = 2000;
            options.MCMC.PTNUTS.maxdepth             = [6,3,3,3,3,3,3,3,3,3];
            options.MCMC.theta0              = repmat(optimum,1,options.MCMC.PTNUTS.nTemps);

        otherwise
            warning('please use PT, NUTS or PTNUTS as sampler')
    end 


    starttic=tic;
    start=cputime;
    parameters = getParameterSamples(par, logP, options);
    parameters.toc = toc(starttic);
    parameters.cputime=cputime-start;
    parameters.randomseed=k;
    parameters.sensitivityCalculation = sensitivityCalculation;
    save([resultFolderName '/JAKSTAT_' sampler '_' sensitivityCalculation '_rng_'  k_str '_niter_' niteration '.mat'],...
        'parameters','options');

end
