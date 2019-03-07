function parameters=mainJAKSTAT(k_str,sampler,sampoptions,niteration)
    dbstop if error
    k=str2num(k_str);
    rng(k);
ym=[0, 5.6060, 3.7180, 9.1150, 9.4730, 13.407, 11.413, 16.251, 13.674, 20.759, 18.196, 23.998, 20.412, 27.290, 23.860, 29.852, 26.679, 29.470, 30.557, 28.206, 33.831, 30.339, 34.782, 33.380, 38.658, 34.839, 43.928, 37.430, 42.167, 40.797, 39.639, 41.158, 40.553, 44.181, 41.428, 46.425, 42.381, 48.313, 44.116, 46.649, 45.211, 47.126, 45.804, 48.210, 45.160, 50.637, 47.403, 53.320, 48.217, 52.555, 48.301, 54.721, 51.334, 50.908, 52.631, 53.549, 50.689, 51.011, 56.259, 53.338, 59.252, 53.006, 60.551, 54.324, 61.331, 55.480, 56.693, 56.566, 55.347, 61.868, 56.979, 62.119, 57.655, 63.279, 57.547, 62.750, 58.245, 60.214, 60.739, 59.195, 61.140, 58.735, 62.918, 59.539, 62.872, 60.759, 64.860, 62.182, 64.782, 61.456, 64.754, 60.671, 61.316, 62.345, 60.914, 62.374, 62.997, 60.783, 61.288, 60.789, 63.226, 58.914, 60.298, 58.257, 57.768, 59.015, 57.314, 61.219, 57.960, 60.024, 58.954, 61.306, 59.615, 64.139, 60.609, 66.178, 59.914, 66.168, 59.455, 66.784, 61.759, 65.957, 60.431, 62.504, 61.800, 61.057, 65.056, 61.657, 65.165, 60.781, 64.104, 60.499, 65.724, 61.099, 65.728, 62.024, 65.761, 59.438, 65.520, 60.641, 63.443, 60.343, 62.660, 60.857, 59.425, 58.236, 60.201, 59.923, 58.400, 65.513]';
t=[1.75000, 1.91667, 2.08333, 2.25000, 2.41667, 2.58333, 2.75000, 2.91667, 3.08333, 3.25000, 3.41667, 3.58333, 3.75000, 3.91667, 4.08333, 4.25000, 4.41667, 4.58333, 4.75000, 4.91667, 5.08333, 5.25000, 5.41667, 5.58333, 5.75000, 5.91667, 6.08333, 6.25000, 6.41667, 6.58333, 6.75000, 6.91667, 7.08333, 7.25000, 7.41667, 7.58333, 7.75000, 7.91667, 8.08333, 8.25000, 8.41667, 8.58333, 8.75000, 8.91667, 9.08333, 9.25000, 9.41667, 9.58333, 9.75000, 9.91667, 10.0833, 10.2500, 10.4167, 10.5833, 10.7500, 10.9167, 11.0833, 11.2500, 11.4167, 11.5833, 11.7500, 11.9167, 12.0833, 12.2500, 12.4167, 12.5833, 12.7500, 12.9167, 13.0833, 13.2500, 13.4167, 13.5833, 13.7500, 13.9167, 14.0833, 14.2500, 14.4167, 14.5833, 14.7500, 14.9167, 15.0833, 15.2500, 15.4167, 15.5833, 15.7500, 15.9167, 16.0833, 16.2500, 16.4167, 16.5833, 16.7500, 16.9167, 17.0833, 17.2500, 17.4167, 17.5833, 17.7500, 17.9167, 18.0833, 18.2500, 18.4167, 18.5833, 18.7500, 18.9167, 19.0833, 19.2500, 19.4167, 19.5833, 19.7500, 19.9167, 20.0833, 20.2500, 20.4167, 20.5833, 20.7500, 20.9167, 21.0833, 21.2500, 21.4167, 21.5833, 21.7500, 21.9167, 22.0833, 22.2500, 22.4167, 22.5833, 22.7500, 22.9167, 23.0833, 23.2500, 23.4167, 23.5833, 23.7500, 23.9167, 24.0833, 24.2500, 24.4167, 24.5833, 24.7500, 24.9167, 25.0833, 25.2500, 25.4167, 25.5833, 25.7500, 25.9167, 26.0833, 26.2500, 26.4167, 26.5833]';
    [resultFolderName,~,~]=fileparts(which('mainMRNA.m'));
    options                          = PestoOptions();
    options.mode                     = 'text';
    options.obj_type                       = 'log-posterior';
    options.MCMC                     = PestoSamplingOptions();
    options.MCMC.objOutNumber        = 1;
    options.MCMC.nIterations         = str2num(niteration);
    options.MCMC.mode                = 'text';
    options.MCMC.debug               = true;

    % objective function
    objectiveFunction = @(theta) logLikelihoodTransfection(theta, t, ym);

    % Set required sampling options.MCMC for Parallel Tempering
      parameters.min    = [-2; -5; -5; -5; -2];
      parameters.max    = [log10(max(t)); 5; 5; 5; 2];
      parameters.number = 5;
      parameters.name   = {'log_{10}(t_0)';'log_{10}(k_{TL}*m_0)';...
         'log_{10}(\beta)';'log_{10}(\delta)';'log_{10}(\sigma)'};


    switch sampler
        case 'PT'
        switch sampoptions
            case 'adapted'
                options.MCMC.samplingAlgorithm   = 'PT';
                options.MCMC.PT.nTemps           = 30;
                options.MCMC.PT.exponentT        = 1000;
                options.MCMC.PT.maxT             = 2000;
                options.MCMC.PT.alpha            = 0.51;
                options.MCMC.PT.temperatureNu    = 1e4;
                options.MCMC.PT.memoryLength     = 1;
                options.MCMC.PT.regFactor        = 1e-8;
                options.MCMC.PT.temperatureEta   = 10;
            case 'default'
                options.MCMC.samplingAlgorithm   = 'PT';
                options.MCMC.PT.nTemps           = 30;

            otherwise
                warning('please use default or adapted as sampler option')
            end
            options.MCMC.sigma0                 = 1e0*diag(ones(1,parameters.number));
            options.MCMC.theta0 = repmat([0.15,1.08,-2.8,-0.9,0.4],options.MCMC.PT.nTemps,1);
            options.MCMC.theta0(1:2:end) = ...
        repmat([0.15,1.08,-0.9,-2.8,0.4],floor(options.MCMC.PT.nTemps/2),1);
            options.MCMC.theta0                 = options.MCMC.theta0';

        case 'NUTS'
            options.MCMC.samplingAlgorithm      = 'NUTS';
            options.MCMC.NUTS.nIterations_adapt = 2000;        
            if mod(k,2)==0
                options.MCMC.theta0             = [0.2,1.08,-2.5,-0.7,0.4]';
            else
                options.MCMC.theta0             = [0.2,1.08,-0.7,-2.5,0.4]';
            end

        case 'PTNUTS'
        switch sampoptions
            case 'adapted'
                options.MCMC.samplingAlgorithm           = 'PTNUTS';
                options.MCMC.PTNUTS.nTemps               = 10;
                options.MCMC.PTNUTS.exponentT            = 1000;
                options.MCMC.PTNUTS.maxT                 = 2000;
                options.MCMC.PTNUTS.temperatureNu        = 1e4;
                options.MCMC.PTNUTS.temperatureEta       = 10;
                options.MCMC.PTNUTS.nIterations_adapt    = 2000;
                options.MCMC.PTNUTS.maxdepth             = [10,6,4,3,3,3,3,3,3,3];
                options.MCMC.theta0              = repmat(optimum,1,options.MCMC.PTNUTS.nTemps);
            case 'default'
                options.MCMC.samplingAlgorithm          = 'PTNUTS';
                options.MCMC.PTNUTS.nIterations_adapt   = 2000;
                options.MCMC.PTNUTS.maxdepth            = 10;
                options.MCMC.PTNUTS.nTemps              = 10;

        otherwise
            warning('please use default or adapted as sampler option')
        end
        options.MCMC.theta0 = repmat([0.15,1.08,-2.8,-0.9,0.4],options.MCMC.PTNUTS.nTemps,1);
        options.MCMC.theta0(1:2:end) = ...
            repmat([0.15,1.08,-0.9,-2.8,0.4],floor(options.MCMC.PTNUTS.nTemps/2),1);
        options.MCMC.theta0 = options.MCMC.theta0';
    end 


    starttic=tic;
    start=cputime;
    parameters = getParameterSamples(parameters, objectiveFunction, options);
    parameters.toc = toc(starttic);
    parameters.cputime=cputime-start;
    parameters.randomseed=k;
    save([resultFolderName '/MRNA_' sampler '_' sampoptions '_rng_'  k_str '_niter_' niteration '.mat'],...
        'parameters','options');

end
