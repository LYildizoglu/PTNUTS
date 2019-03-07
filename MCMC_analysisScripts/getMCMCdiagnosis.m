function [ parameters ] = getMCMCdiagnosis(parameters,firstIterations)

    if ~isfield(parameters,'S')
        error(['Data of MCMC sampling not found!'])

    end
    
    if isfield(parameters,'S_full')
        parameters.S.par = parameters.S_full.par;
    end

    if length(size(parameters.S.par)) == 3
        parameters.burn_in = BurnInBySequentialGeweke(permute(parameters.S.par, [2 1 3]),2);
    else
        parameters.burn_in = BurnInBySequentialGeweke(parameters.S.par',2);
    end
    
    parameters.S_full.par = parameters.S.par;
    
    if parameters.burn_in == 0
        disp([' zero index found, please check MCMCdiagnosis plot'])
        parameters.S.par = parameters.S.par(:,[firstIterations:end]);
        parameters.S.logPost = parameters.S.logPost([firstIterations:end]);
        disp([' cutting ' num2str(firstIterations) 'first iterations!' ])
    else
        parameters.S.par = parameters.S.par(:,[parameters.burn_in:end]);
        parameters.S.logPost = parameters.S.logPost([parameters.burn_in:end]);
    end        
end