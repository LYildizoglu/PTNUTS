% Calculates the burn-in of mcmc chains.
%
% Input:
% - chain: The mcmc chain
% - zscore: The Geweke test threshold. Default 2.
%
% Output:
% - burn_in: Iteration where the first and the last fraction of the chain
%            do not differ significantly regarding Geweke test -> Burn-In
%
% Written by Benjamin Ballnus (2016)


function burn_in = BurnInBySequentialGeweke( chain, zscore )

[nsimu,npar]=size(chain);

n  = 20;
e  = fix(5*nsimu/5);
l  = fix(e/n);
ii = 1:l:e;

z = zeros(length(ii),npar);
for i=1:length(ii)
  z(i,:) = gewekeTest(chain(ii(i):end,:));
end

% Sort z-score for Bonferroni-Holm inverse to sorting p-values
max_z = max(abs(z(:,1:npar))');
[~,idxs] = sort(max_z,'descend');
alpha2 = zscore*ones(1,length(idxs));
for i = 1:length(max_z)
  alpha2(idxs(i)) = alpha2(idxs(i)) / (length(ii) - find(idxs == i) + 1);
end
if ~isempty(find(alpha2 > max_z,1)*l); burn_in = (find(alpha2 > max_z,1)-1)*l; else burn_in = nsimu; end; 

 
