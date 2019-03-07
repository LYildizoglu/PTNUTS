% Performs a multivariate, pairwise Brooks-Gelman-Rubin & Geweke test for
% checking similarity of mcmc runs. Groups together runs, which pass both
% tests.
%
% Input:
% - path: Path where raw mcmc results with names '<label>1.mat', '<label>2.mat',... are
%         stored
% - ids: Array of integers corresponding to the evaluated mcmc data sets
%        '<label>1.mat', '<label>2.mat',... -> ids = [1,2,...].
% - label: Unique identifier for the files. The code e.g. loads '<path>/<label>1.mat'
% - dDense: Thin the results by dDense. If 0 take the whole chains, but
%           shorten the longer one to the size of the shorter one. Ideally
%           the mcmc results were already saved thinned so that loading
%           takes less time. In that case set dDense = 0;
% - BGR_threshold: The smaller the more conservative. Default 0.05. Must be
%                  positive.
% - Geweke_threshold: The smaller the more conservative. Default 1.05. Must be
%                  greater 1.
%
% Output:
% - idxs: Groups as indices
% - z: Matrix with Brooks-Gelman-Rubin test values
% - sfac: Matrix with Geweke test values
%
% Remark:
% This implementation uses a brute force appraoch, meaning you might need a
% lot of RAM. The benefit is that files do not have to be loaded multiple
% times (which from our experience speeds up the process). This method can
% be easily rewritten to a version not storing all runs into memory.
%
% Written by Benjamin Ballnus (2016)


function [idxs,z,sfac] = SimilarityGroupingBrute(path,ids,label,dDense,BGR_threshold,Geweke_threshold)
   
   % Init
   sfac = nan(length(ids),length(ids));
   H = nan(length(ids),length(ids));
   z = nan(length(ids),length(ids));
   J = nan(length(ids),length(ids));
   flag = zeros(length(ids),length(ids));
   
   % Load data brute force
   data = {};
   for i = 1:length(ids)
      try
         file1 = load([path filesep label num2str(ids(i))]);
         data{i} = file1.data;
         clear file1;
      catch e
         data{i} = [];
         disp(e.message);
         flag(i,:) = -1;
      end
   end
   
   for i = 1:length(ids)
      try
         data1 = data{i};
         for j = i+1:length(ids)
            try
               if (flag(i,j) == 0)
                  % Run 2
                  data2 = data{j};
                  
                  % The data sets need to have the same lengths for the GR test, we ensure this
                  % either by thinning onto dDense equvidistant sample points or if
                  % dDense == 0, we shorten the longer singal to the length of the
                  % shorter one
                  if dDense > 0
                     n_x = max(size(data1));
                     n_y = max(size(data2));
                     x_start = 1;
                     y_start = 1;
                     x_idxs = interp1(1:n_x,1:n_x,x_start + (n_x-x_start)/dDense * [0:(dDense-1)],'nearest');
                     y_idxs = interp1(1:n_y,1:n_y,y_start + (n_y-y_start)/dDense * [0:(dDense-1)],'nearest');
                  else
                     n_x = max(size(data1)); n_y = max(size(data2)); n = min(n_x,n_y);
                     x_start = n_x-n+1; y_start = n_y-n+1;
                     x_idxs = x_start:n_x;
                     y_idxs = y_start:n_y;
                  end
                  if (isempty(data1) || isempty(data2))
                     disp(['data empty!  ' num2str(i) '   ' num2str(j)])
                     sfac(i,j) = inf;
                     z(i,j) = inf;
                  else
                     disp(['Success: Loaded data sets ' num2str(i) ' and ' num2str(j) '.'])
                     sfac(i,j) = gelmanRubinBrooksTest(data1(:,x_idxs),data2(:,y_idxs));
                     z(i,j) = max(abs(gewekeTest([data1(:,x_idxs),data2(:,y_idxs)]',0.5,0.5)));
                  end
                  H(i,j) = sfac(i,j) < BGR_threshold;
                  J(i,j) = z(i,j) < Geweke_threshold;
               end
            catch e
               disp(e.message);
               for errIdx = 1:length(e.stack)
                  disp([e.stack(errIdx).file '  ' e.stack(errIdx).name '  ' num2str(e.stack(errIdx).line)]);
               end
               flag(i,j) = -1;
               sfac(i,j) = inf;
               z(i,j) = inf;
            end
            if (flag(i,j) < 0)
               H(i,j) = 0;
               J(i,j) = 0;
            end
         end
      catch e
         disp(e.message);
         for errIdx = 1:length(e.stack)
            disp([e.stack(errIdx).file '  ' e.stack(errIdx).name '  ' num2str(e.stack(errIdx).line)]);
         end
         flag(i,:) = -1;
         H(i,:) = 0;
         J(i,:) = 0;
      end
   end
   
   
   % Merge both test outcomes into one binary value
   L = J .* H;
   
   % Extract indices of data sets which where identified as similar.
   idxs = {};
   for i = 1:length(ids)-1
      for j = i+1:length(ids)
         if L(i,j) > 0
            create_new_group_flag = 1;
            for k = 1:length(idxs)
               if sum(i == idxs{k}) || sum(j == idxs{k})
                  idxs{k} = [idxs{k},i,j];
                  create_new_group_flag = 0;
               end
            end
            if create_new_group_flag == 1
               idxs{end+1} = [i,j];
            end
         end
      end
   end
   
   % Remove non unique indices from groups
   for k = 1:length(idxs)
      idxs{k} = unique(idxs{k});
   end
   
   % Merge groups with an non empty intersection
   for k = 1:length(idxs)
      for l = k+1:length(idxs)
         if ~isempty(intersect(idxs{k},idxs{l}))
            idxs{l} = unique([idxs{k},idxs{l}]);
            idxs{k} = [];
         end
      end
   end
   idxs(cellfun('isempty',idxs)) = [];
   
   % Indices without group are getting their own group
   for i = 1:length(ids)
      has_no_group_flag = 1;
      for k = 1:length(idxs)
         if sum(i==idxs{k}) > 0
            has_no_group_flag = 0;
         end
      end
      if has_no_group_flag
         idxs{end+1} = i;
      end
   end
end




