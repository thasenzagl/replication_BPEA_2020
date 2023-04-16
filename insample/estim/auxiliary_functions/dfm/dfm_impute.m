function the_X = dfm_impute(X_raw, r_imp, tol_imp)

    % Impute missing panel data using DFM

    miss_ser = any(isnan(X_raw),1); % Series with missing data
    disp('Pct. missing obs.');
    disp(mean(isnan(X_raw(:))));
    disp('Pct. series with missing obs.');
    disp(mean(miss_ser));

    stand = @(x) (x-nanmean(x))./nanstd(x); % Standardization function
    X_raw = stand(X_raw); % Standardize initial data
    
    the_X = X_raw(:,~miss_ser);
    
    if sum(miss_ser) == 0
        return; % No need for imputation
    end
    
    the_err = Inf;
    iter = 0;
    
    disp('Imputing missing data iteratively...');

    while the_err > tol_imp % Iterate until convergence
    
        the_F = pca(the_X,r_imp); % Factor estimate given current data

        % Impute missing data using factor
        the_X_new = X_raw;
        for j=find(miss_ser==1) % For every series with missing obs.
            the_miss = isnan(X_raw(:,j)); % Missing obs. for series
            the_beta = the_F(~the_miss,:)\X_raw(~the_miss,j); % Regression on non-missing data
            the_X_new(the_miss,j) = the_F(the_miss,:)*the_beta; % Impute missing data for series
        end

        if iter>0
            the_err = max(abs(the_X_new(:)-the_X(:))); % Max change in imputation
            fprintf('%4d: %10.8e\n', iter, the_err);
        end
        the_X = the_X_new;
        iter = iter+1;
    
    end
    
    the_X = stand(the_X); % Standardize again
    
    disp('Done.')

end