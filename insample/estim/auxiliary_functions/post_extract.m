function [samples_mix_nowarmup, samples_all] = post_extract(samples_struct, param, warmup, thin)

    numchain = length(samples_struct);

    % Extract samples of specific parameter
    the_dim = size(samples_struct(1).(param));
    samples_all = nan([numchain the_dim]);
    for c=1:numchain
        samples_all(c,:,:) = samples_struct(c).(param);
    end
    
    % Discard warm-up and mix samples from all chains
    samples_mix_nowarmup = reshape(samples_all(:,warmup/thin+1:end,:), [], the_dim(2));

end