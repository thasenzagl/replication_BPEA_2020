function [A, Sigma] = sim_post_niw_diffuse(Ahat, Sigmahat_inv_chol, XpX, T)

    % Simulate from maximally diffuse normal-inverse-Wishart posterior
    % See Uhlig (JME 2005), p. 410
    
    Sigma_inv = wishrnd([],T,Sigmahat_inv_chol/sqrt(T));
    Sigma = inv(Sigma_inv);
    A_vec = mvnrnd(Ahat(:), kron(Sigma,inv(XpX)));
    A = reshape(A_vec, size(Ahat));

end