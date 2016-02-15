function [F] = free_energy_variational_wrapper(hyperparams)
params = load('passedparams.mat');
[~,~,M,~,free_energy,~,m_winner,~] = variational_microstates_smooth(params.Y,params.K,params.draw,params.alpha,hyperparams(1),params.MULTI,params.G,params.max_nits,params.kfolds,hyperparams(2),params.learn_rate_decay,params.verbose);
%F = mean(free_energy(m_winner,:,2),2);
[~,Z_est] = max(M,[],1);
run_length_encoding_length = description_length(Z_est);
F=run_length_encoding_length;%-mutualinfo(params.Z_true,Z_est)