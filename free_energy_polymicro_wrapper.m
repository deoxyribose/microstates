%function [F, outparams] = free_energy_polymicro_wrapper(job_id,optparams)
function [F, outparams] = free_energy_polymicro_wrapper(optparams)
    params = load('simulparams.mat');
%     gamma00=log(optparams.gamma00)
%     gamma11=log(optparams.gamma11)
    gamma00=optparams.gamma00;
    gamma11=optparams.gamma11;
     %X0=optparams.X0;
    %[F, outparams] = free_energy_polymicro(params.Y,params.K,params.draw,params.w,params.x,params.M,params.beta,gamma00,1-gamma00,gamma11,1-gamma11,X0);
    [F, outparams] = free_energy_polymicro(params.Y,params.K,params.draw,params.w,params.x,params.M,params.beta,gamma00,1-gamma00,gamma11,1-gamma11,params.x);
end