params.n_iterations = 190;
params.n_init_samples = 10;
params.crit_name = 'cEI';
params.surr_name = 'sStudentTProcessNIG';
params.noise = 1e-6;
params.kernel_name = 'kMaternARD5';
params.kernel_hp_mean = [1];
params.kernel_hp_std = [10];
params.verbose_level = 1;
params.log_filename = 'matbopt.log';

disp('Continuous optimization');
fun = 'wrapper'; n = 2;

tic;
[rbfsigmaandalpha, meanGMD] = bayesoptcont(fun,n,params,[10^-6,10^-6],[10^4, 10^4])
toc;
disp('Press INTRO');
pause;
