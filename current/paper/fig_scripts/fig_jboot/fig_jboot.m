jobdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm/155354_sfm/';

topdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm/';


d = dir(jobdir);

jf = [];
j0 = [];
dj = [];
times = [];
for k = 1:length(d)-2
  simdir = [jobdir d(k+2).name '/'];
%   try
    load([simdir 'sims.mat'])
    if length(sims) >= 3
      load([simdir 'eqs.mat'])
      times = [times; str2num(d(k+2).name)];
      j0 = [j0; eqs{2}.jpar];
      jf = [jf; eqs{end}.jpar];
      dj = [dj; (eqs{end}.jpar - eqs{2}.jpar)];
    end
%   catch
%   end
end
psin = eqs{end}.psibar;

plot(psin, dj)
xlim([0.8 1]);
xlabel('psi_N')
ylabel('Jf - J0')

figure
hold on
plot(psin, j0, 'r')
plot(psin, jf, 'g')



jf_pk = mean(max(jf(:,52:end)'));
j0_pk = mean(max(j0(:,52:end)'));

jf_pk / j0_pk



























