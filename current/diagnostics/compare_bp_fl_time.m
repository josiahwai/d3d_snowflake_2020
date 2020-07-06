ccc
root = '/u/jwai/d3d_snowflake_2020/current/';

% ========
% SETTINGS
% ========
shot = 155352;

% =======================
% LOAD EQ AND FIT BPROBES
% =======================
load('all_shots_times.mat')
idiv = find(all_shots_times(:,1) == shot);
times = sort(unique(all_shots_times(idiv,2)));

efit_dir = ['/p/omfit/users/jwai/projects/snowflake_efits/EFITtime/OUTPUTS/' ...
  num2str(shot) '/gEQDSK/'];

load(['bp_fl_data' num2str(shot)])
struct_to_ws(bp_fl_data.template);  

chisq0 = [];
chisqf = [];
expmpi_all = [];
expmpi_all_smooth = [];
times(times==4258)=[];

for iTime = 1:length(times)
  try
  time_ms = times(iTime);
  
  shotdir_sfm = [root '/sfmodel/jobs/sfm/' num2str(shot) '_sfm/'];
  shotdir_sfp = [root '/sfmodel/jobs/sfp_constrain_sp/' num2str(shot) '_sfp_sp/'];  
  fn_sfm = [shotdir_sfm num2str(time_ms) '/eqs.mat'];
  fn_sfp = [shotdir_sfp num2str(time_ms) '/eqs.mat'];
  
  if isfile(fn_sfm)
    load(fn_sfm);
  else 
    load(fn_sfp);
  end  

  % eq0 = read_eq(shot, time_ms/1000, efit_dir);
  eq0 = eqs{2};
  eqf = eqs{end};
  
  
  bp_fl = getfield(bp_fl_data, ['t' num2str(time_ms)]);
  struct_to_ws(bp_fl);
  
  expmpi_all(iTime,:) = expmpi;  
  
  fwtmp2 = fwtmp2 / max(fwtmp2);
%   fwtmp2 = FWTMP2;
    
  Bp0 = calc_mag_probe(eq0, XMP2, YMP2, AMP2);
  Bpf = calc_mag_probe(eqf, XMP2, YMP2, AMP2);
  
  chisq0(iTime,:) = (fwtmp2 .* (expmpi - Bp0)).^2;  % fitting w
  chisqf(iTime,:) = (fwtmp2 .* (expmpi - Bpf)).^2;
  catch
  end
end


% ========
% PLOT IT
% ========
figure
hold on
plot(times, sum(chisq0,2), 'r')
plot(times, sum(chisqf,2), 'b')
set(gcf,'position',[740 387 654 313])


for k = 1:76
  expmpi_all_smooth(:,k) = smooth(expmpi_all(:,k), 10);
end
idiv = 52;% :52;
ndiv = length(idiv);

figure
for k = 1:ndiv
  subplot(ndiv,1,k)
  hold on
  y0 = expmpi_all_smooth(5,idiv(k));
  
  plot(times, expmpi_all(:,idiv(k)) - y0, 'r')
  plot(times, expmpi_all_smooth(:,idiv(k)) - y0, 'b')
  
  mean(abs(expmpi_all(10:end,idiv(k)) - expmpi_all_smooth(10:end,idiv(k))))
  title(['Probe ' num2str(idiv(k))])
  ylabel('B field [T]')
  xlim([3000 5000])
  
  if k==1, legend('"Raw"', 'smoothed'); end
end
xlabel('time [ms]')



% figure(1)
% chisq_smooth = [];  
% for iTime = 1:length(times)
%   chisq_smooth(iTime,:) = (fwtmp2 .* (expmpi_all(iTime,:) - expmpi_all_smooth(iTime,:))).^2; 
% end
% plot(times, sum(chisq_smooth, 2), 'k')
































