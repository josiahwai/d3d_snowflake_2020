% Finds the thermal diffusivity chi and scrape off layer power
% decay lengths lambdaq for the inboard and outboard divertor. Takes an
% average over the whole shot. Params are found from the eich profile

% lambdaq in [cm], chi in Wb^2/s

ccc
shotlist = 155344:2:155352;
saveit = 1;
plotit = 0;
  
root = '/u/jwai/d3d_snowflake_2020/current/';
load('d3d_obj_mks_struct_6565.mat')
warning('off', 'all')

for iShot = 1:length(shotlist)
  try
  shot = shotlist(iShot)
  
  coeff = [];
  valid_times = [];
  
  % heat flux data was stored in two places, load and merge
  load(['qperp_std_' num2str(shot) '.mat']) 
  t0 = t; 
  qperp0 = qperp;
  
  load(['qperp_' num2str(shot) '.mat'])
  qperp = [qperp0; qperp];
  t = double([t0 t]);
  
  % loop over times
  for iTime = 1:length(t)
    try
      time_ms = t(iTime)
      q = qperp(iTime,:)';
      
      % load eq and calculate flux expansion
      eqdir = [root 'inputs/eqs/cake/' num2str(shot)];
      eq = read_eq(shot, time_ms/1000, eqdir);
      
      iX = find(s > 120 & s < 145);
      snowPlus = 1;
      if max(q(iX)) > 3 * median(q), snowPlus = 0; end
      
      ef = eich_fitter_geo(s, q, eq, tok_data_struct, plotit);
      
      % sanity checks on the fit
      if (ef.gofi.rsquare > 0.7) && (ef.gofo.rsquare > 0.7) && ...
          (ef.fito.lambdaQ > .03) && (ef.fito.q0 < 1000)
        
        coeff = [coeff; ef.fiti.lambdaQ ef.fito.lambdaQ ef.chi_i ef.chi_o snowPlus];
        valid_times = [valid_times; time_ms];
        
      else
        warning('on', 'all')
        warning('Fit quality low. Skipping time slice...')
        warning('off', 'all')
      end
      
    catch
      warning('on', 'all')
      warning(['Could not fit for ' num2str(shot) ': ' num2str(time_ms) 'ms'])
      warning('off', 'all')
    end
  end
  
  avg_coeff = mean(coeff);
  
  [lambdaq_i, lambdaq_o, chi_i, chi_o] = unpack(avg_coeff);
  
  sol_params = struct('lambdaq_i', lambdaq_i, 'lambdaq_o', lambdaq_o, ...
    'chi_i', chi_i, 'chi_o', chi_o, 'coeff', coeff, 'times', valid_times);
  
  if saveit
    savedir = '/u/jwai/d3d_snowflake_2020/current/eich/geo/sol_params/';
    fn = [savedir 'sol_params2_' num2str(shot) '.mat'];
    save(fn, 'sol_params');
  end
  
  
  if plotit    
    sf = 1:length(coeff);
    sfp = find(coeff(:,5));
    sfm = setdiff(sf,sfp);
    
    figure
    subplot(1,2,1)
    bar(coeff(sfp,2))
    title(['SFP: LambdaQ: ' num2str(shot)])
    ylim([0 0.4])
    
    subplot(1,2,2)
    bar(coeff(sfm,2),'r')
    title(['SFM: LambdaQ: ' num2str(shot)])
    ylim([0 0.4])
    
    set(gcf, 'position', [671 640 577 238])
  end
  catch
  end
end
  
  
  
  
  
  
  
  
  
  
  
