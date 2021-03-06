clear all; clc; close all;

% ===========
% SETTINGS
% ===========
saveit = 1;

% sfdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm/';
% save_fn = 'strike_sfm.mat';
% expect_snowplus = 0;

% sfdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfp_constrain_sp/';
% save_fn = 'strike_sfp.mat';
% expect_snowplus = 1;

sfdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm_large_lambdaq/';
save_fn = 'strike_sfm_large_lambdaq.mat';
expect_snowplus = 0;


% ====================
% FIND STRIKE PT DIFFS
% ====================
load('d3d_obj_mks_struct_6565.mat')
warning('off', 'all')


sfmdir_info = dir(sfdir);
sfmdir_info(1:2) = [];
isim = 0;
strike.err0 = [];
strike.errf = [];

for ishot = 1:length(sfmdir_info)
  
  shot = str2num(sfmdir_info(ishot).name(1:6))
  
  shotdir = [sfdir sfmdir_info(ishot).name '/'];
  shotdir_info = dir(shotdir);
  shotdir_info(1:2) = [];
  
  for itime = 1:length(shotdir_info)
    try
      time = str2num(shotdir_info(itime).name);
      
      jobdir = [shotdir shotdir_info(itime).name '/'];
      load( [jobdir 'sims.mat'])
      load( [jobdir 'eqs.mat'])
      load( [jobdir 'xps.mat'])
      
      
      valid_sim = 1;
      if length(eqs) < 3  || any(isnan(xps{end})) || ...
          any(isnan(eqs{end}.jpar))  || ~isreal(xps{end}) ||...
          ~isreal(eqs{end}.jpar) || ~isreal(eqs{2}.jpar) || ...
          length(sims) < 3
        
        valid_sim = 0;
      end
      
      if valid_sim
        snowf = analyzeSnowflake( eqs{end});
        if snowf.snowPlus ~= expect_snowplus
          valid_sim = 0;
        end
      end
      
      
      if valid_sim
        isim = isim + 1;
        isim
        
        ef = eich_fitter( sims{1}.sir*100, sims{1}.qir, eqs{1}, tok_data_struct);
        
        snow0 = ef.snow;
        dsp0 = [];
        dspf = [];
        
        % snow minus
        if ~snowf.snowPlus
          
          dsp0(1) = norm([ef.rsp(1) - snow0.rSPP(1); ef.zsp(1) - snow0.zSPP(1)]);
          dsp0(2) = norm([ef.rsp(2) - snow0.rSPP(2); ef.zsp(2) - snow0.zSPP(2)]);
          dsp0(3) = norm([ef.rsp(3) - snow0.rSPS(end); ef.zsp(3) - snow0.zSPS(end)]);
          
          dspf(1) = norm([ef.rsp(1) - snowf.rSPP(1); ef.zsp(1) - snowf.zSPP(1)]);
          dspf(2) = norm([ef.rsp(2) - snowf.rSPP(2); ef.zsp(2) - snowf.zSPP(2)]);
          dspf(3) = norm([ef.rsp(3) - snowf.rSPS(end); ef.zsp(3) - snowf.zSPS(end)]);
          
          % power fraction
          strike.pf_true(isim) = sims{1}.qirmaxN(3) / sum(sims{1}.qirmaxN(2:3));
          strike.pf_0(isim) = sims{1}.qmaxN(3) / sum(sims{1}.qmaxN(2:3));
          strike.pf_f(isim) = sims{end}.qmaxN(3) / sum(sims{end}.qmaxN(2:3));
          
          % snowplus
        else
          
          dsp0(1) = norm([ef.rsp(1) - snow0.rSPP(1); ef.zsp(1) - snow0.zSPP(1)]);
          dsp0(2) = norm([ef.rsp(end) - snow0.rSPP(end); ef.zsp(end) - snow0.zSPP(end)]);
          
          dspf(1) = norm([ef.rsp(1) - snowf.rSPP(1); ef.zsp(1) - snowf.zSPP(1)]);
          dspf(2) = norm([ef.rsp(end) - snowf.rSPP(end); ef.zsp(end) - snowf.zSPP(end)]);
        end
        
        strike.err0(isim,:) = dsp0;
        strike.errf(isim,:) = dspf;
        strike.shot(isim) = shot;
        strike.time(isim) = time;
        
        
      end
    catch
    end
  end
end


if saveit
  save(['/u/jwai/d3d_snowflake_2020/current/paper/analysis/strike_pt_errs/' save_fn], 'strike')
end




































