load('/u/jwai/d3d_snowflake_2020/current/dev/3594/eqs.mat')
load('/u/jwai/d3d_snowflake_2020/current/dev/3594/sims.mat')
load('d3d_obj_mks_struct_6565.mat')

ef = eich_fitter_dev( sims{1}.sir*100, sims{1}.qir, eqs{1}, tok_data_struct, 1);

snow0 = ef.snow;
% snowf = analyzeSnowflake( eqs{end});

% snow minus
dsp0(1) = norm([ef.rsp(1) - snow0.rSPP(1); ef.zsp(1) - snow0.zSPP(1)]);
dsp0(2) = norm([ef.rsp(2) - snow0.rSPP(2); ef.zsp(2) - snow0.zSPP(2)]);
dsp0(3) = norm([ef.rsp(3) - snow0.rSPS(end); ef.zsp(3) - snow0.zSPS(end)]);

dspf(1) = norm([ef.rsp(1) - snowf.rSPP(1); ef.zsp(1) - snowf.zSPP(1)]);
dspf(2) = norm([ef.rsp(2) - snowf.rSPP(2); ef.zsp(2) - snowf.zSPP(2)]);
dspf(3) = norm([ef.rsp(3) - snowf.rSPS(end); ef.zsp(3) - snowf.zSPS(end)]);

plot_eq(eqs{end})
scatter(snowf.rSPP, snowf.zSPP, 'b', 'filled')
scatter(snowf.rSPS, snowf.zSPS, 'b', 'filled') 
scatter(ef.rsp, ef.zsp, 'r', 'filled')         % target


dsp0
dspf




