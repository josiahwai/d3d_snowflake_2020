load('/u/jwai/d3d_snowflake_2020/current/inputs/edge_current/155355/jdata_155355.mat')
close 

k = find(times==3900);
figure
hold on
plot(psi_n,jtot(k,:),'g','linewidth',2)
% plot(psi_n,j_efit(k,:))

load('cake_eq_155355_3900.mat')
plot(eq.psibar,eq.jpar,'b','linewidth',2)

load('/u/jwai/d3d_snowflake_2020/current/ml/train/jobs/155355/3900/452/eq.mat')
plot(eq.psibar,eq.jpar,'r','linewidth',2)





xlim([0.8304 1.0922])
ylim([-1.0482e+05 1.2952e+06])


labels = {'CAKE','constrained_{original x-pts}', 'constrained_{best x-pts}'};
co  = {[0 1 0], [0 0 1], [1 0 0]};


mylegend(labels, [],[], co);
ylabel('j')
xlabel('Psi_N')
title('155358: 3900ms')