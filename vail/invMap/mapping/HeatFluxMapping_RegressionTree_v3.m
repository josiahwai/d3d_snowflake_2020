
% REGRESSION TREE
% Predictors: s_SP1, s_SP2, s_SP3, q0_SP1, q0_SP2, q0_SP3

% clear all

plotit = 1;
saveit = 0;

%...............................................................................
% Gather the heat flux simulation data

HeatFluxMapping_GatherData

neq = size(spData_ALL,1);

%...............................................................................
% Partition data for training and testing

pTrain = 0.99;

tf = false(neq,1);
tf(1:round(pTrain*neq)) = true;
tf = tf(randperm(neq));

spData_Train = spData_ALL( tf,:);
spData_Test  = spData_ALL(~tf,:);

xpData_Train = xpData_ALL( tf,:);
xpData_Test  = xpData_ALL(~tf,:);

% Determine which of the training and testing data are HFS/LFS

idxHFS_Train = find(xpData_Train(:,1) > xpData_Train(:,2));
idxLFS_Train = setdiff(1:size(xpData_Train,1), idxHFS_Train);

idxHFS_Test  = find(xpData_Test(:,1) > xpData_Test(:,2));
idxLFS_Test  = setdiff(1:size(xpData_Test,1), idxHFS_Test);

%...............................................................................
% Fit model to the training data

tbl_rxP = table(spData_Train(:,7),  spData_Train(:,8),  spData_Train(:,9),  ...
                spData_Train(:,1),  spData_Train(:,2),  spData_Train(:,3),  ...
                xpData_Train(:,1), 'VariableNames', {'s_SP1', 's_SP2',      ...
                's_SP3', 'q0_SP1', 'q0_SP2', 'q0_SP3', 'rxP'});

tbl_rxS = table(spData_Train(:,7),  spData_Train(:,8),  spData_Train(:,9),  ...
                spData_Train(:,1),  spData_Train(:,2),  spData_Train(:,3),  ...
                xpData_Train(:,2), 'VariableNames', {'s_SP1', 's_SP2',      ...
                's_SP3', 'q0_SP1', 'q0_SP2', 'q0_SP3', 'rxS'});           
       
tbl_zxP = table(spData_Train(:,7),  spData_Train(:,8),  spData_Train(:,9),  ...
                spData_Train(:,1),  spData_Train(:,2),  spData_Train(:,3),  ...
                xpData_Train(:,3), 'VariableNames', {'s_SP1', 's_SP2',      ...
                's_SP3', 'q0_SP1', 'q0_SP2', 'q0_SP3', 'zxP'});
                 
tbl_zxS = table(spData_Train(:,7),  spData_Train(:,8),  spData_Train(:,9),  ...
                spData_Train(:,1),  spData_Train(:,2),  spData_Train(:,3),  ...
                xpData_Train(:,4), 'VariableNames', {'s_SP1', 's_SP2',      ...
                's_SP3', 'q0_SP1', 'q0_SP2', 'q0_SP3', 'zxS'});      

rtree_rxP = fitrtree(tbl_rxP, 'rxP');
rtree_rxS = fitrtree(tbl_rxS, 'rxS');
rtree_zxP = fitrtree(tbl_zxP, 'zxP');
rtree_zxS = fitrtree(tbl_zxS, 'zxS');

rxP_PRED_TRAIN = predict(rtree_rxP, spData_Train(:,[7:9 1:3]));
rxS_PRED_TRAIN = predict(rtree_rxS, spData_Train(:,[7:9 1:3]));
zxP_PRED_TRAIN = predict(rtree_zxP, spData_Train(:,[7:9 1:3]));
zxS_PRED_TRAIN = predict(rtree_zxS, spData_Train(:,[7:9 1:3]));

rxP_TRUE_TRAIN = xpData_Train(:,1);
rxS_TRUE_TRAIN = xpData_Train(:,2);
zxP_TRUE_TRAIN = xpData_Train(:,3);
zxS_TRUE_TRAIN = xpData_Train(:,4);

error_mean_rxP_TRAIN = mean(abs((rxP_TRUE_TRAIN - rxP_PRED_TRAIN)));
error_mean_rxS_TRAIN = mean(abs((rxS_TRUE_TRAIN - rxS_PRED_TRAIN)));
error_mean_zxP_TRAIN = mean(abs((zxP_TRUE_TRAIN - zxP_PRED_TRAIN)));
error_mean_zxS_TRAIN = mean(abs((zxS_TRUE_TRAIN - zxS_PRED_TRAIN)));

if saveit
    save('error_mean_rxP_TRAIN.mat', 'error_mean_rxP_TRAIN')
    save('error_mean_rxS_TRAIN.mat', 'error_mean_rxS_TRAIN')
    save('error_mean_zxP_TRAIN.mat', 'error_mean_zxP_TRAIN')
    save('error_mean_zxS_TRAIN.mat', 'error_mean_zxS_TRAIN')
end

error_stdv_rxP_TRAIN = std(abs((rxP_TRUE_TRAIN - rxP_PRED_TRAIN)));
error_stdv_rxS_TRAIN = std(abs((rxS_TRUE_TRAIN - rxS_PRED_TRAIN)));
error_stdv_zxP_TRAIN = std(abs((zxP_TRUE_TRAIN - zxP_PRED_TRAIN)));
error_stdv_zxS_TRAIN = std(abs((zxS_TRUE_TRAIN - zxS_PRED_TRAIN)));

if saveit
    save('error_stdv_rxP_TRAIN.mat', 'error_stdv_rxP_TRAIN')
    save('error_stdv_rxS_TRAIN.mat', 'error_stdv_rxS_TRAIN')
    save('error_stdv_zxP_TRAIN.mat', 'error_stdv_zxP_TRAIN')
    save('error_stdv_zxS_TRAIN.mat', 'error_stdv_zxS_TRAIN')
end

%...............................................................................
% Test model using the testing data

rxP_PRED_TEST = predict(rtree_rxP, spData_Test(:,[7:9 1:3]));
rxS_PRED_TEST = predict(rtree_rxS, spData_Test(:,[7:9 1:3]));
zxP_PRED_TEST = predict(rtree_zxP, spData_Test(:,[7:9 1:3]));
zxS_PRED_TEST = predict(rtree_zxS, spData_Test(:,[7:9 1:3]));

rxP_TRUE_TEST = xpData_Test(:,1);
rxS_TRUE_TEST = xpData_Test(:,2);
zxP_TRUE_TEST = xpData_Test(:,3);
zxS_TRUE_TEST = xpData_Test(:,4);

error_mean_rxP_TEST = mean(abs((rxP_TRUE_TEST - rxP_PRED_TEST)));
error_mean_rxS_TEST = mean(abs((rxS_TRUE_TEST - rxS_PRED_TEST)));
error_mean_zxP_TEST = mean(abs((zxP_TRUE_TEST - zxP_PRED_TEST)));
error_mean_zxS_TEST = mean(abs((zxS_TRUE_TEST - zxS_PRED_TEST)));

if saveit
    save('error_mean_rxP_TEST.mat', 'error_mean_rxP_TEST')
    save('error_mean_rxS_TEST.mat', 'error_mean_rxS_TEST')
    save('error_mean_zxP_TEST.mat', 'error_mean_zxP_TEST')
    save('error_mean_zxS_TEST.mat', 'error_mean_zxS_TEST')
end

error_stdv_rxP_TEST = std(abs((rxP_TRUE_TEST - rxP_PRED_TEST)));
error_stdv_rxS_TEST = std(abs((rxS_TRUE_TEST - rxS_PRED_TEST)));
error_stdv_zxP_TEST = std(abs((zxP_TRUE_TEST - zxP_PRED_TEST)));
error_stdv_zxS_TEST = std(abs((zxS_TRUE_TEST - zxS_PRED_TEST)));

if saveit
    save('error_stdv_rxP_TEST.mat', 'error_stdv_rxP_TEST')
    save('error_stdv_rxS_TEST.mat', 'error_stdv_rxS_TEST')
    save('error_stdv_zxP_TEST.mat', 'error_stdv_zxP_TEST')
    save('error_stdv_zxS_TEST.mat', 'error_stdv_zxS_TEST')
end

%...............................................................................
% Plot the training data

if plotit

figure()

set(gcf, 'Position', [200 200 800 800])

subplot(2,2,1)
    plot(linspace(-2,2), linspace(-2,2),'--k')
subplot(2,2,2)
    plot(linspace(-2,2), linspace(-2,2), '--k')
subplot(2,2,3)
    plot(linspace(-2,2), linspace(-2,2), '--k')
subplot(2,2,4)
    plot(linspace(-2,2), linspace(-2,2), '--k')

subplot(2,2,1)
    hold on
    plot(rxP_PRED_TRAIN, rxP_TRUE_TRAIN,   'ok')
    hold on
    
    minx = min(rxP_PRED_TRAIN) - 0.02; maxx = max(rxP_PRED_TRAIN) + 0.02;
    miny = min(rxP_TRUE_TRAIN) - 0.02; maxy = max(rxP_TRUE_TRAIN) + 0.02;

    xlim([minx maxx])
    ylim([miny maxy])

    title('rxP [m] TRAIN')
    xlabel('PRED', 'Interpreter', 'None')
    ylabel('TRUE', 'Interpreter', 'None')
    grid on

text(1.01*minx, 0.99*maxy, ['Error = ' num2str(100*error_mean_rxP_TRAIN, ...
    '%5.3f') '\pm' num2str(100*error_stdv_rxP_TRAIN, '%5.3f') ' cm'],    ...
    'FontSize', 10)

subplot(2,2,2)
    hold on
    plot(rxS_PRED_TRAIN, rxS_TRUE_TRAIN,   'ok')
    hold on
    
    minx = min(rxS_PRED_TRAIN) - 0.02; maxx = max(rxS_PRED_TRAIN) + 0.02;
    miny = min(rxS_TRUE_TRAIN) - 0.02; maxy = max(rxS_TRUE_TRAIN) + 0.02;

    xlim([minx maxx])
    ylim([miny maxy])

    title('rxS [m] TRAIN')
    xlabel('PRED', 'Interpreter', 'None')
    ylabel('TRUE', 'Interpreter', 'None')
    grid on

text(1.01*minx, 0.99*maxy, ['Error = ' num2str(100*error_mean_rxS_TRAIN, ...
    '%5.3f') '\pm' num2str(100*error_stdv_rxS_TRAIN, '%5.3f') ' cm'],    ...
    'FontSize', 10)
  
subplot(2,2,3)
    hold on
    plot(zxP_PRED_TRAIN, zxP_TRUE_TRAIN,   'ok')
    hold on
    
    minx = min(zxP_PRED_TRAIN) - 0.02; maxx = max(zxP_PRED_TRAIN) + 0.02;
    miny = min(zxP_TRUE_TRAIN) - 0.02; maxy = max(zxP_TRUE_TRAIN) + 0.02;

    xlim([minx maxx])
    ylim([miny maxy])

    title('zxP [m] TRAIN')
    xlabel('PRED', 'Interpreter', 'None')
    ylabel('TRUE', 'Interpreter', 'None')
    grid on

text(0.99*minx, 1.01*maxy, ['Error = ' num2str(100*error_mean_zxP_TRAIN, ...
    '%5.3f') '\pm' num2str(100*error_stdv_zxP_TRAIN, '%5.3f') ' cm'],    ...
    'FontSize', 10)

subplot(2,2,4)
    hold on
    plot(zxS_PRED_TRAIN, zxS_TRUE_TRAIN,   'ok')
    hold on
    
    minx = min(zxS_PRED_TRAIN) - 0.02; maxx = max(zxS_PRED_TRAIN) + 0.02;
    miny = min(zxS_TRUE_TRAIN) - 0.02; maxy = max(zxS_TRUE_TRAIN) + 0.02;

    xlim([minx maxx])
    ylim([miny maxy])

    xlim([min(zxS_PRED_TRAIN)-0.02 max(zxS_PRED_TRAIN)+0.02])
    ylim([min(zxS_TRUE_TRAIN)-0.02 max(zxS_TRUE_TRAIN)+0.02])

    title('zxS [m] TRAIN')
    xlabel('PRED', 'Interpreter', 'None')
    ylabel('TRUE', 'Interpreter', 'None')
    grid on
    
text(0.99*minx, 1.01*maxy, ['Error = ' num2str(100*error_mean_zxS_TRAIN, ...
    '%5.3f') '\pm' num2str(100*error_stdv_zxS_TRAIN, '%5.3f') ' cm'],    ...
    'FontSize', 10)    
    
if saveit
    set(gcf, 'PaperPositionMode', 'auto')
    print -depsc2 HeatFluxMapping_RegressionTree_v1_Training.eps
end

end

%...............................................................................
% Plot the testing data

if plotit

figure()

set(gcf, 'Position', [200 200 800 800])

subplot(2,2,1)
    plot(linspace(-2,2), linspace(-2,2),'--k')
subplot(2,2,2)
    plot(linspace(-2,2), linspace(-2,2), '--k')
subplot(2,2,3)
    plot(linspace(-2,2), linspace(-2,2), '--k')
subplot(2,2,4)
    plot(linspace(-2,2), linspace(-2,2), '--k')

subplot(2,2,1)
    hold on
    plot(rxP_PRED_TEST, rxP_TRUE_TEST,   'ok')
    hold on
    
    minx = min(rxP_PRED_TEST) - 0.02; maxx = max(rxP_PRED_TEST) + 0.02;
    miny = min(rxP_TRUE_TEST) - 0.02; maxy = max(rxP_TRUE_TEST) + 0.02;

    xlim([minx maxx])
    ylim([miny maxy])

    title('rxP [m] TEST')
    xlabel('PRED', 'Interpreter', 'None')
    ylabel('TRUE', 'Interpreter', 'None')
    grid on

text(1.01*minx, 0.99*maxy, ['Error = ' num2str(100*error_mean_rxP_TEST, ...
    '%5.3f') '\pm' num2str(100*error_stdv_rxP_TEST, '%5.3f') ' cm'],    ...
    'FontSize', 10)

subplot(2,2,2)
    hold on
    plot(rxS_PRED_TEST, rxS_TRUE_TEST,   'ok')
    hold on
    
    minx = min(rxS_PRED_TEST) - 0.02; maxx = max(rxS_PRED_TEST) + 0.02;
    miny = min(rxS_TRUE_TEST) - 0.02; maxy = max(rxS_TRUE_TEST) + 0.02;

    xlim([minx maxx])
    ylim([miny maxy])

    title('rxS [m] TEST')
    xlabel('PRED', 'Interpreter', 'None')
    ylabel('TRUE', 'Interpreter', 'None')
    grid on

text(1.01*minx, 0.99*maxy, ['Error = ' num2str(100*error_mean_rxS_TEST, ...
    '%5.3f') '\pm' num2str(100*error_stdv_rxS_TEST, '%5.3f') ' cm'],    ...
    'FontSize', 10)
  
subplot(2,2,3)
    hold on
    plot(zxP_PRED_TEST, zxP_TRUE_TEST,   'ok')
    hold on
    
    minx = min(zxP_PRED_TEST) - 0.02; maxx = max(zxP_PRED_TEST) + 0.02;
    miny = min(zxP_TRUE_TEST) - 0.02; maxy = max(zxP_TRUE_TEST) + 0.02;

    xlim([minx maxx])
    ylim([miny maxy])

    title('zxP [m] TEST')
    xlabel('PRED', 'Interpreter', 'None')
    ylabel('TRUE', 'Interpreter', 'None')
    grid on

text(0.99*minx, 1.01*maxy, ['Error = ' num2str(100*error_mean_zxP_TEST, ...
    '%5.3f') '\pm' num2str(100*error_stdv_zxP_TEST, '%5.3f') ' cm'],    ...
    'FontSize', 10)

subplot(2,2,4)
    hold on
    plot(zxS_PRED_TEST, zxS_TRUE_TEST,   'ok')
    hold on
    
    minx = min(zxS_PRED_TEST) - 0.02; maxx = max(zxS_PRED_TEST) + 0.02;
    miny = min(zxS_TRUE_TEST) - 0.02; maxy = max(zxS_TRUE_TEST) + 0.02;

    xlim([minx maxx])
    ylim([miny maxy])

    xlim([min(zxS_PRED_TEST)-0.02 max(zxS_PRED_TEST)+0.02])
    ylim([min(zxS_TRUE_TEST)-0.02 max(zxS_TRUE_TEST)+0.02])

    title('zxS [m] TEST')
    xlabel('PRED', 'Interpreter', 'None')
    ylabel('TRUE', 'Interpreter', 'None')
    grid on
    
text(0.99*minx, 1.01*maxy, ['Error = ' num2str(100*error_mean_zxS_TEST, ...
    '%5.3f') '\pm' num2str(100*error_stdv_zxS_TEST, '%5.3f') ' cm'],    ...
    'FontSize', 10)    
    
if saveit
    set(gcf, 'PaperPositionMode', 'auto')
    print -depsc2 HeatFluxMapping_RegressionTree_v1_Testing.eps
end

end