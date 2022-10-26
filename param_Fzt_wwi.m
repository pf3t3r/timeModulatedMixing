clc; clear; close all;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 22]);
set(0,'defaultAxesFontSize',18);

% Add subfolders to path
addpath(genpath(pwd));

load E_wwi;
% load('E_hil.mat','zq','N2_int_3D_overTime'); % didn't save zq in WWI so including from here
load('Mixing_parameterization_fields.mat');
% load N2v2.mat;
load N_is.mat;
load IC3.mat;
mooring = "IC3";

%% Compute the vertical structure of eps_wwi

% There is no topography component to the vertical structure this time. It
% is purely a function of the buoyancy frequency.

% depthsIC3 = zq;
depthsIC3 = zmid_IC3(:,1);
N2 = N2_IC3;


%%
%% Integrate N^2
dz_IC3 = cat(1,10*ones(5,1),25*ones(6,1),50*ones(30,1));
dz_M1 = cat(1,10*ones(5,1),25*ones(6,1),50*ones(31,1));

for i=1:length(time)
    N2_int_overTime_IC3(i) = sum(N2_IC3(:,i).*dz_IC3,'omitnan');
    N2_int_overTime_M1(i) = sum(N2_M1(:,i).*dz_M1,'omitnan');
end

N2_int_3D_overTime_IC3 = repmat(N2_int_overTime_IC3,1,1,length(zmid_IC3(:,1)));
N2_int_3D_overTime_IC3 = squeeze(N2_int_3D_overTime_IC3(1,:,:))';

N2_int_3D_overTime_M1 = repmat(N2_int_overTime_M1,1,1,length(zmid_M1(:,1)));
N2_int_3D_overTime_M1 = squeeze(N2_int_3D_overTime_M1(1,:,:))';

%%

if isfile('verticalStructureWwi.mat')
    disp('Time-dependent vertical structure already calculated');
    load verticalStructureWwi.mat;
else
    disp('Calculating vertical structure...');
    verticalStructureIC3Wwi = NaN(length(depthsIC3),length(time));
    for t = 1:length(time)
        for k=1:length(depthsIC3)
            verticalStructureIC3Wwi(k,t) = (N2(k,t)./N2_int_3D_overTime_IC3(k,t));
        end
    if mod(t,1000)==0
        disp(t);
    end
    end
    disp('Vertical structure calculated');
    save Matfiles/verticalStructureWwi.mat verticalStructureIC3Wwi;
end

%% Set up 2D-grid for contourf

time2 = time';
timeNum = datenum(time2);
[timeNum2,depths2] = meshgrid(timeNum,depthsIC3);

%% Histogram: F(z,t)

% Use this as a justification for binning scheme.

ax1 = figure;
histogram(verticalStructureIC3Wwi);
title('F_{wwi}(z,t): before binning');

% savefig('figures/main/param/' + mooring + '_hist_verticalStructureWwi');
exportgraphics(ax1,'figures/main/param/' + mooring + '_hist_verticalStructureWwi.png');
 
%% Bin the vertical structure of WWI

timeEnd = length(time);
edges = [0:1e-4:0.9e-3 1e-3:2e-3:7e-3];
NCbar = length(edges) - 1;
binnedEpsilon = discretize(verticalStructureIC3Wwi,edges);

ax2 = figure;
contourf(timeNum2(:,1:timeEnd)',depths2(:,1:timeEnd)',binnedEpsilon(:,1:timeEnd)','LineColor','none');
datetick('x','yyyy mmm','keeplimits');
cmap = flipud(cbrewer2('div','Spectral',NCbar));
colormap(cmap);
% 'Ticks',[0,1/14,2/14,3/14,4/14,5/14,6/14,7/14,8/14,9/14,10/14,11/14,12/14,13/14,1],
c = colorbar('YTick',1:14,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14))});
c.Label.String = 'F(z,t) [m^{-1}]';
xlabel('time');
ylabel('depth [m]');
title('F_{wwi}(z,t): binned values');
% savefig('figures/main/param/' + mooring + '_verticalStructureWWIbinned');
exportgraphics(ax2,'figures/main/param/' + mooring + '_verticalStructureWWIbinned.png');

%% Vertical Structure of WWI for IC3: not binned

ax2a = figure;
contourf(timeNum2,depths2,log10(verticalStructureIC3Wwi),'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
ylabel('Depth [m]');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
ylim([-inf -100]);
c.Label.String = 'log_{10}(F) [log m^{-1}]';
% title('IC3: Vertical Structure due to WWI');
exportgraphics(ax2a,'figures/main/param/' + mooring + '_verticalStructureWwi.png');

%% Zoomed-in contourf

depthToPlotFrom = 33;
t1 = 9000;
t2 = 10000;
t3 = 30000;
t4 = 31000;

ax2a = figure;
sgtitle('F(z,t): zoomed-in');
subplot(1,2,1)
contourf(timeNum2(depthToPlotFrom:end,t1:t2)',depths2(depthToPlotFrom:end,t1:t2)',binnedEpsilon(depthToPlotFrom:end,t1:t2)','LineColor','none');
datetick('x','yyyy mmm','keeplimits');
cmap = flipud(cbrewer2('div','Spectral',NCbar));
colormap(cmap);
c = colorbar('YTick',1:14,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14))});
xlabel('time');
ylabel('depth [m]');

subplot(1,2,2)
contourf(timeNum2(depthToPlotFrom:end,t3:t4)',depths2(depthToPlotFrom:end,t3:t4)',binnedEpsilon(depthToPlotFrom:end,t3:t4)','LineColor','none');
datetick('x','yyyy mmm','keeplimits');
cmap = flipud(cbrewer2('div','Spectral',NCbar));
colormap(cmap);
c = colorbar('YTick',1:14,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14))});
xlabel('time');
ylabel('depth [m]');

clear depthToPlotFrom t1 t2 t3 t4;

% savefig('figures/main/param/' + mooring + '_verticalStructureFunctionZoom');
exportgraphics(ax2a,'figures/main/param/' + mooring + '_verticalStructureFunctionZoom.png');

%% Histogram after binning
ax3 = figure;
histogram(binnedEpsilon);
title('F_{wwi}(z,t): after binning');

% savefig('figures/main/param/' + mooring + '_hist_verticalStructure_binned');
exportgraphics(ax3,'figures/main/param/' + mooring + '_hist_verticalStructure_binned.png');
%% Kurtosis of F(z,t)

% kurtF = kurtosis(verticalStructureIC3Hil,0);
% 
% figure
% plot(time,kurtF);
% xlabel('time');
% ylabel('K_F (Kurtosis-value of F(z,t)');
% 
% 
% lessProneToOutliers = 0;
% moreProneToOutliers = 0;
% 
% for i=1:length(time)
%     if kurtF(i) < 3
%         lessProneToOutliers = lessProneToOutliers + 1;
%     end
%     if kurtF(i) > 3
%         moreProneToOutliers = moreProneToOutliers + 1;
%     end
% end
% 
% propProne = 100*moreProneToOutliers./length(time);

%% Plot of Vertical Structure components T1 and T2
% 
% ax4 = figure;
% 
% subplot(2,1,1)
% contourf(timeNum2',depths2',T1','LineColor','none');
% datetick('x','yyyy mmm','keeplimits');
% c = colorbar;
% c.Label.String = 'F(z,t) [m^{-1}]';
% colormap(flipud(cbrewer2('Spectral',NCbar)));
% title('Vertical Structure F(z,t): T1');
% 
% subplot(2,1,2)
% contourf(timeNum2',depths2',T2','LineColor','none');
% datetick('x','yyyy mmm','keeplimits');
% c = colorbar;
% c.Label.String = 'F(z,t) [m^{-1}]';
% colormap(flipud(cbrewer2('Spectral',NCbar)));
% title('Vertical Structure F(z,t): T2');
% 
% % subplot(3,1,3)
% % contourf(timeNum2',depths2',sqrt(N2)','LineColor','none');
% % datetick('x','yyyy mmm','keeplimits');
% % c = colorbar;
% % c.Label.String = 'N(z,t) [s{-1}]';
% % colormap(flipud(cbrewer2('Spectral',NCbar)));
% % title('Buoyancy Frequency');
% 
% savefig('figures/main/param/' + mooring + '_verticalStructureT1_T2');
% exportgraphics(ax4,'figures/main/param/' + mooring + '_verticalStructureT1_T2.png');
% 
% % [rho_Nv,p_Nv] = corr(N2(1,:),verticalStructureIC3Hil(1,:));
% % 
% % % disp(corrN2vert);

%% Save parameters for next file

save Matfiles/Fzt_wwi.mat verticalStructureIC3Wwi;