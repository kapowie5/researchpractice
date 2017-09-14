%% Contour map of conductivity   For Fluoroscence
% close all
% clear all 
%% Import the data
% [~, ~, raw] = xlsread('contour info for map.xlsx','Sheet2','A1:AP42');
[~, ~, raw] = xlsread('contour info for map.xlsx','Sheet3','A1:CE83');
raw = cell2mat(raw);
x = raw(2:end,1);
y = raw(1,2:end);
alphaSpace = raw(2:end,2:end);

%% 'contour info for map.xlsx'
fsize = 20;
figure;contourf(x,y,alphaSpace)
% colormap(hot)
xlabel('x axis (mm)','Fontsize',fsize,'FontWeight','bold')
ylabel('y axis (mm)','Fontsize',fsize,'FontWeight','bold')
c = colorbar('Location','northoutside','Fontsize',fsize,'FontWeight','bold');%,'AxisLocation','in');
c.Label.String = '\alpha (mm^2/s)';
set(gca,'Fontsize',fsize,'FontWeight','bold')

yIndex=20;
xIndex=4; %4, 24, 40


save('Spatial Alpha.mat','x','y','yIndex','xIndex','alphaSpace')

a_matlab_3D_evo