%%----------------HEADER---------------------------%%
%Author:           Boris Segret
%Version & Date:
%                  V1.2 06-05-2016 (dd/mm/yyyy)
%                  - forked from plotEME.m v1.1
%CL=0
%Version & Date: (initial)
%
% Reads & plots a scenario. To access the data execute the line "global" after running
% All inputs are hard-coded below.

% function myData = plotTrajEph()
clear;
% global cubesat ceresData jupiterData earthData marsData saturnData;
au2km = 149597870.7; % km (source obspm.fr, 29/01/2014)
scnPath = '../ifod/inputs/Trajectories/Cas_Y/';
[NbL0, T0, xyz, vel] = readTraj([scnPath 'Y-line_vts.xyzv'], 1./au2km);
[NbL1, T1, lt1, lg1, ds1] = readEphem([scnPath 'Y-line_P001_vts.eph'], 1.);
[NbL2, T2, lt2, lg2, ds2] = readEphem([scnPath 'Y-line_P002_vts.eph'], 1.);
[NbL3, T3, lt3, lg3, ds3] = readEphem([scnPath 'Y-line_P003_vts.eph'], 1.);
[NbL4, T4, lt4, lg4, ds4] = readEphem([scnPath 'Y-line_P004_vts.eph'], 1.);

t0c=cubesat(1,1);
t0p=ceresData(1,1)-t0c;

clrs(1,1:3)=[0 0 1];
clrs(2,1:3)=[1 0 0];
clrs(3,1:3)=[0 1 0];
clrs(4,1:3)=[.5 0 1];
clrs(5,1:3)=[0 1 1];
clrs(6,1:3)=[1 1 0];

figure(10);
plot3(cubesat(:,2),cubesat(:,3),cubesat(:,4),'b');
hold on;
plot3(ceresData(:,2),ceresData(:,3),ceresData(:,4),'k');
plot3(jupiterData(:,2),jupiterData(:,3),jupiterData(:,4),'k');
plot3(earthData(:,2),earthData(:,3),earthData(:,4),'k');
plot3(marsData(:,2),marsData(:,3),marsData(:,4),'k');

figure(11); clf;
% mettre en AU, echelle egales (log? fonction loglog(X1,Y1)), légendes => puis montrer la
% différence entre 0dv et 1m/s dv "expected"
plot(cubesat(:,2),cubesat(:,3),'.k');
hold on; axis equal; minx=0.;maxx=0.;miny=0.;maxy=0.;
plot(earthData(:,2),earthData(:,3), 'color', clrs(1,1:3));
minx = min([minx earthData(:,2)']); maxx = max([maxx earthData(:,2)']); miny = min([miny earthData(:,3)']); maxy = max([maxy earthData(:,3)']);
plot(marsData(:,2),marsData(:,3), 'color', clrs(2,1:3));
minx = min([minx marsData(:,2)']); maxx = max([maxx marsData(:,2)']); miny = min([miny marsData(:,3)']); maxy = max([maxy marsData(:,3)']);
plot(ceresData(:,2),ceresData(:,3), 'color', clrs(3,1:3));
minx = min([minx ceresData(:,2)']); maxx = max([maxx ceresData(:,2)']); miny = min([miny ceresData(:,3)']); maxy = max([maxy ceresData(:,3)']);
%plot(jupiterData(:,2),jupiterData(:,3), 'color', clrs(4,1:3));
%minx = min([minx jupiterData(:,2)']); maxx = max([maxx jupiterData(:,2)']); miny = min([miny jupiterData(:,3)']); maxy = max([maxy jupiterData(:,3)']);
%plot(saturnData(:,2),saturnData(:,3), 'color', clrs(5,1:3));
%minx = min([minx saturnData(:,2)']); maxx = max([maxx saturnData(:,2)']); miny = min([miny saturnData(:,3)']); maxy = max([maxy saturnData(:,3)']);
legend('CubeSat', 'Earth','Mars','Ceres','Jupiter','Saturn');
axis([floor(minx)-0.5 floor(maxx+1) floor(miny)-0.5 floor(maxy)+1]);
myData=1;
% end