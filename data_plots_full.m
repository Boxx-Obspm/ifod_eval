%%----------------HEADER---------------------------%%
%Author:           Boris Segret
%Version & Date:   V2.1 11-04-2016 (dd/mm/yyyy)
%                  - specific plots for quick analysis of results
%                  - adapted to output files produced by data_extraction.m in version v2.1
%                  REMAINING ISSUE: some figures are not compatible with multiple foreground objects
%                  (forked from data_plots.m)
%CL=2
%                  V1 11-09-2015 (dd/mm/yyyy) Tristan Mallet
%
%
% This program plots the critical outputs from the dataExtraction file.
%
% 1. Inputs: fname = name of the result file to read (in the workspace)
%
% 2. Outputs: some plots
%

fname = outputs;

plot_file=fopen(fname,'rt');
l=' ';
while 1
  l=fgetl(plot_file);
  if strfind(l,'SCENARIO')>0
   ttl=l;
   scn=ttl(strfind(ttl,':')+1:length(ttl));
  end
  if strfind(l,'META_STOP')>0
    break;
  end;
end;

data=fscanf(plot_file, '%g', [53 inf]); data=data';
fclose(plot_file);
tt = -data(1,1)+data(:,1)+data(:,2)/86400.;
mm = median(abs(data));


figure(1); clf; % (transversal & longitudinal shifts)
%----------------------------------------------------
subplot(9,3,[1  3]); axis([-1 1 -1 1]);
text(0,0,ttl,'FontWeight','bold', 'horizontalalignment','center'); set(gca, 'Visible', 'off');
subplot(9,3,[4 13]); hold on; plot(tt, data(:,7), 'ob', tt, data(:,3), 'xk');
idy=find(abs(data(:,7))<=2*mm(7)); mn1=mean(data(idy,7)); sg1 = max([1 std(data(idy,7))]); ylim([mn1-3*sg1 mn1+3*sg1]);
ylabel('shift to Reference (km)'); title('Transversal shift');
legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
subplot(9,3,[16 25]); hold on; plot(tt, data(:,7)-data(:,3), '-r');
idy=find(abs(data(:,7)-data(:,3))<=mm(3)+mm(7)); mn1=mean(data(idy,7)-data(idy,3)); sg1 = max([.001 std(data(idy,7)-data(idy,3))]); ylim([mn1-3*sg1 mn1+3*sg1]);
xlabel('time (days)'); ylabel('shift difference (km)');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

subplot(9,3,[5 14]); hold on; plot(tt, data(:,8), 'ob', tt, data(:,4), 'xk');
idy=find(abs(data(:,8))<=2*mm(8)); mn2=mean(data(idy,8)); sg2 = max([1 std(data(idy,8))]); ylim([mn2-3*sg2 mn2+3*sg2]);
ylabel('shift to Reference (km)'); title('Longitudinal shift');
legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
subplot(9,3,[17 26]); hold on; plot(tt, data(:,8)-data(:,4), '-r');
idy=find(abs(data(:,8)-data(:,4))<=median(abs(data(:,8)-data(:,4)))); mn2=mean(data(idy,8)-data(idy,4)); sg2 = max([.001 std(data(idy,8)-data(idy,4))]);
ylim([mn2-3*sg2 mn2+3*sg2]);
xlabel('time (days)'); ylabel('shift difference (km)');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

subplot(9,3,[6 15]); hold on; plot(data(:,7), data(:,8), 'ob', data(:,3), data(:,4), 'xk');
mn=mean(data(:,3)); mn1=(mn1+mn)/2.; sg1=max([sg1 abs(mn1-mn)]); xlim([mn1-sg1 mn1+sg1]);
mn=mean(data(:,4)); mn2=(mn2+mn)/2.; sg2=max([sg2 abs(mn2-mn)]); ylim([mn2-sg2 mn2+sg2]);
xlabel('Transversal shift (km)'); ylabel('Longitudinal shift (km)'); title('(Longit.,Transv.) shifts');
legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
subplot(9,3,[18 27]); hold on; plot(data(:,7)-data(:,3), data(:,8)-data(:,4), '*r');
mn1=mean(data(:,7)-data(:,3)); sg1=max([.001 std(data(:,7)-data(:,3))]); xlim([mn1-3*sg1 mn1+3*sg1]);
mn2=mean(data(:,8)-data(:,4)); sg2=max([.001 std(data(:,8)-data(:,4))]); ylim([mn2-3*sg2 mn2+3*sg2]);
xlabel('Transversal diff. (km)'); ylabel('Longitudinal diff (km)');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

figure(2); clf; % (dX, Vx; dY,Vy; dZ,Vz)
%-------------------------------------------
subplot(7,2,[ 1  2]); axis([-1 1 -1 1]);
text(0,0,ttl,'FontWeight','bold', 'horizontalalignment','center'); set(gca, 'Visible', 'off');
%
subplot(7,2,[ 3  5]); hold on; plot(tt, data(:,16), 'ob', tt, data(:,35), 'xk');
idy=find(abs(data(:,16))<=2*mm(16)); mn1=mean(data(idy,16)); sg1 = max([1 std(data(idy,16))]); ylim([mn1-3*sg1 mn1+3*sg1]);
%xlabel('time (days)');
ylabel('dX (km)');    legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
%
subplot(7,2,[ 7  9]); hold on; plot(tt, data(:,17), 'ob', tt, data(:,36), 'xk');
idy=find(abs(data(:,17))<=2*mm(17)); mn1=mean(data(idy,17)); sg1 = max([1 std(data(idy,17))]); ylim([mn1-3*sg1 mn1+3*sg1]);
%xlabel('time (days)');
ylabel('dY (km)');    legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
%
subplot(7,2,[11 13]); hold on; plot(tt, data(:,18), 'ob', tt, data(:,37), 'xk');
idy=find(abs(data(:,18))<=2*mm(18)); mn1=mean(data(idy,18)); sg1 = max([1 std(data(idy,18))]); ylim([mn1-3*sg1 mn1+3*sg1]);
xlabel('time (days)'); ylabel('dZ (km)');    legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
%
subplot(7,2,[ 4  6]); hold on; plot(tt, data(:,32), 'ob', tt, data(:,51), 'xk');
idy=find(abs(data(:,32))<=2*mm(32)); mn1=mean(data(idy,32)); sg1 = max([.0001 std(data(idy,32))]); ylim([mn1-3*sg1 mn1+3*sg1]);
%xlabel('time (days)');
ylabel('dVx (km/s)'); legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
%
subplot(7,2,[ 8 10]); hold on; plot(tt, data(:,33), 'ob', tt, data(:,52), 'xk');
idy=find(abs(data(:,33))<=2*mm(33)); mn1=mean(data(idy,33)); sg1 = max([.0001 std(data(idy,33))]); ylim([mn1-3*sg1 mn1+3*sg1]);
%xlabel('time (days)');
ylabel('dVy (km/s)'); legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
%
subplot(7,2,[12 14]); hold on; plot(tt, data(:,34), 'ob', tt, data(:,53), 'xk');
idy=find(abs(data(:,34))<=2*mm(34)); mn1=mean(data(idy,34)); sg1 = max([.0001 std(data(idy,34))]); ylim([mn1-3*sg1 mn1+3*sg1]);
xlabel('time (days)'); ylabel('dVz (km/s)'); legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

figure(3); clf; % (dr), shift of distance to the foreground object
%-----------------------------------------------------------------
hold on; plot(tt, data(:,28), '-b', tt, data(:,47), '-k');
idy=find(abs(data(:,28))<=2*mm(28)); mn1=mean(data(idy,28)); sg1 = max([1 std(data(idy,28))]); ylim([mn1-3*sg1 mn1+3*sg1]);
xlabel('time (days)'); ylabel('dr (km)');    legend('reconstructed', 'expected');
title([scn ' - Distance shift to foreground body'], 'FontWeight','bold');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

figure(4); clf; % uniformity of speed
%------------------------------------
subplot(5,2,[1  2]); axis([-1 1 -1 1]);
text(0,0,[scn ' - Uniformity of speed during OD'], 'FontWeight','bold', 'horizontalalignment','center');
set(gca, 'Visible', 'off');
subplot(5,2,[3  9]); hold on; plot(tt, data(:,9)/3600., '-b'); %, tt, data(:,9), '-k');
xlabel('time (days)'); ylabel('Deviation of speed (deg)');
%legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
subplot(5,2,[4 10]); hold on; plot(tt, data(:,10), '-b'); %, tt, data(:,10), '-k');
xlabel('time (days)'); ylabel('Increasing of speed (mm/s)');
%legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

figure(5); clf; % (transversal & longitudinal shifts)
%----------------------------------------------------
subplot(10,3,[1  3]); axis([-1 1 -1 1]);
text(0,0,[ttl ' - Direction diff.of foreground body'],'FontWeight','bold', 'horizontalalignment','center');
set(gca, 'Visible', 'off');
subplot(10,3,[4 28]); hold on; plot(tt, data(:,5), '-b'); %, tt, data(:,3), '-k');
xlabel('time (days)'); ylabel('Latitude diff. (arcsec)'); set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
subplot(10,3,[5 29]); hold on; plot(tt, data(:,6), '-b'); %, tt, data(:,4), '-k');
xlabel('time (days)'); ylabel('Longitude diff. (arcsec)'); set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
subplot(10,3,[6 30]); hold on; plot(data(:,6), data(:,5), 'ob'); %, data(:,3), data(:,4), 'ok');
xlabel('Longitude diff.(arcsec)'); ylabel('Latitude diff.(arcsec)'); set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

%- figure savings to PNG
ff=[scn '_shifts.png']; print(1,ff,'-dpng');
ff=[scn '_dXdV.png']; print(2,ff,'-dpng');
ff=[scn '_dr.png']; print(3,ff,'-dpng');
ff=[scn '_uniV.png']; print(4,ff,'-dpng');
ff=[scn '_lglt.png']; print(5,ff,'-dpng');
%- figure savings to SVG (may crash, thus at the end)
%ff=[scn '_shifts.svg']; print(1,ff,'-dsvg');
%ff=[scn '_dXdV.svg']; print(2,ff,'-dsvg');
%ff=[scn '_dr.svg']; print(3,ff,'-dsvg');
%ff=[scn '_uniV.svg']; print(4,ff,'-dsvg');
%ff=[scn '_lglt.svg']; print(5,ff,'-dsvg');

