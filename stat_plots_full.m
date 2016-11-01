%%----------------HEADER---------------------------%%
%Author:           Boris Segret
%Version & Date:   V3.2 16-08-2016 (dd/mm/yyyy)
%                  - adaptations to N=5 measurements
%                  REMAINING ISSUE: some figures are not compatible with multiple foreground objects
%CL=0
%Version & Date:
%                  V3.1 07-05-2016
%                  V3.0 27-04-2016
%                  (forked from data_plots_full.m, v2.1)
%                  V2.1 11-04-2016 Boris Segret
%                  - specific plots for quick analysis of results
%                  - adapted to output files produced by data_extraction.m in version v2.1
%                  (forked from data_plots.m)
%                  V1 11-09-2015 Tristan Mallet
%
%
% This program plots the critical outputs from the dataExtraction file.
%
% 1. Inputs: fname = name of the result file to read (in the workspace)
%
% 2. Outputs: some plots
%

fname = outputs;
ix0=16; % rank before the first X data
lx=26;  % length of the X-vector
exx=ix0+3*lx; % length of data

pf=fopen(fname,'rt');
while not(feof(pf))
  l=fgetl(pf);
  if not(isempty(l))
    if not(isempty(strfind(l,'META_STOP')))
      l=fgetl(pf);
      data = sscanf(l, '%f', [1 exx]);
      break;
    end;
    if strfind(l,'SCENARIO')>0
      ttl=l;
      scn=ttl(strfind(ttl,':')+1:length(ttl));
    end
  end;
end;
while not(feof(pf));
  l=fgetl(pf);
  if not(isempty(l))
    data = [data; sscanf(l, '%f', [1 exx])];
  end;
end;

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
subplot(7,2,[ 3  5]); hold on;
%plot(tt, data(:,16), 'ob', tt, data(:,35), 'xk');
plot(tt, data(:,ix0+lx+1), '-k');
errorbar(tt, data(:,ix0+1), data(:,ix0+2*lx+1), 'ob');
idy=find(abs(data(:,ix0+1))<=2*mm(ix0+1));
mn1=mean(data(idy,ix0+1));
sg1 = max([1 std(data(idy,ix0+1))]);
% ylim([mn1-3*sg1 mn1+3*sg1]);
%xlabel('time (days)');
ylabel('dX (km)');  %  legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
%
subplot(7,2,[ 7  9]); hold on;
%plot(tt, data(:,17), 'ob', tt, data(:,36), 'xk');
plot(tt, data(:,ix0+lx+2), '-k');
errorbar(tt, data(:,ix0+2), data(:,ix0+2*lx+2), 'ob');
idy=find(abs(data(:,ix0+2))<=2*mm(ix0+2));
mn1=mean(data(idy,ix0+2));
sg1 = max([1 std(data(idy,ix0+2))]);
% ylim([mn1-3*sg1 mn1+3*sg1]);
%xlabel('time (days)');
ylabel('dY (km)');  %  legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
%
subplot(7,2,[11 13]); hold on;
% plot(tt, data(:,18), 'ob', tt, data(:,37), 'xk');
plot(tt, data(:,ix0+lx+3), '-k');
errorbar(tt, data(:,ix0+3), data(:,ix0+2*lx+3), 'ob');
idy=find(abs(data(:,ix0+3))<=2*mm(ix0+3));
mn1=mean(data(idy,ix0+3));
sg1 = max([1 std(data(idy,ix0+3))]); % ylim([mn1-3*sg1 mn1+3*sg1]);
xlabel('time (days)'); ylabel('dZ (km)');  % legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
%
subplot(7,2,[ 4  6]); hold on; 
% plot(tt, data(:,32), 'ob', tt, data(:,51), 'xk');
plot(tt, data(:,ix0+2*lx-5), '-k');
errorbar(tt, data(:,ix0+lx-5), data(:,ix0+3*lx-5), 'ob');
idy=find(abs(data(:,ix0+lx-5))<=2*mm(ix0+lx-5));
mn1=mean(data(idy,ix0+lx-5));
sg1 = max([.0001 std(data(idy,ix0+lx-5))]);
% ylim([mn1-3*sg1 mn1+3*sg1]);
%xlabel('time (days)');
ylabel('dVx (km/s)'); %legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
%
subplot(7,2,[ 8 10]); hold on;
% plot(tt, data(:,33), 'ob', tt, data(:,52), 'xk');
plot(tt, data(:,ix0+2*lx-4), '-k');
errorbar(tt, data(:,ix0+lx-4), data(:,ix0+3*lx-4), 'ob');
idy=find(abs(data(:,ix0+lx-4))<=2*mm(ix0+lx-4));
mn1=mean(data(idy,ix0+lx-4));
sg1 = max([.0001 std(data(idy,ix0+lx-4))]); % ylim([mn1-3*sg1 mn1+3*sg1]);
%xlabel('time (days)');
ylabel('dVy (km/s)'); % legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
%
subplot(7,2,[12 14]); hold on;
% plot(tt, data(:,34), 'ob', tt, data(:,53), 'xk');
plot(tt, data(:,ix0+2*lx-3), '-k');
errorbar(tt, data(:,ix0+lx-3), data(:,ix0+3*lx-3), 'ob');
idy=find(abs(data(:,ix0+lx-3))<=2*mm(34));
mn1=mean(data(idy,ix0+lx-3));
sg1 = max([.0001 std(data(idy,ix0+lx-3))]); % ylim([mn1-3*sg1 mn1+3*sg1]);
xlabel('time (days)'); ylabel('dVz (km/s)'); % legend('reconstructed', 'expected');
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

figure(3); clf; % (dr), shift of distance to the foreground object
%-----------------------------------------------------------------
hold on; plot(tt, data(:,ix0+16), '-b', tt, data(:,ix0+lx+16), '-k');
idy=find(abs(data(:,ix0+16))<=2*mm(ix0+16)); 
mn1=mean(data(idy,ix0+16));
sg1 = max([1 std(data(idy,ix0+16))]);
ylim([mn1-3*sg1 mn1+3*sg1]);
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

