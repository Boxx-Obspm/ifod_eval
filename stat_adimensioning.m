%%----------------HEADER---------------------------%%
%Author:           Boris Segret
% Version & Date:
%                  V1.1, 05-06-2017
%                  - fixes & settings
%                  V1.0, 28-05-2017, from stat_processing in v1.3 26-05-2017
%
% CL=2 (v1.0)
%
% The expected shifts and residual errors are expressed in dimensionmless values.
% The optical error is used to set the scale factor over the distance of the 3rd
% foreground body.
%
% I/
%    <result_file> & <optical error>
%    needed trajectories: <foreground object's>, <ref> and <actual> trajectories
% O/
%    fig.400: transversal errors without and with KF, adimensioned
%    fig.410: longitudinal errors without and with KF, adimensioned
%    fig.411: dist.from ref (km), scale factor (km), ang.drift of Foreground body 3
%

clear;
% outputs_bin = '../cas_DAV/DAV_01as,50MCx8KF,t_bin'; opterr=0.1;
outputs_bin = '../cas_EME/outs/E+rgECMJE_01as,450MCx192KF,1393809_bin'; opterr=0.1;
% outputs_bin = '../cas_Y/outs/Y0_41324_01as,400MCx192KF,1393470_bin'; opterr=0.1;
% outputs_bin = '../cas_Y/outs/YYv4_01as,400MCx192KF,1420098_bin'; opterr=0.1;
fg_body_Nb3 = '../cas_EME/Mars_imcce.xva';
actual_traj = '../cas_EME/MI/58122+SOI_v6.4_jdv-1_312_vts.xyzv';
ref_traj    = '../cas_EME/58122+SOI_v6.4_jdv000_312_vts.xyzv';

T0N=datenum([2000 1 1 0 0 0]);
T0JD=2451544.5;

% reads foreground body's and actual trajectory
addpath('../ifod');
[NbLT0, TimeList0, data] = readTraj(ref_traj);
coord0 = data(:,1:3); vel0   = data(:,4:6);
clear data; fprintf('Reading ref.traj''s data: %s, %.1f...%.1f = %i records\n', ...
  ref_traj, TimeList0(1), TimeList0(NbLT0), NbLT0);
[NbLT1, TimeList1, data] = readTraj(actual_traj);
coord1 = data(:,1:3); vel1   = data(:,4:6);
clear data; fprintf('Reading act.traj''s data: %s, %.1f...%.1f = %i records\n', ...
  actual_traj, TimeList1(1), TimeList1(NbLT1), NbLT1);
[nbLgn, timeSteps, data] = readTraj(fg_body_Nb3);
coord3 = data(:,1:3); %vel3   = data(:,4:6);
clear data; fprintf('Body #3 trajecto''s data: %s, %.1f...%.1f = %i records\n', ...
  fg_body_Nb3, timeSteps(1), timeSteps(nbLgn), nbLgn);
clear data;

% initizes the binary outputs with the number of written columns
fw = fopen(outputs_bin,'r');

lastData = false; ntStep = 0;
while not(lastData)
% fwrite(fw, [obstime(ik+nKF-2) obstime(ik+nKF)], 'double');
T = fread(fw, 2, 'double'); % will provide a 1x2 double array
%                             (time for KF-measured and for KF-predicted)
t1=T(1)-T0JD+T0N; % epoch of the last 3D geometric solution (3rd point)
t2=T(2)-T0JD+T0N; % epoch of the last KF-anticipated solution (5th point)
dtKF=(t2-t1)/2.;
if T(2)>2400000.5
    T(1)=T(1)-2400000.5;
    T(2)=T(2)-2400000.5;
end

% fwrite(fw, [Nobs nKF nbCycles (ik+nKF==length(obstime))], 'uint32');
Z = fread(fw, 4, 'uint32'); % will provide a 1x3 uint32 array
% fwrite(fw, ...
%     [rex rrme rrkf lKg ldP vkf dtrm' dlgm' dtrk' dlgk' mmkf' mmtk' mmlk'], ...
%     'double');
Nobs=Z(1);
nKF=Z(2); iKF=floor(nKF/3); jKF=floor(nKF*2/3);
nbCycles=Z(3);
% nbPts=Z(2)-Z(1)+1; % (nbPts=nKF-Nobs+1;)
nbPts=nKF;
lastData = (Z(4)==1);

if (ntStep==0)
%   t0=t2;
  t0=t2-nKF*Nobs*dtKF;
  nstats=1; % nb of values in the statistical post-processing
  nvKF=11;  % nb of stored values per KF-step
  ep = double(zeros(100,1));
  moyKF=double(zeros(100,nvKF)); moyiKF=double(zeros(100,nvKF)); moyjKF=double(zeros(100,nvKF));
  stdKF=double(zeros(100,nvKF)); stdiKF=double(zeros(100,nvKF)); stdjKF=double(zeros(100,nvKF));
  moyMd=double(zeros(100,5));
  stdMd=double(zeros(100,5));
  drift=double(zeros(100,6));
  stats=double(zeros(100,nstats));
  ntStep=1;
else
  if (mod(ntStep,100)==0)
      % memory re-allocation every 100 steps
      Y=ep;
      ep=double(zeros(ntStep+100, 1));
      ep(1:ntStep)=Y; clear Y;
      Y=moyKF;
      moyKF=double(zeros(ntStep+100, nvKF));
      moyKF(1:ntStep,:)=Y; clear Y;
      Y=moyiKF;
      moyiKF=double(zeros(ntStep+100, nvKF));
      moyiKF(1:ntStep,:)=Y; clear Y;
      Y=moyjKF;
      moyjKF=double(zeros(ntStep+100, nvKF));
      moyjKF(1:ntStep,:)=Y; clear Y;
      Y=stdKF;
      stdKF=double(zeros(ntStep+100, nvKF));
      stdKF(1:ntStep,:)=Y; clear Y;
      Y=stdiKF;
      stdiKF=double(zeros(ntStep+100, nvKF));
      stdiKF(1:ntStep,:)=Y; clear Y;
      Y=stdjKF;
      stdjKF=double(zeros(ntStep+100, nvKF));
      stdjKF(1:ntStep,:)=Y; clear Y;
      Y=moyMd;
      moyMd=double(zeros(ntStep+100, 5));
      moyMd(1:ntStep,:)=Y; clear Y;
      Y=stdMd;
      stdMd=double(zeros(ntStep+100, 5));
      stdMd(1:ntStep,:)=Y; clear Y;
      Y=drift;
      drift=double(zeros(ntStep+100, 6));
      drift(1:ntStep,:)=Y; clear Y;
      Y=stats;
      stats=double(zeros(ntStep+100, nstats));
      stats(1:ntStep,:)=Y; clear Y;
  end
  ntStep = ntStep+1;
end

ep(ntStep)=t2-t0;
stf=datestr(t2,'yyyy-mm-dd HH:MM:SS');
dtKF=(t2-t1)/2.;
sdt=['~' num2str(24*Nobs*dtKF,'%4.1f') 'h/step'];

solKF = double(zeros(nbCycles, nvKF));
soliKF = double(zeros(nbCycles, nvKF));
soljKF = double(zeros(nbCycles, nvKF));
solMd = double(zeros(nbCycles*nbPts, 5));
method = 'linear';
timeB3 = double([T(2)-2*dtKF-(nKF-1)*(Nobs*dtKF):(Nobs*dtKF):T(2)-2*dtKF T(2)]); % nKF+1 values
vec3B3 = interp1(timeSteps, coord3, timeB3);
vec0B3 = vec3B3 - interp1(TimeList0, coord0, timeB3);
% norm0B3 = sqrt(vec0B3(:,1).^2 + vec0B3(:,2).^2 + vec0B3(:,3).^2);
vec1B3 = vec3B3 - interp1(TimeList1, coord1, timeB3);
norm1B3 = sqrt(vec1B3(:,1).^2 + vec1B3(:,2).^2 + vec1B3(:,3).^2);
scalef = (opterr/3600.).*(pi()/180.).*norm1B3;
drift(ntStep,1:3)= vec1B3(nKF+1,:)-vec0B3(nKF+1,:);
drift(ntStep,5)  = sqrt(drift(ntStep,1).^2 + drift(ntStep,2).^2 + drift(ntStep,3).^2);
drift(ntStep,4)  = asind(drift(ntStep,5)./norm1B3(nKF+1)); % parallaxis in degrees
% drift(ntStep,5)  = drift(ntStep,5)./scalef(nKF+1); % later, it should be transv. vs. long. shift
drift(ntStep,6)  = scalef(nKF+1); % (at the moment) meters
drift(ntStep,1:3)= drift(ntStep,1:3)./scalef(nKF+1);
clear vec0B3 vec1B3 vec3B3 norm1B3 timeB3;

for nC=1:nbCycles
  % il faudrait cycler
  rawDATA = fread(fw, (3*3+10)*(nbPts+1), 'double');
  % will read the detailed value for 1 Monte-Carlo cycle,
  % to be reshaped into (nbPts+1) lines with:
  DATA = reshape(rawDATA, (3*3+10), nbPts+1)'; % (nbPts+1) lines x (3*3+10) columns
  
  rex  = DATA(:,1:3)./repmat(scalef, [1 3]); % Expected results
  rrme = DATA(:,4:6)./repmat(scalef, [1 3]); % "Measured" results -rex (3D-OD)
  rrkf = DATA(:,7:9)./repmat(scalef, [1 3]); % "Kalman Filtered" results -rex
%   lKg  = DATA(:,10);
%   ldP  = DATA(:,11);
%   vkf  = DATA(:,12);  % "Kalman Filtered" velocity
  dtrm = (DATA(:,13)./scalef)'; % transverse residual of "measured" results
  dlgm = (DATA(:,14)./scalef)'; % longitudin residual of "measured" results
  dtrk = (DATA(:,15)./scalef)'; % transverse residual of "K-Filter" results
  dlgk = (DATA(:,16)./scalef)'; % longitudin residual of "K-Filter" results
  mmkf = (DATA(:,17)./scalef)'; % max of the remaining norm residuals in KF-results
  mmtk = (DATA(:,18)./scalef)'; % max of the remaining transverse residuals in KF-results
  mmlk = (DATA(:,19)./scalef)'; % max of the remaining longitudin residuals in KF-results
  
  nbCol=1;
  solKF(nC, nbCol:nbCol+2) = rex(nbPts+1,1:3); %1-3 EXPECTED
  soliKF(nC, nbCol:nbCol+2) = rex(iKF,1:3); %1-3 EXPECTED
  soljKF(nC, nbCol:nbCol+2) = rex(jKF,1:3); %1-3 EXPECTED
  nbCol=nbCol+3;
  solKF(nC, nbCol:nbCol+2) = rrkf(nbPts+1,1:3); %4-6 KF solution
  soliKF(nC, nbCol:nbCol+2) = rrkf(iKF,1:3); %4-6 KF solution
  soljKF(nC, nbCol:nbCol+2) = rrkf(jKF,1:3); %4-6 KF solution
  nbCol=nbCol+3;
  solKF(nC, nbCol) = dtrk(nbPts+1); %7 Transverse error with KF solution
  soliKF(nC, nbCol) = dtrk(iKF); %7 Transverse error with KF solution
  soljKF(nC, nbCol) = dtrk(jKF); %7 Transverse error with KF solution
  nbCol=nbCol+1;
  solKF(nC, nbCol) = dlgk(nbPts+1); %8 Longitudinal error with KF solution
  soliKF(nC, nbCol) = dlgk(iKF); %8 Longitudinal error with KF solution
  soljKF(nC, nbCol) = dlgk(jKF); %8 Longitudinal error with KF solution
  nbCol=nbCol+1;
  solKF(nC, nbCol) = mmkf(nbPts+1); nbCol=nbCol+1; %9
  solKF(nC, nbCol) = mmtk(nbPts+1); nbCol=nbCol+1; %10
  solKF(nC, nbCol) = mmlk(nbPts+1); nbCol=nbCol+1; %11

  solMd(1+(nC-1)*nbPts:nC*nbPts, 1:3) = rrme(1:nbPts, 1:3);
  solMd(1+(nC-1)*nbPts:nC*nbPts, 4) = dtrm(1:nbPts);
  solMd(1+(nC-1)*nbPts:nC*nbPts, 5) = dlgm(1:nbPts);
end

% statistics for ntStep:
moyKF(ntStep,:) = mean(solKF);
stdKF(ntStep,:) =  std(solKF);
moyiKF(ntStep,:) = mean(soliKF);
stdiKF(ntStep,:) =  std(soliKF);
moyjKF(ntStep,:) = mean(soljKF);
stdjKF(ntStep,:) =  std(soljKF);
moyMd(ntStep,:) = mean(solMd);
stdMd(ntStep,:) =  std(solMd);
if (mod(ntStep,20)==0) fprintf('(%i)',ntStep); else fprintf('.'); end;
if (mod(ntStep,100)==0) fprintf('\n'); end;

end
fprintf('\n');
fclose(fw);

%% ------------------------------
% TRANSVERSAL
figure(400); clf; ylim([ -10 10]); xlim([0 t2-t0]); hold on; set(gca, 'LooseInset', [0 0 0 0]);
% figure(400);
% xlim([0 t2-t0]);
errorbar(ep(1:ntStep), moyMd(1:ntStep,4), stdMd(1:ntStep,4), '-r');
errorbar(ep(1:ntStep), moyKF(1:ntStep,7), stdKF(1:ntStep,7), '-k');
% plot(ep(1:ntStep), moyMd(1:ntStep,4), '-r');
% plot(ep(1:ntStep), moyKF(1:ntStep,7), '-k');
xlabel('time (days)'); ylabel('Residual (in units of scale factor)'); 
title(['Transversal errors, ' datestr(t0,'yyyy-mm-dd') '..' datestr(t2,'yyyy-mm-dd')]);
legend('without KF', [num2str(nKF) '-step KF']);
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

% LONGITUDINAL
figure(410); clf; ylim([ -10 10]); xlim([0 t2-t0]); hold on; set(gca, 'LooseInset', [0 0 0 0]);
% figure(410);
% xlim([0 t2-t0]);
errorbar(ep(1:ntStep), moyMd(1:ntStep,5), stdMd(1:ntStep,5), '-r');
errorbar(ep(1:ntStep), moyKF(1:ntStep,8), stdKF(1:ntStep,8), '-k');
% plot(ep(1:ntStep), moyMd(1:ntStep,5), '-r');
% plot(ep(1:ntStep), moyKF(1:ntStep,8), '-k');
xlabel('time (days)'); ylabel('Residual (in units of scale factor)');
title(['Longitudinal errors, ' datestr(t0,'yyyy-mm-dd') '..' datestr(t2,'yyyy-mm-dd')]);
legend('without KF', [num2str(nKF) '-step KF']);
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

% DRIFT from Ref
figure(411); clf; %set(gca, 'LooseInset', [0.09 0.09 0.07 0.07]);
[AX,dst,ang] = plotyy(ep(1:ntStep), drift(1:ntStep,5), ep(1:ntStep), drift(1:ntStep,4), ...
    'semilogy', 'semilogy'); hold on;
% help "Using Multiple X- and Y-Axes"
set(AX(1),'YGrid','on'); set(AX(1),'YMinorGrid','off');
set(dst, 'Color', [0 0 1]);  set(AX(1),'YLim',[1 1E5], 'YColor', [0 0 1]); 
set(get(AX(1),'Ylabel'),'String','(kilometers)');
set(ang, 'Color', [0 .5 0]); set(AX(2),'YLim',[0.0001 10], 'YColor', [0 .5 0], 'XColor', 'm'); 
set(get(AX(2),'Ylabel'),'String','(degrees)');
ax3 = line(ep(1:ntStep), drift(1:ntStep,6),'LineStyle','--','Color','b','Parent',AX(1));
legend('Distance drift from Ref.', 'Scale factor (km)', 'Angular drift to FGB', 'Location', 'NorthWest');
xlabel('time (days)');
set(gca, 'XColor', 'm');
title(['Drift wrt Ref, ' datestr(t0,'yyyy-mm-dd') '..' datestr(t2,'yyyy-mm-dd')]);

% HISTOGRAMS
smootime = ep(1):(ep(ntStep)-ep(1))/1000.:ep(ntStep);
figure(412); clf; hold on;
smoothed = interp1(ep(1:ntStep),stdKF(1:ntStep,7),smootime,method);
% hist(stdKF(1:ntStep,7), [0:0.1:10]); xlim([0 2]);
hist(smoothed, [0:0.05:10]); xlim([0 3]);
h = findobj(gca,'Type','patch'); set(h,'FaceColor','g','EdgeColor','w');
smoothed = interp1(ep(1:ntStep),stdKF(1:ntStep,8),smootime,method);
% hist(stdKF(1:ntStep,8), [0:0.1:10]); xlim([0 2]);
hist(smoothed, [0:0.05:10]); xlim([0 3]);
set(gca, 'ytick', []);
legend('transversely', 'longitudinaly');
title('Hist.of standard deviations');
xlabel('Std.Dev in units of scale factor');
set(gca, 'XColor', 'm');

% figure(120);
% %
% subplot(3,1,1);
% % xlim([0 t2-t0]);
% errorbar(t2-t0, moyMd(ntStep,1), stdMd(ntStep,1), '-r');
% errorbar(t2-t0, moyKF(ntStep,4), stdKF(ntStep,4), '-k');
% % plot(ep(1:ntStep), moyMd(1:ntStep,1), '-r');
% % plot(ep(1:ntStep), moyKF(1:ntStep,4), '-k');
% xlabel('time (days)'); ylabel('X residual (sc.fact.)'); 
% set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
% %
% subplot(3,1,2);
% % xlim([0 t2-t0]);
% errorbar(t2-t0, moyMd(ntStep,2), stdMd(ntStep,2), '-r');
% errorbar(t2-t0, moyKF(ntStep,5), stdKF(ntStep,5), '-k');
% % plot(ep(1:ntStep), moyMd(1:ntStep,2), '-r');
% % plot(ep(1:ntStep), moyKF(1:ntStep,5), '-k');
% xlabel('time (days)'); ylabel('Y residual (sc.fact.)'); 
% set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
% %
% subplot(3,1,3);
% % xlim([0 t2-t0]);
% errorbar(t2-t0, moyMd(ntStep,3), stdMd(ntStep,3), '-r');
% errorbar(t2-t0, moyKF(ntStep,6), stdKF(ntStep,6), '-k');
% % plot(ep(1:ntStep), moyMd(1:ntStep,3), '-r');
% % plot(ep(1:ntStep), moyKF(1:ntStep,6), '-k');
% xlabel('time (days)'); ylabel('Z residual (sc.fact.)'); 
% set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

figure(400); ylim([ -10 10]); %figure(201); ylim([ 0 200]);
figure(410); ylim([ -12 12]); % figure(101); ylim([ -400 200]);
% figure(120);
% subplot(3,1,1); ylim([ -10  10]);
% subplot(3,1,2); ylim([ -10  10]);
% subplot(3,1,3); ylim([ -10  10]);
