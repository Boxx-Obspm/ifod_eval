%%----------------HEADER---------------------------%%
%Author:           Boris Segret
% Version & Date:
%                  V1.1, 25-04-2017
%                  - additional plots
%                  until V1   23-03-2017
%
% CL=2 (v1.1)
%
%
% This produces basic graphs of the reconstructed OD with error bars,
% as stored in the binary file of results of a given scenario.
% 
% Option: it can also plot the bahavior of each Kalman filter step
%
% I/
%    <result_file>
%    <option "graphs">
% O/
%    <shift errors with error bars in transverse and longitudinal directions>
%    <shift errors with error bars in X, Y, Z, axis>
%    (with option "graphs") <shifts over the Kalman filter>


clear;
% outputs_bin = '../cas_EME/outs/E0_ECMJE_01as,450MCx192KF,1393471_bin';
% outputs_bin = '../cas_Y/outs/Y0_41324_01as,400MCx192KF,1393470_bin';
% outputs_bin = '../cas_Y/outs/Y+Y_E-1,400x200h,P41324_fixed_bin';
outputs_bin = '../ifod_tests/outs/Y0_41324v4_01as,400MCx192KF,t_bin'; % non-fonctionnel
graphs = true;
graphs = false;

T0N=datenum([2000 1 1 0 0 0]);
T0JD=2451544.5;

figure(200); clf; xlim([0  200]); hold on;
figure(101); clf;
subplot(1,2,1);
ylim([ -600  600]); xlim([0  200]); hold on;
% subplot(1,2,2); ylim([-2500 500]); xlim([ -1  200]); hold on;
subplot(1,2,2); ylim([-600 600]); xlim([0  200]); hold on;
figure(120); clf;
subplot(3,1,1); ylim([ -300  300]); xlim([0  200]); hold on;
subplot(3,1,2); ylim([ -300  300]); xlim([0  200]); hold on;
subplot(3,1,3); ylim([ -300  300]); xlim([0  200]); hold on;

% initizes the binary outputs with the number of written columns
fw = fopen(outputs_bin,'r');

lastData = false; ntStep = 0;
while not(lastData)
% fwrite(fw, [obstime(ik+nKF-2) obstime(ik+nKF)], 'double');
T = fread(fw, 2, 'double'); % will provide a 1x2 double array
%                             (time for KF-measured and for KF-predicted)
t1=T(1)-T0JD+T0N;
t2=T(2)-T0JD+T0N;

% fwrite(fw, [Nobs nKF nbCycles (ik+nKF==length(obstime))], 'uint32');
Z = fread(fw, 4, 'uint32'); % will provide a 1x3 uint32 array
% fwrite(fw, ...
%     [rex rrme rrkf lKg ldP vkf dtrm' dlgm' dtrk' dlgk' mmkf' mmtk' mmlk'], ...
%     'double');
Nobs=Z(1);
nKF=Z(2); dtKF=12./1440.; % dtKF, in decimal days, not stored!!!
nbCycles=Z(3);
% nbPts=Z(2)-Z(1)+1; % (nbPts=nKF-Nobs+1;)
nbPts=nKF;
lastData = (Z(4)==1);

if (ntStep==0)
%   t0=t2;
  t0=t2-nKF*Nobs*dtKF;
  nstats=1; % nb of values in the statistical post-processing
  nvKF=19;  % nb of stored values per KF-step
  ep = double(zeros(100,1));
  moyKF=double(zeros(100,nvKF));
  stdKF=double(zeros(100,nvKF));
  moyMd=double(zeros(100,5));
  stdMd=double(zeros(100,5));
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
      moyKF(1:ntStep, 1:nvKF)=Y; clear Y;
      Y=stdKF;
      stdKF=double(zeros(ntStep+100, nvKF));
      stdKF(1:ntStep, 1:nvKF)=Y; clear Y;
      Y=moyMd;
      moyMd=double(zeros(ntStep+100, 2));
      moyMd(1:ntStep, 1:5)=Y; clear Y;
      Y=stdMd;
      stdMd=double(zeros(ntStep+100, 2));
      stdMd(1:ntStep, 1:5)=Y; clear Y;
      Y=stats;
      stats=double(zeros(ntStep+100, nstats));
      stats(1:ntStep, 1:nstats)=Y; clear Y;
  end
  ntStep = ntStep+1;
end

ep(ntStep)=t2-t0;
stf=datestr(t2,'yyyy-mm-dd HH:MM:SS');
dtKF=24*Nobs*((t2-t1)/2);
sdt=['~' num2str(dtKF,'%4.1f') 'h/step'];

solKF = double(zeros(nbCycles, nvKF));
solMd = double(zeros(nbCycles*nbPts, 5));

for nC=1:nbCycles
  % il faudrait cycler
  rawDATA = fread(fw, (3*3+10)*(nbPts+1), 'double');
  % will read the detailed value for 1 Monte-Carlo cycle,
  % to be reshaped into (nbPts+1) lines with:
  DATA = reshape(rawDATA, (3*3+10), nbPts+1)'; % (nbPts+1) lines x (3*3+10) columns
  
  nbCol=1;
  rex  = DATA(:,nbCol:nbCol+2); nbCol=nbCol+3;
  rrme = DATA(:,nbCol:nbCol+2); nbCol=nbCol+3;
  rrkf = DATA(:,nbCol:nbCol+2); nbCol=nbCol+3;
  lKg  = DATA(:,nbCol); nbCol=nbCol+1;
  ldP  = DATA(:,nbCol); nbCol=nbCol+1;
  vkf  = DATA(:,nbCol); nbCol=nbCol+1;
  dtrm = DATA(:,nbCol)'; nbCol=nbCol+1;
  dlgm = DATA(:,nbCol)'; nbCol=nbCol+1;
  dtrk = DATA(:,nbCol)'; nbCol=nbCol+1;
  dlgk = DATA(:,nbCol)'; nbCol=nbCol+1;
  mmkf = DATA(:,nbCol)'; nbCol=nbCol+1;
  mmtk = DATA(:,nbCol)'; nbCol=nbCol+1;
  mmlk = DATA(:,nbCol)'; nbCol=nbCol+1;
  
  if (graphs)
    % chronogrammes
    figure(102); clf;
    subplot(4,3, [1 3]);
    plot(rrkf(:,1), 'r-'); hold on;
    plot(rrkf(:,2), 'g-');
    plot(rrkf(:,3), 'b-');
    plot(rrme(:,1), 'rx:');
    plot(rrme(:,2), 'gx:');
    plot(rrme(:,3), 'bx:');
%     ylim([min(min(rrme)) max(max(rrme))]);
   ylim([-5*mmkf(floor(nbPts*3/4)) 5*mmkf(floor(nbPts*3/4))]);
    title(['Residuals dX, dY, dZ (km) ' sdt ', ' stf]);
    legend('dX (-KF, ..3D-OD)', 'dY (-KF, ..3D-OD)', ...
        'dZ (-KF, ..3D-OD)', 'Location', 'SouthWest');
    plot(mmkf, 'k:'); plot(-mmkf, 'k:');
    subplot(4,3, [4 6]);
    plot(dtrk, 'r-'); hold on;
    plot(dtrm, 'gx:');
%     ylim([-2*std(dtrm) 2*std(dtrm)]);
    ylim([-5*mmtk(floor(nbPts*3/4)) 5*mmtk(floor(nbPts*3/4))]);
    plot(mmtk, 'r:'); plot(-mmtk, 'r:');
    title(['Residual in transversal shift (km) ' sdt ', ' stf]);
    subplot(4,3, [7 9]);
    plot(dlgk, 'r-'); hold on;
    plot(dlgm, 'gx:');
%     ylim([-2*std(dlgm) 2*std(dlgm)]);
    ylim([-5*mmlk(floor(nbPts*3/4)) 5*mmlk(floor(nbPts*3/4))]);
    plot(mmlk, 'r:'); plot(-mmlk, 'r:');
    title(['Residual in longitudinal shift (km) ' sdt ', ' stf]);

    subplot(4,3, 10);
    plot(rrme(:,1),rrme(:,2),'gx:'); hold on;
    plot(rrkf(:,1),rrkf(:,2),'r-');
    plot(rrkf(nbPts,1),rrkf(nbPts,2),'ks');
    axis([min(rrme(:,1)) max(rrme(:,1)) min(rrme(:,2)) max(rrme(:,2))]);
%     axis equal;
%     legend('KF', '3D-OD', 'final KF');
    title(['dY=f(dX), km']);
    subplot(4,3, 11);    
    plot(rrme(:,1),rrme(:,3),'gx:'); hold on;
    plot(rrkf(:,1),rrkf(:,3),'r-');
    plot(rrkf(nbPts,1),rrkf(nbPts,3),'ks');
    axis([min(rrme(:,1)) max(rrme(:,1)) min(rrme(:,3)) max(rrme(:,3))]);
%     axis equal; 
%     legend('KF', '3D-OD', 'final KF');
    title(['dZ=f(dX), km']);
    subplot(4,3, 12);
    plot(rrme(:,2),rrme(:,3),'gx:'); hold on;
    plot(rrkf(:,2),rrkf(:,3),'r-');
    plot(rrkf(nbPts,2),rrkf(nbPts,3),'ks');
    axis([min(rrme(:,2)) max(rrme(:,2)) min(rrme(:,3)) max(rrme(:,3))]);
%     axis equal; 
%     legend('KF', '3D-OD', 'final KF');
    title(['dZ=f(dY), km']);

    figure(103); clf;
    subplot(2,1,1);
	% 
	% Plutot visualiser les evolutions d'erreur en transveral et logitudinal p.r.vitesse
	% de plus, la position finale obtenue est √† une date passee, il faudrait comparer
    % la pr√©diction √† la date en cours avec l'attendu
	%
    semilogy(lKg, 'b-'); hold on;
    semilogy(ldP, 'b:');
    ylim([1e-5 1e5]);
    title(['Kalman filtering, ' sdt ', ' stf]);
    legend('Kalman gain, det(K''*K)', 'det(P)');
    
    subplot(2,1,2);
	% 
	% on dirait un delta-V et non un Velocity
	%
    semilogy(vkf, 'kx-'); hold on;
    legend('Velocity (km/s)');
  end
  nbCol=1;
  solKF(nC, nbCol:nbCol+2) = rex(nbPts+1,1:3);  nbCol=nbCol+3; %1-3 EXPECTED
  solKF(nC, nbCol:nbCol+2) = rrme(nbPts+1,1:3); nbCol=nbCol+3; %4-6: Áa n'a pas de sens!
  solKF(nC, nbCol:nbCol+2) = rrkf(nbPts+1,1:3); nbCol=nbCol+3; %7-9 KF solution
  solKF(nC, nbCol) = lKg(nbPts+1); nbCol=nbCol+1; %10
  solKF(nC, nbCol) = ldP(nbPts+1); nbCol=nbCol+1; %11
  solKF(nC, nbCol) = vkf(nbPts+1);  nbCol=nbCol+1; %12
  solKF(nC, nbCol) = dtrm(nbPts+1); nbCol=nbCol+1; %13: Áa n'a pas de sens!
  solKF(nC, nbCol) = dlgm(nbPts+1); nbCol=nbCol+1; %14: Áa n'a pas de sens!
  solKF(nC, nbCol) = dtrk(nbPts+1); nbCol=nbCol+1; %15
  solKF(nC, nbCol) = dlgk(nbPts+1); nbCol=nbCol+1; %16
  solKF(nC, nbCol) = mmkf(nbPts+1); nbCol=nbCol+1; %17
  solKF(nC, nbCol) = mmtk(nbPts+1); nbCol=nbCol+1; %18
  solKF(nC, nbCol) = mmlk(nbPts+1); nbCol=nbCol+1; %19
  solMd(1+(nC-1)*nbPts:nC*nbPts, 1) = rrme(1:nbPts,1);
  solMd(1+(nC-1)*nbPts:nC*nbPts, 2) = rrme(1:nbPts,2);
  solMd(1+(nC-1)*nbPts:nC*nbPts, 3) = rrme(1:nbPts,3);
  solMd(1+(nC-1)*nbPts:nC*nbPts, 4) = dtrm(1:nbPts);
  solMd(1+(nC-1)*nbPts:nC*nbPts, 5) = dlgm(1:nbPts);
end

% statistics for ntStep:
moyKF(ntStep,:) = mean(solKF);
stdKF(ntStep,:) =  std(solKF);
moyMd(ntStep,:) = mean(solMd);
stdMd(ntStep,:) =  std(solMd);

if (mod(ntStep,2)==0)
    figure(200);
%
errorbar(t2-t0, moyMd(ntStep,4), stdMd(ntStep,4), 'r');
errorbar(t2-t0, moyKF(ntStep,15), stdKF(ntStep,15), 'k');
xlabel('time (days)'); ylabel('Shift error (km)'); 
title(['Transversal errors, ' datestr(t0,'yyyy-mm-dd') '..' datestr(t2,'yyyy-mm-dd')]);
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
end 
figure(101);
%
subplot(1,2,1); hold on;
errorbar(t2-t0, moyMd(ntStep,4), stdMd(ntStep,4), 'r');
errorbar(t2-t0, moyKF(ntStep,15), stdKF(ntStep,15), 'k');
xlabel('time (days)'); ylabel('Shift error (km)'); 
title(['Transversal errors, ' datestr(t0,'yyyy-mm-dd') '..' datestr(t2,'yyyy-mm-dd')]);
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

subplot(1,2,2); hold on;
errorbar(t2-t0, moyMd(ntStep,5), stdMd(ntStep,5), 'r');
errorbar(t2-t0, moyKF(ntStep,16), stdKF(ntStep,16), 'k');
xlabel('time (days)'); ylabel('Shift error (km)');
title(['Longitudinal errors, ' datestr(t0,'yyyy-mm-dd') '..' datestr(t2,'yyyy-mm-dd')]);
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

figure(120);
%
subplot(3,1,1); hold on;
errorbar(t2-t0, moyMd(ntStep,1), stdMd(ntStep,1), 'r');
errorbar(t2-t0, moyKF(ntStep,7), stdKF(ntStep,7), 'k');
xlabel('time (days)'); ylabel('X shift error (km)'); 
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
%
subplot(3,1,2); hold on;
errorbar(t2-t0, moyMd(ntStep,2), stdMd(ntStep,2), 'r');
errorbar(t2-t0, moyKF(ntStep,8), stdKF(ntStep,8), 'k');
xlabel('time (days)'); ylabel('Y shift error (km)'); 
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
%
subplot(3,1,3); hold on;
errorbar(t2-t0, moyMd(ntStep,3), stdMd(ntStep,3), 'r');
errorbar(t2-t0, moyKF(ntStep,9), stdKF(ntStep,9), 'k');
xlabel('time (days)'); ylabel('Z shift error (km)'); 
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

end
fclose(fw);

figure(101);

subplot(1,2,1);
xlim([0 t2-t0]);
plot(ep(1:ntStep), moyMd(1:ntStep,4), '-r');
plot(ep(1:ntStep), moyKF(1:ntStep,15), '-k');

subplot(1,2,2);
xlim([0 t2-t0]);
plot(ep(1:ntStep), moyMd(1:ntStep,5), '-r');
plot(ep(1:ntStep), moyKF(1:ntStep,16), '-k');

figure(120);
subplot(3,1,1);
xlim([0 t2-t0]);
plot(ep(1:ntStep), moyMd(1:ntStep,1), '-r');
plot(ep(1:ntStep), moyKF(1:ntStep,7), '-k');
subplot(3,1,2);
xlim([0 t2-t0]);
plot(ep(1:ntStep), moyMd(1:ntStep,2), '-r');
plot(ep(1:ntStep), moyKF(1:ntStep,8), '-k');
subplot(3,1,3);
xlim([0 t2-t0]);
plot(ep(1:ntStep), moyMd(1:ntStep,3), '-r');
plot(ep(1:ntStep), moyKF(1:ntStep,9), '-k');

figure(101); subplot(1,2,1); xlim([0  200]); ylim([ -150 150]); 
figure(101); subplot(1,2,2); xlim([0  200]); ylim([ -300 300]); 
figure(120); subplot(3,1,1); xlim([0  200]); ylim([ -300 300]);
figure(120); subplot(3,1,2); xlim([0  200]); ylim([ -200 200]);
figure(120); subplot(3,1,3); xlim([0  200]); ylim([ -150 150]);
