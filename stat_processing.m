clear;
% outputs_bin = '../cas_EME/outs/E2M-58122Y_E-1,050h,a-h';
% outputs_bin = '../cas_Y/outs/Y+Y_out_E-1asec,MC0-8';
outputs_bin = '../ifod_tests/outs/test_E-1,300x50h,tests_bin';
graphs = true;
% graphs = false;

T0N=datenum([2000 1 1 0 0 0]);
T0JD=2451544.5;

figure(101); clf;
subplot(2,1,1); ylim([ -600  600]); hold on;
subplot(2,1,2); ylim([-1000 1000]); hold on;

% initizes the binary outputs with the number of written columns
fw = fopen(outputs_bin,'r');

lastData = false; ntStep = 0;
while not(lastData)
% fwrite(fw, [obstime(ik+nKF-2) obstime(ik+nKF)], 'double');
T = fread(fw, 2, 'double'); % will provide a 1x2 double array
%                             (time for KF-measured and for KF-predicted)
t1=T(1)-T0JD+T0N;
t2=T(2)-T0JD+T0N;

if (ntStep==0)
  t0=t2;
  nstats=1; % nb of values in the statistical post-processing
  nvKF=19;  % nb of stored values per KF-step
  moyKF=double(zeros(100,nvKF));
  stdKF=double(zeros(100,nvKF));
  moyMd=double(zeros(100,2));
  stdMd=double(zeros(100,2));
  stats=double(zeros(100,nstats));
  ntStep=1;
else
  if (mod(ntStep,100)==0)
      % memory re-allocation every 100 steps
      Y=moyKF;
      moyKF=double(zeros(ntStep+100, nvKF));
      moyKF(1:ntStep, 1:nvKF)=Y; clear Y;
      Y=stdKF;
      stdKF=double(zeros(ntStep+100, nvKF));
      stdKF(1:ntStep, 1:nvKF)=Y; clear Y;
      Y=moyMd;
      moyMd=double(zeros(ntStep+100, 2));
      moyMd(1:ntStep, 1:2)=Y; clear Y;
      Y=stdMd;
      stdMd=double(zeros(ntStep+100, 2));
      stdMd(1:ntStep, 1:2)=Y; clear Y;
      Y=stats;
      stats=double(zeros(ntStep+100, nstats));
      stats(1:ntStep, 1:nstats)=Y; clear Y;
  end
  ntStep = ntStep+1;
end

stf=datestr(t2,'yyyy-mm-dd HH:MM:SS');
dtKF=24*(t2-t1)/2;
sdt=['~' num2str(dtKF,'%4.1f') 'h/step'];
% fwrite(fw, [Nobs nKF nbCycles (ik+nKF==length(obstime))], 'uint32');
Z = fread(fw, 4, 'uint32'); % will provide a 1x3 uint32 array
% fwrite(fw, ...
%     [rex rrme rrkf lKg ldP vkf dtrm' dlgm' dtrk' dlgk' mmkf' mmtk' mmlk'], ...
%     'double');
nbPts=Z(2)-Z(1)+1; % (nbPts=nKF-Nobs+1;)
Nobs=Z(1);
nKF=Z(2);
nbCycles=Z(3);
lastData = (Z(4)==1);

solKF = double(zeros(nbCycles, nvKF));
dtlKF = double(zeros(nbCycles*nbPts, 2));

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
	% de plus, la position finale obtenue est à une date passee, il faudrait comparer
    % la prédiction à la date en cours avec l'attendu
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
  solKF(nC, nbCol:nbCol+2) = rex(nbPts+1,1:3);  nbCol=nbCol+3; %1-3
  solKF(nC, nbCol:nbCol+2) = rrme(nbPts+1,1:3); nbCol=nbCol+3; %4-6
  solKF(nC, nbCol:nbCol+2) = rrkf(nbPts+1,1:3); nbCol=nbCol+3; %7-9
  solKF(nC, nbCol) = lKg(nbPts+1); nbCol=nbCol+1; %10
  solKF(nC, nbCol) = ldP(nbPts+1); nbCol=nbCol+1; %11
  solKF(nC, nbCol) = vkf(nbPts+1);  nbCol=nbCol+1; %12
  solKF(nC, nbCol) = dtrm(nbPts+1); nbCol=nbCol+1; %13
  solKF(nC, nbCol) = dlgm(nbPts+1); nbCol=nbCol+1; %14
  solKF(nC, nbCol) = dtrk(nbPts+1); nbCol=nbCol+1; %15
  solKF(nC, nbCol) = dlgk(nbPts+1); nbCol=nbCol+1; %16
  solKF(nC, nbCol) = mmkf(nbPts+1); nbCol=nbCol+1; %17
  solKF(nC, nbCol) = mmtk(nbPts+1); nbCol=nbCol+1; %18
  solKF(nC, nbCol) = mmlk(nbPts+1); nbCol=nbCol+1; %19
  dtlKF(1+(nC-1)*nbPts:nC*nbPts, 1) = dtrm(1:nbPts);
  dtlKF(1+(nC-1)*nbPts:nC*nbPts, 2) = dlgm(1:nbPts);
end

% statistics for ntStep:
moyKF(ntStep,:) = mean(solKF);
stdKF(ntStep,:) =  std(solKF);
moyMd(ntStep,:) = mean(dtlKF);
stdMd(ntStep,:) =  std(dtlKF);

figure(101);
%
subplot(2,1,1); hold on;
% plot(t2-t0, moyKF(ntStep,15), 'or');
errorbar(t2-t0, moyMd(ntStep,1), stdMd(ntStep,1), 'r');
errorbar(t2-t0, moyKF(ntStep,15), stdKF(ntStep,15), 'k');
xlabel('time (days)'); ylabel('Shift error (km)'); 
title(['Transversal errors, ' datestr(t0,'yyyy-mm-dd') '..' datestr(t2,'yyyy-mm-dd')]);
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
%
subplot(2,1,2); hold on;
% plot(t2-t0, moyKF(ntStep,16), 'or');
errorbar(t2-t0, moyMd(ntStep,2), stdMd(ntStep,2), 'r');
errorbar(t2-t0, moyKF(ntStep,16), stdKF(ntStep,16), 'k');
xlabel('time (days)'); ylabel('Shift error (km)');
title(['Longitudinal errors, ' datestr(t0,'yyyy-mm-dd') '..' datestr(t2,'yyyy-mm-dd')]);
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
end
subplot(2,1,1); xlim([0 t2-t0]);
subplot(2,1,2); xlim([0 t2-t0]);
fclose(fw);
