%%----------------HEADER---------------------------%%
%Author:           Boris Segret
% Version & Date:
%                  V1.3, 26-05-2017
%                  - various plot of KF results after n iterations (n<nKF)
%                  V1.2, 21-05-2017
%                  - (bug) dtKF, in dec.days, not stored, then reconstructed
%                  V1.1, 25-04-2017
%                  V1.0   23-03-2017
%
% CL=2 (v1.3)
%
% I/
%    <result_file>
%    <option "graphs">
% O/
%    fig.100: longitudinal errors without and with KF
%    fig.101: longitudinal errors with raw, iKF, jKF, nKF iterations
%    fig.102: (with graph option) KF detailed behavior at each MC-cycle
%    fig.200: transversal errors without and with KF
%    fig.201: transversal errors with raw, iKF, jKF, nKF iterations
%    fig.120: X-, Y-, Z-axis errors without and with KF


clear;
% outputs_bin = '../cas_DAV/DAV_01as,50MCx8KF,t_bin';
outputs_bin = '../cas_EME/outs/Eb+YECMJE_01as,400MCx192KF';
% outputs_bin = '../cas_Y/outs/Y0_41324_01as,400MCx192KF,1393470_bin';
% outputs_bin = '../cas_Y/outs/YYv4_01as,400MCx192KF,1420098_bin';
% outputs_bin = '../ifod_tests/outs/Y0_41324v4_01as,400MCx192KF,t_bin'; % non-fonctionnel
graphs = true;
graphs = false;

T0N=datenum([2000 1 1 0 0 0]);
T0JD=2451544.5;

figure(200); clf; ylim([ -300 300]); xlim([0 200]); hold on; set(gca, 'LooseInset', [0 0 0 0]);
figure(100); clf; ylim([ -300 300]); xlim([0 200]); hold on; set(gca, 'LooseInset', [0 0 0 0]);
figure(201); clf; ylim([ -300 300]); xlim([0 200]); hold on; set(gca, 'LooseInset', [0 0 0 0]);
figure(101); clf; ylim([ -300 300]); xlim([0 200]); hold on; set(gca, 'LooseInset', [0 0 0 0]);

figure(120); clf;
subplot(3,1,1); ylim([ -300  300]); xlim([0  200]); hold on;
subplot(3,1,2); ylim([ -300  300]); xlim([0  200]); hold on;
subplot(3,1,3); ylim([ -300  300]); xlim([0  200]); hold on;
% subplot(3,1,1); ylim([ -1  1]); xlim([0  10]); hold on;
% subplot(3,1,2); ylim([ -1  1]); xlim([0  10]); hold on;
% subplot(3,1,3); ylim([ -1  1]); xlim([0  10]); hold on;

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
      Y=moyiKF;
      moyiKF=double(zeros(ntStep+100, nvKF));
      moyiKF(1:ntStep, 1:nvKF)=Y; clear Y;
      Y=moyjKF;
      moyjKF=double(zeros(ntStep+100, nvKF));
      moyjKF(1:ntStep, 1:nvKF)=Y; clear Y;
      Y=stdKF;
      stdKF=double(zeros(ntStep+100, nvKF));
      stdKF(1:ntStep, 1:nvKF)=Y; clear Y;
      Y=stdiKF;
      stdiKF=double(zeros(ntStep+100, nvKF));
      stdiKF(1:ntStep, 1:nvKF)=Y; clear Y;
      Y=stdjKF;
      stdjKF=double(zeros(ntStep+100, nvKF));
      stdjKF(1:ntStep, 1:nvKF)=Y; clear Y;
      Y=moyMd;
      moyMd=double(zeros(ntStep+100, 5));
      moyMd(1:ntStep, 1:5)=Y; clear Y;
      Y=stdMd;
      stdMd=double(zeros(ntStep+100, 5));
      stdMd(1:ntStep, 1:5)=Y; clear Y;
      Y=stats;
      stats=double(zeros(ntStep+100, nstats));
      stats(1:ntStep, 1:nstats)=Y; clear Y;
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

for nC=1:nbCycles
  % il faudrait cycler
  rawDATA = fread(fw, (3*3+10)*(nbPts+1), 'double');
  % will read the detailed value for 1 Monte-Carlo cycle,
  % to be reshaped into (nbPts+1) lines with:
  DATA = reshape(rawDATA, (3*3+10), nbPts+1)'; % (nbPts+1) lines x (3*3+10) columns
  
  nbCol=1;
  rex  = DATA(:,1:3); % Expected results
  rrme = DATA(:,4:6); % "Measured" results -rex (3D-OD)
  rrkf = DATA(:,7:9); % "Kalman Filtered" results -rex
%   lKg  = DATA(:,10);
%   ldP  = DATA(:,11);
%   vkf  = DATA(:,12);  % "Kalman Filtered" velocity
  dtrm = DATA(:,13)'; % transverse residual of "measured" results
  dlgm = DATA(:,14)'; % longitudin residual of "measured" results
  dtrk = DATA(:,15)'; % transverse residual of "K-Filter" results
  dlgk = DATA(:,16)'; % longitudin residual of "K-Filter" results
  mmkf = DATA(:,17)'; % max of the remaining norm residuals in KF-results
  mmtk = DATA(:,18)'; % max of the remaining transverse residuals in KF-results
  mmlk = DATA(:,19)'; % max of the remaining longitudin residuals in KF-results
  
  if (graphs)
    % chronogrammes
    figure(102); clf;
    subplot(3,3, [1 3]);
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
    subplot(3,3, [4 6]);
    plot(dtrk, 'r-'); hold on;
    plot(dtrm, 'gx:');
%     ylim([-2*std(dtrm) 2*std(dtrm)]);
    ylim([-5*mmtk(floor(nbPts*3/4)) 5*mmtk(floor(nbPts*3/4))]);
    plot(mmtk, 'r:'); plot(-mmtk, 'r:');
    title(['Residual in transversal shift (km) ' sdt ', ' stf]);
    subplot(3,3, [7 9]);
    plot(dlgk, 'r-'); hold on;
    plot(dlgm, 'gx:');
%     ylim([-2*std(dlgm) 2*std(dlgm)]);
    ylim([-5*mmlk(floor(nbPts*3/4)) 5*mmlk(floor(nbPts*3/4))]);
    plot(mmlk, 'r:'); plot(-mmlk, 'r:');
    title(['Residual in longitudinal shift (km) ' sdt ', ' stf]);

%     subplot(4,3, 10);
%     plot(rrme(:,1),rrme(:,2),'gx:'); hold on;
%     plot(rrkf(:,1),rrkf(:,2),'r-');
%     plot(rrkf(nbPts,1),rrkf(nbPts,2),'ks');
%     axis([min(rrme(:,1)) max(rrme(:,1)) min(rrme(:,2)) max(rrme(:,2))]);
% %     axis equal;
% %     legend('KF', '3D-OD', 'final KF');
%     title(['dY=f(dX), km']);
%     subplot(4,3, 11);    
%     plot(rrme(:,1),rrme(:,3),'gx:'); hold on;
%     plot(rrkf(:,1),rrkf(:,3),'r-');
%     plot(rrkf(nbPts,1),rrkf(nbPts,3),'ks');
%     axis([min(rrme(:,1)) max(rrme(:,1)) min(rrme(:,3)) max(rrme(:,3))]);
% %     axis equal; 
% %     legend('KF', '3D-OD', 'final KF');
%     title(['dZ=f(dX), km']);
%     subplot(4,3, 12);
%     plot(rrme(:,2),rrme(:,3),'gx:'); hold on;
%     plot(rrkf(:,2),rrkf(:,3),'r-');
%     plot(rrkf(nbPts,2),rrkf(nbPts,3),'ks');
%     axis([min(rrme(:,2)) max(rrme(:,2)) min(rrme(:,3)) max(rrme(:,3))]);
% %     axis equal; 
% %     legend('KF', '3D-OD', 'final KF');
%     title(['dZ=f(dY), km']);

%     figure(103); clf;
%     subplot(2,1,1);
% 	% 
% 	% Plutot visualiser les evolutions d'erreur en transveral et logitudinal p.r.vitesse
% 	% de plus, la position finale obtenue est à une date passee, il faudrait comparer
%     % la prédiction à la date en cours avec l'attendu
% 	%
%     semilogy(lKg, 'b-'); hold on;
%     semilogy(ldP, 'b:');
%     ylim([1e-5 1e5]);
%     title(['Kalman filtering, ' sdt ', ' stf]);
%     legend('Kalman gain, det(K''*K)', 'det(P)');
%     
%     subplot(2,1,2);
% 	% 
% 	% on dirait un delta-V et non un Velocity
% 	%
%     semilogy(vkf, 'kx-'); hold on;
%     legend('Velocity (km/s)');
  end
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
%   solMd(1+(nC-1)*nbPts:nC*nbPts, 1) = rrme(1:nbPts,1);
%   solMd(1+(nC-1)*nbPts:nC*nbPts, 2) = rrme(1:nbPts,2);
%   solMd(1+(nC-1)*nbPts:nC*nbPts, 3) = rrme(1:nbPts,3);
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

if (mod(ntStep,2)==0)
    figure(100);
    errorbar(t2-t0, moyMd(ntStep,5), stdMd(ntStep,5), 'r');
    errorbar(t2-t0, moyKF(ntStep,8), stdKF(ntStep,8), 'k');

    figure(101);
    errorbar(t2-t0, moyMd(ntStep,5), stdMd(ntStep,5), 'r');
    errorbar(t2-t0, moyiKF(ntStep,8), stdiKF(ntStep,8), 'g');
    errorbar(t2-t0, moyjKF(ntStep,8), stdjKF(ntStep,8), 'b');
    errorbar(t2-t0, moyKF(ntStep,8), stdKF(ntStep,8), 'k');

    figure(200);
    errorbar(t2-t0, moyMd(ntStep,4), stdMd(ntStep,4), 'r');
    errorbar(t2-t0, moyKF(ntStep,7), stdKF(ntStep,7), 'k');
    
    figure(201);
    errorbar(t2-t0, moyMd(ntStep,4), stdMd(ntStep,4), 'r');
    errorbar(t2-t0, moyiKF(ntStep,7), stdiKF(ntStep,7), 'g');
    errorbar(t2-t0, moyjKF(ntStep,7), stdjKF(ntStep,7), 'b');
    errorbar(t2-t0, moyKF(ntStep,7), stdKF(ntStep,7), 'k');
    
end 
figure(120);
%
subplot(3,1,1);
errorbar(t2-t0, moyMd(ntStep,1), stdMd(ntStep,1), 'r');
errorbar(t2-t0, moyKF(ntStep,4), stdKF(ntStep,4), 'k');
%
subplot(3,1,2);
errorbar(t2-t0, moyMd(ntStep,2), stdMd(ntStep,2), 'r');
errorbar(t2-t0, moyKF(ntStep,5), stdKF(ntStep,5), 'k');
%
subplot(3,1,3);
errorbar(t2-t0, moyMd(ntStep,3), stdMd(ntStep,3), 'r');
errorbar(t2-t0, moyKF(ntStep,6), stdKF(ntStep,6), 'k');

end
fclose(fw);

figure(100);
xlim([0 t2-t0]);
plot(ep(1:ntStep), moyMd(1:ntStep,5), '-r');
plot(ep(1:ntStep), moyKF(1:ntStep,8), '-k');
xlabel('time (days)'); ylabel('Shift error (km)');
title(['Longitudinal errors, ' datestr(t0,'yyyy-mm-dd') '..' datestr(t2,'yyyy-mm-dd')]);
legend('without KF', [num2str(nKF) '-step KF']);
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

figure(101);
xlim([0 t2-t0]);
plot(ep(1:ntStep), moyMd(1:ntStep,5), '-r');
plot(ep(1:ntStep), moyiKF(1:ntStep,8), '-g');
plot(ep(1:ntStep), moyjKF(1:ntStep,8), '-b');
plot(ep(1:ntStep), moyKF(1:ntStep,8), '-k');
xlabel('time (days)'); ylabel('Shift error (km)');
title(['Longitudinal errors, ' datestr(t0,'yyyy-mm-dd') '..' datestr(t2,'yyyy-mm-dd')]);
legend('without KF', [num2str(iKF) '-step KF'], [num2str(jKF) '-step KF'], [num2str(nKF) '-step KF']);
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

figure(200);
xlim([0 t2-t0]);
plot(ep(1:ntStep), moyMd(1:ntStep,4), '-r');
plot(ep(1:ntStep), moyKF(1:ntStep,7), '-k');
xlabel('time (days)'); ylabel('Shift error (km)'); 
title(['Transversal errors, ' datestr(t0,'yyyy-mm-dd') '..' datestr(t2,'yyyy-mm-dd')]);
legend('without KF', [num2str(nKF) '-step KF']);
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

figure(201);
xlim([0 t2-t0]);
plot(ep(1:ntStep), moyMd(1:ntStep,4), '-r');
plot(ep(1:ntStep), moyiKF(1:ntStep,7), '-g');
plot(ep(1:ntStep), moyjKF(1:ntStep,7), '-b');
plot(ep(1:ntStep), moyKF(1:ntStep,7), '-k');
xlabel('time (days)'); ylabel('Shift error (km)'); 
title(['Transversal errors, ' datestr(t0,'yyyy-mm-dd') '..' datestr(t2,'yyyy-mm-dd')]);
legend('without KF', [num2str(iKF) '-step KF'], [num2str(jKF) '-step KF'], [num2str(nKF) '-step KF']);
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

figure(120);
subplot(3,1,1);
xlim([0 t2-t0]);
plot(ep(1:ntStep), moyMd(1:ntStep,1), '-r');
plot(ep(1:ntStep), moyKF(1:ntStep,4), '-k');
xlabel('time (days)'); ylabel('X shift error (km)'); 
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
subplot(3,1,2);
xlim([0 t2-t0]);
plot(ep(1:ntStep), moyMd(1:ntStep,2), '-r');
plot(ep(1:ntStep), moyKF(1:ntStep,5), '-k');
xlabel('time (days)'); ylabel('Y shift error (km)'); 
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
subplot(3,1,3);
xlim([0 t2-t0]);
plot(ep(1:ntStep), moyMd(1:ntStep,3), '-r');
plot(ep(1:ntStep), moyKF(1:ntStep,6), '-k');
xlabel('time (days)'); ylabel('Z shift error (km)'); 
set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');

figure(200); ylim([-500 500]); figure(201); ylim([-500 500]);
figure(100); ylim([ -500 500]); figure(101); ylim([ -500 500]);
figure(120);
subplot(3,1,1); ylim([ -300  300]);
subplot(3,1,2); ylim([ -300  300]);
subplot(3,1,3); ylim([ -300  300]);
