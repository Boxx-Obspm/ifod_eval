%%----------------HEADER---------------------------%%
%Author:           Boris Segret
% Version & Date:
%                  V1.0, 21-05-2017, from stat_processing in v1.1 25-04-2017
%                  - deep analysis of point #3 in asynchronous triangulation
%                  PREVIOUSLY : stat_processing V1.2, 21-05-2017
%
% CL=1 (v1.0)
%
% This produces graphs to investigate the statistical results of the point
% #3 in the geometric problem inversion. A raw result file is considered
% and the point #3 solutions are converted into a dimensionless frame
% oriented toward the foreground body #3.
% 
% I/
%    <result_file> & <optical error>
%    <foreground object trajectory> & <actual trajectory>
%    <option "graphs">
% O/
%    <shift errors with error bars in transverse and longitudinal directions>
%    <shift errors with error bars in X, Y, Z, axis>
%    (with option "graphs") <shifts over the Kalman filter>


clear;
% actual_traj = '../cas_EME/Y/58122+SOI_v6.4_jdv010_312_vts.xyzv';
actual_traj = '../cas_EME/58122+SOI_v6.4_jdv000_312_vts.xyzv';
fg_body_Nb3 = '../cas_EME/Mars_imcce.xva';
% outputs_bin = '../cas_EME/outs/Eb+YECMJE_01as,400MCx192KF'; opterr=0.1;
outputs_bin = '../cas_EME/outs/E0_ECMJE_01as,450MCx192KF,1393471_bin'; opterr=0.1;
% outputs_bin = '../cas_Y/outs/Y0_41324_01as,400MCx192KF,1393470_bin';
% outputs_bin = '../cas_Y/outs/YYv4_01as,400MCx192KF,1420098_bin';
% outputs_bin = '../ifod_tests/outs/Y0_41324v4_01as,400MCx192KF,t_bin'; % non-fonctionnel
graphs = true;
% graphs = false;

T0N=datenum([2000 1 1 0 0 0]);
T0JD=2451544.5;
T0MJD=51544;

addpath('../ifod');
[NbLT1, TimeList1, data] = readTraj(actual_traj);
coord1 = data(:,1:3); vel1   = data(:,4:6);
clear data; fprintf('Reading act.traj''s data: %s, %.1f...%.1f = %i records\n', ...
  actual_traj, TimeList1(1), TimeList1(NbLT1), NbLT1);
[nbLgn, timeSteps, data] = readTraj(fg_body_Nb3);
coord3 = data(:,1:3); %vel3   = data(:,4:6);
clear data; fprintf('Body #3 trajecto''s data: %s, %.1f...%.1f = %i records\n', ...
  fg_body_Nb3, timeSteps(1), timeSteps(nbLgn), nbLgn);

% initizes the binary outputs with the number of written columns
fw = fopen(outputs_bin,'r');

lastData = false; ntStep = 0;
while not(lastData)
% fwrite(fw, [obstime(ik+nKF-2) obstime(ik+nKF)], 'double');
T = fread(fw, 2, 'double'); % will provide a 1x2 double array
%                             (time for KF-measured and for KF-predicted)
% t1=T(1)-T0MJD+T0N; % epoch of the last 3D geometric solution (3rd point)
% t2=T(2)-T0MJD+T0N; % epoch of the last KF-anticipated solution (5th point)
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
nKF=Z(2);
nbCycles=Z(3);
nbPts=nKF;
lastData = (Z(4)==1);

if (ntStep==0)
  t0=t2-nKF*Nobs*dtKF;
  nvKF =3;  % nb of stored values per KF-step
  ep = double(zeros(100,1));
  moyKF = double(zeros(100,nvKF));
  stdKF = double(zeros(100,nvKF));
  moyMd = double(zeros(100,3));
  stdMd = double(zeros(100,3));
  nstats =9; % nb of values in the statistical post-processing
  stats = double(zeros(100,nstats));
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
      moyMd=double(zeros(ntStep+100, 3));
      moyMd(1:ntStep, :)=Y; clear Y;
      Y=stdMd;
      stdMd=double(zeros(ntStep+100, 3));
      stdMd(1:ntStep, :)=Y; clear Y;
      Y=stats;
      stats=double(zeros(ntStep+100, nstats));
      stats(1:ntStep, 1:nstats)=Y; clear Y;
  end
  ntStep = ntStep+1;
end

ep(ntStep)=t2-t0;
stf=datestr(t2,'yyyy-mm-dd HH:MM:SS');
sdt=['~' num2str(24*Nobs*dtKF,'%4.1f') 'h/step'];

solKF = double(zeros(nbCycles, nvKF));
% solMd = double(zeros(nbCycles*nbPts, 5));
solMd = double(zeros(nbCycles*nbPts, 3));
method = 'linear';
timeB3 = double([T(2)-2*dtKF-(nKF-1)*(Nobs*dtKF):(Nobs*dtKF):T(2)-2*dtKF T(2)]); % nKF+1 values
vectB3 = interp1(timeSteps, coord3, timeB3) - interp1(TimeList1, coord1, timeB3);
%       obsXYZ(i,1:3) = interp1(ttimes, coords, obstime(i), method);
normB3 = sqrt(vectB3(:,1).^2 + vectB3(:,2).^2 + vectB3(:,3).^2);
scalef = (opterr/3600.).*(pi()/180.).*normB3;
xuniB3 = vectB3./repmat(normB3,[1 3]);
velTrj = interp1(TimeList1, vel1, timeB3);
zB3 = cross(xuniB3, velTrj);
zuniB3 = zB3./repmat(sqrt(zB3(:,1).^2 + zB3(:,2).^2 + zB3(:,3).^2), [1 3]);
yuniB3 = cross(zuniB3, xuniB3);
rotfromB3 = [xuniB3 yuniB3 zuniB3]; rottoB3 = double(zeros(nKF,9));
% fprintf('Computing rotation matrices to fg-Body linked fame, timeStep: %s\n', stf);
for i=1:nKF+1
  rotfrB3 = reshape(rotfromB3(i,:),3,3);
  rottoB3(i,:) = reshape(rotfrB3',1,9);
%   if (mod(i,10) ==0) fprintf('.'); end
%   if (mod(i,100)==0) fprintf(':%i\n',i); end
end
% fprintf(':%i\n',i);
clear vectB3 normB3 velTrj zB3 rotfromB3 rotfrB3; clear timeB3;

for nC=1:nbCycles
  % il faudrait cycler
  rawDATA = fread(fw, (3*3+10)*(nbPts+1), 'double');
  % will read the detailed value for 1 Monte-Carlo cycle,
  % to be reshaped into (nbPts+1) lines with:
  DATA = reshape(rawDATA, (3*3+10), nbPts+1)'; % (nbPts+1) lines x (3*3+10) columns
  
  nbCol=1;
  rex  = DATA(:,nbCol:nbCol+2); nbCol=nbCol+3; % Expected results
  rrme = DATA(:,nbCol:nbCol+2); nbCol=nbCol+3; % "Measured" results -rex (3D-OD)
  rrkf = DATA(:,nbCol:nbCol+2); nbCol=nbCol+3; % "Kalman Filtered" results -rex
%   lKg  = DATA(:,nbCol); nbCol=nbCol+1;
%   ldP  = DATA(:,nbCol); nbCol=nbCol+1;
%   vkf  = DATA(:,nbCol); nbCol=nbCol+1;  % "Kalman Filtered" velocity
%   dtrm = DATA(:,nbCol)'; nbCol=nbCol+1; % transverse residual of "measured" results
%   dlgm = DATA(:,nbCol)'; nbCol=nbCol+1; % longitudin residual of "measured" results
%   dtrk = DATA(:,nbCol)'; nbCol=nbCol+1; % transverse residual of "K-Filter" results
%   dlgk = DATA(:,nbCol)'; nbCol=nbCol+1; % longitudin residual of "K-Filter" results
%   mmkf = DATA(:,nbCol)'; nbCol=nbCol+1; % max of the remaining norm residuals in KF-results
%   mmtk = DATA(:,nbCol)'; nbCol=nbCol+1; % max of the remaining transverse residuals in KF-results
%   mmlk = DATA(:,nbCol)'; nbCol=nbCol+1; % max of the remaining longitudin residuals in KF-results
  
%   nbCol=1;
%   solKF(nC, nbCol:nbCol+2) = rex(nbPts+1,1:3);  nbCol=nbCol+3; %1-3 EXPECTED
%   solKF(nC, nbCol:nbCol+2) = rrme(nbPts+1,1:3); nbCol=nbCol+3; %4-6: ça n'a pas de sens!
%   solKF(nC, nbCol:nbCol+2) = rrkf(nbPts+1,1:3); nbCol=nbCol+3; %7-9 KF solution
%   solKF(nC, nbCol) = lKg(nbPts+1); nbCol=nbCol+1; %10
%   solKF(nC, nbCol) = ldP(nbPts+1); nbCol=nbCol+1; %11
%   solKF(nC, nbCol) = vkf(nbPts+1);  nbCol=nbCol+1; %12
%   solKF(nC, nbCol) = dtrm(nbPts+1); nbCol=nbCol+1; %13: ça n'a pas de sens!
%   solKF(nC, nbCol) = dlgm(nbPts+1); nbCol=nbCol+1; %14: ça n'a pas de sens!
%   solKF(nC, nbCol) = dtrk(nbPts+1); nbCol=nbCol+1; %15
%   solKF(nC, nbCol) = dlgk(nbPts+1); nbCol=nbCol+1; %16
%   solKF(nC, nbCol) = mmkf(nbPts+1); nbCol=nbCol+1; %17
%   solKF(nC, nbCol) = mmtk(nbPts+1); nbCol=nbCol+1; %18
%   solKF(nC, nbCol) = mmlk(nbPts+1); nbCol=nbCol+1; %19

%   fprintf('rrme in fg-Body frame, Cycle %i:', nC);
  for i=1:nKF
    rotB3 = reshape(rottoB3(i,:),3,3);
    if rank(rotB3)==3 
      solMd((nC-1)*nKF+i, :) = (rotB3*rrme(i,:)')'./scalef(i);
%   solMd(1+(nC-1)*nbPts:nC*nbPts, 1) = rrme(1:nbPts,1);
%   solMd(1+(nC-1)*nbPts:nC*nbPts, 2) = rrme(1:nbPts,2);
%   solMd(1+(nC-1)*nbPts:nC*nbPts, 3) = rrme(1:nbPts,3);
%       if (mod(i,10) ==0) fprintf('.'); end
%       if (mod(i,100)==0) fprintf('/'); end
    else
      fprintf('0');
    end
  end
  nbCol=1;
  solKF(nC, nbCol:nbCol+2) = (rotB3*rrkf(nKF+1,:)')'./scalef(nKF+1); nbCol=nbCol+3; %7-9 KF solution

end % loop with nC on nbCycles

% statistics for ntStep:
moyKF(ntStep,:) = mean(solKF.*scalef(nKF+1));
stdKF(ntStep,:) =  std(solKF.*scalef(nKF+1));
moyMd(ntStep,:) = mean(solMd);
stdMd(ntStep,:) =  std(solMd);
% stats(:,1) = Corr(X_Body, Y)
% stats(:,2) = Corr(X_Body, Z)
% stats(:,3) = Corr(Y, Z)
% stats(:,4) = min(sigma_X*scalef(iKF))
% stats(:,5) = min(sigma_Y*scalef(iKF))
% stats(:,6) = min(sigma_Z*scalef(iKF))
% stats(:,7) = max(sigma_X*scalef(iKF))
% stats(:,8) = max(sigma_Y*scalef(iKF))
% stats(:,9) = max(sigma_Z*scalef(iKF))
stats(ntStep,1:3) = [corr(solMd(:,1),solMd(:,2)) corr(solMd(:,1),solMd(:,3)) corr(solMd(:,2),solMd(:,3))];
stats(ntStep,4) = min(scalef.*stdMd(ntStep,1));
stats(ntStep,5) = min(scalef.*stdMd(ntStep,2));
stats(ntStep,6) = min(scalef.*stdMd(ntStep,3));
stats(ntStep,7) = max(scalef.*stdMd(ntStep,1));
stats(ntStep,8) = max(scalef.*stdMd(ntStep,2));
stats(ntStep,9) = max(scalef.*stdMd(ntStep,3));

if (graphs)
    % dfltALI = get(0, 'DefaultAxesLooseInset'); % [0.1300    0.1100    0.0950    0.0750]
    % set(gca, 'LooseInset', get(gca,'TightInset'));
    % set(0, 'DefaultAxesLooseInset', [0 0 0 0]);
    figure(301); % dfltTI = get(gca, 'TightInset');
    mh = max([max(abs(solMd(:,1))) max(abs(solMd(:,3)))]);
    mv = max(abs(solMd(:,2)));
    % subplot(3,6,[1 9]); hold off; 
    % zz=(solMd(:,1)<0 & solMd(:,2)<0 & solMd(:,3)>0);
    % plot3(solMd(zi,1), solMd(zi,2), solMd(zi,3),'.b'); hold on;
    % zz=(solMd(:,1)<0);
    % zi=find(zz>0);
    % plot3(zz.*solMd(:,1), zz.*solMd(:,2), zz.*solMd(:,3),'.b'); hold on;

    %     subplot(3,2,2); hold off; plot([1 1 100], [10 0 0], '-'); axis off;
    %     text( 5,9,['Step:' num2str(ntStep,'%4i') ', date: ' stf]);
    %     text(15,7,['sigma_Z = ' num2str(stdMd(ntStep,3),'%.3f')]);
    %     text(15,6,['sigma_Y = ' num2str(stdMd(ntStep,2),'%.3f')]);
    %     text(50,6,['sigma_X = ' num2str(stdMd(ntStep,1),'%.3f')]);
    %     text( 5,4,['Scale: unit = ' num2str(scalef(1),'%.1f') ' ... ' num2str(scalef(nKF), '%.1f km')]);
    subplot(10,10, [6 50]); hold off; 
    % subplot(3,2,2); hold off; 
    plot(solMd(:,1), solMd(:,2),'.g'); hold on; axis equal;
    plot(solKF(:,1), solKF(:,2),'.r');
    plot([moyMd(ntStep,1)-4*stdMd(ntStep,1) moyMd(ntStep,1)+4*stdMd(ntStep,1)], [moyMd(ntStep,2) moyMd(ntStep,2)],'-b', 'LineWidth', 4);
    plot([moyMd(ntStep,1) moyMd(ntStep,1)], [moyMd(ntStep,2)-4*stdMd(ntStep,2) moyMd(ntStep,2)+4*stdMd(ntStep,2)],'-k', 'LineWidth', 4);
%     plot(moyMd(ntStep,1),moyMd(ntStep,2), 'o', 'LineWidth',2, 'MarkerSize',12, 'MarkerEdgeColor','k');
    axis([-4.5*stdMd(ntStep,1) 4.5*stdMd(ntStep,1) -4.5*stdMd(ntStep,2) 4.5*stdMd(ntStep,2)]);
    text(-4*stdMd(ntStep,1), -4*stdMd(ntStep,2), 'LoS x LoS_Y');

    subplot(10,10, [51 95]); hold off; 
    % subplot(3,2,3); hold off; 
    plot(solMd(:,2), solMd(:,3),'.g'); hold on; axis equal;
    plot(solKF(:,2), solKF(:,3),'.r');
    plot([moyMd(ntStep,2)-4*stdMd(ntStep,2) moyMd(ntStep,2)+4*stdMd(ntStep,2)], [moyMd(ntStep,3) moyMd(ntStep,3)],'-k', 'LineWidth', 4);
    plot([moyMd(ntStep,2) moyMd(ntStep,2)], [moyMd(ntStep,3)-4*stdMd(ntStep,3) moyMd(ntStep,3)+4*stdMd(ntStep,3)],'-k', 'LineWidth', 4);
    plot(moyMd(ntStep,2),moyMd(ntStep,3), 'o', 'LineWidth',1, 'MarkerSize',12, 'MarkerFaceColor','b', 'MarkerEdgeColor','k');
    axis([-4.5*stdMd(ntStep,1) 4.5*stdMd(ntStep,1) -4.5*stdMd(ntStep,3) 4.5*stdMd(ntStep,3)]);
%   axis([-mh mh -mv mv]);
    text(-4*stdMd(ntStep,1), -4*stdMd(ntStep,3), 'Perp.to LoS: LoS_Y x LoS_Z');
    
    subplot(10,10, [56 100]); hold off; 
    % subplot(3,2,4); hold off; 
    plot(solMd(:,1), solMd(:,3),'.g'); hold on; axis equal;
    plot(solKF(:,1), solKF(:,3),'.r');
    plot([moyMd(ntStep,1)-4*stdMd(ntStep,1) moyMd(ntStep,1)+4*stdMd(ntStep,1)], [moyMd(ntStep,3) moyMd(ntStep,3)],'-b', 'LineWidth', 4);
    plot([moyMd(ntStep,1) moyMd(ntStep,1)], [moyMd(ntStep,3)-4*stdMd(ntStep,3) moyMd(ntStep,3)+4*stdMd(ntStep,3)],'-k', 'LineWidth', 4);
%     plot(moyMd(ntStep,1),moyMd(ntStep,3), 'o', 'LineWidth',2, 'MarkerSize',12, 'MarkerEdgeColor','k');
    axis([-4.5*stdMd(ntStep,1) 4.5*stdMd(ntStep,1) -4.5*stdMd(ntStep,3) 4.5*stdMd(ntStep,3)]);
    text(-4*stdMd(ntStep,1), -4*stdMd(ntStep,3), 'LoS x LoS_Z');

    nbins=100;
    cc=1;
    subplot(10,10, [34 45]); hold off; 
    % subplot(3,2,6); hold off; 
    hist(solMd(:,cc), nbins); h = findobj(gca,'Type','patch'); set(h,'FaceColor','g','EdgeColor','g');
    hold on;
    plot([moyMd(ntStep,cc)-stdMd(ntStep,cc) moyMd(ntStep,cc)-stdMd(ntStep,cc)], [0 nbCycles*nKF/nbins], '-b', 'LineWidth', 3);
    plot([moyMd(ntStep,cc)+stdMd(ntStep,cc) moyMd(ntStep,cc)+stdMd(ntStep,cc)], [0 nbCycles*nKF/nbins], '-b', 'LineWidth', 3);
    xlim([-3*stdMd(ntStep,cc) 3*stdMd(ntStep,cc)]);
    set(gca,'xtick', 0.01*floor(100.*[moyMd(ntStep,cc)-stdMd(ntStep,cc) moyMd(ntStep,cc)+stdMd(ntStep,cc)]), 'ytick', []); ylabel('LoS');
    cc=2;
    subplot(10,10, [14 25]); hold off; 
    % subplot(3,2,5); hold off; 
    hist(solMd(:,cc), nbins); h = findobj(gca,'Type','patch'); set(h,'FaceColor','g','EdgeColor','g');
    hold on;
    plot([moyMd(ntStep,cc)-stdMd(ntStep,cc) moyMd(ntStep,cc)-stdMd(ntStep,cc)], [0 nbCycles*nKF/nbins], '-b', 'LineWidth', 3);
    plot([moyMd(ntStep,cc)+stdMd(ntStep,cc) moyMd(ntStep,cc)+stdMd(ntStep,cc)], [0 nbCycles*nKF/nbins], '-b', 'LineWidth', 3);
    xlim([-3*stdMd(ntStep,cc) 3*stdMd(ntStep,cc)]);
    set(gca,'xtick', 0.01*floor(100.*[moyMd(ntStep,cc)-stdMd(ntStep,cc) moyMd(ntStep,cc)+stdMd(ntStep,cc)]), 'ytick', []); ylabel('LoS_Y');
    cc=3;
    subplot(10,10, [32 43]); hold off; 
    % subplot(3,2,1); hold off; 
    hist(solMd(:,cc), nbins); h = findobj(gca,'Type','patch'); set(h,'FaceColor','g','EdgeColor','g');
    hold on;
    plot([moyMd(ntStep,cc)-stdMd(ntStep,cc) moyMd(ntStep,cc)-stdMd(ntStep,cc)], [0 nbCycles*nKF/nbins], '-b', 'LineWidth', 3);
    plot([moyMd(ntStep,cc)+stdMd(ntStep,cc) moyMd(ntStep,cc)+stdMd(ntStep,cc)], [0 nbCycles*nKF/nbins], '-b', 'LineWidth', 3);
    xlim([-3*stdMd(ntStep,cc) 3*stdMd(ntStep,cc)]);
    set(gca,'xtick', 0.01*floor(100.*[moyMd(ntStep,cc)-stdMd(ntStep,cc) moyMd(ntStep,cc)+stdMd(ntStep,cc)]), 'ytick', []); ylabel('LoS_Z');

    % set(gca, 'XColor', 'm'); set(gca, 'YColor', 'm');
else
    fprintf(['ntStep %i: corr. %6.3f(X,Y) %6.3f(X,Z) %3.6f(Y,Z), ' ...
        'min-max(km) %.1f..%.1f(X) %.1f..%.1f(Y) %.1f..%.1f(Z)\n'], ...
        ntStep, stats(ntStep,[1:3 4 7 5 8 6 9]));
end


figure(302); clf;
% subplot(2,1,1);
semilogy(ep(1:ntStep),stats(1:ntStep,7),'--b');
hold on;
set(gca, 'LooseInset', [0 0 0 0]); % semilogy(stats(1:ntStep,4),'--r');
semilogy(ep(1:ntStep),stats(1:ntStep,8),'--k'); % semilogy(stats(1:ntStep,5),'--k'); 
% semilogy(ep(1:ntStep),stats(1:ntStep,9),'--k'); % semilogy(stats(1:ntStep,6),'--k');
semilogy(ep(1:ntStep),stdKF(1:ntStep,1),'-r'); % semilogy(stats(1:ntStep,6),'--k');
semilogy(ep(1:ntStep),stdKF(1:ntStep,2),'-k'); % semilogy(stats(1:ntStep,6),'--k');
% semilogy(ep(1:ntStep),stdKF(1:ntStep,3),'-k'); % semilogy(stats(1:ntStep,6),'--k');
legend('LoS', 'LoS_y', 'LoS w.KF', 'LoS_y w.KF');
% legend('LoS', 'LoS_y', 'LoS_Z', 'LoS w.KF', 'LoS_y w.KF', 'LoS_Z w.KF');
ylabel(['Standard deviation (km)']); xlabel('days');
% hold off;

figure(303); clf;
% subplot(2,1,2);
plot(ep(1:ntStep),stats(1:ntStep,1),'-b');
hold on;
set(gca, 'LooseInset', [0 0 0 0]);
plot(ep(1:ntStep),stats(1:ntStep,2),'--b');  
plot(ep(1:ntStep),stats(1:ntStep,3),'-k');
legend('(LoS,LoS_Y)', '(LoS,LoS_Z)', '(LoS_Y,LoS_Z)');
% legend('LoS', 'LoS_y', 'LoS_Z', 'LoS w.KF', 'LoS_y w.KF', 'LoS_Z w.KF');
ylabel(['Correlation factors']); xlabel('days');
% hold off;

end; % loop with ntStep (while not lastData)
fclose(fw);
