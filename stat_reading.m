outputs_o  = '../cas_EME/outs/E0_ECMJE_001as,450MCx192KF,1393475_bin';
% outputs_o  = '../cas_Y/outs/Y+Y41324_E-1,400x200h,P41324(101+),1355334_bin';
% outputs_o  = '../ifod_tests/outs/test_E-1,4x50h,tests_bin';
graphs = true;
graphs = false;

% --------------------------
% Statistics after stacking:
% --------------------------
% For each time-step, the mean and the standard deviation of each
% parameters are computed, plotted and compared with the mean and standard
% deviation of N random extractions of the half of the values of each
% parameter. The dimensionless shifts are plotted.

T0N=datenum([2000 1 1 0 0 0]);
T0JD=2451544.5;
nvKF=19;  % nb of stored values per KF-step
nix = 10;  % nb of extraction for sub-statistitcs

fo = fopen(outputs_o,'r');
fprintf('Statistics: %s\n', outputs_o);
lastData = false;
ntStep = 0;
while not(lastData)
 % fwrite(fw, [obstime(ik+nKF-2) obstime(ik+nKF)], 'double');
 T = fread(fo, 2, 'double'); % will provide a 1x2 double array
 %                             (time for KF-measured and for KF-predicted)

 if (ntStep==0)
   moyKF=double(zeros(100,nvKF));
   stdKF=double(zeros(100,nvKF));
%    moyMD=double(zeros(100,nvKF));
%    stdMD=double(zeros(100,nvKF));
   moyKFi=double(zeros(100,nvKF,nix));
   stdKFi=double(zeros(100,nvKF,nix));
   stda = double(zeros(100,1)); moya = double(zeros(100,1));
   stdb = double(zeros(100,1)); moyb = double(zeros(100,1));
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
%        Y=moyMD;
%        moyMD=double(zeros(ntStep+100, nvKF));
%        moyMD(1:ntStep, 1:nvKF)=Y; clear Y;
%        Y=stdMD;
%        stdMD=double(zeros(ntStep+100, nvKF));
%        stdMD(1:ntStep, 1:nvKF)=Y; clear Y;
       Y=moyKFi;
       moyKFi=double(zeros(ntStep+100, nvKF, nix));
       moyKFi(1:ntStep, 1:nvKF, 1:nix)=Y; clear Y;
       Y=stdKFi;
       stdKFi=double(zeros(ntStep+100, nvKF, nix));
       stdKFi(1:ntStep, 1:nvKF, 1:nix)=Y; clear Y;
       Y=stda;
       stda=double(zeros(ntStep+100, 1));
       stda(1:ntStep)=Y; clear Y;
       Y=moya;
       moya=double(zeros(ntStep+100, 1));
       moya(1:ntStep)=Y; clear Y;
       Y=stdb;
       stdb=double(zeros(ntStep+100, 1));
       stdb(1:ntStep)=Y; clear Y;
       Y=moyb;
       moyb=double(zeros(ntStep+100, 1));
       moyb(1:ntStep)=Y; clear Y;
   end
   ntStep = ntStep+1;
 end
 
% fwrite(fw, [Nobs nKF nbCycles (ik+nKF==length(obstime))], 'uint32');
Z = fread(fo, 4, 'uint32'); % will provide a 1x3 uint32 array
Nobs=Z(1);
nKF=Z(2);
nbPts=nKF;
nbCycles=Z(3);
lastData = (Z(4)==1);

% DATAmd = double(zeros(nbPts*nbCycles, nvKF));
DATAkf = double(zeros(nbCycles, nvKF));

  for nC=1:nbCycles
    % fwrite(fw, ...
    %     [rex rrme rrkf lKg ldP vkf dtrm' dlgm' dtrk' dlgk' mmkf' mmtk' mmlk'], ...
    %     'double');
    rawDATA = fread(fo, nvKF*(nbPts+1), 'double');
    fmtDATA = reshape(rawDATA, nvKF, nbPts+1)'; % (nbPts+1) lines x (3*3+10) columns
%     DATAmd(1+(nC-1)*nbPts:nC*nbPts, 1:nvKF) = fmtDATA(1:nbPts, 1:nvKF);
    DATAkf(nC, 1:nvKF) = fmtDATA(nbPts+1, 1:nvKF);
  end
 %   nbCol=1;
 %   solKF(nC, nbCol:nbCol+2) = rex(nbPts+1,1:3);  nbCol=nbCol+3; %1-3
 %   solKF(nC, nbCol:nbCol+2) = rrme(nbPts+1,1:3); nbCol=nbCol+3; %4-6
 %   solKF(nC, nbCol:nbCol+2) = rrkf(nbPts+1,1:3); nbCol=nbCol+3; %7-9
 %   solKF(nC, nbCol) = lKg(nbPts+1); nbCol=nbCol+1; %10
 %   solKF(nC, nbCol) = ldP(nbPts+1); nbCol=nbCol+1; %11
 %   solKF(nC, nbCol) = vkf(nbPts+1);  nbCol=nbCol+1; %12
 %   solKF(nC, nbCol) = dtrm(nbPts+1); nbCol=nbCol+1; %13
 %   solKF(nC, nbCol) = dlgm(nbPts+1); nbCol=nbCol+1; %14
 %   solKF(nC, nbCol) = dtrk(nbPts+1); nbCol=nbCol+1; %15
 %   solKF(nC, nbCol) = dlgk(nbPts+1); nbCol=nbCol+1; %16
 %   solKF(nC, nbCol) = mmkf(nbPts+1); nbCol=nbCol+1; %17
 %   solKF(nC, nbCol) = mmtk(nbPts+1); nbCol=nbCol+1; %18
 %   solKF(nC, nbCol) = mmlk(nbPts+1); nbCol=nbCol+1; %19
 %   dtlKF(1+(nC-1)*nbPts:nC*nbPts, 1) = dtrm(1:nbPts);
 %   dtlKF(1+(nC-1)*nbPts:nC*nbPts, 2) = dlgm(1:nbPts);  

 % statistics for ntStep:
 moyKF(ntStep,:) = mean(DATAkf);
 stdKF(ntStep,:) =  std(DATAkf);
%  moyMD(ntStep,:) = mean(DATAmd);
%  stdMD(ntStep,:) =  std(DATAmd);
 ix = floor(unifrnd(1,nbCycles,nix,floor(nbCycles/2))); % selection of half of the indices
 for iix=1:nix
   % stats on the the extracted sample iix
   moyKFi(ntStep, 4:nvKF, iix) = mean(DATAkf(ix(iix,:),4:nvKF));
   stdKFi(ntStep, 4:nvKF, iix) =  std(DATAkf(ix(iix,:),4:nvKF));
 end
 if (graphs) figure(104); clf; hold on; end;
%  iKF=11;
 plotXi = double(zeros(1,5*nix));
 plotYi = double(zeros(1,5*nix));
 iiKF=0;
 for iKF=[7 8 9 15 16]
   plotX = (reshape(moyKFi(ntStep,iKF,:),1,nix) -moyKF(ntStep,iKF))./ stdKF(ntStep,iKF);
   plotY = (reshape(stdKFi(ntStep,iKF,:),1,nix) -stdKF(ntStep,iKF))./ stdKF(ntStep,iKF);
   if (graphs) plot(plotX, plotY, 'xb'); end
   plotXi(1+iiKF*nix:(iiKF+1)*nix) = plotX;
   plotYi(1+iiKF*nix:(iiKF+1)*nix) = plotY;
   iiKF=iiKF+1;
 end
 moya(ntStep) = mean(plotXi); moyb(ntStep) = mean(plotYi);
 stda(ntStep) =  std(plotXi); stdb(ntStep) =  std(plotYi);
 if (graphs)
   plot(0,0,'sr');
   plot(moya(ntStep), moyb(ntStep), 'xr');
   plot(moya(ntStep)+3*stda(ntStep)*cos(2.*pi()*[0:100]/100.), ...
       moyb(ntStep)+3*stdb(ntStep)*sin(2.*pi()*[0:100]/100.), '-r');
 end
 if (mod(ntStep,10)==0) fprintf('/'); else fprintf('.'); end
 if (mod(ntStep,100)==0) fprintf(':%i\n',ntStep); end
end
fprintf('\n');
fclose(fo);
% a = mean(moyKF(:, [7 8 9 15 16])); %amax = max(moyKF(:, [7 8 9 15 16]));
% b = mean(stdKF(:, [7 8 9 15 16])); %bmax = max(stdKF(:, [7 8 9 15 16]));
% c = 100*mean(stdb); % accuracy of standard deviations stdKF

fprintf('Statistics:\n- Nb.of time-steps: %i\n', ntStep);
% fprintf('Average results:\n');
% fprintf('==>  %f +/- %f km\n', [a' b']); fprintf('     @ %3.1f%% accuracy on sigma\n', c);
fprintf('- Nb.of Cycles per time-step: %i\n', nbCycles);
fprintf('- Nb.of Kalman Filter steps per Cycle: %i\n', nKF);
fprintf('- Nb.of extracted samples (half the syze of a Cycle, i.e.%i): %i\n', nbCycles/2, nix);
% fprintf('- Maximum shift of mean of samples wrt total: %f km\n', 0.);

figure(102); clf;

subplot(2,3,[1 2]); hold on; xlim([0 ntStep]); ylim([-0.2 0.2]);
errorbar(moya(1:ntStep), stda(1:ntStep),'-b');
plot([0 ntStep],   [1/30. 1/30.], '-r', 'LineWidth', 1.5);
plot([0 ntStep], [-1/30. -1/30.], '-r', 'LineWidth', 1.5);
legend('(Mean.i-Mean)/SIGMA = f(t)', 'Convergence targets'); %title('(Mean.i-Mean)/SIGMA');
% plot(moya(1:ntStep),'-k');
plot(moya(1:ntStep)+0.1,'-r', 'LineWidth', 2.5);
plot(moya(1:ntStep)-0.1,'-r', 'LineWidth', 2.5);

subplot(2,3,[4 5]); hold on; xlim([0 ntStep]); ylim([-0.2 0.2]);
errorbar(moyb(1:ntStep), stdb(1:ntStep),'-b');
plot([0 ntStep],   [1/30. 1/30.], '-r', 'LineWidth', 1.5);
plot([0 ntStep], [-1/30. -1/30.], '-r', 'LineWidth', 1.5);
legend('(SIGMA.i-SIGMA)/SIGMA = f(t)','Convergence targets'); %title('(SIGMA.i-SIGMA)/SIGMA = f(t)');
% plot(moyb(1:ntStep), '-k');
plot(moyb(1:ntStep)+0.1,'-r', 'LineWidth', 2.5);
plot(moyb(1:ntStep)-0.1,'-r', 'LineWidth', 2.5);

subplot(2,3,3); hold on;
hist(log10(stda(1:ntStep)), 50); % title('Log10[|Mean.i-Mean|/SIGMA]');
% h = findobj(gca,'Type','patch'); set(h, 'FaceColor','r','EdgeColor','w')
h = findobj(gca,'Type','patch'); set(h, 'FaceColor','b', 'EdgeColor','b');
plot([-1 -1], [0 ntStep/50], '-r', 'LineWidth', 2.5);
legend('Nb.of occurences', 'Convergence targets (X in Log)');
plot([log10(1/30.) log10(1/30.)], [0 ntStep/50], '-r', 'LineWidth', 3);

subplot(2,3,6); hold on;
hist(log10(stdb(1:ntStep)), 50); % title('Log10[|SIGMA.i-SIGMA|/SIGMA]');
h = findobj(gca,'Type','patch'); set(h, 'FaceColor','b', 'EdgeColor','b');
plot([-1 -1], [0 ntStep/50], '-r', 'LineWidth', 2.5);
legend('Nb.of occurences/SIGMA', 'Convergence targets (X in Log)');
plot([log10(1/30.) log10(1/30.)], [0 ntStep/50], '-r', 'LineWidth', 3);


