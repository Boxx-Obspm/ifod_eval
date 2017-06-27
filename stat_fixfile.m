% Version 1.0, Author B.Segret
% multi-purpose IFOD files troubleshooting
bin1 = '../cas_Y/outs/Y+Y41324_E-1,400x200h,P41324,1355332_bin'; % filename must be provided
bin2 = '../cas_Y/outs/Y+Y41324_E-1,400x200h,P41324(101+),1355334_bin'; % filename must be provided
outf = [bin1 '_fixed'];
n1 = 50; % first n1 time steps from bin1
n2=50; % first n2 time steps from bin2 (required: same nb of cycles in bin1 & bin2)
graphs = true;
graphs = false;

T0N=datenum([2000 1 1 0 0 0]);
T0JD=2451544.5;

fo  = fopen(outf,'w'); fclose(fo);

% initizes the binary outputs with the number of written columns
fi = fopen(bin1,'r');
fprintf('File 1: %s\n', bin1);
ntStep = 0;
lastData = false;
while (not(lastData) & ntStep < n1)
    % fwrite(fw, [obstime(ik+nKF-2) obstime(ik+nKF)], 'double');
    T = fread(fi, 2, 'double'); % will provide a 1x2 double array
    %                             (time for KF-measured and for KF-predicted)
    t2=T(2)-T0JD+T0N;
    stf=datestr(t2,'yyyy-mm-dd HH:MM:SS');

    % fwrite(fw, [Nobs nKF nbCycles (ik+nKF==length(obstime))], 'uint32');
    Z = fread(fi, 4, 'uint32'); % will provide a 1x4 uint32 array
    ntStep = ntStep+1;

    nbPts=Z(2)-Z(1)+1; % (nbPts=nKF-Nobs+1;)
    Nobs=Z(1);
    nKF=Z(2);
    nbCycls1=Z(3);
    lastData = (Z(4)==1);

    fo  = fopen(outf,'a');
    fwrite(fo, T, 'double');
    fwrite(fo, Z, 'uint32');

    % solKF = double(zeros(nbCycles, nvKF));
    % dtlKF = double(zeros(nbCycles*nbPts, 2));

    for nC=1:nbCycls1
      % fwrite(fw, ...
      %     [rex rrme rrkf lKg ldP vkf dtrm' dlgm' dtrk' dlgk' mmkf' mmtk' mmlk'], ...
      %     'double');
      rawDATA = fread(fi, (3*3+10)*(nbPts+1), 'double');
      fwrite(fo, rawDATA, 'double');
    end
    fclose(fo);
    % if (mod(ntStep,100)==1) fprintf('#%04i (%s) ', ntStep, stf); end
    % if (mod(ntStep,10)==0)  fprintf('/'); else fprintf('.'); end
    % if (mod(ntStep,100)==0) fprintf(':#%i = (%s)\n', ntStep, stf); end
    fprintf('%04i: %s, %i %i %i %i\n', ntStep, stf, Z);

end
nf1=ntStep;
fclose(fi);

fi = fopen(bin2,'r');
fprintf('File 2: %s\n', bin2);
ntStep = 0;
lastData = false;
while (not(lastData) & ntStep < n2)
    % fwrite(fw, [obstime(ik+nKF-2) obstime(ik+nKF)], 'double');
    T = fread(fi, 2, 'double'); % will provide a 1x2 double array
    %                             (time for KF-measured and for KF-predicted)
    t2=T(2)-T0JD+T0N;
    stf=datestr(t2,'yyyy-mm-dd HH:MM:SS');

    % fwrite(fw, [Nobs nKF nbCycles (ik+nKF==length(obstime))], 'uint32');
    Z = fread(fi, 4, 'uint32'); % will provide a 1x4 uint32 array
    ntStep = ntStep+1;

    nbPts=Z(2)-Z(1)+1; % (nbPts=nKF-Nobs+1;)
    Nobs=Z(1);
    nKF=Z(2);
    nbCycls1=Z(3);
    if (ntStep==n2) Z(4)=true; end
    lastData = (Z(4)==1);

    fo  = fopen(outf,'a');
    fwrite(fo, T, 'double');
    fwrite(fo, Z, 'uint32');

    % solKF = double(zeros(nbCycles, nvKF));
    % dtlKF = double(zeros(nbCycles*nbPts, 2));

    for nC=1:nbCycls1
      % fwrite(fw, ...
      %     [rex rrme rrkf lKg ldP vkf dtrm' dlgm' dtrk' dlgk' mmkf' mmtk' mmlk'], ...
      %     'double');
      rawDATA = fread(fi, (3*3+10)*(nbPts+1), 'double');
      fwrite(fo, rawDATA, 'double');
    end
    fclose(fo);
    % if (mod(ntStep,100)==1) fprintf('#%04i (%s) ', ntStep, stf); end
    % if (mod(ntStep,10)==0)  fprintf('/'); else fprintf('.'); end
    % if (mod(ntStep,100)==0) fprintf(':#%i = (%s)\n', ntStep, stf); end
    fprintf('%04i: %s, %i %i %i %i\n', ntStep, stf, Z);

end
nf2=ntStep;
fclose(fi);

fprintf('Concatenated:\n - %s\n - %s\n => %s\n', bin1, bin2, outf);
fprintf('=> Nb.of time-steps (nf1+nf2): %i\n', nf1+nf2);
fprintf('=> Nb.of Cycles per time-step (unchanged): %i\n', nbCycls1);
outputs_o = outf;
stat_reading;
