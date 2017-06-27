%%----------------HEADER---------------------------%%
%Author:           Boris Segret
% Version & Date:
%                  V1.0, 28-05-2017, from stat_extracting in v4.0 24-04-2017
% CL=0 (v1.0)
%
% This produces some basic plots to assess the geometry of the scenario.
%
% I/
%    <scenario>
%    <runifod_MCdrivers>
% O/

clear;
path_to_scn = '..\cas_EME\';
fscenario   = 'scnErg_v4';
sigma_obs = 0.1;
% nKF = 8*24; dtKF = 3600./5.; % (dtKF in seconds)
nKF = 8; dtKF = 1800./5.; % (dtKF in seconds)

Nobs = 5;
addpath('../ifod');
data=[];
MJD_0=2400000.5;
SEC_0=86400;
scn=fopen([path_to_scn fscenario],'r'); %we open the scenario file
l=' ';
while 1
    l=fgetl(scn);
    if strfind(l,'META_STOP')>0 %We start reading the file from the META_STOP tag
        break;
    end;
end;
%l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; scenario_num=lgn;
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; outputs = strcat(path_to_scn,lgn);
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; algo=lgn;
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; simLims = double(str2num(lgn));
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; reftrajectory = strcat(path_to_scn,lgn);
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; actual_trajectory = strcat(path_to_scn,lgn);
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; nbofBodies = str2num(lgn);

refEphNameLength(1)=0; refEphemerid=''; BodyName='';
for ii=1:nbofBodies
    l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end;
    refEphNameLength(ii+1)=refEphNameLength(ii)+length(path_to_scn)+length(lgn)-11;
    refEphemerid = [refEphemerid path_to_scn lgn(12:length(lgn))];
    BodyName = [BodyName lgn(1:10)];
end
fclose(scn);

% Extraction of the trajectory data
% > TimeList0(:), TimeList1(:) in decimal JD, time-stamped of the trajectories ("0" ref & "1" actual)
% > NbLT0, NbLT1: Nb of lines in the trajectory files
% > coord0(:), coord1(:): in km
% > vel0(:), vel1(:): in m/s
[NbLT0, TimeList0, data] = readTraj(reftrajectory);
coord0 = data(:,1:3);
vel0   = data(:,4:6);
fprintf('Reading ref.traj''s data: %s, %.1f...%.1f = %i records\n', ...
  reftrajectory, TimeList0(1), TimeList0(NbLT0), NbLT0);
[NbLT1, TimeList1, data] = readTraj(actual_trajectory);
coord1 = data(:,1:3);
vel1   = data(:,4:6);
fprintf('Reading act.traj''s data: %s, %.1f...%.1f = %i records\n', ...
  actual_trajectory, TimeList1(1), TimeList1(NbLT1), NbLT1);

% Extraction of foreground body data, ii being the number of a given body:
% > NbLE(ii) => Nb of records in the ephemeris files
% > TimeListE(ii, 1:NbLE(ii)) => in decimal JD, time-stamped of the ephemerides files
% > ephXYZ(ii, 1:NbLE(ii), 1:3) => X(t),Y(t),Z(t) rectangular coordinates of the body
for ii=1:nbofBodies
  bodyFile = refEphemerid(1+refEphNameLength(ii):refEphNameLength(ii+1));
  [nbLgn, timeSteps, data] = readTraj(bodyFile);
  NbLE(ii) = nbLgn;
  mxLgn = max(NbLE);
  if (ii>1)
    tempoSteps = TimeListE;
    TimeListE = double(zeros(ii,mxLgn));
    TimeListE(1:ii-1, 1:size(tempoSteps,2)) = tempoSteps;
    tempoEph = ephXYZ;
    ephXYZ   = double(zeros(ii,mxLgn,3));
    ephXYZ(1:ii-1, 1:size(tempoEph,2), 1:3) = tempoEph;
  end
  TimeListE(ii,1:nbLgn) = timeSteps;
  ephXYZ(ii,1:nbLgn,1:3) = data(:,1:3);
  fprintf('Reading Body %i''s data: %s, %.1f...%.1f = %i records\n', ...
      ii, bodyFile, TimeListE(ii,1), TimeListE(ii,nbLgn), nbLgn);
end
clear tempoSteps tempoEph data;
T0N=datenum([2000 1 1 0 0 0]);
T0JD=2451544.5;
T0MJD= 51544.;
t_start = max([TimeList0(1) TimeList1(1) TimeListE(:,1)' simLims(1)]);
t_stop  = min([TimeList0(NbLT0) TimeList1(NbLT1) max(TimeListE(:,NbLE(:))) simLims(3)]);
starttxt = datestr(t_start-T0MJD+T0N, 'yyyy-mm-ddTHH:MM:SS');
stoptxt  = datestr(t_stop-T0MJD+T0N, 'yyyy-mm-ddTHH:MM:SS');
nbEpochs = 1+floor((t_stop-t_start)/simLims(2));
fprintf('Ready to run simulation from %.1f to %.1f, %.1fH intervals (%i steps)\n', ...
  t_start, t_stop, simLims(2)*24., nbEpochs);
if (nbEpochs <= 0) fprintf(' === WARNING! Nb.Steps < 0 ===\n'); end;

ii_MAX=length(TimeList1);

%--------- output preparations -----------
% dataExtraction=fopen(outputs,'w');
fprintf('Optical accuracy (arcsec): %11.5f\n', sigma_obs);
fprintf('IFOD-KF processing: %d steps * %.1fH = %.1f hours\n', nKF, dtKF/3600., (nKF)*dtKF/3600.);
fprintf('STARTING = %8.2f (dec.MJD), %s (yyyy-mm-ddTHH:MM:SS)\n', t_start, starttxt);
fprintf('ENDING   = %8.2f (dec.MJD), %s (yyyy-mm-ddTHH:MM:SS)\n', t_stop, stoptxt);
fprintf('STEPS = %i (Nb), %.1f (H)\n', nbEpochs, simLims(2)*24.);
%-----------------------------------------
% Memory allocations
% obstime = double(zeros(1,Nobs*nKF));
obstime = (t_start+(Nobs*nKF+2)*dtKF/86400.):simLims(2):t_stop;
FG_body = double(zeros(length(obstime), 3, nbofBodies));
observd = double(zeros(length(obstime), 3, nbofBodies));
predict = double(zeros(length(obstime), 3, nbofBodies));
drifts  = double(zeros(length(obstime), 4));
optics  = double(zeros(length(obstime), 3, nbofBodies));
method = 'linear';
ref_loc = interp1(TimeList0, coord0, obstime, method);
ref_vel = interp1(TimeList0, vel0, obstime, method);
act_loc = interp1(TimeList1, coord1, obstime, method);
drifts(:,1:3) =  act_loc - ref_loc;
drifts(:,4)   = sqrt(drifts(:,1).^2+drifts(:,2).^2+drifts(:,3).^2);

% ik=0;
for ib = 1:nbofBodies
   FGB(:,1:3) = interp1(TimeListE(ib, 1:NbLE(ib)), ...
       reshape(ephXYZ(ib, 1:NbLE(ib), 1:3), [NbLE(ib) 3]), obstime, method);
   FG_body(:,1:3,ib) = FGB;
   predict(:,1:3,ib) = FGB - ref_loc;
   observd(:,1:3,ib) = FGB - act_loc;
   % distance to the foreground bodies
   optics(:,1,ib) = sqrt(predict(:,1,ib).^2 + predict(:,2,ib).^2 + predict(:,3,ib).^2);
   optics(:,2,ib) = sqrt(observd(:,1,ib).^2 + observd(:,2,ib).^2 + observd(:,3,ib).^2);
   % scale factors wrt the foreground bodies
   optics(:,3,ib) = (sigma_obs./3600.).*(pi()/180.).*optics(:,1,ib);
   % drift angle to the foreground bodies from REFERENCE to ACTUAL location
   optics(:,4,ib) = acosd(dot(...
       predict(:,:,ib)./repmat(optics(:,1,ib),[1 3]), ...
       observd(:,:,ib)./repmat(optics(:,2,ib),[1 3]), ...
       2)); % in degrees
end
clear FGB;

% PLOTS
myColors=[0.5 0   0;...
          0.5 0.5 0;...
          0   0.5 0.5;...
          0   0   0.5;...
          0.4 0.4 0.4];
% myColors=[1 0   0;...
%           1 1 0;...
%           0   1 1;...
%           0   0   1;...
%           0.6 0.6 0.6];
figure(11); clf; %set(gca, 'LooseInset', [0.09 0.09 0.07 0.07]);
figure(12); clf; %set(gca, 'LooseInset', [0.09 0.09 0.07 0.07]);
semilogy(obstime-t_start, drifts(:,4), '--k'); hold on;
for ib = 1:nbofBodies
    % DRIFT ANGLES from Ref
    figure(11);
    semilogy(obstime-t_start, optics(:,4,ib), 'Color', myColors(ib,:)); hold on;

    % DRIFT DISTANCES and SCALE FACTORS
    figure(12);
    semilogy(obstime-t_start, optics(:,3,ib),'-', 'Color', myColors(ib,:)); hold on;
end
figure(11);
legend(['to ' BodyName( 1:10)], ...
       ['to ' BodyName(11:20)], ...
       ['to ' BodyName(21:30)], ...
       ['to ' BodyName(31:40)], ...
       ['to ' BodyName(41:50)], ...
        'Location', 'NorthWest');
ylabel('Ang.drift (degrees)');
xlabel('time (days)');
set(gca, 'XColor', 'm', 'YColor', 'm', 'YGrid', 'on', 'YMinorGrid', 'off');
% title(['Drift angles, ' datestr(t_start,'yyyy-mm-dd') '..' datestr(t_stop,'yyyy-mm-dd')]);
figure(12);
legend('(drift distance)',...
       ['sc.f.to ' BodyName( 1:10)], ...
       ['sc.f.to ' BodyName(11:20)], ...
       ['sc.f.to ' BodyName(21:30)], ...
       ['sc.f.to ' BodyName(31:40)], ...
       ['sc.f.to ' BodyName(41:50)], ...
        'Location', 'SouthWest');
ylabel('(kilometers)');
xlabel('time (days)');
set(gca, 'XColor', 'm', 'YColor', 'm', 'YGrid', 'on', 'YMinorGrid', 'off');
title(['Scale factors, ' datestr(T0N+t_start-T0MJD,'yyyy-mm-dd') '..' ...
    datestr(T0N+t_stop-T0MJD,'yyyy-mm-dd')]);

% figure(1); plot3(ref_loc(:,1)-FG_body(:,1,2), ref_loc(:,2)-FG_body(:,2,2), ref_loc(:,3)-FG_body(:,3,2));
