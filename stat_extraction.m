%%----------------HEADER---------------------------%%
%Author:           Boris Segret
version = '4.0';
% Version & Date:
%                  V4.0, 24-04-2017
%                  - Trajectories are provided in new .xva format, with gravity data
%                  - Ephemerides are provided in rectangular absolute coordinates,
%                  - -> only 1 set, the expected and actual directions are computed
%                  - -> sampling will be still based on the optical accuracy
%                  - scenario format changed
%                  - Nb of time-steps must not be equal in input files anymore
% CL=2 (v3.3.4)
%                  V3.3.4, 23-03-2017
%                  - all bugs fixed
%                  - adaptation to 5 new observations for each KF step
%                  V3.3 03-11-2016 (dd-mm-yyyy) V3.3.2 19-02-2017 (see CL below)
%                  V3.3 03-11-2016 (dd-mm-yyyy) V3.3.1 27-01-2017
%                  - Rebuilt according to ../ifod/oneOD to introduce a Kalman Filter
%                  - Output formats changed: only KF result in "a priori" last prediction
%                  - Input format changed: no more time sampling at scenario level
%                  V3.2 26-08-2016 (dd-mm-yyyy)
%                  - Output formats are changed to focus on M3 (M1..M5)
%                  V3.1 06-08-2016
%                  - Adaptation to N=5 measurements
%                  V3.0 15-05-2016
%                  - Monte-Carlo series introduced to study sensitivity of the algorithm
%                  -> time-sampling forced constant and driven from runifod.m
%                  -> observation accuracy driven from runifod.m
%                  -> sensitivity is studied on the transversal and longitudinal error distributions
%                  - written from "data_extraction.m" v2.1
%                  - REMAINING ISSUE: some output colums are not adapted to several foreground bodies
%                  - REMAINING ISSUE: the outputs are not adapted to monte carlo statistics
%                  V2.1 10-04-2016, Boris Segret
%                  - inputs & outputs in VTS format only
%                  - separeted "in-flight" and "test" modules, here only "test"
%                  - multiple foreground bodies
%                  - rectified header in the results files & added a verion number
%                  V1.2 03-03-2016, Boris Segret
%                  - *no* multiple calls of reference_trajectory.m
%                  - new "scenario" format, and computes only in an interval
%                  - minor adapations
%                  - works with a call from "ifod_tests"
%                  until V1   11-09-2015, Tristan Mallet
%
%
% This produces a comparison *without* error bars at every time step between the
% computed and the expected results. It also produces a fully detailed RAW data file
% for further statistical survey at each KF step or cycle step.
%
% The assumption is to input a reference trajectory, ephemerides of a number
% of foreground objects, and a slightly shifted trajectory wrt the reference.
%
% Various files drive the computation, they are loaded by runifod main.
%
% I/
%    <runifod_scenario>
%    <runifod_MCdrivers>
%               see User Manual
% O/
%    <simplified results> IFOD-computed vs. expected results at successive time steps
%    <raw data> detailed IFOD-results at each Kalman filter loops and cylce loops


addpath('../ifod');

data=[];
MJD_0=2400000.5;
SEC_0=86400;
scn=fopen([datapath fscenario],'r'); %we open the scenario file
l=' ';
while 1
    l=fgetl(scn);
    if strfind(l,'META_STOP')>0 %We start reading the file from the META_STOP tag
        break;
    end;
end;
%l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; scenario_num=lgn;
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; outputs = strcat(datapath,lgn);
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; algo=lgn;
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; simLims = double(str2num(lgn));
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; reftrajectory = strcat(datapath,lgn);
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; actual_trajectory = strcat(datapath,lgn);
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; nbofBodies = str2num(lgn);

refEphNameLength(1)=0; refEphemerid='';
for ii=1:nbofBodies
    l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end;
    refEphNameLength(ii+1)=refEphNameLength(ii)+length(datapath)+length(lgn)-11;
    refEphemerid = [refEphemerid datapath lgn(12:length(lgn))];
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

outputs = [outputs pfix];

%--------- output preparations -----------
dataExtraction=fopen(outputs,'w');
fprintf(dataExtraction,'NAV_Results Version : %s\nGenerated by BIRDY NAV TEAM\nDate : %s\nalgo : %s\n\n', version, datestr(now), algo);
fprintf(dataExtraction,'OBJECT_NAME : BIRDY\nID_NAME : BIRDY\nSCENARIO : %s\n', outputs);
fprintf(dataExtraction,'REFERENCE TRAJECTORY : %s\n', reftrajectory);
fprintf(dataExtraction,'ACTUAL TRAJECTORY : %s\n\n', actual_trajectory);
fprintf(dataExtraction,'Optical accuracy (arcsec): %11.5f\n', sigma_obs);
fprintf(dataExtraction,'IFOD-KF processing: %d steps * %.1fH = %.1f hours\n', nKF, dtKF/3600., (nKF)*dtKF/3600.);
fprintf(dataExtraction,'\nForeground bodies:\n');
for ii=1:nbofBodies
  fprintf(dataExtraction,'%2i, ref. : %s\n', ii, refEphemerid(1+refEphNameLength(ii):refEphNameLength(ii+1)));
  %   fprintf(dataExtraction,'    actual %s\n', actEphemerid(1+actEphNameLength(ii):actEphNameLength(ii+1)));
end
fprintf(dataExtraction,'\nMETA_START');
fprintf(dataExtraction,'\n\nSTARTING = %8.2f (dec.MJD), %s (yyyy-mm-ddTHH:MM:SS)\n', t_start, starttxt);
fprintf(dataExtraction,'ENDING   = %8.2f (dec.MJD), %s (yyyy-mm-ddTHH:MM:SS)\n', t_stop, stoptxt);
fprintf(dataExtraction,'STEPS = %i (Nb), %.1f (H)\n', nbEpochs, simLims(2)*24.);
fprintf(dataExtraction,'COLUMN #01 : Day of the date (in MJD)\n');
fprintf(dataExtraction,'COLUMN #02 : Seconds in the day (in seconds)\n');
fprintf(dataExtraction,'COLUMN #03 : Computed Transversal shift (km) of OD wrt Reference\n');
fprintf(dataExtraction,'COLUMN #04 : Expected Transversal shift (km) of Actual Traj wrt Reference\n');
fprintf(dataExtraction,'COLUMN #05 : Computed Longitudinal shift (km) of OD wrt Reference\n');
fprintf(dataExtraction,'COLUMN #06 : Expected Longitudinal shift (km) of Actual Traj wrt Reference\n');
fprintf(dataExtraction,'COLUMN #07 : Colinearity of speed (in arcsec) during OD\n');
fprintf(dataExtraction,'COLUMN #08 : Uniformity of speed (in m/s) during OD\n');
fprintf(dataExtraction,'COLUMN #09 : Nb.of Monte-Carlo selections\n');
fprintf(dataExtraction,'COLUMN #10 : performance (in CPU milliseconds)\n');
fprintf(dataExtraction,'COLUMN #11 to #13 : Computed Vector X(10:12) if Nobs=4, X(13:15) if Nobs=5 [(dx,dy,dz) in km]\n');
fprintf(dataExtraction,'COLUMN #14 to #16 : Expected Vector X-expected [(dx,dy,dz) in km]\n\n');
fprintf(dataExtraction, '   1            2           3           4           5           6       7         8      9        10           11          12          13          14          15          16\n');
fprintf(dataExtraction,'META_STOP\n');
fclose(dataExtraction);
%-----------------------------------------

% Construction of the observations and predicted observations
%------------------------------------------------------------
% Memory allocations
obstime = double(zeros(1,Nobs*nKF));
obsbody = 1+mod([0:Nobs*nKF-1], nbofBodies); % obsbody = [nbofBody] between 1 and nbofBodies
observd = double(zeros(Nobs*nKF,2)); % observd = [out_lat1 out_long1];
predict = double(zeros(Nobs*nKF,3)); % predict = [out_lat0 out_long0 out_distance0]
obsXYZ  = double(zeros(Nobs*nKF,3)); % rectangular coordinates of the foreground bodies at obstime

% Feeding the IFOD-KF (Kalman Filtered Orbit Determination)
ik=0;
iDebug=10; debug;
tic;
for tt = t_start:simLims(2):t_stop
obstime = (tt+(0:Nobs*nKF-1)*dtKF/86400.);
if obstime(length(obstime)) <= t_stop
  % nbTries=0;
  % ----- HERE IS A STILL UNEXPLAINED BUG! -----
%   if (scnRealistic)
%       % -long1/-long0 avec les donnees realistes (bug!)
%       observd = extractObs(obstime, obsbody, NbLE1, TimeListE1, lat1, -long1);
%       predict = prepareObs(obstime, obsbody, NbLE0, TimeListE0, lat0, -long0, dist0);
%   else
      % +long1/+long0 avec les donnees simulees (bug!)
%       observd = extractObs(obstime, obsbody, NbLE1, TimeListE1, lat1, long1);
%       predict = prepareObs(obstime, obsbody, NbLE0, TimeListE0, lat0, long0, dist0);
    method = 'linear';
    for i=1:Nobs*nKF
      ttimes = TimeListE(obsbody(i),1:NbLE(obsbody(i)));
      coords = reshape(ephXYZ(obsbody(i), 1:NbLE(obsbody(i)), 1:3), [NbLE(obsbody(i)) 3]);
      obsXYZ(i,1:3) = interp1(ttimes, coords, obstime(i), method);
      % fprintf('%11.1f %11.1f %11.1f\n', obsXYZ(i,1:3));
    end
    fromXYZ0 = interp1(TimeList0, coord0, obstime, method);
    vv = obsXYZ - fromXYZ0;
    rr3 = sqrt(vv(:,1).^2+vv(:,2).^2+vv(:,3).^2);
    rr2 = sqrt(vv(:,1).^2+vv(:,2).^2);
    llat = asind(vv(:,3)./rr3); % why in degrees ?!
    llng = asind(vv(:,2)./rr2); % why in degrees ?!
    mirror=(vv(:,1)<0); llng=llng.*(1.-2.*mirror) + mirror*180.; % case of 90° < longitude < 270°
    predict = [llat llng rr3];
    
    fromXYZ1 = interp1(TimeList1, coord1, obstime, method);
    vv = obsXYZ - fromXYZ1;
    rr3 = sqrt(vv(:,1).^2+vv(:,2).^2+vv(:,3).^2);
    rr2 = sqrt(vv(:,1).^2+vv(:,2).^2);
    llat = asind(vv(:,3)./rr3); % why in degrees ?!
    llng = asind(vv(:,2)./rr2); % why in degrees ?!
    mirror=(vv(:,1)<0); llng=llng.*(1.-2.*mirror) + mirror*180.; % case of 90° < longitude < 270°
    observd = [llat llng];
%   end
  clear ttimes coords llat llng mirror;
  clear vv rr3 rr2;
  % clear fromXYZ0 fromXYZ1;
  % -----------------------------------------------
  
  err_obs = normrnd(0., sigma_obs/3600., [nbCycles Nobs*nKF 2]);
  err_obs(:,:,2) = err_obs(:,:,2)./repmat(cosd(observd(:,1)'), nbCycles,1);

%   refState = [ interp1(TimeList0, coord0(:,1), obstime(3), 'linear'); ...
%        interp1(TimeList0, coord0(:,2), obstime(3), 'linear'); ...
%        interp1(TimeList0, coord0(:,3), obstime(3), 'linear'); ...
%        interp1(TimeList0, vel0(:,1),   obstime(3), 'linear'); ...
%        interp1(TimeList0, vel0(:,2),   obstime(3), 'linear'); ...
%        interp1(TimeList0, vel0(:,3),   obstime(3), 'linear') ];
  refState = interp1(TimeList0, [coord0 vel0], obstime(3), method);
  
  iDebug=11; debug; % initizes the monitoring file before each MC simu
  metaX = double(zeros(nbCycles, 3));
  tci0=toc;
  cumul_CPU = 0.;
  for nm=1:nbCycles
%     nbTries = nbTries +1;
    measurd = observd + [err_obs(nm,:,1)' err_obs(nm,:,2)'];
    refLoc = double([0 0 0]);
    iDebug=1; debug;
    %------------------
    % KALMAN FILTER, steps Nobs to nKF
    
      oneOD;
%       tcf0=toc; fprintf('Step %i, MC cycle #%i (ifod_kf): %5.2f ms\n', ii,nm,(tcf0-tci0)*1000.);
%       tci0=toc;

    % returned X vector provides the last and best determination
    % returned Xexp is the associated X-expected (same epochs)
    %------------------
    iDebug=4; debug;
    metaX(nm,:) = X(3*Nobs-2:3*Nobs)';
  end
%   elapsed_time = cumul_CPU./nbTries; % then ok to divide by nbTries
  elapsed_time = cumul_CPU./nbCycles;
  ik=ik+nKF;

  % --- Expected solution: at last epoch(5|4) in the KF (last update, i.e."a priori")
                 
  % Extraction of the date in MJD (integer) and seconds:
  day=floor(epochs(Nobs));
  sec=(epochs(Nobs)-day)*SEC_0;

  % Transversal and Longitudinal shifts in the trajectories as computed and as expected
  % unitvvector was already computed in iDebug==3 
  unitvvector = interp1(TimeList0, vel0, epochs(Nobs), method); % last epoch in KF!
%   unitvvector(1) = interp1(TimeList0, vel0(:,1), epochs(Nobs), 'linear'); % last epoch in KF!
%   unitvvector(2) = interp1(TimeList0, vel0(:,2), epochs(Nobs), 'linear'); % last epoch in KF!
%   unitvvector(3) = interp1(TimeList0, vel0(:,3), epochs(Nobs), 'linear'); % last epoch in KF!
  unitvvector = unitvvector./norm(unitvvector);
  %==> as expected:
  trans_traj = norm(cross(Xexp(3*Nobs-2:3*Nobs)', unitvvector));
  longi_traj = dot(Xexp(3*Nobs-2:3*Nobs), unitvvector);
  %==> as computed:
%   Xobs = mean(metaX,1);
%   trans_err = norm(cross(Xobs, unitvvector));
%   longi_err = dot(Xobs', unitvvector);
%   trans_err = norm(cross(X(7:9)', unitvvector));
%   longi_err = dot(X(7:9)', unitvvector);
  trans_err = mean(sqrt(metaX(:,1).^2+metaX(:,2).^2+metaX(:,3).^2));
  longimeta = dot(metaX(:,1:3)', repmat(unitvvector, [size(metaX,1) 1])');
  longi_err = mean(longimeta);

  % Colinearity and Intensity of the velocity during OD
  vIni(1)   = interp1(TimeList0, vel0(:,1), epochs(1), 'linear');
  vIni(2)   = interp1(TimeList0, vel0(:,2), epochs(1), 'linear');
  vIni(3)   = interp1(TimeList0, vel0(:,3), epochs(1), 'linear');
  vFinal(1) = interp1(TimeList0, vel0(:,1), epochs(Nobs), 'linear');
  vFinal(2) = interp1(TimeList0, vel0(:,2), epochs(Nobs), 'linear');
  vFinal(3) = interp1(TimeList0, vel0(:,3), epochs(Nobs), 'linear');
  speed_angle = (180.*3600./pi)*acos( dot(vIni, vFinal) / (norm(vIni)*norm(vFinal)) ); % asec
  norme=1000.*(norm(vFinal)-norm(vIni)); % m/s
	
  % CPU time
  CPU=elapsed_time;

  % concatenation
  data = [day, sec, trans_err, trans_traj, longi_err, longi_traj, ...
          speed_angle, norme, nbCycles, CPU, ...
          mean(metaX,1), ...
          Xexp(3*Nobs-2:3*Nobs)'];
  iDebug=12; debug;
  %--------- output writing ----------------
  dataExtraction=fopen(outputs,'a');
  fprintf(dataExtraction, ['%0.3d %12.4f %11.3f %11.3f %11.3f %11.3f' ...
                        ' %7.1f %9.2f %5d %13.3e' ...
                        ' %11.3f %11.3f %11.3f' ...
                        ' %11.3f %11.3f %11.3f\n'], data');
  fclose(dataExtraction);
  %-----------------------------------------
else
  fprintf('Job ended: Last step was not compatible with a full KF process\n');
end
end
