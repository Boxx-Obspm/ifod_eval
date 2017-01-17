%%----------------HEADER---------------------------%%
%Author:           Boris Segret
version = '3.3.x';
% Version & Date:
%                  V3.3 03-11-2016 (dd-mm-yyyy)
%                  - Rebuilt according to ../ifod/oneOD to introduce a Kalman Filter
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
% CL=2 (v3.2)
%
%
% This produces a comparison *with* error bars at every time step between the
% computed and the expected results. The assumption is to input a reference
% trajectory with its ephemerides for a number of foreground objects, and a
% slightly shifted trajectory with its ephemerides for the same foreground
% objects. All results are produced in the current workspace for easy access
% afterwards. A "scenario" file drives the computation.
%
% I/
%    <scenario> scenario file (in VTS format) path must be given in datapath below
%               see User Manual to set the scenario file
%    datapath : path for the "scenario" file, it *must* be provided in the workspace
% O/
%    <results> comparison between computed and expected results (VTS format,
%              prefixed as requested in <scenario>)

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
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; scenario_num=lgn;
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; reftrajectory = strcat(datapath,lgn);
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; nbofBodies = str2num(lgn);

% refEphNameLength(i): [path+filename] length of the ephemeris of the i-th foreground body for the reference trajectory
% refEphemerid: [path+filename] are put one after another
% actEphNameLength(i): [path+filename] length of the ephemeris of the i-th foreground body for the actual trajectory
% actEphemerid: [path+filename] are put one after another
refEphNameLength(1)=0; refEphemerid='';
for ii=1:nbofBodies
    l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end;
    refEphNameLength(ii+1)=refEphNameLength(ii)+length(datapath)+length(lgn)-11;
    refEphemerid = [refEphemerid datapath lgn(12:length(lgn))];
end
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; actual_trajectory = strcat(datapath,lgn);
actEphNameLength(1)=0; actEphemerid='';
for ii=1:nbofBodies
    l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end;
    actEphNameLength(ii+1)=actEphNameLength(ii)+length(datapath)+length(lgn)-11;
    actEphemerid = [actEphemerid datapath lgn(12:length(lgn))];
end
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; algo=lgn;
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; sampling = strcat(datapath,lgn);
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; outputs = strcat(datapath,lgn);
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; simLims = uint32(str2num(lgn));
fclose(scn);

%Extraction of the trajectory data, i being the number of a given foreground body
% > TimeList0(:), TimeList1(:) in decimal JD, time-stamped of the trajectories ("0" ref & "1" actual)
% > NbLT0, NbLT1: Nb of lines in the trajectory files
% > NbLE0(i), NbLE1(i): Nb of lines in the ephemeris files
% > TimeListE0(i,:), TimeListE1(i,:) in decimal JD, time-stamped of the ephemeris files ("0" ref & "1" actual)
% > lat0(i,:), long0(i,:), lat1(i,:), long1(i,:) in DEGREES
% > dist0(i,:), dist1(i,:) in km
% > coord0(:), coord1(:): in km
% > vel0(:), vel1(:): in m/s

[NbLT0, TimeList0, coord0, vel0] = readTraj(reftrajectory);
[NbLT1, TimeList1, coord1, vel1] = readTraj(actual_trajectory);
for ii=1:nbofBodies
  [NbLE0(ii), TimeListE0(ii,:), lat0(ii,:), long0(ii,:), dist0(ii,:)] = readEphem(refEphemerid(1+refEphNameLength(ii):refEphNameLength(ii+1)));
  [NbLE1(ii), TimeListE1(ii,:), lat1(ii,:), long1(ii,:), dist1(ii,:)] = readEphem(actEphemerid(1+actEphNameLength(ii):actEphNameLength(ii+1)));
end

ii_MAX=length(TimeList1);

% (prepare for STEP 2)
% alternative:
% (dtConst in hours)
dtKF = 86400.*1/24; % in seconds (measurement of all planets once an hour, with 5 planets)
dtConst = (nKF)*dtKF/3600.;
fprintf('WARNING: dynamic time step (not used: %s)\n', sampling);
fprintf('Time-sampling is nKF*dtKF = %i hours\n', dtConst);
% fprintf('WARNING: file %s is not read\n', sampling);
% fprintf('Time-sampling is forced constant for the full scenario: %i hours\n', dtConst);
% timeStep=fopen(sampling, 'r'); %we open the file timeStep-i
% 	l=' ';
% 	while 1
% 		l=fgetl(timeStep);
% 		if strfind(l,'META_STOP')>0
% 			break;
% 		end;
% 	end;
% [sampling_data,nb] = fscanf(timeStep,'%g %g %g %g %g', [5 inf]);
% sampling_data = sampling_data';
% fprintf('Changes in timesteps:\n');
% T=[sampling_data(:,1),sampling_data(:,2)]; fprintf('%15i %15i\n', T'); % we pick the sampling data and stock it in T and dt
% dt=[sampling_data(:,3),sampling_data(:,4),sampling_data(:,5)]; fprintf('%5i %5i %5i\n', dt');
% % Not tolerant to empty lines !!!
% fclose(timeStep);

outputs = [outputs pfix];

%--------- output preparations -----------
dataExtraction=fopen(outputs,'w');
fprintf(dataExtraction,'NAV_Results Version : %s\nGenerated by BIRDY NAV TEAM\nDate : %s\nalgo : %s\n\n', version, datestr(now), algo);
fprintf(dataExtraction,'OBJECT_NAME : BIRDY\nID_NAME : BIRDY\nSCENARIO : %s\nTRAJECTORY_NAME : %s\n\n',scenario_num, actual_trajectory);
fprintf(dataExtraction,'Observation accuracy (arcsec): %11.5f\n', sigma_obs);
fprintf(dataExtraction,'Time-Sampling between measurements (hours): %d\n', dtConst);
fprintf(dataExtraction,'\nForeground bodies:\n');
for ii=1:nbofBodies
  fprintf(dataExtraction,'%2i, ref. : %s\n', ii, refEphemerid(1+refEphNameLength(ii):refEphNameLength(ii+1)));
  fprintf(dataExtraction,'    actual %s\n', actEphemerid(1+actEphNameLength(ii):actEphNameLength(ii+1)));
end
fprintf(dataExtraction,'\nMETA_START');
fprintf(dataExtraction,'\n\nCOLUMN #01 : Day of the date (in MJD)\n');
fprintf(dataExtraction,'COLUMN #02 : Seconds in the day (in seconds)\n');
fprintf(dataExtraction,'COLUMN #03 : Expected Transversal shift (km) of Actual Traj wrt Reference\n');
fprintf(dataExtraction,'COLUMN #04 : Computed Transversal shift (km) of OD wrt Reference\n');
fprintf(dataExtraction,'COLUMN #05 : Std.dev. Transversal shift (km) of OD wrt Reference\n');
fprintf(dataExtraction,'COLUMN #06 : Expected Longitudinal shift (km) of Actual Traj wrt Reference\n');
fprintf(dataExtraction,'COLUMN #07 : Computed Longitudinal shift (km) of OD wrt Reference\n');
fprintf(dataExtraction,'COLUMN #08 : Std.dev. Longitudinal shift (km) of OD wrt Reference\n');
fprintf(dataExtraction,'COLUMN #09 : Colinearity of speed (in arcsec) during OD\n');
fprintf(dataExtraction,'COLUMN #10 : Uniformity of speed (in mm/s) during OD\n');
fprintf(dataExtraction,'COLUMN #11 : Nb.of tries\n');
fprintf(dataExtraction,'COLUMN #12 : date 2 (in seconds)\n');
fprintf(dataExtraction,'COLUMN #13 : date 3 (in seconds)\n');
fprintf(dataExtraction,'COLUMN #14 : date 4 (in seconds)\n');
fprintf(dataExtraction,'COLUMN #15 : date 5 (in seconds)\n');
fprintf(dataExtraction,'COLUMN #16 : performance (in CPU milliseconds)\n');
fprintf(dataExtraction,'COLUMN #17 to #42 : Vector X [(dx,dy,dz)_i,dr_i in km, dvx,dvy,dvz in km/s, dax,day,daz in km/s2]\n');
fprintf(dataExtraction,'COLUMN #43 to #68 : Vector X-expected (same units like X)\n\n');
fprintf(dataExtraction,'COLUMN #69 to #94 : Vector std(X) (same units like X)\n\n');
fprintf(dataExtraction, '   1            2           3           4           5           6           7           8       9        10     11        12          13          14          15         16   // X...    17          18          19          20          21          22          23          24          25          26          27          28          29          30          31          32          33          34          35          36       37          38          39        40          41           42    // X-exp   43          44          45          46          47          48          49          50          51          52          53          54          55          56          57          58          59          60          61          62       63          64          65        66           67          68    // std(X)   69          70          71          72          73          74          75          76          77          78          79          80          81          82          83          84          85          86          87          88       89          90          91        92          93          94   ///\n');
fprintf(dataExtraction,'META_STOP\n');
fclose(dataExtraction);
%-----------------------------------------
% construction of the observations and predicted observations

t = TimeList1(simLims(1));
t_fin = TimeList1(min([NbLT1 simLims(3)]));
% simLims(2) not used anymore
obstime = [];
for ii = simLims(1):simLims(2):min([NbLT1 simLims(3)])
  obstime = [obstime (TimeList1(simLims(1)+(ii-1))+(0:nKF-1)*dtKF/86400.)];
end
% for ii=1:nKF:length(obstime) fprintf('%12.5f ',obstime(ii:ii+nKF-1)); fprintf('\n'); end
obsbody = 1+mod([0:nKF-1], nbofBodies); % obsbody = [nbofBody] between 1 and nbofBodies
observd = double(zeros(nKF,2)); % observd = [out_lat1 out_long1];
predict = double(zeros(nKF,3)); % predict = [out_lat0 out_long0 out_distance0]

%-----------------------------------------
% feeding the KFd-OD (Kalman Filtered Orbit Determination)

ik=0;
iDebug=10; debug;
for ii = simLims(1):simLims(2):min([NbLT1 simLims(3)])
  post_sigma  = 0.1;
  delta_sigma = 1.;
  nbTries=0;
  metaX = double(zeros(1,4*Nobs+6));
  cumul_CPU = 0.;
  %notStabilized = true;
  notStabilized = false;
  % ----- pas de convergence d'erreur pour le moment
% while or(delta_sigma/post_sigma > 0.01, notStabilized)
% while or(delta_sigma/post_sigma > 0.10, notStabilized)
% the Monte-Carlo analysis is run until sigma has been stabilized at 3% accuracy
% notStabilized = (delta_sigma/post_sigma > 0.01); % previous value
  %
  % ----- HERE IS A STILL UNEXEPLAINED BUG! -----
  if (scnRealistic)
      % -long1/-long0 avec les donnees realistes (bug!)
      observd = extractObs(obstime(ik+1:ik+nKF), obsbody, NbLE1, TimeListE1, lat1, -long1);
      predict = prepareObs(obstime(ik+1:ik+nKF), obsbody, NbLE0, TimeListE0, lat0, -long0, dist0);
  else
      % +long1/+long0 avec les donnees simulees (bug!)
      observd = extractObs(obstime(ik+1:ik+nKF), obsbody, NbLE1, TimeListE1, lat1, long1);
      predict = prepareObs(obstime(ik+1:ik+nKF), obsbody, NbLE0, TimeListE0, lat0, long0, dist0);
  end
  % -----------------------------------------------
  
  err_obs = normrnd(0., sigma_obs/3600., [nbCycles nKF 2]);
  err_obs(:,:,2) = err_obs(:,:,2)./repmat(cosd(observd(:,1)'), nbCycles,1);

  epochs(1:Nobs-1) = obstime(ik+1:ik+Nobs-1);
  refState = [ interp1(TimeList0, coord0(:,1), obstime(ik+3), 'linear'); ...
       interp1(TimeList0, coord0(:,2), obstime(ik+3), 'linear'); ...
       interp1(TimeList0, coord0(:,3), obstime(ik+3), 'linear'); ...
       interp1(TimeList0, vel0(:,1),   obstime(ik+3), 'linear'); ...
       interp1(TimeList0, vel0(:,2),   obstime(ik+3), 'linear'); ...
       interp1(TimeList0, vel0(:,3),   obstime(ik+3), 'linear') ];
  iDebug=11; debug;
  for nm=1:nbCycles
    nbTries = nbTries +1;
    measurd = observd + [err_obs(nm,:,1)' err_obs(nm,:,2)'];
    refLoc = [0;0;0];
    iDebug=1; debug;
    %------------------
    % KALMAN FILTER, steps Nobs to nKF
%     for ij = Nobs:nKF
%       epochs(1:Nobs)   = obstime(ik+ij-Nobs+1:ik+ij);
%       measur(1:Nobs,:) = measurd(ij-Nobs+1:ij,:);
%       predic(1:Nobs,:) = predict(ij-Nobs+1:ij,:);
%       bodies(1:Nobs)   = obsbody(ij-Nobs+1:ij);
      oneOD;
%     end
    % last X vector provides the best determination
    %------------------
    iDebug=4; debug;
    if (nbTries==1)
      metaX = X';
    else
      metaX = [metaX; X'];
    end
    cumul_CPU = cumul_CPU + elapsed_time;
    pre_sigma   = post_sigma;
    post_sigma  = sum(std(metaX(:,7:9)));
    delta_sigma = abs(post_sigma - pre_sigma);
  end
  ik=ik+nKF;

% % % % %         %**********
% % % %           % needed in oneOD at the moment for delta-acceleration estimate
% % % % %           Xexp = expectedOD (TimeList0, NbLE0, TimeListE0, dist0, coord0, vel0, ...
% % % % %                      epochs, nbofBodies, ...
% % % % %                      TimeList1, NbLE1, TimeListE1, dist1,coord1,vel1);
% % % % %         %**********
  X = mean(metaX,1);
  stdX = std(metaX,0,1);
  elapsed_time = cumul_CPU./nbTries;
  % --- expected solution:
                 
  % Extraction of the date in MJD (integer) and seconds:
  day=fix(TimeList1(ii)-MJD_0);
  sec=mod(TimeList1(ii)-MJD_0,1)*SEC_0;
	
  % Transversal and Longitudinal shifts in the trajectories as computed and as expected
  % unitvvector = unit_speed_vector(ii,vel1);
  unitvvector(1) = interp1(TimeList1, vel0(:,1), epochs(3), 'linear');
  unitvvector(2) = interp1(TimeList1, vel0(:,2), epochs(3), 'linear');
  unitvvector(3) = interp1(TimeList1, vel0(:,3), epochs(3), 'linear');
  unitvvector = unitvvector./norm(unitvvector);
  %==> as expected:
  trans_traj = norm(cross(Xexp(7:9)', unitvvector));
  longi_traj = dot(Xexp(7:9), unitvvector);
  %==> as computed:
  transmeta = cross(metaX(:,7:9), repmat(unitvvector, [size(metaX,1) 1]));
  trans_std = std(sqrt(transmeta(:,1).^2+transmeta(:,2).^2+transmeta(:,3).^2));
  longimeta = dot(metaX(:,7:9)', repmat(unitvvector, [size(metaX,1) 1])');
  longi_std = std(longimeta);
  trans_err = norm(cross(X(7:9), unitvvector));
  longi_err = dot(X(7:9)', unitvvector);
  
  % Colinearity and Intensity of the velocity during OD
  % with Nobs=5 the criterion to be watched for is "dacc"
  vIni(1)   = interp1(TimeList0, vel0(:,1), epochs(1), 'linear');
  vIni(2)   = interp1(TimeList0, vel0(:,2), epochs(1), 'linear');
  vIni(3)   = interp1(TimeList0, vel0(:,3), epochs(1), 'linear');
  vFinal(1) = interp1(TimeList0, vel0(:,1), epochs(Nobs), 'linear');
  vFinal(2) = interp1(TimeList0, vel0(:,2), epochs(Nobs), 'linear');
  vFinal(3) = interp1(TimeList0, vel0(:,3), epochs(Nobs), 'linear');
  speed_angle = (180.*3600./pi)*acos( dot(vIni, vFinal) / (norm(vIni)*norm(vFinal)) );
  norme=1000.*(norm(vFinal)-norm(vIni)); % m/s
	
  % 3 next dates that were considered with the current epoch during OD
  nxdates=epochs(2:Nobs)-MJD_0;

  % CPU time
  CPU=elapsed_time;

  % concatenation
  data = [day, sec, trans_traj, trans_err, trans_std, longi_traj, longi_err, longi_std, ...
          speed_angle, norme, nbTries, nxdates, CPU, ...
          X, ...
          Xexp', ...
          stdX];
  iDebug=12; debug;
  %--------- output writing ----------------
  dataExtraction=fopen(outputs,'a');
  fprintf(dataExtraction, ['%0.3d %12.4f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f' ...
                        ' %7.1f %9.2f %5d %11.4f %11.4f %11.4f %11.4f %13.3e' ...
                        ' %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f' ...
                        ' %11.3f %11.3f %11.3f %11.3f %11.3f %11.6f %11.6f %11.6f %11.8f %11.8f %11.8f' ...
                        ' %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f' ...
                        ' %11.3f %11.3f %11.3f %11.3f %11.3f %11.6f %11.6f %11.6f %11.8f %11.8f %11.8f' ...
                        ' %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f' ...
                        ' %11.3f %11.3f %11.3f %11.3f %11.3f %11.6f %11.6f %11.6f %11.8f %11.8f %11.8f\n'], data');
  fclose(dataExtraction);
  %-----------------------------------------
%   fprintf(' : %d => %d tries\n', ii, nbTries);
end
