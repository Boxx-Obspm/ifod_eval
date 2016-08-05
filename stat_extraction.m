%%----------------HEADER---------------------------%%
%Author:           Boris Segret
version = '3.0';
% Version & Date:
%                  V3.0 15-05-2016 (dd-mm-yyyy)
%                  - Monte-Carlo series introduced to study sensitivity of the algorithm
%                  -> time-sampling forced constant and driven from runifod.m
%                  -> observation accuracy driven from runifod.m
%                  -> sensitivity is studied on the transversal and longitudinal error distributions
%                  - written from "data_extraction.m" v2.1
%                  - REMAINING ISSUE: some output colums are not adapted to several foreground bodies
%                  - REMAINING ISSUE: the outputs are not adapted to monte carlo statistics
% Version & Date:
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
% CL=0
%
%
% This produces a comparison *with* error bars at every time step between the
% computed and the expected results. The assumption is to input a reference
% trajectory with its ephemerides for a number of foreground objetcs, and a
% slightly shifted trajectory with  its ephemerides for the same foreground
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
%scn=fopen(strcat(datapath, 'scenario'),'r'); %we open the scenario file
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
% actEphNameLength(i): [path+filename] length of the ephemeris of the i-th foreground body for the actual trajectoru
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
% > vel0(:), vel(:): in m/s

[NbLT0, TimeList0, coord0, vel0] = readTraj(reftrajectory);
[NbLT1, TimeList1, coord1, vel1] = readTraj(actual_trajectory);
for ii=1:nbofBodies
  [NbLE0(ii), TimeListE0(ii,:), lat0(ii,:), long0(ii,:), dist0(ii,:)] = readEphem(refEphemerid(1+refEphNameLength(ii):refEphNameLength(ii+1)));
  [NbLE1(ii), TimeListE1(ii,:), lat1(ii,:), long1(ii,:), dist1(ii,:)] = readEphem(actEphemerid(1+actEphNameLength(ii):actEphNameLength(ii+1)));
end

ii_MAX=length(TimeList1);

fprintf('WARNING: file %s is not read\n', sampling);
fprintf('Time-sampling is forced constant for the full scenario: %i hours\n', dtConst);
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
T=[floor(TimeList1(1)-2400000.5), 0; floor(TimeList1(ii_MAX)-2400000.5)+1, 0];
dt=[dtConst, dtConst, dtConst; 0, 0, 0];
outputs = [outputs pfix];

%--------- output preparations -----------
dataExtraction=fopen(outputs,'w');
fprintf(dataExtraction,'NAV_Results Version : %s\nGenerated by BIRDY NAV TEAM\nDate : %s\nalgo : %s\n\n', version, datestr(now), algo);
fprintf(dataExtraction,'OBJECT_NAME : BIRDY\nID_NAME : BIRDY\nSCENARIO : %s\nTRAJECTORY_NAME : %s\n\n',scenario_num, actual_trajectory);
fprintf(dataExtraction,'META_START\n\n');
fprintf(dataExtraction,'Observation accuracy (arcsec): %11.5f\n', sigma_obs);
fprintf(dataExtraction,'Time-Sampling between measurements (hours): %d\n', dtConst);
%fprintf(dataExtraction,'Sampling of observations are changed at the following dates :\n');
%fprintf(dataExtraction,'%15d%15f\n',T');
%fprintf(dataExtraction,'\nSamplings (hours after previous obs) :\n');
%fprintf(dataExtraction,'%15d%15d%15d\n',dt');
fprintf(dataExtraction,'\n\nCOLUMN #01 : Day of the date (in MJD)\n');
fprintf(dataExtraction,'COLUMN #02 : Seconds in the day (in seconds)\n');
fprintf(dataExtraction,'COLUMN #03 : Transversal difference between the reference and the actual trajectories (in km)\n');
fprintf(dataExtraction,'COLUMN #04 : Longitudinal difference between the reference and the actual trajectories (in km)\n');
fprintf(dataExtraction,'COLUMN #05 : Latitude diff.of a foreground body between the ref.trajectory point of view and the actual trajectory one (in arcsec)\n');
fprintf(dataExtraction,'COLUMN #06 : Longitude diff.of a foreground body between the ref.trajectory point of view and the actual trajectory one(in arcsec)\n');
fprintf(dataExtraction,'COLUMN #07 : Transversal Error (in km) of OD wrt actual trajectory\n');
fprintf(dataExtraction,'COLUMN #08 : Longitudinal Error (in km) of OD wrt actual trajectory\n');
fprintf(dataExtraction,'COLUMN #09 : Colinearity of speed (in arcsec) during OD\n');
fprintf(dataExtraction,'COLUMN #10 : Uniformity of speed (in mm/s) during OD\n');
fprintf(dataExtraction,'COLUMN #11 : Nb.of tries\n');
fprintf(dataExtraction,'COLUMN #12 : date 2 (in seconds)\n');
fprintf(dataExtraction,'COLUMN #13 : date 3 (in seconds)\n');
fprintf(dataExtraction,'COLUMN #14 : date 4 (in seconds)\n');
fprintf(dataExtraction,'COLUMN #15 : performance (in CPU milliseconds)\n');
fprintf(dataExtraction,'COLUMN #16 to #34 : Vector X [(dx,dy,dz)_i,dr_i in km, dvx,dvy,dvz in km/s]\n');
fprintf(dataExtraction,'COLUMN #35 to #53 : Vector X-expected (same units like X)\n\n');
fprintf(dataExtraction,'COLUMN #54 to #72 : Vector std(X) (same units like X)\n\n');
fprintf(dataExtraction, '   1            2           3           4           5           6           7           8       9        10     11         12          13          14              15          16          17          18          19          20          21          22          23          24          25          26          27          28          29          30          31          32          33          34          35          36          37          38          39          40          41          42          43          44          45          46          47          48          49          50          51          52          53\n');
fprintf(dataExtraction,'META_STOP\n');
fclose(dataExtraction);
%-----------------------------------------


observd = double(zeros(Nobs,2)); % observd = [out_lat1 out_long1];
predict = double(zeros(Nobs,3)); % predict = [out_lat0 out_long0 out_distance0]
for ii = max([1 simLims(1)]) : simLims(2) : min([NbLT1 simLims(3)])
  epochs  = slctEpochs(Nobs, TimeList1(ii), T, dt);
  % -long1/-long0 avec les donnees reelles (bug!)
  observd = extractObs(epochs, nbofBodies, NbLE1, TimeListE1, lat1, -long1);
  predict = prepareObs(epochs, nbofBodies, NbLE0, TimeListE0, lat0, -long0, dist0);
  % +long1/+long0 avec les donnees mod�les (bug!)
%   observd = extractObs(epochs, nbofBodies, NbLE1, TimeListE1, lat1, long1);
%   predict = prepareObs(epochs, nbofBodies, NbLE0, TimeListE0, lat0, long0, dist0);
  Xexp = expectedOD (TimeList0, NbLE0, TimeListE0, dist0, coord0, vel0, ...
                     epochs, nbofBodies, ...
                     TimeList1, NbLE1, TimeListE1, dist1,coord1,vel1);
  % Transversal and Longitudinal shift in the trajectories as computed and as expected
  %==> as expected:
  unitvvector = unit_speed_vector(ii,vel1);
  trans_traj=norm(cross(Xexp(1:3), unitvvector));
  long_traj=dot(Xexp(1:3), unitvvector);

  post_sigma  = 0.1;
  delta_sigma = 1.;
  nbTries=1; nbCycle=500;
  metaX = double(zeros(1,19));
  cumul_CPU = 0.;
  %notStabilized = true;
  notStabilized = false;
  %while or(delta_sigma/post_sigma > 0.03, notStabilized)
  while or(delta_sigma/post_sigma > 0.10, notStabilized)
    % the Monte-Carlo analysis is run until sigma has been stabilized at 1% accuracy
    %notStabilized = (delta_sigma/post_sigma > 0.01); % previous value
    err_obs = normrnd(0., sigma_obs/3600., [nbCycle 2*Nobs]);
    err_obs = [err_obs(:,1) err_obs(:,2)./cosd(observd(1,1)) ...
               err_obs(:,3) err_obs(:,4)./cosd(observd(2,1)) ...
               err_obs(:,5) err_obs(:,6)./cosd(observd(3,1)) ...
               err_obs(:,7) err_obs(:,8)./cosd(observd(4,1)) ];
    for nm=1:nbCycle
        nbTries = nbTries +1;
        measured = observd + [err_obs(nm,1:2); err_obs(nm,3:4); err_obs(nm,5:6); err_obs(nm,7:8)];
        [X,A,B,elapsed_time] = computeSolution(epochs, measured, predict, algo);
        if (nbTries==1)
            metaX = X';
        else
            metaX = [metaX; X'];
        end
        cumul_CPU = cumul_CPU + elapsed_time;
    end
    pre_sigma   = post_sigma;
    post_sigma  = sum(std(metaX(:,1:12)));
    delta_sigma = abs(post_sigma - pre_sigma);
  end
  X = mean(metaX);
  stdX = std(metaX);
  elapsed_time = cumul_CPU./nbTries;
  
  % Extraction of the day in MJD (integer):
	day=fix(TimeList1(ii)-MJD_0);
  % Extraction of the seconds
	sec=mod(TimeList1(ii)-MJD_0,1)*SEC_0;
	
  % Transversal and Longitudinal shift in the trajectories as computed and as expected
  %==> as computed:
	trans_err=norm(cross(X(1:3)', unitvvector));
	long_err=dot(X(1:3)', unitvvector);

	% Latitude and Longitude differences of seeing the 1st foreground bodies
  % ... to be adapted in the use of several foreground bodies
  l0   = interp1(TimeListE0(1,:), lat0(1,:),  epochs(1), 'linear');
  L0   = interp1(TimeListE0(1,:), long0(1,:), epochs(1), 'linear');
	lat_angle=3600.*(lat1(ii)-l0);
	long_angle=3600.*(long1(ii)-L0);

	% Colinearity and Intensity of the velocity during OD
  vIni(1)   = interp1(TimeList0, vel0(:,1), epochs(1), 'linear');
  vIni(2)   = interp1(TimeList0, vel0(:,2), epochs(1), 'linear');
  vIni(3)   = interp1(TimeList0, vel0(:,3), epochs(1), 'linear');
  vFinal(1) = interp1(TimeList0, vel0(:,1), epochs(Nobs), 'linear');
  vFinal(2) = interp1(TimeList0, vel0(:,2), epochs(Nobs), 'linear');
  vFinal(3) = interp1(TimeList0, vel0(:,3), epochs(Nobs), 'linear');
  speed_angle = (180.*3600./pi)*acos( dot(vIni, vFinal) / (norm(vIni)*norm(vFinal)) );
	norme=1000.*(norm(vFinal)-norm(vIni)); % m/s
	
	% Rank of the A matrix
	A_rank=rank(A);
	
	% 3 next dates that were considered with the current epoch furing th OD
	nxdates=epochs(2:4)-MJD_0;
	
	% CPU time
	CPU=elapsed_time;
  
  % concatenation
  data = [day, sec, trans_traj, long_traj, lat_angle, long_angle, trans_err, long_err, speed_angle, norme, ...
          nbTries, nxdates, CPU, ...
          X, ...
          Xexp', ...
          stdX];
  %--------- output writing ----------------
  dataExtraction=fopen(outputs,'a');
  fprintf(dataExtraction, ['%0.3d %12.4f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %7.1f %9.2f' ...
                        ' %5d %11.4f %11.4f %11.4f %15.3e' ...
                        ' %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f' ...
                        ' %11.3f %11.3f %11.3f %11.3f %11.6f %11.6f %11.6f' ...
                        ' %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f' ...
                        ' %11.3f %11.3f %11.3f %11.3f %11.6f %11.6f %11.6f' ...
                        ' %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f' ...
                        ' %11.3f %11.3f %11.3f %11.3f %11.6f %11.6f %11.6f\n'], data');
  fclose(dataExtraction);
  %-----------------------------------------
  fprintf('%d => %d tries\n', ii, nbTries);
end
