%%----------------HEADER---------------------------%%
%Author:           Boris Segret
version = '2.0';
%Version & Date:
%                  Vxxx dd-mm-2016 (dd-mm-yyyy)
%                  - create a "prodObs" module (Produce Observable)
%                  - call the solution based on Observables only
%                  - multiple Ephemeride files
%                  - (to be done) RECTIFY header in the results files
%Version & Date:
%                  V1.2 03-03-2016, Boris Segret
%                  - *no* multiple calls of reference_trajectory.m
%                  - new "scenario" format, and computes only in an interval
%                  - minor adapations
%                  - works with a call from "ifod_tests"
%                  until V1   11-09-2015, Tristan Mallet
%CL=2
%
%
% This produces a timestep-by-timestep comparison between the actual trajectory,
% slightly different from a reference trajectory, and the trajectory that is
% reconstructed by the ifod on the solely basis of the (lat,long) of foreground
% bodies.
%
% I/
%    <scenario> scenario file (in VTS format) path must be given in datapath below
%               see User Manual to set the scenario file
%    datapath : path for the "scenario" file, it *must* be provided in the workspace
% O/
%    <results> comparison between actual and reconstructed trajectories
%              (VTS format, prefixed as requested in <scenario>)

% Procedures and function due to be called from the path of "ifod"
addpath('../ifod');
% reference_trajectory.m
% slctEpochs.m
% extractObs.m
% prepareObs.m
% computeSolution.m
% 
%-----------------------------------------


data=[];
MJD_0=2400000.5;
SEC_0=86400;
%for n=0 : 5 scenario=fopen(strcat('./Inputs/Scenarios/scenario-',num2str(n)),'r')
%scenario=fopen('./Inputs/Scenarios/scenario-4','r'); #we open the scenario file
scn=fopen(strcat(datapath, 'scenario'),'r'); %we open the scenario file
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
refEphNameLength(1)=0; refEphemerid='';
for ii=1:nbofBodies
    % works with the last listed body only, yet
    l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end;
    %refephemerid = strcat(datapath, lgn(12:length(lgn)));
    refEphNameLength(ii+1)=refEphNameLength(ii)+length(datapath)+length(lgn)-11;
    refEphemerid = [refEphemerid datapath lgn(12:length(lgn))];
end
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; actual_trajectory = strcat(datapath,lgn);
actEphNameLength(1)=0; actEphemerid='';
for ii=1:nbofBodies
    % works with the last listed body only, yet
    l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end;
    %trajectory_name_ephjup = strcat(datapath, lgn(12:length(lgn)));
    actEphNameLength(ii+1)=actEphNameLength(ii)+length(datapath)+length(lgn)-11;
    actEphemerid = [actEphemerid datapath lgn(12:length(lgn))];
end
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; algo=lgn;
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; sampling = strcat(datapath,lgn);
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; outputs = strcat(datapath,lgn);
l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end; simLims = uint32(str2num(lgn));
fclose(scn);

%Extraction of the trajectory data
% outputs in same units like trajectories and ephemerides files
% > TimeList0, TimeList1: in decimal Julian days
% > lat0, long0, lat1, long1: in DEGREES
% > distance0, distance1: in km
% > coordinates0, coordinates1: in km
% > velocity0, velocity1: in m/s

%[TimeList0,lat0,long0,distance0,coordinates0,velocity0]=reference_trajectory(reftrajectory, refephemerid);
[NbLT0, TimeList0, coord0, vel0] = readTraj(reftrajectory);
%[TimeList1,lat1,long1,distance1,coordinates1,velocity1]=actual_trajectory(trajectory_name,trajectory_name_ephjup);
[NbLT1, TimeList1, coord1, vel1] = readTraj(actual_trajectory);
for ii=1:nbofBodies
  [NbLE0(ii), TimeListE0(ii,:), lat0(ii,:), long0(ii,:), dist0(ii,:)] = readEphem(refEphemerid(1+refEphNameLength(ii):refEphNameLength(ii+1)));
  [NbLE1(ii), TimeListE1(ii,:), lat1(ii,:), long1(ii,:), dist1(ii,:)] = readEphem(actEphemerid(1+actEphNameLength(ii):actEphNameLength(ii+1)));
end

ii_MAX=length(TimeList1);

timeStep=fopen(sampling, 'r'); %we open the file timeStep-i
	l=' ';
	while 1
		l=fgetl(timeStep);
		if strfind(l,'META_STOP')>0
			break;
		end;
	end;
[sampling_data,nb] = fscanf(timeStep,'%g %g %g %g %g', [5 inf]);
sampling_data = sampling_data';
fprintf('Changes in timesteps:\n');
T=[sampling_data(:,1),sampling_data(:,2)]; fprintf('%15i %15i\n', T'); % we pick the sampling data and stock it in T and dt
dt=[sampling_data(:,3),sampling_data(:,4),sampling_data(:,5)]; fprintf('%5i %5i %5i\n', dt');
% Not tolerant to empty lines !!!
fclose(timeStep);

%--------- output preparations -----------
dataExtraction=fopen(outputs,'w');
fprintf(dataExtraction,'NAV_Results Version : %s\nGenerated by BIRDY NAV TEAM\nDate : %s\nalgo : %s\n\n', version, datestr(now), algo);
fprintf(dataExtraction,'OBJECT_NAME : BIRDY\nID_NAME : BIRDY\nSCENARIO : %s\nTRAJECTORY_NAME : %s\n\n',scenario_num, actual_trajectory);
fprintf(dataExtraction,'META_START\n\n');
fprintf(dataExtraction,'Sampling of observations are changed at the following dates :\n');
fprintf(dataExtraction,'%15d%15f\n',T');
fprintf(dataExtraction,'\nSamplings (hours after previous obs) :\n');
fprintf(dataExtraction,'%15d%15d%15d\n',dt');
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
fprintf(dataExtraction,'COLUMN #11 : Rank of A\n');
fprintf(dataExtraction,'COLUMN #12 : date 2 (in seconds)\n');
fprintf(dataExtraction,'COLUMN #13 : date 3 (in seconds)\n');
fprintf(dataExtraction,'COLUMN #14 : date 4 (in seconds)\n');
fprintf(dataExtraction,'COLUMN #15 : performance (in CPU milliseconds)\n');
fprintf(dataExtraction,'COLUMN #16 to #34 : Vector X [(dx,dy,dz)_i,dr_i in km, dvx,dvy,dvz in km/s]\n');
fprintf(dataExtraction,'COLUMN #35 to #53 : Vector X-expected (same units like X)\n\n');
fprintf(dataExtraction, '   1            2           3           4           5           6           7           8      9        10   11         12          13          14              15          16          17          18          19          20          21          22          23          24          25          26          27          28          29          30          31          32          33          34          35          36          37          38          39          40          41          42          43          44          45          46          47          48          49          50          51          52          53\n');
fprintf(dataExtraction,'META_STOP\n');
fclose(dataExtraction);
%-----------------------------------------


Nobs = 4;
observd = double(zeros(Nobs,2)); % observd = [out_lat1 out_long1];
predict = double(zeros(Nobs,3)); % predict = [out_lat0 out_long0 out_distance0]
for ii = max([1 simLims(1)]) : simLims(2) : min([NbLT1 simLims(3)])
  epochs  = slctEpochs(Nobs, TimeList1(ii), T, dt);
  observd = extractObs(epochs, nbofBodies, NbLE1, TimeListE1, lat1, -long1);
  predict = prepareObs(epochs, nbofBodies, NbLE0, TimeListE0, lat0, -long0, dist0);
  [X,A,B,elapsed_time] = computeSolution(epochs, observd, predict, algo);
  
  Xexp = expectedOD (TimeList0, NbLE0, TimeListE0, dist0, coord0, vel0, ...
                     epochs, nbofBodies, ...
                     TimeList1, NbLE1, TimeListE1, dist1,coord1,vel1);

    % Extraction of the day in MJD (integer):
	day=fix(TimeList1(ii)-MJD_0);
	
	% Extraction of the seconds
	sec=mod(TimeList1(ii)-MJD_0,1)*SEC_0;
	
	% Extraction of the transversal difference between the reference and the actual trajectories
    %% ------- ne PAS utiliser ii avec coord0 ! -----
%     X0(1)   = interp1(TimeList0, coord0(:,1), epochs(1), 'linear');
%     X0(2)   = interp1(TimeList0, coord0(:,2), epochs(1), 'linear');
%     X0(3)   = interp1(TimeList0, coord0(:,3), epochs(1), 'linear');
  unitvvector = unit_speed_vector(ii,vel1);
	trans_traj=norm(cross(Xexp(1:3), unitvvector));
	% Extraction of the longitunal difference between the reference and the actual trajectories
	long_traj=dot(Xexp(1:3), unitvvector);
	% Extraction of transversal error of OD wrt actual trajectory :
	trans_err=norm(cross(X(1:3), unitvvector));
	% Extraction of longitudinal error of OD wrt actual trajectory :
	long_err=dot(X(1:3), unitvvector);

	% Extraction of the latitude and the longitude differences of seeing Jupiter from the reference and actual trajectories
    l0   = interp1(TimeListE0(1,:), lat0(1,:),  epochs(1), 'linear');
    L0   = interp1(TimeListE0(1,:), long0(1,:), epochs(1), 'linear');
	lat_angle=3600.*(lat1(ii)-l0);
	long_angle=3600.*(long1(ii)-L0);

	% Extraction of colinearity of the speed during OD
    vIni(1)   = interp1(TimeList0, vel0(:,1), epochs(1), 'linear');
    vIni(2)   = interp1(TimeList0, vel0(:,2), epochs(1), 'linear');
    vIni(3)   = interp1(TimeList0, vel0(:,3), epochs(1), 'linear');
    vFinal(1) = interp1(TimeList0, vel0(:,1), epochs(Nobs), 'linear');
    vFinal(2) = interp1(TimeList0, vel0(:,2), epochs(Nobs), 'linear');
    vFinal(3) = interp1(TimeList0, vel0(:,3), epochs(Nobs), 'linear');
    %speed_angle=atan2(norm(cross(velocity1(ii,:), velocity1(ii+timeStep(1)+timeStep(2)+timeStep(3),:))), dot(velocity1(ii,:), velocity1(ii+timeStep(1)+timeStep(2)+timeStep(3),:)))*180/pi;
    %% ------- ne PAS utiliser ii avec coord0 ! -----
    speed_angle = (180.*3600./pi)*acos( dot(vIni, vFinal) / (norm(vIni)*norm(vFinal)) );
	% Extraction of uniformity of the speed during OD (in mm/s from m/s=> factor 1000.)
	norme=1000.*(norm(vFinal)-norm(vIni)); % m/s
	
	% Extraction of the rank of the A matrix
	A_rank=rank(A);
	
	% Extraction of the 3 next dates of the pictures
	%nxdates=[TimeList1(ii+timeStep(1))-MJD_0,TimeList1(ii+timeStep(1)+timeStep(2))-MJD_0,TimeList1(ii+timeStep(1)+timeStep(2)+timeStep(3))-MJD_0];
	nxdates=epochs(2:4)-MJD_0;
	
	% Extraction of the CPU time
	CPU=elapsed_time;
  
  %data=vertcat(data, [day,sec,trans_traj,long_traj,lat_angle,long_angle,trans_err,long_err,speed_angle,norme,A_rank,nxdates,CPU,X']); %we create the matrix using concatenation
  data=[day,sec,trans_traj,long_traj,lat_angle,long_angle,trans_err,long_err,speed_angle,norme,A_rank,nxdates,CPU,X',Xexp']; %we create the matrix using concatenation
  %--------- output writing ----------------
  dataExtraction=fopen(outputs,'a');
  fprintf(dataExtraction, ['%0.3d %12.4f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %6.1f %9.2f' ...
                        ' %3d %11.4f %11.4f %11.4f %15.3e' ...
                        ' %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f' ...
                        ' %11.3f %11.3f %11.3f %11.3f %11.6f %11.6f %11.6f' ...
                        ' %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f' ...
                        ' %11.3f %11.3f %11.3f %11.3f %11.6f %11.6f %11.6f\n'], data');
  fclose(dataExtraction);
  %-----------------------------------------

end

