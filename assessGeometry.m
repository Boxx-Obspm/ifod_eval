%%----------------HEADER---------------------------%%
%Author:           Boris Segret
% Version & Date:
%                  V2.1, 10-07-2017
%                  - dynamic optical error is implemented BUT
%                  - 0.1arcsec is not realistic (Earth at 86UA!!)
%                  V2.0, 06-07-2017
%                  - compare ref.&act.Traj. in helio-/planeto-centric
%                  - dynamic optical errors not implemented (depending on
%                  apparent radius over time)
%                  - light-distances still not included (bug!)
%                  V1.0, 28-05-2017, fork from stat_extracting in v4.0 24-04-2017
%                  - Scale factors of foreground bodies (FGB)
%                  - Angular drift of FGB from Ref.Traj. to Act.Traj.
% CL=2 (v2.1)
%
% Plots various geometrical assesments of an IFOD scenario.
% With option "WRITE", writes the results in a VTS-formated file "<scn>_vts.txt"
%
% I/
%    <scenario>
%    <optical errors>
% O/ Plots:
%    - fig.11: drifts of Act.from Ref.Traj.
%    - fig.12: - Angular drifts of FGB from Ref.to Act.Traj.
%              - Scale factors from FGB
%    - fig.13: Radius of spherical-equivalent error for each FGB
%    - fig.14: Error contributions along each FGB's LoS
%    - fig.15+...: Elong.& Longit.errors on each FGB's LoS

clear;
% path_to_scn = '..\cas_EME\'; fscenario   = 'scnEYv5'; % sigma_obs = 0.1;
% path_to_scn = '..\cas_EOs\prepFiles\'; fscenario = 'scn_GTO_wrt_Earth_BAPO'; %sigma_obs = 0.1;
% path_to_scn = '..\cas_EOs\prepFiles\'; fscenario = 'scn_HPO_wrt_Earth_350orbits'; %sigma_obs = 0.1;
path_to_scn = '..\cas_EOs\prepFiles\'; fscenario = 'scn_LEO_wrt_Earth_250orbits'; %sigma_obs = 0.1;

% nKF = 8*24; dtKF = 3600./5.; % (dtKF in seconds)
% nKF = 8; dtKF = 1800./5.; % (dtKF in seconds)

% Nobs = 5;
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
sigmas = double(zeros(1,nbofBodies)); % input errors on each FGB
radii = double(zeros(1,nbofBodies)); % input radii of each FGB
for ii=1:nbofBodies
    % [1:10, Bodyname][11:20, sigma_in][21:..., path&fileName]
    % ex.: 'SUN       180       ./Sun_EME2000_wrt_Earth.xva'
    l='#'; while l=='#' lgn=fgetl(scn); if length(lgn)>0 l=lgn(1); end; end;
    BodyName = [BodyName lgn(1:10)];
    sigmas(ii) = (sscanf(lgn, '%*s %f %*f ')./3600).*pi()./180.; % radians
    % radii(ii)  = sscanf(lgn, '%*s %*f %f '); % km
    radii(ii)  = (1./400).*sscanf(lgn, '%*s %*f %f '); % km, this is to neutralize the criterium on angular radius
    refEphNameLength(ii+1)=refEphNameLength(ii)+length(path_to_scn)+length(lgn)-30;
    refEphemerid = [refEphemerid path_to_scn lgn(31:length(lgn))];
end
fclose(scn);

T0N=datenum([2000 1 1 0 0 0]);
T0JD=2451544.5;
T0MJD= 51544.;

%% Extraction of the trajectory data
% > TimeList0(:), TimeList1(:) in decimal JD, time-stamped of the trajectories ("0" ref & "1" actual)
% > NbLT0, NbLT1: Nb of lines in the trajectory files
% > coord0(:), coord1(:): in km
% > vel0(:), vel1(:): in km/s
[NbLT0, TimeList0, data] = readTraj(reftrajectory, simLims(1), simLims(3));
coord0 = data(:,1:3);
vel0   = data(:,4:6);
fprintf('Reading Ref.Traj''s data: %s, %.1f...%.1f = %i records\n', ...
  reftrajectory, TimeList0(1), TimeList0(NbLT0), NbLT0);
[NbLT1, TimeList1, data] = readTraj(actual_trajectory, simLims(1), simLims(3));
coord1 = data(:,1:3);
vel1   = data(:,4:6);
fprintf('Reading Act.Traj''s data: %s, %.1f...%.1f = %i records\n', ...
  actual_trajectory, TimeList1(1), TimeList1(NbLT1), NbLT1);

%--------- output preparations -----------
t_start = max([TimeList0(1) TimeList1(1) simLims(1)]);
t_stop  = min([TimeList0(NbLT0) TimeList1(NbLT1) simLims(3)]);
starttxt = datestr(t_start-T0MJD+T0N, 'yyyy-mm-ddTHH:MM:SS');
stoptxt  = datestr(t_stop-T0MJD+T0N, 'yyyy-mm-ddTHH:MM:SS');
nbEpochs = 1+floor((t_stop-t_start)/simLims(2));
fprintf('\n*******\nRef. & Act. Trajectories Geometry\n');
fprintf('STARTING = %8.2f (dec.MJD), %s (yyyy-mm-ddTHH:MM:SS)\n', t_start, starttxt);
fprintf('ENDING   = %8.2f (dec.MJD), %s (yyyy-mm-ddTHH:MM:SS)\n', t_stop, stoptxt);
fprintf('STEPS = %i (Nb), %.1f (H)\n', nbEpochs, simLims(2)*24.);
%-----------------------------------------
% Memory allocations
% obstime = double(zeros(1,Nobs*nKF));
% obstime = (t_start+(Nobs*nKF+2)*dtKF/86400.):simLims(2):t_stop;
obstime = t_start:simLims(2):t_stop;

drifts  = double(zeros(length(obstime), 9));
method = 'cubic';
ref_loc = interp1(TimeList0, coord0, obstime, method);
vel0obs = interp1(TimeList0, vel0, obstime, method);
act_loc = interp1(TimeList1, coord1, obstime, method);
vel1obs = interp1(TimeList1, vel1, obstime, method);
% clear coord0 coord1;

v0_unit = vel0obs./repmat(sqrt(vel0obs(:,1).^2+vel0obs(:,2).^2+vel0obs(:,3).^2),[1 3]);
v1_unit = vel1obs./repmat(sqrt(vel1obs(:,1).^2+vel1obs(:,2).^2+vel1obs(:,3).^2),[1 3]);
t0_unit = cross(vel0obs, repmat([0 0 1], [length(obstime) 1])); % transv.axis, perp.to (z,v0) plane
n01_obs = cross(vel0obs, vel1obs); % normal axis, perp.to (v0,v1) plane

drifts(:,1:3) = act_loc - ref_loc;
drifts(:,4)   = sqrt(drifts(:,1).^2+drifts(:,2).^2+drifts(:,3).^2); % dr
drifts(:,6)   = dot(drifts(:,1:3), v0_unit, 2); % longitudinal drift
drifts(:,5)   = sqrt(drifts(:,4).^2-drifts(:,6).^2); % transverse drift
drifts(:,7)   = sqrt(vel1obs(:,1).^2+vel1obs(:,2).^2+vel1obs(:,1).^2)- ...
               sqrt(vel0obs(:,1).^2+vel0obs(:,2).^2+vel0obs(:,1).^2); %delta-V
drifts(:,8)   = 3600.*acosd(dot(v0_unit, v1_unit,2)); %delta-V direction (arcsec)
drifts(:,9)   = real(90+acosd(dot( ...
               t0_unit./repmat(sqrt(t0_unit(:,1).^2+t0_unit(:,2).^2+t0_unit(:,3).^2), [1 3]), ...
               n01_obs./repmat(sqrt(n01_obs(:,1).^2+n01_obs(:,2).^2+n01_obs(:,3).^2), [1 3]), ...
               2))); % rotation angle perpendicularly to trajectory, wrt. z axis

% PLOTS
%-------
% T_A & T_R
figure(1); clf;
subplot(2,2,1);
hold off; plot(ref_loc(:,1)-ref_loc(1,1),ref_loc(:,3)-ref_loc(1,3), '-b'); hold on;
plot(act_loc(:,1)-ref_loc(1,1), act_loc(:,3)-ref_loc(1,3), '--r');
plot(0.,0.,'+k');
set(gca, 'XColor', 'm', 'YColor', 'm', 'ZColor', 'm');% axis equal;
xlabel('X (km)'); ylabel('Z (km)');
title('T_R & T_A in (X,Z) Projection', 'Color', 'm');
subplot(2,2,3);
hold off; plot(ref_loc(:,1)-ref_loc(1,1),ref_loc(:,2)-ref_loc(1,2), '-b'); hold on;
plot(act_loc(:,1)-ref_loc(1,1), act_loc(:,2)-ref_loc(1,2), '--r');
plot(0.,0.,'+k');
set(gca, 'XColor', 'm', 'YColor', 'm', 'ZColor', 'm'); axis equal;
xlabel('X (km)'); ylabel('Y (km)');
title('T_R & T_A in (X,Y) Projection', 'Color', 'm');
subplot(2,2,2);
hold off; plot(ref_loc(:,2)-ref_loc(1,2),ref_loc(:,3)-ref_loc(1,3), '-b'); hold on;
plot(act_loc(:,2)-ref_loc(1,2), act_loc(:,3)-ref_loc(1,3), '--r');
plot(0.,0.,'+k');
set(gca, 'XColor', 'm', 'YColor', 'm', 'ZColor', 'm');% axis equal;
xlabel('Y (km)'); ylabel('Z (km)');
title('T_R & T_A in (Y,Z) Projection', 'Color', 'm');
subplot(2,2,4);
hold off;
plot3(ref_loc(:,1)-ref_loc(1,1),ref_loc(:,2)-ref_loc(1,2),ref_loc(:,3)-ref_loc(1,3), '-b');
hold on;
plot3(act_loc(:,1)-ref_loc(1,1),act_loc(:,2)-ref_loc(1,2),act_loc(:,3)-ref_loc(1,3), '--r');
plot(0.,0.,'+k');
set(gca, 'XColor', 'm', 'YColor', 'm', 'ZColor', 'm'); axis equal;
set(gca, 'XGrid', 'on', 'XMinorGrid', 'off');
set(gca, 'YGrid', 'on', 'YMinorGrid', 'off');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D- T_R & T_A (km)', 'Color', 'm');

% Geometry of T_A wrt T_R
figure(11); clf;

subplot(2,2,1);
hold off; plot(drifts(:,6),drifts(:,5), '-b');
ylabel('OFF axis (km)');
xlabel('<-- delay (km) advance -->');
set(gca, 'XColor', 'm', 'YColor', 'm', 'ZColor', 'm');
set(gca, 'XGrid', 'on', 'XMinorGrid', 'off');
set(gca, 'YGrid', 'on', 'YMinorGrid', 'off');
set(gca, 'ZGrid', 'on', 'ZMinorGrid', 'off');
title('Transv./Longit.Drift wrt T_R', 'Color', 'm');

subplot(2,2,2);
hold off; plot3(drifts(:,1),drifts(:,2),drifts(:,3), '-r'); hold on;
plot3(drifts(:,1),drifts(:,2),zeros(length(obstime),1), '-m');
for obs = 1:length(obstime)
    plot3([drifts(obs,1); drifts(obs,1)], [drifts(obs,2); drifts(obs,2)], [0; drifts(obs,3)], '-g');
end
plot3([min(1.1*[0. min(drifts(:,1))]); max([0. max(drifts(:,1))])], [0.;0.], [0.;0.], '-m');
plot3([0.;0.], 1.1*[min([0. min(drifts(:,2))]); max([0. max(drifts(:,2))])], [0.;0.], '-m');
plot3([0.;0.], [0.;0.], 1.1*[min([0. min(drifts(:,3))]); max([0. max(drifts(:,3))])], '-k');
title('3D-Drift wrt Ref. over time (km)', 'Color', 'm');
set(gca, 'XColor', 'm', 'YColor', 'm', 'ZColor', 'k');
set(gca, 'XGrid', 'off', 'XMinorGrid', 'off');
set(gca, 'YGrid', 'off', 'YMinorGrid', 'off');
set(gca, 'ZGrid', 'on', 'ZMinorGrid', 'off');

subplot(2,2,3);
hold off;
[AX, dVnrm, dVdir] = plotyy(obstime(:)-obstime(1), 1000.*drifts(:,7), ...
    obstime(:)-obstime(1), drifts(:,8));
hold on;
% help "Using Multiple X- and Y-Axes"
set(AX(1),'YGrid','on'); set(AX(1),'YMinorGrid','off');
set(dVnrm, 'Color', [0 0 1]);  set(AX(1), 'YColor', [0 0 1]); 
set(get(AX(1),'Ylabel'),'String','( V_A - V_R ), m/s');
set(dVdir, 'Color', [0 .5 0]); set(AX(2), 'YColor', [0 .5 0], 'XColor', 'm'); 
set(get(AX(2),'Ylabel'),'String', 'Ang( V_R , V_A ), arcsec');
xlabel('time (days)');
set(gca, 'XColor', 'm');

subplot(2,2,4);
hold off;
% plot(cosd(drifts(:,9)).*(obstime(:)-obstime(1)+1), sind(drifts(:,9)).*(obstime(:)-obstime(1)+1), '-r');
plot(cosd(drifts(:,9)).*log(obstime(:)-obstime(1)+1), ...
    sind(drifts(:,9)).*log(obstime(:)-obstime(1)+1), ...
    '-r', 'LineWidth', 2);
hold on;
durations = log((obstime(length(obstime))-obstime(1)+1)).*[1:2:12]/10;
for i=1:length(durations)-1
    plot(cosd(360.*[0:0.01:1]).*durations(i), sind(360.*[0:0.01:1]).*durations(i), ':k');
end
taxis = durations'*[cosd(-135) sind(-135)]; tlabels = num2str(floor(-1+exp(durations))');
plot(0,0,'ok');plot(0,0,'xk'); text(0,0,' V_R');
plot([taxis(length(taxis),1) 0.3*taxis(1,1)], [taxis(length(taxis),2) 0.3*taxis(1,2)], '-m');
text(taxis(1,1), taxis(1,2), tlabels(1,:), 'Color', 'm');
text(taxis(2,1), taxis(2,2), tlabels(2,:), 'Color', 'm');
text(taxis(6,1), taxis(6,2), tlabels(6,:), 'Color', 'm');
text(taxis(6,1)*1.2, taxis(6,2)*1.2, 'days', 'Color', 'm');
set(gca, 'xtick', []); set(gca, 'ytick', []);
title('Angle of V_A wrt V_R', 'Color', 'm');

clear vel0 vel0obs v0_unit t0_unit vel1 vel1obs v1_unit n01_obs;
% clear drifts;
clear obstime ref_loc act_loc;
           
%% Extraction of foreground body data, ii being the number of a given body:
% > NbLE(ii) => Nb of records in the ephemeris files
% > TimeListE(ii, 1:NbLE(ii)) => in decimal JD, time-stamped of the ephemerides files
% > ephXYZ(ii, 1:NbLE(ii), 1:3) => X(t),Y(t),Z(t) rectangular coordinates of the body
for ii=1:nbofBodies
  bodyFile = refEphemerid(1+refEphNameLength(ii):refEphNameLength(ii+1));
  [nbLgn, timeSteps, data] = readTraj(bodyFile, simLims(1), simLims(3));
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
  fprintf('Reading FGB %i''s data: %s, %.1f...%.1f = %i records\n', ...
      ii, bodyFile, TimeListE(ii,1), TimeListE(ii,nbLgn), nbLgn);
end
% clear tempoSteps tempoEph data TimeListE;

t_start = max([TimeList0(1) TimeList1(1) TimeListE(:,1)' simLims(1)]);
t_stop  = min([TimeList0(NbLT0) TimeList1(NbLT1) max(TimeListE(:,NbLE(:))) simLims(3)]);
starttxt = datestr(t_start-T0MJD+T0N, 'yyyy-mm-ddTHH:MM:SS');
stoptxt  = datestr(t_stop-T0MJD+T0N, 'yyyy-mm-ddTHH:MM:SS');
nbEpochs = 1+floor((t_stop-t_start)/simLims(2));
fprintf('Ready to run simulation from %.1f to %.1f, %.1fH intervals (%i steps)\n', ...
  t_start, t_stop, simLims(2)*24., nbEpochs);
if (nbEpochs <= 0) fprintf(' === WARNING! Nb.Steps < 0 ===\n'); end;
% clear TimeList0 TimeList1;

%--------- output preparations -----------
fprintf('\nOptical geometry with Foreground bodies\n');
% fprintf('Optical accuracy (arcsec): %11.5f\n', sigma_obs);
fprintf('Optical accuracy & Max radius: %5.1f arcsec %11.3f km\n', [sigmas'.*180*3600/pi() radii']');
fprintf('STARTING = %8.2f (dec.MJD), %s (yyyy-mm-ddTHH:MM:SS)\n', t_start, starttxt);
fprintf('ENDING   = %8.2f (dec.MJD), %s (yyyy-mm-ddTHH:MM:SS)\n', t_stop, stoptxt);
fprintf('STEPS = %i (Nb), %.1f (H)\n', nbEpochs, simLims(2)*24.);
%-----------------------------------------
% Memory allocations
% obstime = double(zeros(1,Nobs*nKF));
% obstime = (t_start+(Nobs*nKF+2)*dtKF/86400.):simLims(2):t_stop;
obstime = t_start:simLims(2):t_stop;
ref_loc = interp1(TimeList0, coord0, obstime, method);
act_loc = interp1(TimeList1, coord1, obstime, method);

FG_body = double(zeros(length(obstime), 3, nbofBodies));
observd = double(zeros(length(obstime), 3, nbofBodies));
predict = double(zeros(length(obstime), 3, nbofBodies));
optics  = double(zeros(length(obstime), 3, nbofBodies));
elongs  = double(zeros(length(obstime), nbofBodies, nbofBodies));
transv  = double(zeros(length(obstime), nbofBodies, nbofBodies));
allTransv  = double(zeros(length(obstime), nbofBodies, nbofBodies));
topTransv  = double(zeros(length(obstime), nbofBodies, nbofBodies));
topSphere  = double(zeros(length(obstime), nbofBodies, nbofBodies));
idxTransv  = double(zeros(length(obstime), nbofBodies, nbofBodies));
idxSphere  = double(zeros(length(obstime), nbofBodies, nbofBodies));
sphrError   = double(zeros(length(obstime), nbofBodies, nbofBodies));
           
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
%    optics(:,3,ib) = (sigma_obs./3600.).*(pi()/180.).*optics(:,1,ib);
   % sigmas already in radians
   use2ndScaleFactor = repmat(sigmas(ib), [length(obstime) 1]) <= ...
       radii(ib)./sqrt(optics(:,1,ib).^2 + radii(ib).^2);
   optics(:,3,ib) = ...
       (~use2ndScaleFactor).*sigmas(ib).*optics(:,1,ib) + ...
       (use2ndScaleFactor).*radii(ib)./sqrt(1+(radii(ib)./optics(:,1,ib)).^2);
   % drift angle to the foreground bodies from REFERENCE to ACTUAL location
   optics(:,4,ib) = acosd(dot(...
       predict(:,:,ib)./repmat(optics(:,1,ib),[1 3]), ...
       observd(:,:,ib)./repmat(optics(:,2,ib),[1 3]), ...
       2)); % in degrees
end
% clear FGB ephXYZ TimeListE NbLE;

for ib = 1:nbofBodies
for jb = 1:nbofBodies
   pscalr = dot(...
       predict(:,1:3,ib)./repmat(optics(:,1,ib),[1 3]), ...
       predict(:,1:3,jb)./repmat(optics(:,1,jb),[1 3]), ...
       2); % unitized scalar product, i.e.cos(angle)
   elongs(:,ib,jb) = acosd(pscalr); % in degrees
   % contribution of FGB jb in longitudinal axis of FGB ib, in km
   transv(:,ib,jb) = abs(sqrt(1.-pscalr.^2).*optics(:,3,jb));
%    % contribution of FGB jb in longitudinal axis of FGB ib, in Nb.of units
%    % of ib-scale factor
%    transv(:,ib,jb) = abs(pscalr.*optics(:,3,jb)./optics(:,3,ib));
end
end
% clear pscalr;

for ib = 1:nbofBodies
   % what are the best contributors in transv.
    for jb = 1:nbofBodies
        sphrError(:,ib,jb) = ((optics(:,3,ib).^2).*transv(:,ib,jb).*3./(4.*pi())).^(1/3.); % en km^3
    end
   [topSphere(:,ib,:) idxSphere(:,ib,:)] = sort(sphrError(:,ib,:),3);
   [topTransv(:,ib,:) idxTransv(:,ib,:)] = sort(transv(:,ib,:),3);
end
    
% OPTICAL GEOMETRY PLOTS
%-------
% myColors=[0.5 0   0;...
%           0.5 0.5 0;...
%           0   0.5 0.5;...
%           0   0   0.5;...
%           0.4 0.4 0.4];
myColors = lines;
strLegend = [];
for ib=1:nbofBodies
    strLegend = [strLegend; BodyName((ib-1)*10+1:ib*10)];
end

figure(12); clf; %set(gca, 'LooseInset', [0.09 0.09 0.07 0.07]);
% Optical geometry with FBG
% -------------------------
subplot(2,1,1); hold off;
subplot(2,1,2); hold off;
for ib = 1:nbofBodies
    % DRIFT ANGLES from Ref
    subplot(2,1,1);
    yvalues = real(reshape(optics(:,4,ib), [1 length(obstime)]));
    semilogy(obstime-t_start, yvalues, 'Color', myColors(ib,:)); hold on;
    % DRIFT DISTANCES and SCALE FACTORS
    subplot(2,1,2);
    yvalues= real(reshape(optics(:,3,ib), [1 length(obstime)]));
    semilogy(obstime-t_start, yvalues, '-', 'Color', myColors(ib,:)); hold on;
end
subplot(2,1,1);
legend([repmat('to ', [size(strLegend,1) 1]) strLegend], 'Location', 'Best');
ylabel('Ang.drift (degrees)');
xlabel('time (days)');
set(gca, 'XColor', 'm', 'YColor', 'm', 'YGrid', 'on', 'YMinorGrid', 'off');
% title(['Drift angles, ' datestr(t_start,'yyyy-mm-dd') '..' datestr(t_stop,'yyyy-mm-dd')]);
subplot(2,1,2);
legend([repmat('sc.f.to ', [size(strLegend,1) 1]) strLegend], 'Location', 'Best');
ylabel('(kilometers)');
xlabel('time (days)');
set(gca, 'XColor', 'm', 'YColor', 'm', 'YGrid', 'on', 'YMinorGrid', 'off');
title(['Scale factors, ' datestr(T0N+t_start-T0MJD,'yyyy-mm-dd') '..' ...
    datestr(T0N+t_stop-T0MJD,'yyyy-mm-dd')]);

% Synthesis: smallest Sphere-equivalents
% -------------------------------------
figure(13); clf; %set(gca, 'LooseInset', [0.09 0.09 0.07 0.07]);
for ib = 1:nbofBodies
    subplot(nbofBodies,2,2*ib-1);
    for rk =4:-1:2
      yvalues = real(reshape(-idxSphere(:,ib,rk), [1 length(obstime)]));
      plot(obstime(:)-obstime(1), yvalues, 'LineWidth', 1.5, 'Color', myColors(rk,:)); hold on;
    end
    for rk =5:nbofBodies
      yvalues = real(reshape(-idxSphere(:,ib,rk), [1 length(obstime)]));
      plot(obstime(:)-obstime(1), yvalues, 'Color', myColors(rk,:)); hold on;
    end
    set(gca, 'ytick', []); ylim([-nbofBodies-1 0]);
    for jb =1:nbofBodies
      if (jb ~= ib)
        text((obstime(length(obstime))-obstime(1))*1.005, -jb, BodyName((jb-1)*10+1:jb*10));
      end
    end
    xlim([0 obstime(length(obstime))-obstime(1)]); % xlabel('days');
    text((obstime(length(obstime))-obstime(1))*1.005, 0, ...
        BodyName((ib-1)*10+1:ib*10), 'EdgeColor', 'k');
end
legend('3rd', '2nd', '1st  best FGB', 'Location', 'Best');
for ib = 1:nbofBodies
    subplot(nbofBodies,2,2*ib);
    for rk =4:-1:2
      semilogy(obstime(:)-obstime(1), ...
          max([repmat(1E-3,[length(obstime) 1]) topSphere(:,ib,rk)],[],2), ...
          'LineWidth', 1.5, 'Color', myColors(rk,:)); hold on;
    end
    for rk =5:nbofBodies
      semilogy(obstime(:)-obstime(1), ...
          max([repmat(1E-3,[length(obstime) 1]) topSphere(:,ib,rk)],[],2), ...
          'Color', myColors(rk,:)); hold on;
    end
    ylabel('(km)'); %set(gca, 'ytick', [0.01 0.1 1 10 100]); ylim([1E-3 1E3]);
    set(gca, 'YGrid', 'on', 'YMinorGrid', 'on');
end
subplot(nbofBodies,2,1);
text((obstime(length(obstime))-obstime(1))*1.5, 1.2, ...
    'Radius of spherical-equivalent error (km)');

% Synthesis: smallest Errors on FGB's LoS
% ---------------------------------------
figure(14); clf; %set(gca, 'LooseInset', [0.09 0.09 0.07 0.07]);
for ib = 1:nbofBodies
    subplot(nbofBodies,2,2*ib-1);
    for rk =4:-1:2
      yvalues = real(reshape(-idxTransv(:,ib,rk), [1 length(obstime)]));
      plot(obstime(:)-obstime(1), yvalues, 'LineWidth', 1.5, 'Color', myColors(rk,:)); hold on;
    end
    for rk =5:nbofBodies
      yvalues = real(reshape(-idxTransv(:,ib,rk), [1 length(obstime)]));
      plot(obstime(:)-obstime(1), yvalues, 'Color', myColors(rk,:)); hold on;
    end
    set(gca, 'ytick', []); ylim([-nbofBodies-1 0]);
    for jb =1:nbofBodies
      if (jb ~= ib)
        text((obstime(length(obstime))-obstime(1))*1.005, -jb, BodyName((jb-1)*10+1:jb*10));
      end
    end
    xlim([0 obstime(length(obstime))-obstime(1)]); % xlabel('days');
    text((obstime(length(obstime))-obstime(1))*1.005, 0, ...
        BodyName((ib-1)*10+1:ib*10), 'EdgeColor', 'k');
end
legend('3rd', '2nd', '1st  best FGB', 'Location', 'Best');
for ib = 1:nbofBodies
    subplot(nbofBodies,2,2*ib);
    for rk =4:-1:2
      semilogy(obstime(:)-obstime(1), ...
          max([repmat(1E-3,[length(obstime) 1]) topTransv(:,ib,rk)],[],2), ...
          'LineWidth', 1.5, 'Color', myColors(rk,:)); hold on;
    end
    for rk =5:nbofBodies
      semilogy(obstime(:)-obstime(1), ...
          max([repmat(1E-3,[length(obstime) 1]) topTransv(:,ib,rk)],[],2), ...
          'Color', myColors(rk,:)); hold on;
    end
    ylabel('(km)'); %set(gca, 'ytick', [0.01 0.1 1 10 100]); ylim([1E-3 1E3]);
    set(gca, 'YGrid', 'on', 'YMinorGrid', 'on');
end
subplot(nbofBodies,2,1);
text((obstime(length(obstime))-obstime(1))*1.5, 1.2, ...
    'Error contributions along FGB''s LoS (km)');

% Elongations and Scale factors
% -----------------------------
for ib = 1:nbofBodies
figure(15+ib-1); clf; %set(gca, 'LooseInset', [0.09 0.09 0.07 0.07]);
    for jb = 1:nbofBodies
        subplot(1,4,1);
        yvalues = real(reshape(elongs(:,ib,jb), [1 length(obstime)]));
        plot(obstime(:)-obstime(1), yvalues, 'Color', myColors(jb,:)); hold on;
        subplot(1,4,[2 4]);
        semilogy(obstime(:)-obstime(1), ...
            max([repmat(1E-3,[length(obstime) 1]) transv(:,ib,jb)],[],2), ...
            'Color', myColors(jb,:)); hold on;
    end
    subplot(1,4,1);
    title(['Elongations wrt FGB ' BodyName((ib-1)*10+1:ib*10)], 'Color', 'm');
    ylabel('degrees');
    set(gca, 'ytick', [0. 45. 90. 135. 180.]);
    xlabel('time (days)');
    set(gca, 'XColor', 'm', 'YColor', 'm', 'YGrid', 'on', 'YMinorGrid', 'off');
    subplot(1,4,[2 4]);
    title(['Longit.errors on LoS to FGB ' BodyName((ib-1)*10+1:ib*10)], 'Color', 'm');
    ylabel('kilometers');
%     ylabel(['Nb.of scale factors for ' BodyName((ib-1)*10+1:ib*10)]);
    xlabel('time (days)');
    set(gca, 'XColor', 'm', 'YColor', 'm', 'YGrid', 'on', 'YMinorGrid', 'off');
    legend(strLegend, 'Location', 'Best');
end


