%%----------------HEADER---------------------------%%
%Author:           Boris Segret
%Version & Date:
%                  V3.0 7-05-2017 (dd-mm-yyyy)
%                  - new I/O format
%                  V2.4 17-03-2017 (dd-mm-yyyy)
%                  - minor - and dirty - adaptations to extract the acceleration
%                  V2.3 08-08-2016 (dd-mm-yyyy)
%                  - adaptations to N = 5 measurement (X19 to X26)
%                  V2.2 31-03-2016 (dd-mm-yyyy)
%                  - use direct epochs values instead of indices
%                  (agregates former "test_interpolation" and "Calculate_Xexpected")
%                  V2.1 05-02-2016 (dd-mm-yyyy), Boris Segret
%                  - *no* call to reference_trajectory.m
%                  - *no* changes of the inputs
%                  until V2 11-09-2015 Oussema SLEIMI & Tristan Mallet
% CL=1
%
% Ideally the ifod should return X19 as given here from the knownledge of both
% the ref. & the actual trajectories (for tests only) interpolated at the dates of observations
%
% 1. Input:
%     <ref_trajectory>
%     epochs    = 4-vector of dates of observations
%     <actual_trajectory>
% 2. Outputs:
%     X26: 26-vector of the expected result for the vector of unknowns (km and km/s)

function X26 = expectedOD(TL0, vXYZ0, XYZ0, TL1, vXYZ1, XYZ1, tt, whatXYZ, method)
% vXYZ(epochs) were sill not computed
% but XYZ are already set for epochs (5x1 vector), like whatXYZ

% function X26 = expectedOD(TimeList0, NbLE0, TimeListE0, dist0, XYZ0, vXYZ0, ...
%                      epochs, bodies, ...
%                      TimeList1, NbLE1, TimeListE1, dist1, XYZ1, vXYZ1)
% tci=toc;
    
% % Interpolation for the reference trajectory
% out_coordinates0(:,1) = interp1(TimeList0, XYZ0(:,1), tt, 'linear');
% out_coordinates0(:,2) = interp1(TimeList0, XYZ0(:,2), tt, 'linear');
% out_coordinates0(:,3) = interp1(TimeList0, XYZ0(:,3), tt, 'linear');
% out_velocity0(:,1)    = interp1(TimeList0, vXYZ0(:,1),   tt, 'linear');
% out_velocity0(:,2)    = interp1(TimeList0, vXYZ0(:,2),   tt, 'linear');
% out_velocity0(:,3)    = interp1(TimeList0, vXYZ0(:,3),   tt, 'linear');
% % Interpolation for the actual trajectory
% out_coordinates1(:,1) = interp1(TimeList1, XYZ1(:,1), tt, 'linear');
% out_coordinates1(:,2) = interp1(TimeList1, XYZ1(:,2), tt, 'linear');
% out_coordinates1(:,3) = interp1(TimeList1, XYZ1(:,3), tt, 'linear');
% out_velocity1(:,1)    = interp1(TimeList1, vXYZ1(:,1),   tt, 'linear');
% out_velocity1(:,2)    = interp1(TimeList1, vXYZ1(:,2),   tt, 'linear');
% out_velocity1(:,3)    = interp1(TimeList1, vXYZ1(:,3),   tt, 'linear');

%Dvectr = out_coordinates0 - out_coordinates1; 
% Dxyz = out_coordinates1 - out_coordinates0; % km
% Dvelocity = (out_velocity1 - out_velocity0);  % km/s
% drho = out_distance1 - out_distance0;

Dxyz = XYZ1 - XYZ0;
vel0 = interp1(TL0, vXYZ0, tt, method);
vel1 = interp1(TL1, vXYZ1, tt, method);
Dvelocity =  mean(vel1 - vel0,1);
vv = whatXYZ - XYZ0; rr0 = sqrt(vv(:,1).^2+vv(:,2).^2+vv(:,3).^2);
vv = whatXYZ - XYZ1; rr1 = sqrt(vv(:,1).^2+vv(:,2).^2+vv(:,3).^2);
drho = rr1-rr0;

% mean difference of acceleration
Nobs=4; % ok, not very nice, I know! (Nobs = 4 is a minimum)
acc0(1:Nobs-1,1:3)=(vel0(2:Nobs,1:3)-vel0(1:Nobs-1,1:3))./...
    repmat((86400.*(tt(2:Nobs)-tt(1:Nobs-1)))',1,3);
acc1(1:Nobs-1,1:3)=(vel1(2:Nobs,1:3)-vel1(1:Nobs-1,1:3))./...
    repmat((86400.*(tt(2:Nobs)-tt(1:Nobs-1)))',1,3);
dacc(1:3) = mean(acc1-acc0,1);

% for ii=1:length(epochs)
%     i=bodies(ii);
%     out_distance0(ii) = interp1(TimeListE0(i,1:NbLE0(i)), dist0(i,1:NbLE0(i)), tt(ii), 'linear');
%     out_distance1(ii) = interp1(TimeListE1(i,1:NbLE1(i)), dist1(i,1:NbLE1(i)), tt(ii), 'linear');
% end

% tcf=toc; fprintf('expectedOD: %5.2f ms, ', (tcf-tci)*1000.);

X26 =[ Dxyz(1,:)'; Dxyz(2,:)'; Dxyz(3,:)'; Dxyz(4,:)'; Dxyz(5,:)'; ...
    drho;
    Dvelocity';
    dacc'];
    %    Dxyz(1,1); Dxyz(1,2); Dxyz(1,3);...
    %    Dxyz(2,1); Dxyz(2,2); Dxyz(2,3);...
    %    Dxyz(3,1); Dxyz(3,2); Dxyz(3,3);...
    %    Dxyz(4,1); Dxyz(4,2); Dxyz(4,3);...
    %    Dxyz(5,1); Dxyz(5,2); Dxyz(5,3);...
    %    drho(1);...
    %    drho(2);...
    %    drho(3);...
    %    drho(4);...
    %    drho(5);...
    %    Dvelocity(1,1); Dvelocity(1,2); Dvelocity(1,3);...
    %    dacc(1); dacc(2); dacc(3)];
% X19 =[Dvectr(1,1);...
%    Dvectr(1,2);...
%    Dvectr(1,3);...
%    Dvectr(2,1);...
%    Dvectr(2,2);...
%    Dvectr(2,3);...
%    Dvectr(3,1);...
%    Dvectr(3,2);...
%    Dvectr(3,3);...
%    Dvectr(4,1);...
%    Dvectr(4,2);...
%    Dvectr(4,3);...
%    dr(1);...
%    dr(2);...
%    dr(3);...
%    dr(4);...
%    Dvelocity(1,1);...
%    Dvelocity(1,2);...
%    Dvelocity(1,3)];

end
