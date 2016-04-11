%%----------------HEADER---------------------------%%
%Author:           Boris Segret
%Version & Date:
%                  Vxx dd-mm-2016 (dd-mm-yyyy)
%                  - use direct epochs values instead of indices
%                  - agregates former "test_interpolation" and "Calculate_Xexpected"
%                  V2.1 05-02-2016 (dd-mm-yyyy), Boris Segret
%                  - *no* call to reference_trajectory.m
%                  - *no* changes of the inputs
%                  until V2 11-09-2015 Oussema SLEIMI & Tristan Mallet
% CL=1
%
% Ideally the inflightOD shall return X19, given here from the ref. & the actual
% trajectories (for tests only) interpolated at the dates of observations
%
% 1. Input:
%     <ref_trajectory>
%     epochs    = 4-vector of dates of observations
%     <actual_trajectory>
% 2. Outputs:
%     X19: 19-vector of expected result for the vector of unknowns (km and km/s)
%     (note: velocities are input in m/s and output in km/s)

function X19 = expectedOD(TimeList0, NbLE0, TimeListE0, dist0, coord0, vel0, ...
                     epochs, nbofBodies, ...
                     TimeList1, NbLE1, TimeListE1, dist1, coord1, vel1)

% Interpolation for the reference trajectory
out_coordinates0(:,1) = interp1(TimeList0, coord0(:,1), epochs, 'linear');
out_coordinates0(:,2) = interp1(TimeList0, coord0(:,2), epochs, 'linear');
out_coordinates0(:,3) = interp1(TimeList0, coord0(:,3), epochs, 'linear');
out_velocity0(:,1)    = interp1(TimeList0, vel0(:,1),   epochs, 'linear');
out_velocity0(:,2)    = interp1(TimeList0, vel0(:,2),   epochs, 'linear');
out_velocity0(:,3)    = interp1(TimeList0, vel0(:,3),   epochs, 'linear');
% Interpolation for the actual trajectory
out_coordinates1(:,1) = interp1(TimeList1, coord1(:,1), epochs, 'linear');
out_coordinates1(:,2) = interp1(TimeList1, coord1(:,2), epochs, 'linear');
out_coordinates1(:,3) = interp1(TimeList1, coord1(:,3), epochs, 'linear');
out_velocity1(:,1)    = interp1(TimeList1, vel1(:,1),   epochs, 'linear');
out_velocity1(:,2)    = interp1(TimeList1, vel1(:,2),   epochs, 'linear');
out_velocity1(:,3)    = interp1(TimeList1, vel1(:,3),   epochs, 'linear');

for ii=1:length(epochs)
  i=1+mod(ii-1,nbofBodies);
    out_distance0(ii) = interp1(TimeListE0(i,1:NbLE0(i)), dist0(i,1:NbLE0(i)), epochs(ii), 'linear');
    out_distance1(ii) = interp1(TimeListE1(i,1:NbLE1(i)), dist1(i,1:NbLE1(i)), epochs(ii), 'linear');
end

%Dvectr = out_coordinates0 - out_coordinates1; 
Dvectr = out_coordinates1 - out_coordinates0; % km
Dvelocity = (out_velocity1 - out_velocity0);  % km/s
dr = out_distance1 - out_distance0;

X19 =[Dvectr(1,1);...
   Dvectr(1,2);...
   Dvectr(1,3);...
   Dvectr(2,1);...
   Dvectr(2,2);...
   Dvectr(2,3);...
   Dvectr(3,1);...
   Dvectr(3,2);...
   Dvectr(3,3);...
   Dvectr(4,1);...
   Dvectr(4,2);...
   Dvectr(4,3);...
   dr(1);...
   dr(2);...
   dr(3);...
   dr(4);...
   Dvelocity(1,1);...
   Dvelocity(1,2);...
   Dvelocity(1,3)];

end
