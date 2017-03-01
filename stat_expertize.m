figure(104); clf;
% details of the 3D-OD solutions
subplot(2,1,1); hold on;
mo=mean(solMd(:,2));
sd=std(solMd(:,2));
x=mo-3*sd:0.3*sd:mo+3*sd;
hist(solMd(:,2),x); plot([mo mo],[0 20],'r-'); xlim([mo-3*sd mo+3*sd]);
fprintf('3D-measures=%i+/-%i, n=%i ... ', floor(mo), floor(sd), length(solMd(:,2)));
% details of the KF solution => solKF(:,7-9) & solKF(:,15-16)
subplot(2,1,2); hold on;
mo=mean(solKF(:,8));
sd=std(solKF(:,8));
x=mo-3*sd:0.3*sd:mo+3*sd;
hist(solKF(:,8),x); plot([mo mo],[0 20],'r-'); xlim([mo-3*sd mo+3*sd]);
fprintf('KF-solution=%i+/-%i, n=%i\n', floor(mo), floor(sd), length(solKF(:,8)));

%     plot(rrkf(:,1), 'r-');
%     plot(rrkf(:,2), 'g-');
%     plot(rrkf(:,3), 'b-');
%     plot(rrme(:,1), 'rx:');
%     plot(rrme(:,2), 'gx:');
%     plot(rrme(:,3), 'bx:');
% ylim([-1000 1000]);
