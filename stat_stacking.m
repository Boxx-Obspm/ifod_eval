% Version v1.0, author B.Segret
% outputs_bin = ''; % filename must be provided
outputs_i1 = '../cas_EME/outs/Eb+YECMJE_01as,300MCx192KF';
outputs_i2 = '../cas_EME/outs/Eb+YECMJE_01as,100MCx192KF,1391998_bin';

outputs_o  = '../cas_EME/outs/Eb+YECMJE_01as,400MCx192KF'; % aller à 450/500 cycles...
% outputs_o  = '../cas_Y/outs/Y+Y41324_01as,400MCx192KF'; % ok (300 aurait été juste suffisant)
% outputs_o  = '../cas_EME/outs/E2M-58122Y_E-1,200h,a-h'; % sigmas ok (a-d); mean (a-h +1 pending!!)
% outputs_o  = '../cas_EME/outs/E2M-58122Y_E-1,150h,a-h'; % sigma within tolerances
% outputs_o  = '../cas_EME/outs/E2M-58122Y_E-1,100h,0-7'; % sigma within tolerances
% outputs_o  = '../cas_EME/outs/E2M-58122Y_E-1,050h,a-h'; % sigma within tolerances
% outputs_o  = '../cas_EME/outs/E2M-58122Y_1as,100h,a-h'; % sigma within tolerances
% outputs_o  = '../cas_Y/outs/Y+Y_E-1,200h,1-8'; % means & sigma within tolerances
% outputs_o  = '../cas_Y/outs/Y+Y_E-1,150h,0-7'; % means & sigma within tolerances
% outputs_o  = '../cas_Y/outs/Y+Y_1as,100h,0-6'; % means & sigma within tolerances
% outputs_o  = '../cas_Y/outs/Y+Y_E-1,100h,a-g'; % means & sigma within tolerances
% outputs_o  = '../cas_Y/outs/Y+Y_E-1,050h,0-7'; % within tolerances (~1 or 2 mean values above)
% outputs_o  = '../cas_Y/outs/Y+Y_E-1asec,100h,270t,0-8'; % means & sigma within tolerances
graphs = true;
% graphs = false;

% on peut améliorer en n'imposant plus de stacker sur les mêmes dates:
% - tant qu'on n'est pas a la fin d'1 des 2 fichiers...
% - regarder s'il y a des données sur la même date, si oui alors comme déjà
% - sinon alors on recopie d'abord les data de la première des 2 dates
% - puis on boucle en lisant la date suivante qu'on vient de recopier
% (en fait, l'exception c'est quand la date est la même on fait 1 seule
% recopie. La regle generale c'est qu'on en fait 1 par date et par fichier)

T0N=datenum([2000 1 1 0 0 0]);
T0JD=2451544.5;

% initizes the binary outputs with the number of written columns
fi1 = fopen(outputs_i1,'r');
fi2 = fopen(outputs_i2,'r');
fo  = fopen(outputs_o,'w'); fclose(fo);

lastData = false;
ntStep = 0;
while not(lastData)
% fwrite(fw, [obstime(ik+nKF-2) obstime(ik+nKF)], 'double');
Ti1 = fread(fi1, 2, 'double'); % will provide a 1x2 double array
Ti2 = fread(fi2, 2, 'double'); % will provide a 1x2 double array
%                             (time for KF-measured and for KF-predicted)
if max(abs(Ti1-Ti2))>1/86400. fprintf('Time inconsistency\n'); end
T  = Ti1;
t2=T(2)-T0JD+T0N;
stf=datestr(t2,'yyyy-mm-dd HH:MM:SS');

% fwrite(fw, [Nobs nKF nbCycles (ik+nKF==length(obstime))], 'uint32');
Zi1 = fread(fi1, 4, 'uint32'); % will provide a 1x3 uint32 array
Zi2 = fread(fi2, 4, 'uint32'); % will provide a 1x3 uint32 array
if (Zi1(1:2)~=Zi2(1:2)) fprintf('Inconsistency in the Nb.of points!\n'); end
ntStep = ntStep+1; fprintf('Stacking time-step #%i: %s\n', ntStep, stf);

Z = Zi1;
Nobs=Z(1);
nKF=Z(2);
nbPts=nKF;
nbCycls1=Zi1(3);
nbCycls2=Zi2(3);
Z(3)=nbCycls1+nbCycls2;
lastData = (Z(4)==1);

fo  = fopen(outputs_o,'a');
fwrite(fo, T, 'double');
fwrite(fo, Z, 'uint32');

% solKF = double(zeros(nbCycles, nvKF));
% dtlKF = double(zeros(nbCycles*nbPts, 2));
  
for nC=1:nbCycls1
  % fwrite(fw, ...
  %     [rex rrme rrkf lKg ldP vkf dtrm' dlgm' dtrk' dlgk' mmkf' mmtk' mmlk'], ...
  %     'double');
  rawDATA = fread(fi1, (3*3+10)*(nbPts+1), 'double');
  fwrite(fo, rawDATA, 'double');
end
for nC=1:nbCycls2
  rawDATA = fread(fi2, (3*3+10)*(nbPts+1), 'double');
  fwrite(fo, rawDATA, 'double');
end
fclose(fo);

end
fclose(fi1);
fclose(fi2);

fprintf('Stacked:\n - %s\n - %s\n=> %s\n', outputs_i1, outputs_i2, outputs_o);
fprintf('=> Nb.of time-steps (unchanged): %i\n', ntStep);
fprintf('=> Nb.of Cycles per time-step (1+2): %i\n', nbCycls1+nbCycls2);
stat_reading;
