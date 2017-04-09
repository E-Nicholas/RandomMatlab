% nice looking odd symmetric
b = gaussmf(x,[pi/1.8 pi-pi/8]);
f = 1.5; 
wave = sin(x*f+3*pi/4);
figure(1); plot(x,wave); hold on; plot(x,b); set(gca,'Ytick',[0],'Xlim',[-0.05 2*pi+0.05]);
OS_1_gabr = b.*wave;

figure(2)
plot(d,OS_1_gabr); set(gca,'Ytick',[0],'Xlim', [-1.01 1.01]);
                                              %area under lobes of gabor 
OS_1_on = trapz(abs(OS_1_gabr(OS_1_gabr > 0)))   %positive lobe
OS_1_off = trapz(abs(OS_1_gabr(OS_1_gabr < 0)))  %negative lobes

%1d Gabor difference of gaussians
x=linspace(0, 2*pi, 2000);
b = gaussmf(x,[pi/7 pi*1.2]);
z = gaussmf(x,[pi/6 pi/2]);
sf = 1;
wave = sin(x*sf); %+pi to phase shift

zz=b-z;
%plot(d,wave); hold on;
plot(d,b); hold on; plot(d,z); plot(d,zz);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Making 1D gabors: example for method (and creating x and d, filter space and RF space, respectively)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scrsz = get(0,'ScreenSize');

%%%%
RF_size = 2;             % 2 degree SimpleCell RF's in both Q and H&W model
%%%%
points = RF_size * 1000; % number of data points
x=linspace(0, 2*pi, points); % Filter space, from 0 to 2pi with n points

b = gaussmf(x,[pi/2 pi]);    % Create gaussian in this space, gaussmf(f,[sigma u])
f = 1;                       % Set frequency
wave = sin(x*f);        % Create wave in this space, add pi for 180 phase

%We're supposed to be thinking in the dimension of space, so d is RF 'space'
%our units are degrees...
d = linspace(0-RF_size/2,0+RF_size/2,points);   % d is RF space, 0 is center of RF
dd = 1/length(d);            % our resolution
RF = length(d);              % size of RF in dd

%Plot example wave and gabor functions
figure(1)
plot(x,wave,'r'); hold on; plot(x,b); set(gca,'Ytick',[0],'Xlim',[-0.05 2*pi+0.05]);

gabr = b.*wave;       %Multiply wave by gaussian for gabor

% Plot resulting gabor against x and d spaces (equivalent spaces...)
figure('Name','Plotting Gabor against x and d','NumberTitle','off')
set(gcf, 'Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/4 scrsz(4)/2]); 
subplot(2,1,1);
plot(x,gabr); hold on;
set(gca,'Ytick',[0],'Xlim', [-0.05 2*pi+0.05]);
title('Example Gabor')
subplot(2,1,2);
plot(d,gabr);
set(gca,'Ytick',[0],'Xlim', [-1.0 1.0]);


%%

% Simple Cell filters for H&W model (HW model sums filters)

%The classic even on center (EOC) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = gaussmf(x,[pi/1.8 pi]);
f = 1.5; 
wave = sin(x*f+pi);
EOC_gabr = b.*wave;
%%%%%%
figure('Name','Even on Center filter','NumberTitle','off')
set(gcf, 'Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3)/1.5 scrsz(4)/2]); 
subplot(1,2,1)
plot(x,wave,'r','linewidth',2); hold on; plot(x,b,'linewidth',2); set(gca,'Ytick',[0],'Xlim',[-0.05 2*pi+0.05],'Ylim', [-1 1]);
title('Sine wave and gaussian')
subplot(1,2,2)
plot(d,EOC_gabr,'linewidth',2); set(gca,'Ytick',[0 1],'Xlim', [-1.01 1.01], 'Ylim', [-1 1]); xlabel('RF space (degrees');
title('Resulting Product = gabor');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                              % calculating area under lobes...
                                              % of gabor to ensure roughly equal 
EOC_on = trapz(abs(EOC_gabr(EOC_gabr > 0)))   % positive lobe
EOC_off = trapz(abs(EOC_gabr(EOC_gabr < 0)))  % negative lobes

%The classic even off center (EFC) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = gaussmf(x,[pi/1.8 pi]);
f = 1.5; 
wave = sin(x*f);
EFC_gabr = b.*wave;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Even off Center filter','NumberTitle','off')
plot(d,EFC_gabr); set(gca,'Ytick',[0],'Xlim', [-1.01 1.01]); xlabel('RF (degrees)');
title('Even off center (EFC)')

EFC_on = trapz(abs(EFC_gabr(EFC_gabr > 0)))   %positive lobes
EFC_off = trapz(abs(EFC_gabr(EFC_gabr < 0)))  %negative lobe

%The odd symmetric right (OSR) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = gaussmf(x,[pi/4 pi]);
f = 1; 
wave = sin(x*f+pi);
OSR_gabr = b.*wave;
c = 1/max(OSR_gabr);  %These filters get scaled up to so thier max value = 1
OSR_gabr = OSR_gabr*c; %This wasn't neccessary for the even functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','Odd symmetric right-sided','NumberTitle','off')
set(gcf, 'Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3)/1.5 scrsz(4)/2]);
subplot(1,2,1)
plot(x,wave,'r'); hold on; plot(x,b); set(gca,'Ytick',[0],'Xlim',[-0.05 2*pi+0.05]);
title('Sine wave and gaussian')
subplot(1,2,2)
plot(d,OSR_gabr); set(gca,'Ytick',[0],'Xlim', [-1.01 1.01],'Ylim',[-1 1]); xlabel('RF (degrees)');
title('Resulting gabor')

OSR_on = trapz(abs(OSR_gabr(OSR_gabr > 0)))   %positive lobes
OSR_off = trapz(abs(OSR_gabr(OSR_gabr < 0)))  %negative lobe

%Odd symmetric left (OSL) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = gaussmf(x,[pi/4 pi]);
f = 1; 
wave = sin(x*f);
OSL_gabr = b.*wave;
c = 1/max(OSL_gabr);
OSL_gabr = OSL_gabr*c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Odd symmetric left-sided','NumberTitle','off')
plot(d,OSL_gabr); set(gca,'Ytick',[0],'Xlim', [-1.01 1.01]); xlabel('RF (degrees)');
title('Odd symmetric left (OSL)')

OSL_on = trapz(abs(OSL_gabr(OSL_gabr > 0)))   %positive lobes
OSL_off = trapz(abs(OSL_gabr(OSL_gabr < 0)))  %negative lobe

% Plot all H&W RF types%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('name','HW model filters','NumberTitle','off')
subplot(2,2,1)
plot(d,EOC_gabr); hold on; set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1]);
title('Even on center (EOC)');
subplot(2,2,2)
plot(d,EFC_gabr); hold on; set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1]);
title('Even off center (EFC)'); 
subplot(2,2,3)
plot(d,OSR_gabr); hold on; set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1]);
title('Odd symmetric right (OSR)'); xlabel('RF (degrees');
subplot(2,2,4)
plot(d,OSL_gabr); hold on; set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1]);
title('Odd symmetric left (OSL)'); xlabel('RF (degrees');

% Creating simple cell filters for Quadriture model (Q model sums squared responses)

% Even on center shift_1 (-45 deg from center, left) (EOCs1) %%%%%%%%%%%%%
b = gaussmf(x,[pi/1.8 pi-pi/6]);
f = 1.5; 
wave = sin(x*f+pi/4+pi);
EOCs1_gabr = b.*wave;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sm = find(wave(800:1000)==max(wave)); % me ensuring functions are centered
% gm = find(b==max(b));

figure('Name','Even on center -45 degree shift','NumberTitle','off');
set(gcf, 'Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3)/1.5 scrsz(4)/2]);
subplot(1,2,1)
plot(x,wave,'r'); hold on; plot(x,b); set(gca,'Ytick',[0],'Xlim',[-0.05 2*pi+0.05]); title('wave and guassian');
subplot(1,2,2)
plot(d,EOCs1_gabr); set(gca,'Ytick',[0],'Xlim', [-1.01 1.01]); title('Resulting gabor'); xlabel('RF (degrees)');
                                              %area under lobes of gabor 
EOCs1_on = trapz(abs(EOCs1_gabr(EOCs1_gabr > 0)))   %positive lobe
EOCs1_off = trapz(abs(EOCs1_gabr(EOCs1_gabr < 0)))  %negative lobes

% Even on center shift_2 (+45 deg from center, right) (EOCs2) %%%%%%%%%%%%
b = gaussmf(x,[pi/1.8 pi+pi/6]);
f = 1.5; 
wave = sin(x*f-pi/4+pi);
EOCs2_gabr = b.*wave;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Even on center +45 degree shift','NumberTitle','off');
set(gcf, 'Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3)/1.5 scrsz(4)/2]);
subplot(1,2,1)
plot(x,wave,'r'); hold on; plot(x,b); set(gca,'Ytick',[0],'Xlim',[-0.05 2*pi+0.05]); title('wave and guassian');
subplot(1,2,2)
plot(d,EOCs2_gabr); set(gca,'Ytick',[0],'Xlim', [-1.01 1.01]); xlabel('RF (degrees)'); title('Resulting gabor');
                                              %area under lobes of gabor 
EOCs2_on = trapz(abs(EOCs2_gabr(EOCs2_gabr > 0)))   %positive lobe
EOCs2_off = trapz(abs(EOCs2_gabr(EOCs2_gabr < 0)))  %negative lobes

% Even off center, left shift (EFCs1) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = gaussmf(x,[pi/1.8 pi-pi/6]);
f = 1.5; 
wave = sin(x*f+pi/4);
EFCs1_gabr = b.*wave;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Even off center -45 degree shift','NumberTitle','off'); 
set(gcf, 'Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3)/1.5 scrsz(4)/2]);
subplot(1,2,1)
plot(x,wave,'r'); hold on; plot(x,b); set(gca,'Ytick',[0],'Xlim',[-0.05 2*pi+0.05]); title('wave and guassian');
subplot(1,2,2)
plot(d,EFCs1_gabr); set(gca,'Ytick',[0],'Xlim', [-1.01 1.01]); xlabel('RF (degrees)'); title('Resulting gabor');
                                              %area under lobes of gabor 
EFCs1_on = trapz(abs(EFCs1_gabr(EFCs1_gabr > 0)))   %positive lobe
EFCs1_off = trapz(abs(EFCs1_gabr(EFCs1_gabr < 0)))  %negative lobes

% Even off center, right shift (EFCs2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = gaussmf(x,[pi/1.8 pi+pi/6]);
f = 1.5; 
wave = sin(x*f-pi/4);
EFCs2_gabr = b.*wave;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Even off center +45 degree shift','NumberTitle','off');
set(gcf, 'Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3)/1.5 scrsz(4)/2]);
subplot(1,2,1)
plot(x,wave,'r'); hold on; plot(x,b); set(gca,'Ytick',[0],'Xlim',[-0.05 2*pi+0.05]); title('wave and guassian');
subplot(1,2,2)
plot(d,EFCs2_gabr); set(gca,'Ytick',[0],'Xlim', [-1.01 1.01]); xlabel('RF (degrees)'); title('Resulting gabor');
                                              %area under lobes of gabor 
EFCs2_on = trapz(abs(EFCs2_gabr(EFCs2_gabr > 0)))   %positive lobe
EFCs2_off = trapz(abs(EFCs2_gabr(EFCs2_gabr < 0)))  %negative lobes

% plotting Quadriture filters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
figure('name','Quadriture model filters','NumberTitle','off')
set(gcf, 'Position',[scrsz(3)/6 scrsz(4)/6 scrsz(3)/1.5 scrsz(4)/1.5]);
subplot(2,2,1)
plot(d,EOCs1_gabr); hold on; set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1]);
title('Even on center -45 (EOCs1)');
subplot(2,2,2)
plot(d,EOCs2_gabr); hold on; set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1]);
title('Even on center +45 (EOCs2)');
subplot(2,2,3)
plot(d,EFCs1_gabr); hold on; set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1]);
title('Even on center -45 (EFCs2)');
subplot(2,2,4)
plot(d,EFCs2_gabr); hold on; set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1]);
title('Even on center +45 (EOFs2)');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example Simple cell responses to bars generated by filters (SC responses)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stim = zeros(size(d));        % Creating impulse function
pulse_width = 0.1;            % in degrees. Remember RF width is 2 degrees
pw = pulse_width*1000;
pulse_pos = -0.265;               % in degrees, relative to center (-1.0:1.0)
pp = points/2+1+pulse_pos*1000;
stim(pp:pp+pw-1) = 0.5;           %impulse

EOC_response = dot(EOC_gabr,stim); % Calculating response of RF to stimulus with dot product
if EOC_response <0;                % Rectification, responses below zero are zero
    EOC_response = 0;
end

EFC_response = dot(EFC_gabr,stim);
if EFC_response <0;
    EFC_response = 0;
end

OSR_response = dot(OSR_gabr,stim);
if OSR_response <0;
    OSR_response = 0;
end

OSL_response = dot(OSL_gabr,stim);
if OSL_response <0;
    OSL_response = 0;
end

figure('name','Example responses to bar stimuli','NumberTitle','off') % Plotting filters, spatial impulse, and response (scalar)
set(gcf, 'Position',[scrsz(3)/6 scrsz(4)/6 scrsz(3)/1.5 scrsz(4)/1.5]);
subplot(2,2,1)
plot(d,EOC_gabr); hold on; plot(d,stim); set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1]);
text(-0.4,1.1,['EOC Response: ' num2str(EOC_response)],'FontSize',12);
subplot(2,2,3)
plot(d,EFC_gabr); hold on; plot(d,stim); set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1]);
xlabel('RF (deg)');
text(-0.4,1.1,['EFC Response: ' num2str(EFC_response)],'FontSize',12);
subplot(2,2,2)
plot(d,OSR_gabr); hold on; plot(d,stim); set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1]);
text(-0.4,1.1,['OSR Response: ' num2str(OSR_response)],'FontSize',12);
subplot(2,2,4)
plot(d,OSL_gabr); hold on; plot(d,stim); set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1]);
xlabel('RF (deg)');
text(-0.4,1.1,['OSL Response: ' num2str(OSL_response)],'FontSize',12);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Responses of filters to sine stimuli to determine optimal stimulus frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear HW_drift_R HW_drift_r Q_drift_R Q_drift_r

stim = zeros(size(d));
f = [0.1 0.3 0.5 0.8 1 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.8 2 2.2 2.4 3 6];
phase = linspace(-pi,pi,1000);
HW_filters = [EOC_gabr; EFC_gabr; OSR_gabr; OSL_gabr];
Q_filters = [EOCs1_gabr; EOCs2_gabr; EFCs1_gabr; EFCs2_gabr];

for i = 1:length(f)
    for j = 1:length(phase)
        stim = sin(x*f(i)+phase(j));
        for k = 1:max(size(HW_filters(:,1)))             % k is iterating filters
            HW_drift_R(i,j,k) = dot(HW_filters(k,:),stim);   % Q_resp(k,i,j) where...
        end
        for l = 1:max(size(Q_filters(:,1)))             % l is iterating filters
            Q_drift_R(i,j,l) = dot(Q_filters(l,:),stim);   % Q_resp(l,i,j) where...
        end
    end
end
% max response of each filter given by max value returned from dot
% product of filter with sine stimuli of optimal frequency with amplitude = 1
EOC_scale = 1/max(max(HW_drift_R(:,:,1),[],2));
EFC_scale = 1/max(max(HW_drift_R(:,:,2),[],2));
OSR_scale = 1/max(max(HW_drift_R(:,:,3),[],2));
OSL_scale = 1/max(max(HW_drift_R(:,:,4),[],2));
HW_scales = [EOC_scale; EFC_scale; OSR_scale; OSL_scale]; % create matrix for loops

EOCs1_scale = 1/max(max(Q_drift_R(:,:,1),[],2));
EOCs2_scale = 1/max(max(Q_drift_R(:,:,2),[],2));
EFCs1_scale = 1/max(max(Q_drift_R(:,:,3),[],2));
EFCs2_scale = 1/max(max(Q_drift_R(:,:,4),[],2));
Q_scales = [EOCs1_scale; EOCs2_scale; EFCs1_scale; EFCs2_scale];

thresh = 0.2;

for i = 1:length(f)
    for j = 1:length(phase)
        stim = sin(x*f(i)+phase(j));
        for k = 1:max(size(HW_filters(:,1)))             % k is iterating filters 
            HW_drift_r(i,j,k) = HW_scales(k)*dot(HW_filters(k,:),stim);   % Q_resp(i,j,k) where...
            if HW_drift_r(i,j,k) < thresh;                                % i frequency, j is phase, k is fitler, 
                HW_drift_r(i,j,k) = 0;                      % Thresholding, or rectification in 'if' statement
            end
        end
        for l = 1:max(size(Q_filters(:,1)))             % l is iterating filters   
            Q_drift_r(i,j,l) = Q_scales(l).*dot(Q_filters(l,:),stim);   % Q_resp(i,j,k) where...
            if Q_drift_r(i,j,l) < thresh;                               % i is bar position, j is bar contrast, l is fitler, 
                Q_drift_r(i,j,l) = 0;
            end
        end
    end
end
for i = 1:length(f)      % Calculating Q_CC and HW_CC responses
    for j = 1:length(phase)
        CC_HW_drift_r(i,j) = sum(HW_drift_r(i,j,:));
        CC_Q_drift_r(i,j) = sqrt(sum((Q_drift_r(i,j,:).^2)));
    end
end

plot(phase,HW_drift_r(10,:,1)); set(gca,'Ylim',[0 1.1],'Xlim',[-pi pi])

plot(phase,CC_HW_drift_r(17,:)); set(gca,'Ylim',[0 2],'Xlim',[-pi pi])

t = (0:points-1)*dd;
NFFT = 2^nextpow2(points);
Y = fft(CC_HW_drift_r(1,:),NFFT)/points;
fr = 2000/2*linspace(0,1,NFFT/2+1);

plot(fr,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
set(gca,'Xlim',[0 10])
for i = 1:length(f)
Y = fft(CC_Q_drift_r(i,:),NFFT)/points;
ratio(i) = abs(Y(2))/abs(Y(1));
end
plot(f,ratio);

figure(1)              % plot peak responses for each frequency
subplot(2,1,1)         % complex cells
plot(f/2,max(CC_HW_drift_r,[],2))
subplot(2,1,2)
plot(f/2,max(CC_Q_drift_r,[],2))

figure(2)              % phase response curves for HW model filters
subplot(4,1,1)
semilogx(f./2,max(HW_drift_r(:,:,1),[],2))
set(gca,'Xlim',[0.1 3])
subplot(4,1,2)
semilogx(f./2,max(HW_drift_r(:,:,2),[],2))
set(gca,'Xlim',[0.1 3])
subplot(4,1,3)
semilogx(f./2,max(HW_drift_r(:,:,3),[],2))
set(gca,'Xlim',[0.1 3])
subplot(4,1,4)
semilogx(f./2,max(HW_drift_r(:,:,4),[],2))
set(gca,'Xlim',[0.1 3])

figure(3)              % phase response curves for Q model filters
subplot(4,1,1)
semilogx(f./2,max(Q_drift_r(:,:,1),[],2))
set(gca,'Xlim',[0.1 3])
subplot(4,1,2)
semilogx(f./2,max(Q_drift_r(:,:,1),[],2))
set(gca,'Xlim',[0.1 3])
subplot(4,1,3)
semilogx(f./2,max(Q_drift_r(:,:,1),[],2))
set(gca,'Xlim',[0.1 3])
subplot(4,1,4)
semilogx(f./2,max(Q_drift_r(:,:,1),[],2))
set(gca,'Xlim',[0.1 3])

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating filter and 'Complex cell' responses to single bars across RF's

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear HW_resp HW_R CC_HW_r CC_HW_R CC_Q_r CC_Q_R Q_R Q_resp
stim = zeros(size(d)); %Creating impulse function

num_pos = 11; % number of positions to test
pw = points/num_pos;
barposn = round(linspace(min(d),max(d),num_pos).*100)/100;
contrast = [-0.5 0.5]; % bar contrasts being tested
thresh = 0;

HW_filters = [EOC_gabr; EFC_gabr; OSR_gabr; OSL_gabr];
HW_scales = [EOC_scale; EFC_scale; OSR_scale; OSL_scale];

Q_filters = [EOCs1_gabr; EOCs2_gabr];%; EFCs1_gabr; EFCs2_gabr]; [EOC_gabr; OSL_gabr];%
Q_scales = [EOCs1_scale; EOCs2_scale; EFCs1_scale; EFCs2_scale]; %[EOC_scale; OSL_scale];%

for j = 1:length(contrast)                             % j is iterating over contrasts
    for i = 1:num_pos                                  % i is iterating bar positions
        if i == 1                                      % to reset bar position at -1.0 on sucessive loops
            pp = 1;
        end
        stim = zeros(size(d));                         % reset stim function
        stim(pp:pp+pw-1) = contrast(j);                % generate current stim
        stim = stim(1:2000);                           % prune stim when it gets extended by previous
        pp = pp + pw;                                  % step pulse position by pulse width for next loop
        
        for k = 1:max(size(HW_filters(:,1)))             % k is iterating filters
            HW_R(k,i,j)= HW_scales(k)*dot(HW_filters(k,:),stim); % holding on to unrectified responses
            HW_resp(k,i,j) = HW_scales(k)*dot(HW_filters(k,:),stim);   % Q_resp(k,i,j) where...
            if HW_resp(k,i,j) < thresh;                       % k is fitler, i is bar position, j is bar contrast
                HW_resp(k,i,j) = 0;                      % Thresholding, or rectification in 'if' statement
            end
        end
                                                        % now Q model filter responses
        for l = 1:max(size(Q_filters(:,1)))             % l is iterating filters
            Q_R(l,i,j) = Q_scales(l).*dot(Q_filters(l,:),stim);
            Q_resp(l,i,j) = Q_scales(l).*dot(Q_filters(l,:),stim);   % Q_resp(l,i,j) where...
            if Q_resp(l,i,j) < thresh;                       % l is fitler, i is bar position, j is bar contrast
                Q_resp(l,i,j) = 0;
            end
        end
    end
end

for i = 1:num_pos      % Calculating Q_CC and HW_CC responses
    for j = 1:length(contrast)
        CC_HW_R(i,j) = sum(HW_R(:,i,j));  % unrectified, for no great reason...
        CC_HW_r(i,j) = sum(HW_resp(:,i,j));
        CC_Q_R(i,j) = sqrt(sum((Q_R(:,i,j).^2)));
        CC_Q_r(i,j) = sqrt(sum((Q_resp(:,i,j).^2)));
    end
end

%FIGURE - plotting Q filter responses to single bars of positive and negative contrast
% TRY DIFFERENT COMBINATIONS OF FITLERS IN Q MODEL (EOC and OSL)
% Rename to energy model
figure('name','Q filter responses to single bars','NumberTitle','off')
set(gcf,'Position',[scrsz(3)/8 scrsz(4)/4 scrsz(3)/2 scrsz(4)/1.8])
for i = 1:1:max(size(Q_filters(:,1)))
    
    subplot(max(size(Q_filters(:,1))), 2, (i+i-1));
    plot(d,Q_filters(i,:)); hold on; axis equal; apos = get(gca,'Position');
    set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on',...
    'Position', [apos(1)/4 apos(2) apos(3) apos(4)],'GridLineStyle','-')
    if i == 1
        title('Filter');
    end
    
    subplot(max(size(Q_filters(:,1))), 2, i+i);
    bar(barposn,-(Q_resp(i,:,1))); hold on; bar(barposn,Q_resp(i,:,2)); apos = get(gca,'Position');
    set(gca,'Xlim',[min(d)+min(d)/10 max(d)+max(d)/10],'Position', [apos(1)/1.5 apos(2) apos(3)*1.7 apos(4)]);
    if i == 1
        title('Response by Bar Position');
    end
    if i == max(size(HW_filters(:,1)))
        xlabel('Bar position (deg)')
    end    
end

str = num2str(barposn);
tmp = regexp(str,'([^ ]*)','tokens'); barpositions = cat(2,tmp{:});

%FIGURE - Plotting Q model complex cell responses to bars positive and negative contrast
figure('name','Q model Complex Cell Response','NumberTitle','off')
set(gcf,'Position',[scrsz(3)/8 scrsz(4)/3 scrsz(3)/2 scrsz(4)/3])
bar(barposn,-CC_Q_r(:,1)); hold on; bar(barposn,CC_Q_r(:,2)); apos = get(gca,'Position');
title('Complex cell response of Quadriture model');xlabel('Bar position (deg)');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10]); ylabel('Response');

figure(4)
bar(barposn,-CC_Q_R(:,1)); hold on; bar(barposn,CC_Q_R(:,2)); apos = get(gca,'Position');
set(gcf,'Position',[scrsz(3)/8 scrsz(4)/3 scrsz(3)/2 scrsz(4)/3])
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10]); ylabel('Response');

% FIGURE -  %Plotting H&W simple filters and their responses to bars
figure('name','H&W filter responses to single bars','NumberTitle','off')
set(gcf,'Position',[scrsz(3)/8 scrsz(4)/20 scrsz(3)/2 scrsz(4)/1.2])
for i = 1:max(size(HW_filters(:,1)))
    
    subplot(max(size(HW_filters(:,1))), 2, (i+i-1));
    plot(d,HW_filters(i,:)); hold on; axis equal; apos = get(gca,'Position');
    set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],...
        'Position', [0 apos(2) apos(3) apos(4)],'YGrid','on','GridLineStyle','-')
    if i == 1
        title('Filter');
    end
    
    
    subplot(max(size(HW_filters(:,1))), 2, i+i);
    bar(barposn,-(HW_resp(i,:,1))); hold on; bar(barposn,HW_resp(i,:,2)); apos = get(gca,'Position');
    set(gca,'Position', [apos(1)/1.5 apos(2) apos(3)*1.5 apos(4)],'Xlim',[-1.1 1.1])
    if i == 1
        title('Response by Bar Position');
    end
    if i == max(size(HW_filters(:,1)))
        xlabel('Bar position (deg)')
    end
end

% FIGURE - Plotting H&W complex cell responses to bars positive and negative contrast
figure('name','H&W Complex Cell Response','NumberTitle','off')
set(gcf,'Position',[scrsz(3)/8 scrsz(4)/3 scrsz(3)/2 scrsz(4)/3])
bar(barposn,-CC_HW_r(:,1)); hold on; bar(barposn,CC_HW_r(:,2)); apos = get(gca,'Position');
title('Complex cell response of H&W model');xlabel('Bar position (deg)');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10]);
% 
% figure(2)    % This is plot of complex cell response to unrectified response, for no good reason
% bar(bars,-CC_HW_R(:,1)); hold on; bar(bars,CC_HW_R(:,2)); apos = get(gca,'Position');
% set(gca,'Xlim', [0 num_pos+1], 'Position', [apos(1)/1.8 apos(2) apos(3) apos(4)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Responses of H&W & Q filters to two-bar stimuli

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear HW_2b_resp HW_2b_R CC_HW_2b_r CC_HW_2b_R HW_CB_R HW_CB_resp Q_CB_R CC_Q_2b_R
clear Q_CB_resp HW_CC_CB_R HW_CC_CB_r Q_CC_CB_R Q_CC_CB_r Q_2b_resp Q_2b_R CC_Q_2b_r

stim = zeros(size(d)); %Creating impulse function

num_pos = 11; % number of positions to test using odd number so center is tested
pw = points/num_pos;
barposn = round(linspace(min(d),max(d),num_pos).*100)/100;
contrast = [-0.5 0.5]; % bar contrasts being tested
thresh = 0;

HW_filters = [EOC_gabr; EFC_gabr; OSR_gabr; OSL_gabr];
HW_scales = [EOC_scale; EFC_scale; OSR_scale; OSL_scale];

Q_filters = [EOCs1_gabr; EOCs2_gabr];%; EFCs1_gabr; EFCs2_gabr]; [EOC_gabr; OSL_gabr];%
Q_scales = [EOCs1_scale; EOCs2_scale]; % EFCs1_scale EFCs2_scale]; [EOC_scale; OSL_scale];

cond_bar = 5; % given as index position of bar
cp = points/2+(barposn(cond_bar)*1000+1);

% calculating response to conditioning bar alone, for baseline in two bar interactions
for j = 1:length(contrast)
    stim = zeros(size(d));
    stim(cp:cp+pw-1) = contrast(j);
    for k = 1:max(size(HW_filters(:,1)))             % k is iterating filters
        HW_CB_R(k,j) = HW_scales(k)*dot(HW_filters(k,:),stim);
        HW_CB_resp(k,j) = HW_scales(k)*dot(HW_filters(k,:),stim);   % CB_resp(k,j) where...
        if HW_CB_resp(k,j) < thresh;                       % k is fitler, j is bar contrast
            HW_CB_resp(k,j) = 0;                      % Thresholding, or rectification
        end
    end
    
    for l = 1:max(size(Q_filters(:,1)))             % k is iterating filters
        Q_CB_R(l,j) = Q_scales(l)*dot(Q_filters(l,:),stim);
        Q_CB_resp(l,j) = Q_scales(l)*dot(Q_filters(l,:),stim);   % CB_resp(k,j) where...
        if Q_CB_resp(l,j) < thresh;                       % k is fitler, j is bar contrast
            Q_CB_resp(l,j) = 0;                      % Thresholding, or rectification
        end
    end
end
for i = 1:length(contrast)                        % Calculating CC response to conditioning bar alone
    HW_CC_CB_R(i) = sum(HW_CB_R(:,i));
    HW_CC_CB_r(i) = sum(HW_CB_resp(:,i));
    Q_CC_CB_R(i) = sum(Q_CB_R(:,i).^2);           % conditioning bar alone Q model
    Q_CC_CB_r(i) = sum(Q_CB_resp(:,i).^2);
end

% now two bar interactions
for p = 1:length(contrast)                                 % p is iterating over contrasts for conditioning bars
    for j = 1:length(contrast)                             % j is iterating over contrasts for test bars
        for i = 1:num_pos                                  % i is iterating bar positions
            if i == 1                                      % to reset bar position at -1.0 on sucessive loops
                pp = 1;
            end
            stim = zeros(size(d));
            stim(cp:cp+pw-1) = contrast(p);
            stim(pp:pp+pw-1) = stim(pp:pp+pw-1) + contrast(j);                % generate current stim
            stim = stim(1:2000);                           % prune stim when it gets extended by previous
            pp = pp + pw;                                  % step pulse position by pulse width for next loop
            
            for k = 1:max(size(HW_filters(:,1)))           % k is iterating filters of HW model
                HW_2b_R(k,i,j,p) = HW_scales(k)*dot(HW_filters(k,:),stim);
                HW_2b_resp(k,i,j,p) = HW_scales(k)*dot(HW_filters(k,:),stim); % Q_resp(k,i,j,p) where...
                if HW_2b_resp(k,i,j,p) < thresh;           % k is fitler, i is bar position, j is test bar contrast
                    HW_2b_resp(k,i,j,p) = 0;               % Thresholding, or rectification statement
                end
            end
            
            for l = 1:max(size(Q_filters(:,1)))             % l is iterating filters of Q model
                Q_2b_R(l,i,j,p) = Q_scales(l)*dot(Q_filters(l,:),stim);
                Q_2b_resp(l,i,j,p) = Q_scales(l)*dot(Q_filters(l,:),stim);   % Q_resp(k,i,j,p) where...
                if Q_2b_resp(l,i,j,p) < thresh;             % k is fitler, i is bar position, j is bar contrast
                    Q_2b_resp(l,i,j,p) = 0;                 % Thresholding, or rectification statement
                end
            end
        end
    end
end

for p = 1:length(contrast)
    for i = 1:num_pos      % Calculating HW_CC and Q_CC responses for 2 bar interactions
        for j = 1:length(contrast)
            CC_HW_2b_R(i,j,p) = sum(HW_2b_R(:,i,j,p));
            CC_HW_2b_r(i,j,p) = sum(HW_2b_resp(:,i,j,p));
            CC_Q_2b_R(i,j,p) = sqrt(sum((Q_2b_R(:,i,j,p).^2)));
            CC_Q_2b_r(i,j,p) = sqrt(sum((Q_2b_resp(:,i,j,p).^2)));
        end
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting Q filters and their responses to two bars
figure('name','H&W filter responses to two bars','NumberTitle','off')
set(gcf,'Position',[scrsz(3)/8 scrsz(4)/20 scrsz(3)/1.5 scrsz(4)/1.2]);
for i = 1:1:max(size(Q_filters(:,1)))
    
    subplot(max(size(Q_filters(:,1))), 3, 3*(i-1)+1);
    plot(d,Q_filters(i,:)); hold on; axis equal; apos = get(gca,'Position');
    set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1], 'Position', [apos(1)/4 apos(2) apos(3) apos(4)])
    if i == 1
        title('Filter');
    end
    
    subplot(max(size(Q_filters(:,1))), 3, 2*i+i-1);
    bar(barposn,-(Q_2b_resp(i,:,1,2))); hold on; bar(barposn,Q_2b_resp(i,:,2,2));
    text(barposn(cond_bar),0.25,'*','FontSize',18); apos = get(gca,'Position');
    set(gca,'Xlim',[min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim', [-0.4 0.4],'Position', [apos(1)/1.35 apos(2) apos(3)*1.4 apos(4)]);
        if i == 1
        title('Response with bright conditioning bar');
        end
        if i == max(size(Q_filters(:,1)))
            xlabel('Bar position (deg)')
        end
        
    subplot(max(size(Q_filters(:,1))), 3, i+2*i);
    bar(barposn,-(Q_2b_resp(i,:,1,1))); hold on; bar(barposn,Q_2b_resp(i,:,2,1)); 
    text(barposn(cond_bar),-0.25,'*','FontSize',18); apos = get(gca,'Position');
    set(gca,'Xlim',[min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim', [-0.4 0.4],'Position', [apos(1)/1.05 apos(2) apos(3)*1.4 apos(4)]);
    if i == 1
        title('Response with dark conditioning bar');
    end
    if i == max(size(Q_filters(:,1)))
        xlabel('Bar position (deg)')
    end
end

%%%Quadriture model 'Complex Cell' response to two bars
figure('name','Q filter responses to two bars','NumberTitle','off') %Q CC response to two bars
set(gcf,'Position',[scrsz(3)/8 scrsz(4)/8 scrsz(3)/2 scrsz(4)/1.5])
subplot(2,1,1)
bar(barposn,-CC_Q_2b_r(:,1,2)); hold on; bar(barposn,CC_Q_2b_r(:,2,2)); apos = get(gca,'Position');
text(barposn(cond_bar),0.45,'*','FontSize',18); title('Bright Conditioning Bar');
set(gca,'Xlim',[min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim', [-0.6 0.6], 'Position', [apos(1) apos(2) apos(3) apos(4)]);
ylabel('Response')
subplot(2,1,2)
bar(barposn,-CC_Q_2b_r(:,1,1)); hold on; bar(barposn,CC_Q_2b_r(:,2,1)); apos = get(gca,'Position');
text(barposn(cond_bar),-0.45,'*','FontSize',18); title('Dark Conditioning Bar'); xlabel('Bar position (deg)');
ylabel('Response')
set(gca,'Xlim',[min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim', [-0.6 0.6], 'Position', [apos(1) apos(2) apos(3) apos(4)]);
% 
% subplot(3,1,3); bar(barposn,-CC_Q_2b_R(:,1,2)); hold on; bar(barposn,CC_Q_2b_R(:,2,2)); apos = get(gca,'Position');
% text(barposn(cond_bar),0.45,'*','FontSize',18);
% set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim', [-.6 .6], 'Position', [apos(1) apos(2) apos(3) apos(4)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting H&W filters and their responses to two bars
figure('name','H&W filter responses to two bars','NumberTitle','off')
set(gcf,'Position',[scrsz(3)/8 scrsz(4)/20 scrsz(3)/1.5 scrsz(4)/1.2]);
for i = 1:1:max(size(HW_filters(:,1)))
    
    subplot(max(size(HW_filters(:,1))), 3, 3*(i-1)+1);
    plot(d,HW_filters(i,:)); hold on; axis equal; apos = get(gca,'Position');
    set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1], 'Position', [apos(1)/4 apos(2) apos(3) apos(4)])
    if i == 1
        title('Filter');
    end
 
    subplot(max(size(HW_filters(:,1))), 3, 2*i+i-1);
    bar(barposn,-(HW_2b_resp(i,:,1,2))); hold on; bar(barposn,HW_2b_resp(i,:,2,2)); 
    text(barposn(cond_bar),0.25,'*','FontSize',18); apos = get(gca,'Position');
    set(gca,'Xlim',[min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim', [-0.4 0.4],'Position', [apos(1)/1.5 apos(2) apos(3)*1.4 apos(4)]);
    if i == 1
        title('Response with bright conditioning bar');
    end
    if i == max(size(HW_filters(:,1)))
        xlabel('Bar position (deg)')
    end
    
    subplot(max(size(HW_filters(:,1))), 3, i+2*i);
    bar(barposn,-(HW_2b_resp(i,:,1,1))); hold on; bar(barposn,HW_2b_resp(i,:,2,1)); 
    text(barposn(cond_bar),-0.25,'*','FontSize',18); apos = get(gca,'Position');
    set(gca,'Xlim',[min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim', [-0.4 0.4],'Position', [apos(1)/1.1 apos(2) apos(3)*1.4 apos(4)]);
    if i == 1
        title('Response with dark conditioning bar');
    end
    if i == max(size(HW_filters(:,1)))
        xlabel('Bar position (deg)')
    end
end

figure('name','H&W filter responses to two bars','NumberTitle','off') %H&W CC response to two bars
set(gcf,'Position',[scrsz(3)/8 scrsz(4)/8 scrsz(3)/2 scrsz(4)/1.5])
subplot(2,1,1)
bar(barposn,-CC_HW_2b_r(:,1,2)); hold on; bar(barposn,CC_HW_2b_r(:,2,2)); apos = get(gca,'Position');
text(barposn(cond_bar),0.45,'*','FontSize',18); title('Bright Conditioning Bar');
ylabel('Response')
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim', [-0.6 0.6], 'Position', [apos(1) apos(2) apos(3) apos(4)]);
subplot(2,1,2)
bar(barposn,-CC_HW_2b_r(:,1,1)); hold on; bar(barposn,CC_HW_2b_r(:,2,1)); apos = get(gca,'Position');
text(barposn(cond_bar),-0.45,'*','FontSize',18); title('Dark Conditioning Bar'); xlabel('Bar position (deg)');
ylabel('Response')
set(gca,'Xlim',[min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim', [-0.6 0.6], 'Position', [apos(1) apos(2) apos(3) apos(4)]);
% 
% subplot(3,1,3); bar(barposn,-CC_HW_2b_R(:,1,2)); hold on; bar(barposn,CC_HW_2b_R(:,2,2)); apos = get(gca,'Position');
% text(barposn(cond_bar),0.45,'*','FontSize',18);
% set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim', [-.05 .05], 'Position', [apos(1) apos(2) apos(3) apos(4)]);

% Highlighting underlying filters with relative responses
%Evaluate response code for 2 bars with conditioning bar position 5
figure('name','Highlighting the EOCs1 filter','NumberTitle','off')
set(gcf,'Position',[scrsz(3)/16 scrsz(4)/8 scrsz(3)/1.2 scrsz(4)/1.5])
subplot(2,4,1)
bar(barposn,CC_Q_r(:,2)); hold on; lim = get(gca,'Ylim');
title('Response to single bars');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10]); ylabel('Response');
subplot(2,4,5)
bar(barposn,CC_Q_r(:,1)); xlabel('Bar position (deg)');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',lim);
subplot(2,4,2)
bar(barposn,((CC_Q_2b_r(:,2,2)-Q_CC_CB_r(2))/Q_CC_CB_r(2))); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),max(abs(lim))*.8,'*','FontSize',18); 
subplot(2,4,6)
bar(barposn,((CC_Q_2b_r(:,1,1)-Q_CC_CB_r(1))/Q_CC_CB_r(1))); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),-max(abs(lim))*.8,'*','FontSize',18); 
subplot(2,4,3)
bar(barposn,((CC_Q_2b_r(:,2,1)-Q_CC_CB_r(1))/Q_CC_CB_r(1))); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),-max(abs(lim))*.8,'*','FontSize',18); 
subplot(2,4,7)
bar(barposn,((CC_Q_2b_r(:,1,2)-Q_CC_CB_r(2))/Q_CC_CB_r(2))); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),max(abs(lim))*.8,'*','FontSize',18); 
subplot(2,4,4)
plot(d,EOCs1_gabr); hold on; plot(d,EOCs2_gabr); axis equal; apos = get(gca,'Position');
text(barposn(cond_bar),0.45,'*','FontSize',18);
set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on','GridLineStyle','-')
subplot(2,4,8)
plot(d,EOCs2_gabr); hold on; axis equal; apos = get(gca,'Position');
text(barposn(cond_bar),-0.45,'*','FontSize',18);
set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on','GridLineStyle','-')

% Evaluate with bar position 5
figure('name','Highlighting the EOC and EFC filters','NumberTitle','off')
set(gcf,'Position',[scrsz(3)/16 scrsz(4)/8 scrsz(3)/1.2 scrsz(4)/1.5])
subplot(2,4,1)
bar(barposn,CC_HW_r(:,2)); hold on; apos = get(gca,'Position'); lim = get(gca,'Ylim');
title('Response to single bars');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',lim); ylabel('Response');
subplot(2,4,5)
bar(barposn,CC_HW_r(:,1)); xlabel('Bar position (deg)');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10]);
subplot(2,4,2)
bar(barposn,((CC_HW_2b_r(:,2,2)-HW_CC_CB_r(2))/HW_CC_CB_r(2))); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),max(abs(lim))*.8,'*','FontSize',18); 
subplot(2,4,6)
bar(barposn,((CC_HW_2b_r(:,1,1)-HW_CC_CB_r(1))/HW_CC_CB_r(1))); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),-max(abs(lim))*.8,'*','FontSize',18); 
subplot(2,4,3)
bar(barposn,((CC_HW_2b_r(:,2,1)-HW_CC_CB_r(1))/HW_CC_CB_r(1))); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),-max(abs(lim))*.8,'*','FontSize',18); 
subplot(2,4,7)
bar(barposn,((CC_HW_2b_r(:,1,2)-HW_CC_CB_r(2))/HW_CC_CB_r(2))); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),max(abs(lim))*.8,'*','FontSize',18); 
subplot(2,4,4)
plot(d,EOC_gabr); hold on; axis equal; apos = get(gca,'Position');
text(barposn(cond_bar),0.45,'*','FontSize',18);
set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on','GridLineStyle','-')
subplot(2,4,8)
plot(d,EFC_gabr); hold on; axis equal; apos = get(gca,'Position');
text(barposn(cond_bar),-0.45,'*','FontSize',18);
set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on','GridLineStyle','-')


% Evaluate with bar position 4
figure('name','Highlighting the OSL filter','NumberTitle','off')
set(gcf,'Position',[scrsz(3)/16 scrsz(4)/8 scrsz(3)/1.2 scrsz(4)/1.5])
subplot(2,4,1)
bar(barposn,CC_HW_r(:,2)); hold on; apos = get(gca,'Position'); lim = get(gca,'Ylim');
title('Response to single bars');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',lim); ylabel('Response');
subplot(2,4,5)
bar(barposn,CC_HW_r(:,1)); xlabel('Bar position (deg)');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10]);
subplot(2,4,2)
bar(barposn,((CC_HW_2b_r(:,2,2)-HW_CC_CB_r(2))/HW_CC_CB_r(2))); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),max(abs(lim))*.8,'*','FontSize',18); 
subplot(2,4,6)
bar(barposn,((CC_HW_2b_r(:,1,1)-HW_CC_CB_r(1))/HW_CC_CB_r(1))); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),-max(abs(lim))*.8,'*','FontSize',18); 
subplot(2,4,3)
bar(barposn,((CC_HW_2b_r(:,2,1)-HW_CC_CB_r(1))/HW_CC_CB_r(1))); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),-max(abs(lim))*.8,'*','FontSize',18); 
subplot(2,4,7)
bar(barposn,((CC_HW_2b_r(:,1,2)-HW_CC_CB_r(2))/HW_CC_CB_r(2))); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),max(abs(lim))*.8,'*','FontSize',18); 
subplot(2,4,4)
plot(d,OSL_gabr); hold on; plot(d,EOC_gabr); plot(d,OSR_gabr),axis equal; apos = get(gca,'Position');
text(barposn(cond_bar),0.45,'*','FontSize',18);
set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on','GridLineStyle','-')
subplot(2,4,8)
plot(d,OSL_gabr); hold on; axis equal; apos = get(gca,'Position');
text(barposn(cond_bar),-0.45,'*','FontSize',18);
set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on','GridLineStyle','-')


figure(3)
set(gcf,'Position', [101    62   560   609])
for i = 1:1:max(size(HW_filters(:,1)))
    
    subplot(max(size(HW_filters(:,1))), 2, i+i-1);
    plot(d,HW_filters(i,:)); hold on; axis equal; apos = get(gca,'Position');
    set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1])
    
    subplot(max(size(HW_filters(:,1))), 2, 2);
    bar(barposn,((CC_HW_2b_r(:,2,2)-HW_CC_CB_r(2))/HW_CC_CB_r(2))); hold on; apos = get(gca,'Position');
    text(barposn(cond_bar),0.45,'*','FontSize',18);
    set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10]);
end

figure(3)
set(gcf,'Position', [101    62   560   609])
for i = 1:1:max(size(Q_filters(:,1)))
    
    subplot(2, max(size(Q_filters(:,1))), i);
    plot(d,Q_filters(i,:)); hold on; axis equal; apos = get(gca,'Position');
    set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1])
    
    subplot(2, max(size(Q_filters(:,1))), 3);
    bar(bars,((CC_Q_2b_r(:,2)-Q_CC_CB_r(2))/Q_CC_CB_r(2))); hold on; apos = get(gca,'Position');
    set(gca,'Xlim', [0 num_pos+1]);
end


figure(3)
bar(barposn,(CC_Q_2b_r(:,2)-Q_CC_CB_r(2))); hold on; apos = get(gca,'Position');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10], 'Position', [apos(1)/1.8 apos(2) apos(3) apos(4)]);


%The ringing odd symmetric right (ORR) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = gaussmf(x,[pi/4 pi]);
f = 2; 
wave = sin(x*f+pi);
OSR_gabr = b.*wave;
c = 1/max(OSR_gabr);  %These filters get scaled up to so thier max value = 1
OSR_gabr = OSR_gabr*c; %This wasn't neccessary for the even functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','Odd symmetric right-sided','NumberTitle','off')
set(gcf, 'Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3)/1.5 scrsz(4)/2]);
subplot(1,2,1)
plot(x,wave,'r'); hold on; plot(x,b); set(gca,'Ytick',[0],'Xlim',[-0.05 2*pi+0.05]);
title('Sine wave and gaussian')
subplot(1,2,2)
plot(d,OSR_gabr); set(gca,'Ytick',[0],'Xlim', [-1.01 1.01],'Ylim',[-1 1]); xlabel('RF (degrees)');
title('Resulting gabor')

OSR_on = trapz(abs(OSR_gabr(OSR_gabr > 0)))   %positive lobes
OSR_off = trapz(abs(OSR_gabr(OSR_gabr < 0)))  %negative lobe


% plotting overlyaing filters for pres fig


figure('Name','Quadrature Filters','NumberTitle','off')
plot(d,EOCs1_gabr,'r','linewidth',2); hold on; plot(d,EOCs2_gabr,'b','linewidth',2);
set(gca,'Ytick',[0 1],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on'); xlabel('RF space (degrees');
title('Qaudrature filters');

figure('Name','H&W Filters','NumberTitle','off')
plot(d,EOC_gabr,'r','linewidth',2); hold on; plot(d,OSLn_gabr,'b','linewidth',2); plot(d,OSRn_gabr,'g','linewidth',2); plot(d,EFC_gabr,'c','linewidth',2)
set(gca,'Ytick',[0 1],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on'); xlabel('RF space (degrees');
title('H&W filters');

%plotting random stims for STC

figure('Name','Gaussian bar stims','NumberTitle','off')
subplot(4,1,1)
area(stim); box off; set(gca,'Xtick',[],'Ytick',[],'Ylim', [-1 1])
title('Stims');

subplot(4,1,2)
area(stim); box off; set(gca,'Xtick',[],'Ytick',[],'Ylim', [-1 1])

subplot(4,1,3)
area(stim); box off; set(gca,'Xtick',[],'Ytick',[],'Ylim', [-1 1])

subplot(4,1,4)
area(stim); box off; set(gca,'Xtick',[],'Ytick',[],'Ylim', [-1 1])