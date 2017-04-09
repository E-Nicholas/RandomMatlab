%% Complex Cells in 1 Dimension - Eric Nicholas - University of Washington

%  This is a 1 dimensional model of complex cell responses in primary visual cortex.
%  This model produces a scalar response produced by a dot product, for each underlying simple
%  filter that contributes to the complex cell response (there is no
%  temporal component). The idea of creating these models was to develop
%  intuition as to why two-bar interactions produce the reponses they do in complex
%  cells. To this end I tested two methods for generating 'Complex Cell'
%  responses, a quadrature model with two filters spatially offset by 90 degrees, and a 
%  more classic, simple, Hubel and Wiesel-type model that sums the responses of
%  multiple spatially offset half wave rectified filters.
% 
%  This code was developed at the University of Washington as a project in
%  a computational neuroscience course run by Fred Reike, Wyeth Bair, and
%  Adrienne Fairhall.
%
%  The code here is broken up into modules. The first
%  being the creation of the simple filters, the second is finding the max
%  responses of the simple fitlers and testing their spatial frequency
%  tuning. The third cell produces the responses of the model to single and
%  two bar stimuli. The fourth cell is PCA analysis of the complex cell
%  responses. The remaining code from there is plotting functions.
%  Importantly, there is a block at the end of the first cell where one can
%  set a few free model parameters. Namely, bar contrasts to test, 
%  number of bar positions across the RF to test (for line weighting stims),
%  the postion of the conditioning bar in two line interactions,
%  and the rectification threshold for H&W filters. I apologize for not
%  creating loops to handle and show responses for all possible conditioning 
%  bar locations...
%  
%  I'm posting this code years after writing and developing it. I remember
%  it being a fun project and I hope to come back to this in the future. If
%  you're interested in learning more or using this code feel free to
%  contact me.

%%
%Some setup
%%%%
scrsz = get(0,'ScreenSize');
RF_size = 2;             % 2 degree SimpleCell RF's in both Q and H&W model
%%%%
points = RF_size * 1000; % number of points (this is really too many, but I wanted things to look smooth)
x = linspace(0, 2*pi, points); % Filter space, from 0 to 2pi with n points for generating waves
%We're supposed to be thinking in the dimension of space, so d is RF 'space'
%our units are degrees...
d = linspace(0-RF_size/2,0+RF_size/2,points);   % d is RF space, 0 is center of RF
dd = 1/length(d);            % our resolution
RF = length(d);              % size of RF in dd



% Simple Cell filters for H&W model (HW model sums filters)

%The classic even on center (EOC) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = gaussmf(x,[pi/1.8 pi]); % Generating the gaussian
f = 1.5;                    % setting frequency
wave = sin(x*f+pi);         % Generating the wave
EOC_gabr = b.*wave;         % multiplying the two for gabor
                                              % calculating area under lobes...
                                              % of gabor to ensure roughly equal 
EOC_on = trapz(abs(EOC_gabr(EOC_gabr > 0)));   % positive lobe
EOC_off = trapz(abs(EOC_gabr(EOC_gabr < 0)));  % negative lobes. Checked every filter but only showing this here.
%The classic even off center (EFC) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = gaussmf(x,[pi/1.8 pi]);
f = 1.5; 
wave = sin(x*f);
EFC_gabr = b.*wave;
%The odd symmetric right (OSR) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = gaussmf(x,[pi/4 pi]);
f = 1; 
wave = sin(x*f+pi);
OSR_gabr = b.*wave;
c = 1/max(OSR_gabr);  %These odd filters get scaled up to so their max value = 1
OSR_gabr = OSR_gabr*c; %This wasn't neccessary for the even functions
%Odd symmetric left (OSL) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = gaussmf(x,[pi/4 pi]);
f = 1; 
wave = sin(x*f);
OSL_gabr = b.*wave;
c = 1/max(OSL_gabr);
OSL_gabr = OSL_gabr*c;
%The new OSL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = gaussmf(x,[(pi/1.2)/2 pi]);
f = 1.5; 
wave = sin(x*f+3*pi/2);
OSLn_gabr = b.*wave;
c = 1/max(OSLn_gabr);  
OSLn_gabr = OSLn_gabr*c; 
%The new OSR  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = gaussmf(x,[(pi/1.2)/2 pi]);
f = 1.5; 
wave = sin(x*f+pi/2);
OSRn_gabr = b.*wave;
c = 1/max(OSRn_gabr);  
OSRn_gabr = OSRn_gabr*c; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creating simple cell filters for Energy model (Q-model sums squared responses)

% Even on center shift_1 (-45 deg from center, left) (EOCs1) %%%%%%%%%%%%%
b = gaussmf(x,[pi/1.8 pi-pi/6]);
f = 1.5; 
wave = sin(x*f+pi/4+pi);
EOCs1_gabr = b.*wave;

% Even on center shift_2 (+45 deg from center, right) (EOCs2) %%%%%%%%%%%%
b = gaussmf(x,[pi/1.8 pi+pi/6]);
f = 1.5; 
wave = sin(x*f-pi/4+pi);
EOCs2_gabr = b.*wave;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Free Model parameters

HW_filters = [EOC_gabr; EFC_gabr; OSRn_gabr; OSLn_gabr];
Q_filters = [EOCs1_gabr; EOCs2_gabr];  
                                                              
num_pos = 16; % number of bar positions to test
thresh = 0.0; % Response threshold for H&W filters (0:1) (response must be above to be recorded as response)
contrast = [-0.5 0.5]; % bar contrasts being tested

cond_bar = 7; % given as index position of conditioning bar
%Initially I had this set to take degrees, but this was got confusing when
%I started to play wit the numbers of bars to test..

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pw = points/num_pos;
phase = linspace(-pi,pi,1000);
barposn = round(linspace(min(d),max(d),num_pos).*100)/100;
cp = points/2+(barposn(cond_bar)*1000+1);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Responses of filters to sine stimuli to determine optimal stimulus frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stim = zeros(size(d));
f = [0.1 0.3 0.5 0.8 1 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.8 2 2.2 2.4 3 6]; %set these by hand and never went back... :/

for i = 1:length(f)  % Getting max response for scales
    for j = 1:length(phase)
        stim = sin(x*f(i)+phase(j));
        for k = 1:max(size(HW_filters(:,1)))             % k is iterating filters
            HW_drift_R(i,j,k) = dot(HW_filters(k,:),stim);   % Q_r(k,i,j) where...
        end
        for l = 1:max(size(Q_filters(:,1)))             % l is iterating filters
            Q_drift_R(i,j,l) = dot(Q_filters(l,:),stim);   % Q_r(l,i,j) where...
        end
    end
end
% max response of each filter given by max value returned from dot
% product of filter with sine stimuli of optimal frequency with amplitude = 1
EOC_scale = 1/max(max(HW_drift_R(:,:,1),[],2)); %Setting scale factors for each filter to normalize responses
EFC_scale = 1/max(max(HW_drift_R(:,:,2),[],2));
OSR_scale = 1/max(max(HW_drift_R(:,:,3),[],2));
OSL_scale = 1/max(max(HW_drift_R(:,:,4),[],2));

HW_scales = [EOC_scale; EFC_scale; OSR_scale; OSL_scale]; % create matrix for loops

EOCs1_scale = 1/max(max(Q_drift_R(:,:,1),[],2));
EOCs2_scale = 1/max(max(Q_drift_R(:,:,2),[],2));

Q_scales = [EOCs1_scale; EOCs2_scale]; 

for i = 1:length(f) % Testing frequency tuning
    for j = 1:length(phase)
        stim = sin(x*f(i)+phase(j));
        for k = 1:max(size(HW_filters(:,1)))             % k is iterating filters
            HW_drift_r(i,j,k) = HW_scales(k)*dot(HW_filters(k,:),stim);   % Q_r(i,j,k) where...
            if HW_drift_r(i,j,k) < thresh;                                % i frequency, j is phase, k is fitler,
                HW_drift_r(i,j,k) = 0;                      % Thresholding, or rectification in 'if' statement
            end
            for l = 1:max(size(Q_filters(:,1)))             % l is iterating filters
                Q_drift_norm_R(i,j,l) = Q_scales(l)*dot(Q_filters(l,:),stim);   % Q_r(l,i,j) where...
            end
        end
    end
end
for i = 1:length(f)      % Calculating Q_CC and HW_CC responses
    for j = 1:length(phase)
        CC_HW_drift_r(i,j) = sum(HW_drift_r(i,j,:)); %summing rectified responses
        CC_Q_drift_R(i,j) = sqrt(sum((Q_drift_norm_R(i,j,:).^2))); %summing squared unrectified responses
    end
end
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating filter and 'Complex cell' responses to single bars across RF's

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%single bars
clear HW_r HW_R CC_HW_r CC_HW_R CC_Q_r CC_Q_R Q_R Q_r

stim = zeros(size(d)); %Reseting

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
            HW_r(k,i,j) = HW_scales(k)*dot(HW_filters(k,:),stim);   % Q_r(k,i,j) where...
            if HW_r(k,i,j) < thresh;                       % k is fitler, i is bar position, j is bar contrast
                HW_r(k,i,j) = 0;                      % Thresholding, or rectification in 'if' statement
            end
        end
                                                        % now Q model filter responses
        for l = 1:max(size(Q_filters(:,1)))             % l is iterating filters
            Q_R(l,i,j) = Q_scales(l).*dot(Q_filters(l,:),stim);  % Q_r(l,i,j) where...
        end                                                      % l is fitler, i is bar position, j is bar contrast
    end
end
for i = 1:num_pos      % Calculating Q_CC and HW_CC responses
    for j = 1:length(contrast)
        CC_HW_R(i,j) = sum(HW_R(:,i,j));  % unrectified, for no great reason...
        CC_HW_r(i,j) = sum(HW_r(:,i,j));
        CC_Q_R(i,j) = sqrt(sum((Q_R(:,i,j).^2)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Responses of H&W & Q filters to two-bar stimuli

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear HW_2b_r HW_2b_R CC_HW_2b_r CC_HW_2b_R HW_CB_R HW_CB_resp Q_CB_R CC_Q_2b_R
clear Q_CB_resp HW_CC_CB_R HW_CC_CB_r Q_CC_CB_R Q_CC_CB_r Q_2b_r Q_2b_R CC_Q_2b_r

stim = zeros(size(d)); %Reseting impulse function

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
        Q_CB_R(l,j) = Q_scales(l)*dot(Q_filters(l,:),stim);  % CB_resp(k,j) where k is fitler, j is bar contrast
    end
end
for i = 1:length(contrast)                        % Calculating CC response to conditioning bar alone
    HW_CC_CB_R(i) = sum(HW_CB_R(:,i));
    HW_CC_CB_r(i) = sum(HW_CB_resp(:,i));
    Q_CC_CB_R(i) = sqrt(sum(Q_CB_R(:,i).^2));           % conditioning bar alone Q model
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
                HW_2b_r(k,i,j,p) = HW_scales(k)*dot(HW_filters(k,:),stim); % Q_r(k,i,j,p) where...
                if HW_2b_r(k,i,j,p) < thresh;           % k is fitler, i is bar position, j is test bar contrast
                    HW_2b_r(k,i,j,p) = 0;               % Thresholding, or rectification statement
                end
            end
            
            for l = 1:max(size(Q_filters(:,1)))             % l is iterating filters of Q model
                Q_2b_R(l,i,j,p) = Q_scales(l)*dot(Q_filters(l,:),stim);  % Q_r(k,i,j,p) where k is fitler, 
            end                                                          % i is bar position, j is bar contrast
        end 
    end 
end

for p = 1:length(contrast)
    for i = 1:num_pos      % Calculating HW_CC and Q_CC responses for 2 bar interactions
        for j = 1:length(contrast)
            CC_HW_2b_R(i,j,p) = sum(HW_2b_R(:,i,j,p));
            CC_HW_2b_r(i,j,p) = sum(HW_2b_r(:,i,j,p));
            CC_Q_2b_R(i,j,p) = sqrt(sum((Q_2b_R(:,i,j,p).^2)));
        end
    end
end

%% STC code

num_bar = 16;
trials = 2000;
for i = 1:trials % Generating random stims
    for j = 1:num_bar
        bars(i,j) = normrnd(0, 0.341);
    end
end
stim = zeros(size(d));
pw = points/num_bar;

for i = 1:trials % Simple filter responses to random stims
    stim = zeros(size(d));
    pp = 1;
    
    for j = 1:num_bar %Generating stim to tests
        stim(pp:pp+pw-1) = bars(i,j);
        stim = stim(1:2000);                           % prune stim when it gets extended by previous
        pp = pp + pw;
    end
    
    for k = 1:max(size(HW_filters(:,1)))             % Getting H&W filter responses to stim
        HW_g_R(k,i)= HW_scales(k)*dot(HW_filters(k,:),stim); % holding on to unrectified responses
        HW_g_r(k,i) = HW_scales(k)*dot(HW_filters(k,:),stim);   
        if HW_g_r(k,i) < thresh;                       % k is fitler, i 
            HW_g_r(k,i) = 0;                      % Thresholding
        end
    end
    for l = 1:max(size(Q_filters(:,1)))            % Getting H&W filter responses to stim 
        Q_g_R(l,i) = Q_scales(l)*dot(Q_filters(l,:),stim);
    end
end

for i = 1:max(size(HW_filters(:,1))) %Getting distributions
    HW_g_mean(i) = mean(HW_g_r(i,:));
    HW_g_sd(i) = std(HW_g_r(i,:));
end
for i = 1:max(size(Q_filters(:,1))) %Getting distributions
    Q_g_mean(i) = mean(Q_g_R(i,:));
    Q_g_sd(i) = std(Q_g_R(i,:));
end

for i = 1:trials      % Calculating Q_CC and HW_CC responses
        CC_HW_g_R(i) = sum(HW_g_R(:,i));  % unrectified, for no great reason...
        CC_HW_g_r(i) = sum(HW_g_r(:,i));
        CC_Q_g_R(i) = sqrt(sum(Q_g_R(:,i).^2));
end
cc_HW_g_mean = mean(CC_HW_g_r);
cc_HW_g_sd = std(CC_HW_g_r);
cc_HW_spikes = find(CC_HW_g_r > 3*cc_HW_g_sd); %Getting spike triggering stims
cc_HW_g_STA = mean(bars(cc_HW_spikes,:));

cc_Q_g_sd = std(CC_Q_g_R);
cc_Q_spikes = find(CC_Q_g_R > 3*cc_Q_g_sd); %Getting spike triggering stims
cc_Q_g_STA = mean(bars(cc_Q_spikes,:));

for i = 1:size(bars(cc_HW_spikes,:),1)
    cc_HW_STA_matrix(i,:) = mean(bars(cc_HW_spikes,:));
end
for i = 1:size(bars(cc_Q_spikes,:),1)
    cc_Q_STA_matrix(i,:) = mean(bars(cc_Q_spikes,:));
end

% Now subtract mean stimulus from spike triggering stims to get "zero mean"
% ensemble. Using Rust, et al. zero mean equation.

% HW_STC_ensemble = minus(bar(cc_HW_spikes,:),cc_HW_STA_matrix);
% Q_STC_ensemble = minus(bar(cc_Q_spikes,:),cc_Q_STA_matrix);
HW_STC_ensemble = minus(bars(cc_HW_spikes,:),(cc_HW_STA_matrix.*bars(cc_HW_spikes,:)).*cc_HW_STA_matrix);
Q_STC_ensemble = minus(bars(cc_Q_spikes,:),(cc_Q_STA_matrix.*bars(cc_Q_spikes,:)).*cc_Q_STA_matrix);

%Get covariance matrix
CC_HW_g_STC = cov(HW_STC_ensemble);
CC_Q_g_STC = cov(Q_STC_ensemble);

%Eigen decompostion
[HW_EigVec, HW_EigVal] = eig(CC_HW_g_STC);
[Q_EigVec, Q_EigVal] = eig(CC_Q_g_STC);

%Show H&W components
for j = 1:num_bar
    pp = 1;
    EigV = zeros(size(d));
    for i = 1:num_bar
        EigV(pp:pp+pw-1) = HW_EigVec(i,j);
        EigV = EigV(1:2000);
        pp = pp + pw;
    end
    figure('name',strcat('HW-Model component:', num2str(j)),'NumberTitle','off')
    area(d,EigV);
    title(strcat('HW-Model component:', num2str(j)))
    set(gca,'Xtick',[],'Ytick',[])
end
%Show Q components
for j = 1:num_bar
    pp = 1;
    EigV = zeros(size(d));
    for i = 1:num_bar
        EigV(pp:pp+pw-1) = Q_EigVec(i,j);
        EigV = EigV(1:2000);
        pp = pp + pw;
    end
    figure('name',strcat('Q-Model component:', num2str(j)),'NumberTitle','off')
    area(d,EigV);
    title(strcat('Q-Model component:', num2str(j)))
    set(gca,'Xtick',[],'Ytick',[])
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting functions
scrsz = get(0,'ScreenSize');

% Energy Model Filters
figure('name','Energy model filters','NumberTitle','off')
subplot(2,1,1)
plot(d,EOCs1_gabr); hold on; set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on');
title('Even on center -45 (EOCs1)');
subplot(2,1,2)
plot(d,EOCs2_gabr); hold on; set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on');
title('Even on center +45 (EOCs2)');

% H&W Model Filters
figure('name','HW model filters','NumberTitle','off')
subplot(2,2,1)
plot(d,EOC_gabr); hold on; set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on');
title('Even on center (EOC)');
subplot(2,2,2)
plot(d,EFC_gabr); hold on; set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on');
title('Even off center (EFC)'); 
subplot(2,2,3)
plot(d,OSRn_gabr); hold on; set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on');
title('Odd symmetric right (OSR)'); xlabel('RF (degrees');
subplot(2,2,4)
plot(d,OSLn_gabr); hold on; set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on');
title('Odd symmetric left (OSL)'); xlabel('RF (degrees');

% Spatial frequency tuning of filters
% frequency tuning for Energy model filters
figure('name','Spatial Frequency Tuning, Q filters','NumberTitle','off')
subplot(2,1,1)
semilogx(f./2,max(Q_drift_norm_R(:,:,1),[],2))
set(gca,'Xlim',[0.1 3],'Ylim',[0 1]); title('Spatial Freq tuning: Q Filters');
subplot(2,1,2)
semilogx(f./2,max(Q_drift_norm_R(:,:,2),[],2))
set(gca,'Xlim',[0.1 3],'Ylim',[0 1]); xlabel('Frequency (cpd)');ylabel('Response');

% frequency tuning for HW model filters
figure('name','Spatial Frequency Tuning, HW filters','NumberTitle','off')
subplot(4,1,1)
semilogx(f./2,max(HW_drift_r(:,:,1),[],2))
set(gca,'Xlim',[0.1 3]); title('Spatial Freq tuning: H&W Filters');
subplot(4,1,2)
semilogx(f./2,max(HW_drift_r(:,:,2),[],2))
set(gca,'Xlim',[0.1 3]); 
subplot(4,1,3)
semilogx(f./2,max(HW_drift_r(:,:,3),[],2))
set(gca,'Xlim',[0.1 3])
subplot(4,1,4)
semilogx(f./2,max(HW_drift_r(:,:,4),[],2))
set(gca,'Xlim',[0.1 3]); xlabel('Frequency (cpd)');ylabel('Response');

%% Responses to single bar stimuli

%FIGURE 1 - plotting Q filter responses to single bars of positive and negative contrast
figure('name','Q filter responses to single bars','NumberTitle','off')
set(gcf,'Position',[scrsz(3)/8 scrsz(4)/4 scrsz(3)/2 scrsz(4)/1.8])
for i = 1:max(size(Q_filters(:,1)))
    
    subplot(max(size(Q_filters(:,1))), 2, (i+i-1));
    plot(d,Q_filters(i,:)); hold on; axis equal; apos = get(gca,'Position');
    set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on',...
    'Position', [0 apos(2) apos(3) apos(4)],'GridLineStyle','-')
    if i == 1
        title('Filter');
    end
    
    subplot(max(size(Q_filters(:,1))), 2, i+i);
    bar(barposn,Q_R(i,:,1),'k'); hold on; bar(barposn,Q_R(i,:,2),'w'), apos = get(gca,'Position'); ylabel('Response');
    set(gca,'Xlim',[min(d)+min(d)/10 max(d)+max(d)/10],'Position', [apos(1)/1.5 apos(2) apos(3)*1.7 apos(4)]);
    if i == 1
        title('Response by Bar Position');
    end
    if i == max(size(HW_filters(:,1)))
        xlabel('Bar position (deg)')
    end    
end

%FIGURE 2- Plotting Energy model complex cell responses to bars positive and negative contrast
figure('name','Q model Complex Cell Response','NumberTitle','off')
set(gcf,'Position',[scrsz(3)/8 scrsz(4)/3 scrsz(3)/2 scrsz(4)/3])
bar(barposn,-CC_Q_R(:,1),'k'); hold on; bar(barposn,CC_Q_R(:,2),'w'); apos = get(gca,'Position');
title('Complex cell response of Quadrature model');xlabel('Bar position (deg)');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim',[-0.2 0.2]); ylabel('Response');

%%
%%% H&W filters
% FIGURE 3 -  %Plotting H&W simple filters and their responses to single bars
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
    bar(barposn,-(HW_r(i,:,1)),'k'); hold on; bar(barposn,HW_r(i,:,2),'w'); apos = get(gca,'Position');
    set(gca,'Position', [apos(1)/1.5 apos(2) apos(3)*1.5 apos(4)],'Xlim',[-1.1 1.1]); ylabel('Response');
    if i == 1
        title('Response by Bar Position');
    end
    if i == max(size(HW_filters(:,1)))
        xlabel('Bar position (deg)')
    end
end

% FIGURE 4 - Plotting H&W complex cell responses to bars positive and negative contrast
figure('name','H&W Complex Cell Response','NumberTitle','off')
set(gcf,'Position',[scrsz(3)/8 scrsz(4)/3 scrsz(3)/2 scrsz(4)/3])
bar(barposn,-CC_HW_r(:,1),'k'); hold on; bar(barposn,CC_HW_r(:,2),'w'); apos = get(gca,'Position');
title('Complex cell response of H&W model');xlabel('Bar position (deg)');ylabel('Response');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10]);

%% Two bar stimui responses

%IMPORTANT: evaluate this plotting section after placing the conditioning bar at
%position 7 (of 16) to generate the figure I will be presenting in class.
%FIGURE
%%%Energy model 'Complex Cell' response to two bars
figure('name','Q filter responses to two bars','NumberTitle','off') %Q CC response to two bars
set(gcf,'Position',[scrsz(3)/8 scrsz(4)/8 scrsz(3)/2 scrsz(4)/1.5])
subplot(2,1,1)
bar(barposn,-CC_Q_2b_R(:,1,2),'k'); hold on; bar(barposn,CC_Q_2b_R(:,2,2),'w'); apos = get(gca,'Position');
text(barposn(cond_bar),0.45,'*','FontSize',18); title('Bright Conditioning Bar');
set(gca,'Xlim',[min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim', [-0.6 0.6], 'Position', [apos(1) apos(2) apos(3) apos(4)]);
ylabel('Response')
subplot(2,1,2)
bar(barposn,-CC_Q_2b_R(:,1,1),'k'); hold on; bar(barposn,CC_Q_2b_R(:,2,1),'w'); apos = get(gca,'Position');
text(barposn(cond_bar),-0.45,'*','FontSize',18); title('Dark Conditioning Bar'); xlabel('Bar position (deg)');
ylabel('Response')
set(gca,'Xlim',[min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim', [-0.6 0.6], 'Position', [apos(1) apos(2) apos(3) apos(4)]);

% Highlighting underlying filters with relative responses
%Evaluate with conditioning bar in position 7, 16 bars total. Energy model filters
figure('name','Highlighting the EOC filters of the Q model','NumberTitle','off')
set(gcf,'Position',[scrsz(3)/16 scrsz(4)/8 scrsz(3)/1.2 scrsz(4)/1.5])
subplot(2,4,1)
bar(barposn,CC_Q_R(:,2),'w'); hold on; lim = get(gca,'Ylim');
title('Response to single bars');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10]); ylabel('Response');
subplot(2,4,5)
bar(barposn,CC_Q_R(:,1),'k'); xlabel('Bar position (deg)');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',lim);
subplot(2,4,2)
bar(barposn,((CC_Q_2b_R(:,2,2)-Q_CC_CB_R(2))/Q_CC_CB_R(2)),'w'); hold on;plot(d,EOCs1_gabr); lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),max(abs(lim))*.8,'*','FontSize',18); 
subplot(2,4,6)
bar(barposn,((CC_Q_2b_R(:,1,1)-Q_CC_CB_R(1))/Q_CC_CB_R(1)),'k'); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),-max(abs(lim))*.8,'*','FontSize',18); 
subplot(2,4,3)
bar(barposn,((CC_Q_2b_R(:,2,1)-Q_CC_CB_R(1))/Q_CC_CB_R(1)),'w'); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),-max(abs(lim))*.8,'*','FontSize',18); 
subplot(2,4,7)
bar(barposn,((CC_Q_2b_R(:,1,2)-Q_CC_CB_R(2))/Q_CC_CB_R(2)),'k'); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),max(abs(lim))*.8,'*','FontSize',18); 
subplot(2,4,4)
plot(d,EOCs1_gabr); hold on; axis equal; apos = get(gca,'Position');
text(barposn(cond_bar),0.45,'*','FontSize',18);
set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on')
subplot(2,4,8)
plot(d,EOCs2_gabr); hold on; axis equal; apos = get(gca,'Position');
text(barposn(cond_bar),-0.45,'*','FontSize',18);
set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting Q filters and their responses to two bars
figure('name','Q filter responses to two bars','NumberTitle','off')
set(gcf,'Position',[scrsz(3)/8 scrsz(4)/20 scrsz(3)/1.5 scrsz(4)/1.2]);
for i = 1:1:max(size(Q_filters(:,1)))
    
    subplot(max(size(Q_filters(:,1))), 3, 3*(i-1)+1);
    plot(d,Q_filters(i,:)); hold on; axis equal; apos = get(gca,'Position');
    set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1], 'Position', [apos(1)/4 apos(2) apos(3) apos(4)],...
    'YGrid','on');
    if i == 1
        title('Filter');
    end
    
    subplot(max(size(Q_filters(:,1))), 3, 2*i+i-1);
    bar(barposn,-(Q_2b_R(i,:,1,2)),'k'); hold on; bar(barposn,Q_2b_R(i,:,2,2),'w');
    text(barposn(cond_bar),0.25,'*','FontSize',18); apos = get(gca,'Position');
    set(gca,'Xlim',[min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim', [-0.4 0.4],'Position', [apos(1)/1.35 apos(2) apos(3)*1.4 apos(4)]);
        if i == 1
        title('Response with bright conditioning bar');
        end
        if i == max(size(Q_filters(:,1)))
            xlabel('Bar position (deg)')
        end
        
    subplot(max(size(Q_filters(:,1))), 3, i+2*i);
    bar(barposn,-(Q_2b_R(i,:,1,1)),'k'); hold on; bar(barposn,Q_2b_R(i,:,2,1),'w'); 
    text(barposn(cond_bar),-0.25,'*','FontSize',18); apos = get(gca,'Position');
    set(gca,'Xlim',[min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim', [-0.4 0.4],'Position', [apos(1)/1.05 apos(2) apos(3)*1.4 apos(4)]);
    if i == 1
        title('Response with dark conditioning bar');
    end
    if i == max(size(Q_filters(:,1)))
        xlabel('Bar position (deg)')
    end
end


%%
%FIGURE
%%% H&W model 'Complex Cell' response to two bars
figure('name','H&W Complex cell responses to two bars','NumberTitle','off') %H&W CC response to two bars
set(gcf,'Position',[scrsz(3)/8 scrsz(4)/8 scrsz(3)/2 scrsz(4)/1.5])
subplot(2,1,1)
bar(barposn,-CC_HW_2b_r(:,1,2),'k'); hold on; bar(barposn,CC_HW_2b_r(:,2,2),'w'); apos = get(gca,'Position');
text(barposn(cond_bar),0.45,'*','FontSize',18); title('Bright Conditioning Bar');
ylabel('Response')
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim', [-0.6 0.6], 'Position', [apos(1) apos(2) apos(3) apos(4)]);
subplot(2,1,2)
bar(barposn,-CC_HW_2b_r(:,1,1),'k'); hold on; bar(barposn,CC_HW_2b_r(:,2,1),'w'); apos = get(gca,'Position');
text(barposn(cond_bar),-0.45,'*','FontSize',18); title('Dark Conditioning Bar'); xlabel('Bar position (deg)');
ylabel('Response')
set(gca,'Xlim',[min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim', [-0.6 0.6], 'Position', [apos(1) apos(2) apos(3) apos(4)]);
% 
% subplot(3,1,3); bar(barposn,-CC_HW_2b_R(:,1,2)); hold on; bar(barposn,CC_HW_2b_R(:,2,2)); apos = get(gca,'Position');
% text(barposn(cond_bar),0.45,'*','FontSize',18);
% set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim', [-.05 .05], 'Position', [apos(1) apos(2) apos(3) apos(4)]);

% Evaluate with bar position 5, or 8 of 16 bars. Odd symmetrics(pos 5), Evens(pos8)
figure('name','Two-Bar Interaction Profiles','NumberTitle','off')
set(gcf,'Position',[scrsz(3)/16 scrsz(4)/8 scrsz(3)/1.2 scrsz(4)/1.5])
subplot(2,4,1)
bar(barposn,CC_HW_r(:,2),'w'); hold on; apos = get(gca,'Position'); lim = get(gca,'Ylim');
title('Response to single bars');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',lim); ylabel('Response');
subplot(2,4,5)
bar(barposn,CC_HW_r(:,1),'k'); xlabel('Bar position (deg)');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10]);
subplot(2,4,2)
bar(barposn,((CC_HW_2b_r(:,2,2)-HW_CC_CB_r(2))/HW_CC_CB_r(2)),'w'); hold on; lim = get(gca,'Ylim'); %plot(d,EOC_gabr);
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),max(abs(lim))*.8,'*','FontSize',22); 
subplot(2,4,7)
bar(barposn,((CC_HW_2b_r(:,1,2)-HW_CC_CB_r(2))/HW_CC_CB_r(2)),'k'); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),max(abs(lim))*.8,'*','FontSize',22);
subplot(2,4,3)
bar(barposn,((CC_HW_2b_r(:,2,1)-HW_CC_CB_r(1))/HW_CC_CB_r(1)),'w'); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),-max(abs(lim))*.8,'*','FontSize',22); 
subplot(2,4,6)
bar(barposn,((CC_HW_2b_r(:,1,1)-HW_CC_CB_r(1))/HW_CC_CB_r(1)),'k'); hold on; lim = get(gca,'Ylim');
set(gca,'Xlim', [min(d)+min(d)/10 max(d)+max(d)/10],'Ylim',[-max(abs(lim)) max(abs(lim))]);
text(barposn(cond_bar),-max(abs(lim))*.8,'*','FontSize',22); 
subplot(2,4,4)
plot(d,EOC_gabr); hold on; axis equal; apos = get(gca,'Position');
text(barposn(cond_bar),0.45,'*','FontSize',22);
set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on')
subplot(2,4,8)
plot(d,OSLn_gabr); hold on; axis equal; apos = get(gca,'Position');
text(barposn(cond_bar),-0.45,'*','FontSize',22);
set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1],'YGrid','on')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting H&W filters and their responses to two bars
figure('name','H&W filter responses to two bars','NumberTitle','off')
set(gcf,'Position',[scrsz(3)/8 scrsz(4)/20 scrsz(3)/1.5 scrsz(4)/1.2]);
for i = 1:1:max(size(HW_filters(:,1)))
    
    subplot(max(size(HW_filters(:,1))), 3, 3*(i-1)+1);
    plot(d,HW_filters(i,:)); hold on; axis equal; apos = get(gca,'Position');
    set(gca,'Ytick',[0],'Xlim', [-1.01 1.01], 'Ylim', [-1 1], 'Position', [apos(1)/4 apos(2) apos(3) apos(4)],...
    'YGrid','on');
    if i == 1
        title('Filter');
    end
 
    subplot(max(size(HW_filters(:,1))), 3, 2*i+i-1);
    bar(barposn,-(HW_2b_r(i,:,1,2)),'k'); hold on; bar(barposn,HW_2b_r(i,:,2,2),'w'); 
    text(barposn(cond_bar),0.25,'*','FontSize',18); apos = get(gca,'Position');
    set(gca,'Xlim',[min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim', [-0.4 0.4],'Position', [apos(1)/1.5 apos(2) apos(3)*1.4 apos(4)]);
    if i == 1
        title('Response with bright conditioning bar');
    end
    if i == max(size(HW_filters(:,1)))
        xlabel('Bar position (deg)')
    end
    
    subplot(max(size(HW_filters(:,1))), 3, i+2*i);
    bar(barposn,-(HW_2b_r(i,:,1,1)),'k'); hold on; bar(barposn,HW_2b_r(i,:,2,1),'w'); 
    text(barposn(cond_bar),-0.25,'*','FontSize',18); apos = get(gca,'Position');
    set(gca,'Xlim',[min(d)+min(d)/10 max(d)+max(d)/10], 'Ylim', [-0.4 0.4],'Position', [apos(1)/1.1 apos(2) apos(3)*1.4 apos(4)]);
    if i == 1
        title('Response with dark conditioning bar');
    end
    if i == max(size(HW_filters(:,1)))
        xlabel('Bar position (deg)')
    end
end
