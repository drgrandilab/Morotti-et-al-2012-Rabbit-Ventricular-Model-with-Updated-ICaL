function output = Morotti_et_al_2012_rabbit()
% This model was developed integrating a modified version of the L-type Ca
% current (LTCC) model by Mahajan et al. (Biophys J. 2008 Jan 15;94(2):392
% -410) into the framework of the Shannon-Bers rabbit ventricular model 
% (Biophys J. 2004 Nov;87(5):3351-71). The LTCC model was extended to 
% reflect more faithfully contributions of Ca- and voltage-dependent 
% inactivation (CDI and VDI) to total inactivation.
% 
% Reference:
% S. Morotti, E. Grandi, A. Summa, K.S. Ginsburg, D.M. Bers. Theoretical
% Study of L-type Ca2+ Current Inactivation Kinetics during Action
% Potential Repolarization and Early Afterdepolarizations. J Physiol. 2012
% Sep 15;590(Pt 18):4465-81.

clear all
close all
clc

%% Initial conditions
% Steady-state conditions obtained with the original and the new models
% simulating current-clamp stimulation at 1-Hz

%load yf_rabbit_Shannon_et_al_1Hz % Shannon et al 2004 (flagMica = 0)
load yf_rabbit_Morotti_et_al_1Hz % new model (flagMica = 1)

y0 = yfinal;

%% Single Run Simulation
tic
tspan = [0; 10e3]; % duration (ms)
HR = 1; % stimulation frequency (Hz)
options = odeset('RelTol',1e-5,'MaxStep',1,'Stats','on'); 
p = [HR]; % Parameter array for passing nondefault conditions
[t,y] = ode15s(@f,tspan,y0,options,p);
toc

%% Final conditions
yfinal = y(end,:);
output = yfinal;

%save yf_rabbit_Shannon_et_al_1Hz yfinal % Shannon et al 2004 (flagMica = 0)
%save yf_rabbit_Morotti_et_al_1Hz yfinal % new model (flagMica = 1)

%% Plot
currents = calcCurrents(t,y,p);

figure(1); set(gcf,'color','w')
subplot(4,1,1); hold on; plot(t,y(:,39)); ylabel('Vm (mV)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,2); hold on; plot(t,currents(:,1)); ylabel('ICa (A/F)')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,3); hold on; plot(t,y(:,38)); ylabel('[Ca]i (mM)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,4); hold on; plot(t,y(:,34)); ylabel('[Na]i (mM)');
set(gca,'box','off','tickdir','out','fontsize',12)
xlabel('Time (ms)')

%% Model of rabbit ventricular myocyte
function output = f(t,y,p,runType)
%% 66 state variables
% y(1) INa - gate m - HH model
% y(2) INa - gate h - HH model
% y(3) INa - gate j - HH model 
% y(4-7) ICa - HH model (NOT USED)
% y(8) Itos - gate x
% y(9) Itos - gate y
% y(10) Itof - gate x
% y(11) Itof - gate y
% y(12) IKr - gate x
% y(13) IKs
% y(14) RyR (R)
% y(15) RyR (O)
% y(16) RyR (I)
% y(17) Na Buffer cleft
% y(18) Na Buffer sl
% y(19) Ca Buffer myofilament
% y(20) Ca Buffer TnCH
% y(21) Mg Buffer TnCH
% y(22) Ca Buffer CaM
% y(23) Ca Buffer Myosin
% y(24) Mg Buffer Myosin
% y(25) Ca Buffer SRB
% y(26) Ca Buffer - low cleft
% y(27) Ca Buffer - low sl
% y(28) Ca Buffer - high cleft
% y(29) Ca Buffer - high sl
% y(30) Ca-Csqn
% y(31) [Ca]sr
% y(32) [Na]cleft
% y(33) [Na]sl
% y(34) [Na]i
% y(35) [K]i (NOT USED)
% y(36) [Ca]cleft
% y(37) [Ca]sl
% y(38) [Ca]i
% y(39) Vm
% y(40) Itos - gate r
% y(41) INaL - gate m
% y(42) INaL - gate h
% y(43-48) states ICa Markov Model (m1-cleft)           ¥
% y(49-54) states ICa Markov Model (m2-cleft)           ¥
% y(55-60) states ICa Markov Model (m1-sl)              ¥
% y(61-66) states ICa Markov Model (m2-sl)              ¥
% ¥ [y(i)=C2; y(ii)=C1; y(iii)=I1Ca; y(iv)=I2Ca; y(v)=I1Ba; y(vi)=I2Ba]

ydot = zeros(size(y));

%% Assign passed-in parameters
p_HR = p(1);

% Stimulation protocol ('pace' or 'step')
protocol = 'pace';
%protocol = 'step';

%% Model Flags
% Flag for choice of ICa model (0 or 1)
flagMica = 1;
% 0: HH model for ICa (original model from Shannon et al 2004)
% 1: Markov model for ICa (and other modifications from Morotti et al 2012)

% Flags for EGTA/BAPTA administration (EGTA or BAPTA with 1)
flag_EGTA = 0;
flag_BAPTA = 0;

% Flag for impaired VDI (with 1)
flag_redVDI = 0;

% Flag for CaM1234 expression (reduced CDI with 1)
flag_CaM1234 = 0;

% Flag for Ba-Ca substitution (0 or 1)
flag_Ba = 0; % ICa with 0 (default), IBa (and other modifications) with 1

% Flag for Ba current model (0 or 1)
flag_7_state = 1; % 7-state model with 1 (default), 5-state otherwise

%% Model Parameters
% Constants
R = 8314; % [J/kmol*K]  
Frdy = 96485; % [C/mol]  
Temp = 310; % [K]
FoRT = Frdy/R/Temp; % [1/mV]
Cmem = 1.3810e-10; % [F] membrane capacitance
Qpow = (Temp-310)/10;

% Cell geometry
cellLength = 100; % cell length [um]
cellRadius = 10.25; % cell radius [um]
junctionLength = 160e-3; % junc length [um]
junctionRadius = 15e-3; % junc radius [um]
distSLcyto = 0.45; % dist. SL to cytosol [um]
distJuncSL = 0.5; % dist. junc to SL [um]
DcaJuncSL = 1.64e-6; % Dca junc to SL [cm^2/sec]
DcaSLcyto = 1.22e-6; % Dca SL to cyto [cm^2/sec]
DnaJuncSL = 1.09e-5; % Dna junc to SL [cm^2/sec]
DnaSLcyto = 1.79e-5; % Dna SL to cyto [cm^2/sec] 
Vcell = pi*cellRadius^2*cellLength*1e-15;    % [L]
Vmyo = 0.65*Vcell; Vsr = 0.035*Vcell; Vsl = 0.02*Vcell; Vjunc = 0.0539*.01*Vcell; 
SAjunc = 20150*pi*2*junctionLength*junctionRadius; % [um^2]
SAsl = pi*2*cellRadius*cellLength; % [um^2]
%J_ca_juncsl = DcaJuncSL*SAjunc/distSLcyto*1e-10;% [L/msec] = 1.1074e-13
%J_ca_slmyo = DcaSLcyto*SAsl/distJuncSL*1e-10;  % [L/msec] = 1.5714e-12
%J_na_juncsl = DnaJuncSL*SAjunc/distSLcyto*1e-10;% [L/msec] = 7.36e-13
%J_na_slmyo = DnaSLcyto*SAsl/distJuncSL*1e-10;  % [L/msec] = 2.3056e-11
%J_ca_juncsl = DcaJuncSL*SAjunc/distJuncSL*1e-10;% [L/msec] = 1.1074e-13
%J_ca_slmyo = DcaSLcyto*SAsl/distSLcyto*1e-10;  % [L/msec] = 1.5714e-12
%J_na_juncsl = DnaJuncSL*SAjunc/distJuncSL*1e-10;% [L/msec] = 7.36e-13
%J_na_slmyo = DnaSLcyto*SAsl/distSLcyto*1e-10;  % [L/msec] = 2.3056e-11
% tau's from c-code, not used here
J_ca_juncsl =1/1.2134e12; % [L/msec] = 8.2413e-13
J_ca_slmyo = 1/2.68510e11; % [L/msec] = 3.2743e-12
J_na_juncsl = 1/(1.6382e12/3*100); % [L/msec] = 6.1043e-13
J_na_slmyo = 1/(1.8308e10/3*100); % [L/msec] = 5.4621e-11

% Fractional currents in compartments
Fjunc = 0.11; Fsl = 1-Fjunc;
Fjunc_CaL = 0.9; Fsl_CaL = 1-Fjunc_CaL;

% Fixed ion concentrations     
Cli = 15;  % Intracellular Cl [mM]
Clo = 150; % Extracellular Cl [mM]
Ko = 5.4;  % Extracellular K [mM]
Nao = 140; % Extracellular Na [mM]
Cao = 1.8; % Extracellular Ca [mM]
Mgi = 0.5; % Intracellular Mg [mM]

% Nernst Potentials
ena_junc = (1/FoRT)*log(Nao/y(32));     % [mV]
ena_sl = (1/FoRT)*log(Nao/y(33));       % [mV]
ek = (1/FoRT)*log(Ko/y(35));	        % [mV]
eca_junc = (1/FoRT/2)*log(Cao/y(36));   % [mV]
eca_sl = (1/FoRT/2)*log(Cao/y(37));     % [mV]
ecl = (1/FoRT)*log(Cli/Clo);            % [mV]

%% Na transport parameters
GNa = 16;
GNaL = 0.0045;
GNaB = 0.297e-3;    % [mS/uF] 
IbarNaK = 1.90719;  % [uA/uF]
KmNaip = 11;        % [mM]
KmKo = 1.5;         % [mM]
Q10NaK = 1.63;  
Q10KmNai = 1.39;

%% K current parameters
pNaK = 0.01833;      
GtoSlow = 0.06;     % [mS/uF]
GtoFast = 0.02;     % [mS/uF] 
gkp = 0.001;

%% Cl current parameters
GClCa = 0.109625;   % [mS/uF]
GClB = 9e-3;        % [mS/uF]
KdClCa = 100e-3;    % [mM]

%% Ca transport parameters
Q10CaL = 1.8;       
% ICaL - HH model
pCa = 5.4e-4;       % [cm/sec]
pK = 2.7e-7;        % [cm/sec]
pNa = 1.5e-8;       % [cm/sec]
% ICaL - Markov model
Zca = 2;
Pca = 100*24.3e-6;  % [cm/s] 10*0.45*5.4e-6
aff = 1;
gammaCai = 0.0341;
gammaCao = 0.341;
Zk = 1;
Pk = 100*12.15e-9;  % [cm/s]
gammaKi = 0.75;
gammaKo = 0.75;
Zna = 1;
Pna = 100*0.675e-9; % [cm/s]
gammaNai = 0.75;
gammaNao = 0.75;

IbarNCX = 1.0*9.0;  % [uA/uF]
KmCai = 3.59e-3;    % [mM]
KmCao = 1.3;        % [mM]
KmNai = 12.29;      % [mM]
KmNao = 87.5;       % [mM]
ksat = 0.27;        % [none]  
nu = 0.35;          % [none]
Kdact = 0.256e-3;   % [mM] 
Q10NCX = 1.57;      % [none]
IbarSLCaP = 0.0673; % IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
KmPCa = 0.5e-3;     % [mM] 
GCaB = 2.513e-4;    % [uA/uF] 
Q10SLCaP = 2.35;    % [none]

Q10SRCaP = 2.6;          % [none]
Vmax_SRCaP = 5.3114e-3;  % [mM/msec] (286 umol/L cytosol/sec)
Kmf = 0.246e-3;          % [mM] default
%Kmf = 0.175e-3;          % [mM]
Kmr = 1.7;               % [mM]L cytosol
hillSRCaP = 1.787;       % [mM]
ks = 25;                 % [1/ms]      
koCa = 10;               % [mM^-2 1/ms]   %default 10   modified 20
kom = 0.06;              % [1/ms]     
kiCa = 0.5;              % [1/mM/ms]
kim = 0.005;             % [1/ms]
ec50SR = 0.45;           % [mM]
G_SRleak = 5.348e-6;     % [1/ms] 

%% Buffering parameters
Bmax_Naj = 7.561;       % [mM] % Bmax_Naj = 3.7; (c-code difference?)  % Na buffering
Bmax_Nasl = 1.65;       % [mM]
koff_na = 1e-3;         % [1/ms]
kon_na = 0.1e-3;        % [1/mM/ms]
Bmax_TnClow = 70e-3;    % [mM]                      % TnC low affinity
koff_tncl = 19.6e-3;    % [1/ms] 
kon_tncl = 32.7;        % [1/mM/ms]
Bmax_TnChigh = 140e-3;  % [mM]                      % TnC high affinity 
koff_tnchca = 0.032e-3; % [1/ms] 
kon_tnchca = 2.37;      % [1/mM/ms]
koff_tnchmg = 3.33e-3;  % [1/ms] 
kon_tnchmg = 3e-3;      % [1/mM/ms]
Bmax_CaM = 24e-3;       % [mM] **? about setting to 0 in c-code**   % CaM buffering
koff_cam = 238e-3;      % [1/ms] 
kon_cam = 34;           % [1/mM/ms]
Bmax_myosin = 140e-3;   % [mM]                      % Myosin buffering
koff_myoca = 0.46e-3;   % [1/ms]
kon_myoca = 13.8;       % [1/mM/ms]
koff_myomg = 0.057e-3;  % [1/ms]
kon_myomg = 0.0157;     % [1/mM/ms]
Bmax_SR = 19*.9e-3;     % [mM] (Bers text says 47e-3) 19e-3
koff_sr = 60e-3;        % [1/ms]
kon_sr = 100;           % [1/mM/ms]
Bmax_SLlowsl = 37.4e-3*Vmyo/Vsl;        % [mM]    % SL buffering
Bmax_SLlowj = 4.6e-3*Vmyo/Vjunc*0.1;    % [mM]    %Fei *0.1!!! junction reduction factor
koff_sll = 1300e-3;     % [1/ms]
kon_sll = 100;          % [1/mM/ms]
Bmax_SLhighsl = 13.4e-3*Vmyo/Vsl;       % [mM] 
Bmax_SLhighj = 1.65e-3*Vmyo/Vjunc*0.1;  % [mM] %Fei *0.1!!! junction reduction factor
koff_slh = 30e-3;       % [1/ms]
kon_slh = 100;          % [1/mM/ms]
Bmax_Csqn = 140e-3*Vmyo/Vsr;            % [mM] % Bmax_Csqn = 2.6;      % Csqn buffering
koff_csqn = 65;         % [1/ms] 
kon_csqn = 100;         % [1/mM/ms] 

%% Parameteres modified in presence of Ba, Ca buffers or mutations
if flag_EGTA == 1,     
    % [Ca]c clamped
    ks = 0; % no SR Ca release
    Vmax_SRCaP = 0; % no SR Ca uptake
    G_SRleak = 0; % no SR Ca leak
end

if flag_BAPTA == 1,  
    % [Ca]c,sl,j clamped to Cai0 [mM]
    Cai0 = 1*1e-4; % 1*(100 nM) % DEFINE HERE FREE [Ca]i 
    y(36)=Cai0; y(37)=Cai0; y(38)=Cai0;
    ks = 0; % no SR Ca release
    Vmax_SRCaP = 0; % no SR Ca uptake
    G_SRleak = 0; % no SR Ca leak
end

if flag_Ba == 1,
	Pca = (1700/4200*25/9*0.85)*100*24.3e-6; % [cm/s]
    aff = 1/8;
    % [Ca]c clamped
    gks_ca_dep = 0; % no Ca-dep component of IKs conductance
    GClCa = 0; % no Ca-dep Cl current
    IbarNCX = 0; % no NCX current
    GCaB = 0; % no background Ca current
    IbarSLCaP = 0; % no sarcolemmal Ca pump current
    ks = 0; % no SR Ca release
    Vmax_SRCaP = 0; % no SR Ca uptake
    G_SRleak = 0; % no SR Ca leak
else
    gks_ca_dep = 1;
end

if flag_7_state == 0,
    aff = 1/800000;
	flag_7_state = 1e-15;
end

if flag_redVDI == 1, % impaired VDI
    redVDI = 0.30;
end

if flag_CaM1234 == 1, % impaired CDI
    alpha_CaM1234 = 0.25; % DEFINE HERE FRACTION OF MODE-2 LTCCS
    aff_m2 = 1/800000;
    flag_7_state_m2 = 1e-15;
else
    alpha_CaM1234 = 0;
    aff_m2 = aff;
    flag_7_state_m2 = flag_7_state;
end

%% I_Na: Na Current
% Fast I_Na
am = 0.32*(y(39)+47.13)/(1-exp(-0.1*(y(39)+47.13)));
bm = 0.08*exp(-y(39)/11);
if y(39) >= -40,
    ah = 0; aj = 0;
    bh = 1/(0.13*(1+exp(-(y(39)+10.66)/11.1)));
    bj = 0.3*exp(-2.535e-7*y(39))/(1+exp(-0.1*(y(39)+32)));
else
    ah = 0.135*exp((80+y(39))/-6.8);
    bh = 3.56*exp(0.079*y(39))+3.1e5*exp(0.35*y(39));
    aj = (-1.2714e5*exp(0.2444*y(39))-3.474e-5*exp(-0.04391*y(39)))*(y(39)+37.78)/(1+exp(0.311*(y(39)+79.23)));
    bj = 0.1212*exp(-0.01052*y(39))/(1+exp(-0.1378*(y(39)+40.14)));
end
ydot(1) = am*(1-y(1))-bm*y(1);
ydot(2) = ah*(1-y(2))-bh*y(2);
ydot(3) = aj*(1-y(3))-bj*y(3);

% Late I_Na
aml = 0.32*(y(39)+47.13)/(1-exp(-0.1*(y(39)+47.13))); % = am
bml = 0.08*exp(-y(39)/11);                            % = bm
hlinf = 1/(1+exp((y(39)+91)/6.1));
tauhl = 600;
ydot(41) = aml*(1-y(41))-bml*y(41);
ydot(42) = (hlinf-y(42))/tauhl;

I_NaL_junc = Fjunc*GNaL*y(41)^3*y(42)*(y(39)-ena_junc);
I_NaL_sl = Fsl*GNaL*y(41)^3*y(42)*(y(39)-ena_sl);
I_NaL = I_NaL_junc+I_NaL_sl;

% Total I_Na
I_Na_junc = Fjunc*GNa*y(1)^3*y(2)*y(3)*(y(39)-ena_junc)+I_NaL_junc ;
I_Na_sl = Fsl*GNa*y(1)^3*y(2)*y(3)*(y(39)-ena_sl)+I_NaL_sl ;
I_Natot = I_Na_junc+I_Na_sl;

%% I_nabk: Na Background Current
I_nabk_junc = Fjunc*GNaB*(y(39)-ena_junc);
I_nabk_sl = Fsl*GNaB*(y(39)-ena_sl);
I_nabk = I_nabk_junc+I_nabk_sl;

%% I_nak: Na/K Pump Current
sigma = (exp(Nao/67.3)-1)/7;
fnak = 1/(1+0.1245*exp(-0.1*y(39)*FoRT)+0.0365*sigma*exp(-y(39)*FoRT));
I_nak_junc = Fjunc*IbarNaK*fnak*Ko /(1+(KmNaip/y(32))^4) /(Ko+KmKo);
I_nak_sl = Fsl*IbarNaK*fnak*Ko /(1+(KmNaip/y(33))^4) /(Ko+KmKo);
I_nak = I_nak_junc+I_nak_sl;

%% I_kr: Rapidly Activating K Current
gkr = 1*0.03*sqrt(Ko/5.4);
xrss = 1/(1+exp(-(y(39)+50)/7.5));
tauxr = 1/(1.38e-3*(y(39)+7)/(1-exp(-0.123*(y(39)+7)))+6.1e-4*(y(39)+10)/(exp(0.145*(y(39)+10))-1));
ydot(12) = (xrss-y(12))/tauxr;
rkr = 1/(1+exp((y(39)+33)/22.4));
I_kr = gkr*y(12)*rkr*(y(39)-ek);

%% I_ks: Slowly Activating K Current
pcaks_junc = -log10(y(36))+3.0; 
pcaks_sl = -log10(y(37))+3.0;  
gks_junc = 0.07*(0.057 + gks_ca_dep*0.19/(1+ exp((-7.2+pcaks_junc)/0.6)));
gks_sl = 0.07*(0.057 + gks_ca_dep*0.19/(1+ exp((-7.2+pcaks_sl)/0.6))); 
eks = (1/FoRT)*log((Ko+pNaK*Nao)/(y(35)+pNaK*y(34)));	
xsss = 1/(1+exp(-(y(39)-1.5)/16.7));
tauxs = 1/(7.19e-5*(y(39)+30)/(1-exp(-0.148*(y(39)+30)))+1.31e-4*(y(39)+30)/(exp(0.0687*(y(39)+30))-1)); 
ydot(13) = (xsss-y(13))/tauxs;
I_ks_junc = Fjunc*gks_junc*y(13)^2*(y(39)-eks);
I_ks_sl = Fsl*gks_sl*y(13)^2*(y(39)-eks);
I_ks = I_ks_junc+I_ks_sl;

%% I_kp: Plateau K current
kp_kp = 1/(1+exp(7.488-y(39)/5.98));
I_kp_junc = Fjunc*gkp*kp_kp*(y(39)-ek);
I_kp_sl = Fsl*gkp*kp_kp*(y(39)-ek);
I_kp = I_kp_junc+I_kp_sl;

%% I_to: Transient Outward K Current (slow and fast components)
xtoss = 1/(1+exp(-(y(39)+3.0)/15));
ytoss = 1/(1+exp((y(39)+33.5)/10));
rtoss = 1/(1+exp((y(39)+33.5)/10));
tauxtos = 9/(1+exp((y(39)+3.0)/15))+0.5;
tauytos = 3e3/(1+exp((y(39)+60.0)/10))+30;
%tauytos = 182/(1+exp((y(39)+33.5)/10))+1;
taurtos = 2.8e3/(1+exp((y(39)+60.0)/10))+220; %Fei changed here!! time-dependent gating variable
%taurtos = 8085/(1+exp((y(39)+33.5)/10))+313;
ydot(8) = (xtoss-y(8))/tauxtos;
ydot(9) = (ytoss-y(9))/tauytos;
ydot(40)= (rtoss-y(40))/taurtos; %Fei changed here!! time-dependent gating variable
I_tos = GtoSlow*y(8)*(y(9)+0.5*y(40))*(y(39)-ek); % [uA/uF]

tauxtof = 3.5*exp(-y(39)*y(39)/30/30)+1.5;
%tauxtof = 3.5*exp(-((y(39)+3)/30)^2)+1.5;
tauytof = 20.0/(1+exp((y(39)+33.5)/10))+20.0;
%tauytof = 20.0/(1+exp((y(39)+33.5)/10))+20.0;
ydot(10) = (xtoss-y(10))/tauxtof;
ydot(11) = (ytoss-y(11))/tauytof;
I_tof = GtoFast*y(10)*y(11)*(y(39)-ek);

I_to = I_tos + I_tof;

%% I_ki: Time-Independent K Current
aki = 1.02/(1+exp(0.2385*(y(39)-ek-59.215)));
bki =(0.49124*exp(0.08032*(y(39)+5.476-ek)) + exp(0.06175*(y(39)-ek-594.31))) /(1 + exp(-0.5143*(y(39)-ek+4.753)));
kiss = aki/(aki+bki);
I_ki = 1*0.9*sqrt(Ko/5.4)*kiss*(y(39)-ek);

%% I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/y(36))*(y(39)-ecl);
I_ClCa_sl = Fsl*GClCa/(1+KdClCa/y(37))*(y(39)-ecl);
I_ClCa = I_ClCa_junc+I_ClCa_sl;

I_Clbk = GClB*(y(39)-ecl);

%% I_Ca: L-type Calcium Current (HH model, flagMica=0)
dss = 1/(1+exp(-(y(39)+14.5)/6.0));
taud = dss*(1-exp(-(y(39)+14.5)/6.0))/(0.035*(y(39)+14.5));
fss = 1/(1+exp((y(39)+35.06)/3.6))+0.6/(1+exp((50-y(39))/20));
tauf = 1/(0.0197*exp( -(0.0337*(y(39)+14.5))^2 )+0.02);
ydot(4) = (dss-y(4))/taud;
ydot(5) = (fss-y(5))/tauf;
ydot(6) = 1.7*y(36)*(1-y(6))-11.9e-3*y(6); % fCa_junc   koff!!!!!!!!
ydot(7) = 1.7*y(37)*(1-y(7))-11.9e-3*y(7); % fCa_sl
fcaCaMSL = 0.1/(1+(0.01/y(37)));
fcaCaj = 0.1/(1+(0.01/y(36)));
fcaCaMSL = 0;
fcaCaj = 0;
%y(6)=0;
%y(7)=0;
ibarca_j = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(36)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibarca_sl = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(37)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibark = pK*(y(39)*Frdy*FoRT)*(0.75*y(35)*exp(y(39)*FoRT)-0.75*Ko) /(exp(y(39)*FoRT)-1);
ibarna_j = pNa*(y(39)*Frdy*FoRT) *(0.75*y(32)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);
ibarna_sl = pNa*(y(39)*Frdy*FoRT) *(0.75*y(33)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);

I_Ca_junc1 = (Fjunc_CaL*ibarca_j*y(4)*y(5)*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.45*1;
I_Ca_sl1 = (Fsl_CaL*ibarca_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*0.45*1;
%I_Ca = I_Ca_junc+I_Ca_sl;
I_CaK1 = (ibark*y(4)*y(5)*(Fjunc_CaL*(fcaCaj+(1-y(6)))+Fsl_CaL*(fcaCaMSL+(1-y(7))))*Q10CaL^Qpow)*0.45*1;
I_CaNa_junc1 = (Fjunc_CaL*ibarna_j*y(4)*y(5)*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.45*1;
I_CaNa_sl1 = (Fsl_CaL*ibarna_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*.45*1;
%I_CaNa = I_CaNa_junc+I_CaNa_sl;
%I_Catot = I_Ca+I_CaK+I_CaNa;

%% I_Ca: L-type Calcium Current (Markov model, flagMica=1)
% LTCC Current - Input: state variables (12-24)
% junc, mode-1
Pc2_LCCj_m1=y(43); Pc1_LCCj_m1=y(44); Pi1Ca_LCCj_m1=y(45);
Pi2Ca_LCCj_m1=y(46); Pi1Ba_LCCj_m1=y(47); Pi2Ba_LCCj_m1=y(48);
% junc, mode-2
Pc2_LCCj_m2=y(49); Pc1_LCCj_m2=y(50); Pi1Ca_LCCj_m2=y(51);
Pi2Ca_LCCj_m2=y(52); Pi1Ba_LCCj_m2=y(53); Pi2Ba_LCCj_m2=y(54);
% SL, mode-1
Pc2_LCCsl_m1=y(55); Pc1_LCCsl_m1=y(56); Pi1Ca_LCCsl_m1=y(57);
Pi2Ca_LCCsl_m1=y(58); Pi1Ba_LCCsl_m1=y(59); Pi2Ba_LCCsl_m1=y(60);
% SL, mode-2
Pc2_LCCsl_m2=y(61); Pc1_LCCsl_m2=y(62); Pi1Ca_LCCsl_m2=y(63);
Pi2Ca_LCCsl_m2=y(64); Pi1Ba_LCCsl_m2=y(65); Pi2Ba_LCCsl_m2=y(66);
% dependent state-variables (PO)
Po_LCCj_m1 =  1-Pc2_LCCj_m1-Pc1_LCCj_m1-Pi1Ca_LCCj_m1-Pi2Ca_LCCj_m1-Pi1Ba_LCCj_m1-Pi2Ba_LCCj_m1;
Po_LCCj_m2 =  1-Pc2_LCCj_m2-Pc1_LCCj_m2-Pi1Ca_LCCj_m2-Pi2Ca_LCCj_m2-Pi1Ba_LCCj_m2-Pi2Ba_LCCj_m2;
Po_LCCsl_m1 = 1-Pc2_LCCsl_m1-Pc1_LCCsl_m1-Pi1Ca_LCCsl_m1-Pi2Ca_LCCsl_m1-Pi1Ba_LCCsl_m1-Pi2Ba_LCCsl_m1;
Po_LCCsl_m2 = 1-Pc2_LCCsl_m2-Pc1_LCCsl_m2-Pi1Ca_LCCsl_m2-Pi2Ca_LCCsl_m2-Pi1Ba_LCCsl_m2-Pi2Ba_LCCsl_m2;

% LTCC Current - Input: intracellular concentrations
cajLCC = y(36); kjLCC = y(35); najLCC = y(32);
caslLCC = y(37); kslLCC = y(35); naslLCC = y(33);

% LTCC Current - Temp-dependent Parameters
ICa_scale = 1*Q10CaL^Qpow;
ICa_speed = 1;

% LTCC Current - Fixed Parameters
cpt = 3.75e-3;     % [mM]
cat = 7.617e-3;    % [mM]
s1o = 0.0182688;   % [1/ms]
k1o = 0.024168;    % [1/ms]
k2o = 0.000103615; % [1/ms]
sp0 = 1.5;
sp1 = 3;           % [ms]
sp2 = 40;          % [mV]
sp3 = 3;           % [mV]
sp4 = 4;           % [mV]
sp5 = 11.32;       % [mV]
sp6 = 15.6;        % [mV]
sp7 = 10;          % [ms]
sp8 = 4954;        % [ms]
sp9 = 78.0329;     % [ms]
sp10 = 0.1;        % [ms]
aR2 = 1;
sR2 = -2;          % [mV]
pR2 = 0.145;       % [1/mV]
aT2 = 1;           % [1/ms]
sT2 = -1000;       % [mV]
pT2 = 0.100;       % [1/mV]
aR1 = 0.09091;     
sR1 = -1000;       % [mV]
pR1 = 0.100;       % [1/mV]
aT1 = 0.30303;     % [1/ms]
sT1 = -1000;       % [mV]
pT1 = 0.100;       % [1/mV]
aRv2 = 0.9;
sRv2 = -29;        % [mV]
pRv2 = 0.135;      % [1/mV]
aTv2 = 500;        % [1/ms]
sTv2 = -25;        % [mV]
pTv2 = 0.050;      % [1/mV]
aRv1 = 0.85;
sRv1 = -180;       % [mV]
pRv1 = 0.090;      % [1/mV]
aTv1 = 270;        % [1/ms]
sTv1 = -180;       % [mV]
pTv1 = 0.100;      % [1/mV]
aTrev1 = 205.12;   % [1/ms]
sTrev1 = -65;      % [mV]
pTrev1 = 0.100;    % [1/mV]
aTrev2 = 7e8;      % [1/ms]
sTrev2 = 60;       % [mV]
pTrev2 = 0.130;    % [1/mV]

if flag_redVDI == 1,
    aRv2 = redVDI*aRv2;
    aRv1 = redVDI*aRv1;
end

% junc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% I_Ca JUNC - mode-1
% Voltage- and Ca-dependent Parameters
fcp=1/(1+(cpt/cajLCC/aff)^3);        % Ca-dep
tca=sp9/(1+(cajLCC*aff/cat)^4)+sp10; % Ca-dep
R2=aR2/(1+exp(-(y(39)-sR2)*pR2));
T2=aT2/(1+exp(-(y(39)-sT2)*pT2));
PT=1-(1/(1+exp(-(y(39)+sp2)/sp3)));
R1=aR1/(1+exp(-(y(39)-sR1)*pR1));
T1=aT1/(1+exp(-(y(39)-sT1)*pT1));
RV=sp7+sp8*exp(y(39)/sp6);
Pr=1-(1/(1+exp(-(y(39)+sp2)/sp4)));
Pq=1+sp0/(1+exp(-(y(39)+sp2)/sp4));
    TCa=Pq*((RV-tca)*Pr+tca);
Ps=1/(1+exp(-(y(39)+sp2)/sp5));
Rv1=aRv1/(1+exp(-(y(39)-sRv1)*pRv1));
Tv1=aTv1/(1+exp(-(y(39)-sTv1)*pTv1));
Rv2=aRv2/(1+exp(-(y(39)-sRv2)*pRv2));
Tv2=aTv2/(1+exp(-(y(39)-sTv2)*pTv2));
Trev1=aTrev1/(1+exp(-(y(39)-sTrev1)*pTrev1));
    Frev1=(1-Rv1)/Rv1*R1/(1-R1);
Trev2=aTrev2/(1+exp(-(y(39)-sTrev2)*pTrev2));
    Frev2=(1-Rv2)/Rv2*R2/(1-R2)*Rv1/(1-Rv1);
% Transition Rates (20 rates)
alphaLCC=ICa_speed*R2/T2;
betaLCC=ICa_speed*(1-R2)/T2;
r1=ICa_speed*R1/T1;
r2=ICa_speed*(1-R1)/T1;
k1=flag_7_state*ICa_speed*k1o*fcp;
k2=ICa_speed*k2o;
k3=ICa_speed*PT/sp1;
%k4
k5=ICa_speed*(1-Ps)/TCa;
k6=flag_7_state*ICa_speed*fcp*Ps/TCa;
s1=flag_7_state*ICa_speed*s1o*fcp;
%s2
k1p=ICa_speed*Rv1/Tv1;
k2p=ICa_speed*(1-Rv1)/Tv1;
k3p=ICa_speed*1/(Trev2*(1+Frev2));
    k4p=Frev2*k3p;                            % REV
k5p=ICa_speed*(1-Rv2)/Tv2;
k6p=ICa_speed*Rv2/Tv2;
s1p=ICa_speed*1/(Trev1*(1+Frev1));
    s2p=Frev1*s1p;                            % REV
    k4=k3*(alphaLCC/betaLCC)*(k1/k2)*(k5/k6); % REV
    s2=s1*(k2/k1)*(r1/r2);                    % REV
    
% State transitions for mode-1 junctional LCCs
dPc2_LCCj_m1 = betaLCC*Pc1_LCCj_m1 + k5*Pi2Ca_LCCj_m1 + k5p*Pi2Ba_LCCj_m1 - (k6+k6p+alphaLCC)*Pc2_LCCj_m1;                      % C2_m1j
dPc1_LCCj_m1 = alphaLCC*Pc2_LCCj_m1 + k2*Pi1Ca_LCCj_m1 + k2p*Pi1Ba_LCCj_m1 + r2*Po_LCCj_m1 - (r1+betaLCC+k1+k1p)*Pc1_LCCj_m1;   % C1_m1j
dPi1Ca_LCCj_m1 = k1*Pc1_LCCj_m1 + k4*Pi2Ca_LCCj_m1 + s1*Po_LCCj_m1 - (k2+k3+s2)*Pi1Ca_LCCj_m1;                              % I1Ca_m1j
dPi2Ca_LCCj_m1 = k3*Pi1Ca_LCCj_m1 + k6*Pc2_LCCj_m1 - (k4+k5)*Pi2Ca_LCCj_m1;                                                 % I2Ca_m1j
dPi1Ba_LCCj_m1 = k1p*Pc1_LCCj_m1 + k4p*Pi2Ba_LCCj_m1 + s1p*Po_LCCj_m1 - (k2p+k3p+s2p)*Pi1Ba_LCCj_m1;                        % I1Ba_m1j
dPi2Ba_LCCj_m1 = k3p*Pi1Ba_LCCj_m1 + k6p*Pc2_LCCj_m1 - (k5p+k4p)*Pi2Ba_LCCj_m1;                                             % I2Ba_m1j

ibarca_jm1 = Zca^2*Pca*Frdy*FoRT*y(39)/(exp(Zca*y(39)*FoRT)-1)*(gammaCai*cajLCC*exp(Zca*y(39)*FoRT)-gammaCao*Cao);
I_Ca_junc_m1 = Fjunc_CaL*(ibarca_jm1*Po_LCCj_m1)*ICa_scale;

ibarna_jm1 = Zna^2*Pna*Frdy*FoRT*y(39)/(exp(Zna*y(39)*FoRT)-1)*(gammaNai*najLCC*exp(Zna*y(39)*FoRT)-gammaNao*Nao);
I_Na_junc_m1 = Fjunc_CaL*(ibarna_jm1*Po_LCCj_m1)*ICa_scale;

ibark_jm1 = Zk^2*Pk*Frdy*FoRT*y(39)/(exp(Zk*y(39)*FoRT)-1)*(gammaKi*kjLCC*exp(Zk*y(39)*FoRT)-gammaKo*Ko);
I_K_junc_m1 = Fjunc_CaL*(ibark_jm1*Po_LCCj_m1)*ICa_scale;

%% I_Ca JUNC - mode-2 (channels expressing CaM1234, impaired CDI)
% Re-define all parameters as mode-2 specific parameters
% Voltage- and Ca-dependent Parameters
fcp=1/(1+(cpt/cajLCC/aff_m2)^3);        % Ca-dep
tca=sp9/(1+(cajLCC*aff_m2/cat)^4)+sp10; % Ca-dep
R2=aR2/(1+exp(-(y(39)-sR2)*pR2));
T2=aT2/(1+exp(-(y(39)-sT2)*pT2));
PT=1-(1/(1+exp(-(y(39)+sp2)/sp3)));
R1=aR1/(1+exp(-(y(39)-sR1)*pR1));
T1=aT1/(1+exp(-(y(39)-sT1)*pT1));
RV=sp7+sp8*exp(y(39)/sp6);
Pr=1-(1/(1+exp(-(y(39)+sp2)/sp4)));
Pq=1+sp0/(1+exp(-(y(39)+sp2)/sp4));
    TCa=Pq*((RV-tca)*Pr+tca);
Ps=1/(1+exp(-(y(39)+sp2)/sp5));
Rv1=aRv1/(1+exp(-(y(39)-sRv1)*pRv1));
Tv1=aTv1/(1+exp(-(y(39)-sTv1)*pTv1));
Rv2=aRv2/(1+exp(-(y(39)-sRv2)*pRv2));
Tv2=aTv2/(1+exp(-(y(39)-sTv2)*pTv2));
Trev1=aTrev1/(1+exp(-(y(39)-sTrev1)*pTrev1));
    Frev1=(1-Rv1)/Rv1*R1/(1-R1);
Trev2=aTrev2/(1+exp(-(y(39)-sTrev2)*pTrev2));
    Frev2=(1-Rv2)/Rv2*R2/(1-R2)*Rv1/(1-Rv1);
% Re-calculate transition rates
alphaLCCm2=ICa_speed*R2/T2;
betaLCCm2=ICa_speed*(1-R2)/T2;
r1m2=ICa_speed*R1/T1;
r2m2=ICa_speed*(1-R1)/T1;
k1m2=flag_7_state_m2*ICa_speed*k1o*fcp;
k2m2=ICa_speed*k2o;
k3m2=ICa_speed*PT/sp1;
%k4
k5m2=ICa_speed*(1-Ps)/TCa;
k6m2=flag_7_state_m2*ICa_speed*fcp*Ps/TCa;
s1m2=flag_7_state_m2*ICa_speed*s1o*fcp;
%s2
k1pm2=ICa_speed*Rv1/Tv1;
k2pm2=ICa_speed*(1-Rv1)/Tv1;
k3pm2=ICa_speed*1/(Trev2*(1+Frev2));
    k4pm2=Frev2*k3pm2;                                        % REV
k5pm2=ICa_speed*(1-Rv2)/Tv2;
k6pm2=ICa_speed*Rv2/Tv2;
s1pm2=ICa_speed*1/(Trev1*(1+Frev1));
    s2pm2=Frev1*s1pm2;                                        % REV
    k4m2=k3m2*(alphaLCCm2/betaLCCm2)*(k1m2/k2m2)*(k5m2/k6m2); % REV
    s2m2=s1m2*(k2m2/k1m2)*(r1m2/r2m2);                        % REV
    
% State transitions for mode-2 junctional LCCs
dPc2_LCCj_m2 = betaLCCm2*Pc1_LCCj_m2 + k5m2*Pi2Ca_LCCj_m2 + k5pm2*Pi2Ba_LCCj_m2 - (k6m2+k6pm2+alphaLCCm2)*Pc2_LCCj_m2;                          % C2_m2j
dPc1_LCCj_m2 = alphaLCCm2*Pc2_LCCj_m2 + k2m2*Pi1Ca_LCCj_m2 + k2pm2*Pi1Ba_LCCj_m2 + r2m2*Po_LCCj_m2 - (r1m2+betaLCCm2+k1m2+k1pm2)*Pc1_LCCj_m2;   % C1_m2j
dPi1Ca_LCCj_m2 = k1m2*Pc1_LCCj_m2 + k4m2*Pi2Ca_LCCj_m2 + s1m2*Po_LCCj_m2 - (k2m2+k3m2+s2m2)*Pi1Ca_LCCj_m2;                                  % I1Ca_m2j
dPi2Ca_LCCj_m2 = k3m2*Pi1Ca_LCCj_m2 + k6m2*Pc2_LCCj_m2 - (k4m2+k5m2)*Pi2Ca_LCCj_m2;                                                         % I2Ca_m2j
dPi1Ba_LCCj_m2 = k1pm2*Pc1_LCCj_m2 + k4pm2*Pi2Ba_LCCj_m2 + s1pm2*Po_LCCj_m2 - (k2pm2+k3pm2+s2pm2)*Pi1Ba_LCCj_m2;                            % I1Ba_m2j
dPi2Ba_LCCj_m2 = k3pm2*Pi1Ba_LCCj_m2 + k6pm2*Pc2_LCCj_m2 - (k5pm2+k4pm2)*Pi2Ba_LCCj_m2;                                                     % I2Ba_m2j

ibarca_jm2 = Zca^2*Pca*Frdy*FoRT*y(39)/(exp(Zca*y(39)*FoRT)-1)*(gammaCai*cajLCC*exp(Zca*y(39)*FoRT)-gammaCao*Cao);
I_Ca_junc_m2 = Fjunc_CaL*(ibarca_jm2*Po_LCCj_m2)*ICa_scale;

ibarna_jm2 = Zna^2*Pna*Frdy*FoRT*y(39)/(exp(Zna*y(39)*FoRT)-1)*(gammaNai*najLCC*exp(Zna*y(39)*FoRT)-gammaNao*Nao);
I_Na_junc_m2 = Fjunc_CaL*(ibarna_jm2*Po_LCCj_m2)*ICa_scale;

ibark_jm2 = Zk^2*Pk*Frdy*FoRT*y(39)/(exp(Zk*y(39)*FoRT)-1)*(gammaKi*kjLCC*exp(Zk*y(39)*FoRT)-gammaKo*Ko);
I_K_junc_m2 = Fjunc_CaL*(ibark_jm2*Po_LCCj_m2)*ICa_scale;

%junc_mode2 = 0; % Sum up total fraction of CKII and PKA-shifted mode 2 channels
junc_mode2 = alpha_CaM1234;

% Total junctional I_Ca
I_Ca_junc2 = (1-junc_mode2)*I_Ca_junc_m1 + junc_mode2*I_Ca_junc_m2;
I_Na_junc2 = (1-junc_mode2)*I_Na_junc_m1 + junc_mode2*I_Na_junc_m2;
I_K_junc2 = (1-junc_mode2)*I_K_junc_m1 + junc_mode2*I_K_junc_m2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% I_Ca SL - mode-1
% Voltage- and Ca-dependent Parameters
fcp=1/(1+(cpt/caslLCC/aff)^3);            % Ca-dep
tca=((sp9/(1+(caslLCC*aff/cat)^4))+sp10); % Ca-dep % /1000
R2=aR2/(1+exp(-(y(39)-sR2)*pR2));
T2=aT2/(1+exp(-(y(39)-sT2)*pT2));
PT=1-(1/(1+exp(-(y(39)+sp2)/sp3)));
R1=aR1/(1+exp(-(y(39)-sR1)*pR1));
T1=aT1/(1+exp(-(y(39)-sT1)*pT1));
RV=sp7+sp8*exp(y(39)/sp6);
Pr=1-(1/(1+exp(-(y(39)+sp2)/sp4)));
Pq=1+sp0/(1+exp(-(y(39)+sp2)/sp4));
    TCa=Pq*((RV-tca)*Pr+tca);
Ps=1/(1+exp(-(y(39)+sp2)/sp5));
Rv1=aRv1/(1+exp(-(y(39)-sRv1)*pRv1));
Tv1=aTv1/(1+exp(-(y(39)-sTv1)*pTv1));
Rv2=aRv2/(1+exp(-(y(39)-sRv2)*pRv2));
Tv2=aTv2/(1+exp(-(y(39)-sTv2)*pTv2));
Trev1=aTrev1/(1+exp(-(y(39)-sTrev1)*pTrev1));
    Frev1=(1-Rv1)/Rv1*R1/(1-R1);
Trev2=aTrev2/(1+exp(-(y(39)-sTrev2)*pTrev2));
    Frev2=(1-Rv2)/Rv2*R2/(1-R2)*Rv1/(1-Rv1);
% Transition Rates (20 rates)
alphaLCC=ICa_speed*R2/T2;
betaLCC=ICa_speed*(1-R2)/T2;
r1=ICa_speed*R1/T1;
r2=ICa_speed*(1-R1)/T1;
k1=flag_7_state*ICa_speed*k1o*fcp;
k2=ICa_speed*k2o;
k3=ICa_speed*PT/sp1;
%k4
k5=ICa_speed*(1-Ps)/TCa;
k6=flag_7_state*ICa_speed*fcp*Ps/TCa;
s1=flag_7_state*ICa_speed*s1o*fcp;
%s2
k1p=ICa_speed*Rv1/Tv1;
k2p=ICa_speed*(1-Rv1)/Tv1;
k3p=ICa_speed*1/(Trev2*(1+Frev2));
    k4p=Frev2*k3p;                            % REV
k5p=ICa_speed*(1-Rv2)/Tv2;
k6p=ICa_speed*Rv2/Tv2;
s1p=ICa_speed*1/(Trev1*(1+Frev1));
    s2p=Frev1*s1p;                            % REV
    k4=k3*(alphaLCC/betaLCC)*(k1/k2)*(k5/k6); % REV
    s2=s1*(k2/k1)*(r1/r2);                    % REV
    
% State transitions for 'mode 1' sarcolemmal LCCs
dPc2_LCCsl_m1 = betaLCC*Pc1_LCCsl_m1 + k5*Pi2Ca_LCCsl_m1 + k5p*Pi2Ba_LCCsl_m1 - (k6+k6p+alphaLCC)*Pc2_LCCsl_m1;                      % C2_m1sl
dPc1_LCCsl_m1 = alphaLCC*Pc2_LCCsl_m1 + k2*Pi1Ca_LCCsl_m1 + k2p*Pi1Ba_LCCsl_m1 + r2*Po_LCCsl_m1 - (r1+betaLCC+k1+k1p)*Pc1_LCCsl_m1;    % C1_m1sl
dPi1Ca_LCCsl_m1 = k1*Pc1_LCCsl_m1 + k4*Pi2Ca_LCCsl_m1 + s1*Po_LCCsl_m1 - (k2+k3+s2)*Pi1Ca_LCCsl_m1;                         % I1Ca_m1sl
dPi2Ca_LCCsl_m1 = k3*Pi1Ca_LCCsl_m1 + k6*Pc2_LCCsl_m1 - (k4+k5)*Pi2Ca_LCCsl_m1;                                               % I2Ca_m1sl
dPi1Ba_LCCsl_m1 = k1p*Pc1_LCCsl_m1 + k4p*Pi2Ba_LCCsl_m1 + s1p*Po_LCCsl_m1 - (k2p+k3p+s2p)*Pi1Ba_LCCsl_m1;                       % I1Ba_m1sl
dPi2Ba_LCCsl_m1 = k3p*Pi1Ba_LCCsl_m1 + k6p*Pc2_LCCsl_m1 - (k5p+k4p)*Pi2Ba_LCCsl_m1;                                               % I2Ba_m1sl

ibarca_slm1 = Zca^2*Pca*Frdy*FoRT*y(39)/(exp(Zca*y(39)*FoRT)-1)*(gammaCai*caslLCC*exp(Zca*y(39)*FoRT)-gammaCao*Cao);
I_Ca_sl_m1 = Fsl_CaL*(ibarca_slm1*Po_LCCsl_m1)*ICa_scale;

ibarna_slm1 = Zna^2*Pna*Frdy*FoRT*y(39)/(exp(Zna*y(39)*FoRT)-1)*(gammaNai*naslLCC*exp(Zna*y(39)*FoRT)-gammaNao*Nao);
I_Na_sl_m1 = Fsl_CaL*(ibarna_slm1*Po_LCCsl_m1)*ICa_scale;

ibark_slm1 = Zk^2*Pk*Frdy*FoRT*y(39)/(exp(Zk*y(39)*FoRT)-1)*(gammaKi*kslLCC*exp(Zk*y(39)*FoRT)-gammaKo*Ko);
I_K_sl_m1 = Fsl_CaL*(ibark_slm1*Po_LCCsl_m1)*ICa_scale;

%% I_Ca SL - mode-2 (channels expressing CaM1234, impaired CDI)
% Re-define all parameters as mode-2 specific parameters
fcp=1/(1+(cpt/caslLCC/aff_m2)^3);        % Ca-dep
tca=sp9/(1+(caslLCC*aff_m2/cat)^4)+sp10; % Ca-dep
R2=aR2/(1+exp(-(y(39)-sR2)*pR2));
T2=aT2/(1+exp(-(y(39)-sT2)*pT2));
PT=1-(1/(1+exp(-(y(39)+sp2)/sp3)));
R1=aR1/(1+exp(-(y(39)-sR1)*pR1));
T1=aT1/(1+exp(-(y(39)-sT1)*pT1));
RV=sp7+sp8*exp(y(39)/sp6);
Pr=1-(1/(1+exp(-(y(39)+sp2)/sp4)));
Pq=1+sp0/(1+exp(-(y(39)+sp2)/sp4));
    TCa=Pq*((RV-tca)*Pr+tca);
Ps=1/(1+exp(-(y(39)+sp2)/sp5));
Rv1=aRv1/(1+exp(-(y(39)-sRv1)*pRv1));
Tv1=aTv1/(1+exp(-(y(39)-sTv1)*pTv1));
Rv2=aRv2/(1+exp(-(y(39)-sRv2)*pRv2));
Tv2=aTv2/(1+exp(-(y(39)-sTv2)*pTv2));
Trev1=aTrev1/(1+exp(-(y(39)-sTrev1)*pTrev1));
    Frev1=(1-Rv1)/Rv1*R1/(1-R1);
Trev2=aTrev2/(1+exp(-(y(39)-sTrev2)*pTrev2));
    Frev2=(1-Rv2)/Rv2*R2/(1-R2)*Rv1/(1-Rv1);
% Re-calculate transition rates
alphaLCCm2=ICa_speed*R2/T2;
betaLCCm2=ICa_speed*(1-R2)/T2;
r1m2=ICa_speed*R1/T1;
r2m2=ICa_speed*(1-R1)/T1;
k1m2=flag_7_state_m2*ICa_speed*k1o*fcp;
k2m2=ICa_speed*k2o;
k3m2=ICa_speed*PT/sp1;
%k4
k5m2=ICa_speed*(1-Ps)/TCa;
k6m2=flag_7_state_m2*ICa_speed*fcp*Ps/TCa;
s1m2=flag_7_state_m2*ICa_speed*s1o*fcp;
%s2
k1pm2=ICa_speed*Rv1/Tv1;
k2pm2=ICa_speed*(1-Rv1)/Tv1;
k3pm2=ICa_speed*1/(Trev2*(1+Frev2));
    k4pm2=Frev2*k3pm2;                                        % REV
k5pm2=ICa_speed*(1-Rv2)/Tv2;
k6pm2=ICa_speed*Rv2/Tv2;
s1pm2=ICa_speed*1/(Trev1*(1+Frev1));
    s2pm2=Frev1*s1pm2;                                        % REV
    k4m2=k3m2*(alphaLCCm2/betaLCCm2)*(k1m2/k2m2)*(k5m2/k6m2); % REV
    s2m2=s1m2*(k2m2/k1m2)*(r1m2/r2m2);                        % REV
    
% State transitions for mode-2 sarcolemmal LCCs
dPc2_LCCsl_m2 = betaLCCm2*Pc1_LCCsl_m2 + k5m2*Pi2Ca_LCCsl_m2 + k5pm2*Pi2Ba_LCCsl_m2 - (k6m2+k6pm2+alphaLCCm2)*Pc2_LCCsl_m2;                      % C2_m2sl
dPc1_LCCsl_m2 = alphaLCCm2*Pc2_LCCsl_m2 + k2m2*Pi1Ca_LCCsl_m2 + k2pm2*Pi1Ba_LCCsl_m2 + r2m2*Po_LCCsl_m2 - (r1m2+betaLCCm2+k1m2+k1pm2)*Pc1_LCCsl_m2;% C1_m2sl
dPi1Ca_LCCsl_m2 = k1m2*Pc1_LCCsl_m2 + k4m2*Pi2Ca_LCCsl_m2 + s1m2*Po_LCCsl_m2 - (k2m2+k3m2+s2m2)*Pi1Ca_LCCsl_m2;                       % I1Ca_m2sl
dPi2Ca_LCCsl_m2 = k3m2*Pi1Ca_LCCsl_m2 + k6m2*Pc2_LCCsl_m2 - (k4m2+k5m2)*Pi2Ca_LCCsl_m2;                                               % I2Ca_m2sl
dPi1Ba_LCCsl_m2 = k1pm2*Pc1_LCCsl_m2 + k4pm2*Pi2Ba_LCCsl_m2 + s1pm2*Po_LCCsl_m2 - (k2pm2+k3pm2+s2pm2)*Pi1Ba_LCCsl_m2;                     % I1Ba_m2sl
dPi2Ba_LCCsl_m2 = k3pm2*Pi1Ba_LCCsl_m2 + k6pm2*Pc2_LCCsl_m2 - (k5pm2+k4pm2)*Pi2Ba_LCCsl_m2;                                               % I2Ba_m2sl

ibarca_slm2 = Zca^2*Pca*Frdy*FoRT*y(39)/(exp(Zca*y(39)*FoRT)-1)*(gammaCai*caslLCC*exp(Zca*y(39)*FoRT)-gammaCao*Cao);
I_Ca_sl_m2 = Fsl_CaL*(ibarca_slm2*Po_LCCsl_m2)*ICa_scale;

ibarna_slm2 = Zna^2*Pna*Frdy*FoRT*y(39)/(exp(Zna*y(39)*FoRT)-1)*(gammaNai*naslLCC*exp(Zna*y(39)*FoRT)-gammaNao*Nao);
I_Na_sl_m2 = Fsl_CaL*(ibarna_slm2*Po_LCCsl_m2)*ICa_scale;

ibark_slm2 = Zk^2*Pk*Frdy*FoRT*y(39)/(exp(Zk*y(39)*FoRT)-1)*(gammaKi*kslLCC*exp(Zk*y(39)*FoRT)-gammaKo*Ko);
I_K_sl_m2 = Fsl_CaL*(ibark_slm2*Po_LCCsl_m2)*ICa_scale;

%sl_mode2 = 0; % Sum up total fraction of CKII and PKA-shifted mode 2 channels
sl_mode2 = alpha_CaM1234;

% Total junctional I_Ca
I_Ca_sl2 = (1-sl_mode2)*I_Ca_sl_m1 + sl_mode2*I_Ca_sl_m2; 
I_Na_sl2 = (1-sl_mode2)*I_Na_sl_m1 + sl_mode2*I_Na_sl_m2;
I_K_sl2 = (1-sl_mode2)*I_K_sl_m1 + sl_mode2*I_K_sl_m2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total K current through LCC
I_CaK2 = I_K_junc2 + I_K_sl2;

% output LTCC module ODE
ydot(43:48)=[dPc2_LCCj_m1 dPc1_LCCj_m1 dPi1Ca_LCCj_m1 dPi2Ca_LCCj_m1 dPi1Ba_LCCj_m1 dPi2Ba_LCCj_m1];
ydot(49:54)=[dPc2_LCCj_m2 dPc1_LCCj_m2 dPi1Ca_LCCj_m2 dPi2Ca_LCCj_m2 dPi1Ba_LCCj_m2 dPi2Ba_LCCj_m2];
ydot(55:60)=[dPc2_LCCsl_m1 dPc1_LCCsl_m1 dPi1Ca_LCCsl_m1 dPi2Ca_LCCsl_m1 dPi1Ba_LCCsl_m1 dPi2Ba_LCCsl_m1];
ydot(61:66)=[dPc2_LCCsl_m2 dPc1_LCCsl_m2 dPi1Ca_LCCsl_m2 dPi2Ca_LCCsl_m2 dPi1Ba_LCCsl_m2 dPi2Ba_LCCsl_m2];

%% Compute total I_Ca (depending on flagMica)
I_Ca_junc = (1-flagMica)*I_Ca_junc1 + flagMica*I_Ca_junc2;
I_Ca_sl = (1-flagMica)*I_Ca_sl1 + flagMica*I_Ca_sl2;
I_Ca = I_Ca_junc+I_Ca_sl; % Total Ca curren throuhgh LCC

I_CaK = (1-flagMica)*(I_CaK1) + flagMica*(I_CaK2); % Total K current through LCC

I_CaNa_junc = (1-flagMica)*(I_CaNa_junc1) + (flagMica)*(I_Na_junc2);
I_CaNa_sl = (1-flagMica)*(I_CaNa_sl1) + (flagMica)*(I_Na_sl2);
I_CaNa = I_CaNa_junc + I_CaNa_sl; % Total Na current through LCC

% Collect all currents through LCC
I_Catot = I_Ca + I_CaK + I_CaNa;

%% I_ncx: Na/Ca Exchanger flux
Ka_junc = 1/(1+(Kdact/y(36))^3);
Ka_sl = 1/(1+(Kdact/y(37))^3);
s1_junc = exp(nu*y(39)*FoRT)*y(32)^3*Cao;
s1_sl = exp(nu*y(39)*FoRT)*y(33)^3*Cao;
s2_junc = exp((nu-1)*y(39)*FoRT)*Nao^3*y(36);
s3_junc = KmCai*Nao^3*(1+(y(32)/KmNai)^3) + KmNao^3*y(36)*(1+y(36)/KmCai)+KmCao*y(32)^3+y(32)^3*Cao+Nao^3*y(36);
s2_sl = exp((nu-1)*y(39)*FoRT)*Nao^3*y(37);
s3_sl = KmCai*Nao^3*(1+(y(33)/KmNai)^3) + KmNao^3*y(37)*(1+y(37)/KmCai)+KmCao*y(33)^3+y(33)^3*Cao+Nao^3*y(37);
I_ncx_junc = Fjunc*IbarNCX*Q10NCX^Qpow*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1+ksat*exp((nu-1)*y(39)*FoRT));
I_ncx_sl = Fsl*IbarNCX*Q10NCX^Qpow*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1+ksat*exp((nu-1)*y(39)*FoRT));
I_ncx = I_ncx_junc+I_ncx_sl;

%% I_pca: Sarcolemmal Ca Pump Current
I_pca_junc = Fjunc*Q10SLCaP^Qpow*IbarSLCaP*y(36)^1.6/(KmPCa^1.6+y(36)^1.6);
I_pca_sl = Fsl*Q10SLCaP^Qpow*IbarSLCaP*y(37)^1.6/(KmPCa^1.6+y(37)^1.6);
I_pca = I_pca_junc+I_pca_sl;

%% I_cabk: Ca Background Current
I_cabk_junc = Fjunc*GCaB*(y(39)-eca_junc);
I_cabk_sl = Fsl*GCaB*(y(39)-eca_sl);
I_cabk = I_cabk_junc+I_cabk_sl;

%% SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
% Vm-dependence modulation (with flagMica=1)
aRyR = 1.5; % [-]
sRyR = 30; % [mV]
pRyR = 0.065; % [1/mV]
Vdep_SRCarel = flagMica * aRyR/(1+exp((y(39)-sRyR)*pRyR)) + (1-flagMica) * 1;

MaxSR = 15; MinSR = 1;
kCaSR = MaxSR - (MaxSR-MinSR)/(1+(ec50SR/y(31))^2.5);
koSRCa = koCa/kCaSR;
kiSRCa = kiCa*kCaSR;
RI = 1-y(14)-y(15)-y(16);
ydot(14) = (kim*RI-kiSRCa*y(36)*y(14))-(koSRCa*y(36)^2*y(14)-kom*y(15));   % R
ydot(15) = (koSRCa*y(36)^2*y(14)-kom*y(15))-(kiSRCa*y(36)*y(15)-kim*y(16));% O
ydot(16) = (kiSRCa*y(36)*y(15)-kim*y(16))-(kom*y(16)-koSRCa*y(36)^2*RI);   % I
J_SRCarel = Vdep_SRCarel*ks*y(15)*(y(31)-y(36)); % [mM/ms]

J_serca = Q10SRCaP^Qpow*Vmax_SRCaP*((y(38)/Kmf)^hillSRCaP-(y(31)/Kmr)^hillSRCaP)...
    /(1+(y(38)/Kmf)^hillSRCaP+(y(31)/Kmr)^hillSRCaP);

J_SRleak = G_SRleak*(y(31)-y(36)); % [mM/ms]

%% Sodium and Calcium Buffering
ydot(17) = kon_na*y(32)*(Bmax_Naj-y(17))-koff_na*y(17);        % NaBj      [mM/ms]
ydot(18) = kon_na*y(33)*(Bmax_Nasl-y(18))-koff_na*y(18);       % NaBsl     [mM/ms]

% Cytosolic Ca Buffers
ydot(19) = kon_tncl*y(38)*(Bmax_TnClow-y(19))-koff_tncl*y(19);            % TnCL      [mM/ms]
ydot(20) = kon_tnchca*y(38)*(Bmax_TnChigh-y(20)-y(21))-koff_tnchca*y(20); % TnCHc     [mM/ms]
ydot(21) = kon_tnchmg*Mgi*(Bmax_TnChigh-y(20)-y(21))-koff_tnchmg*y(21);   % TnCHm     [mM/ms]
ydot(22) = kon_cam*y(38)*(Bmax_CaM-y(22))-koff_cam*y(22);                 % CaM       [mM/ms]
ydot(23) = kon_myoca*y(38)*(Bmax_myosin-y(23)-y(24))-koff_myoca*y(23);    % Myosin_ca [mM/ms]
ydot(24) = kon_myomg*Mgi*(Bmax_myosin-y(23)-y(24))-koff_myomg*y(24);      % Myosin_mg [mM/ms]
ydot(25) = kon_sr*y(38)*(Bmax_SR-y(25))-koff_sr*y(25);                    % SRB       [mM/ms]
%J_CaB_cytosol = sum(ydot(19:25)); % wrong formulation
J_CaB_cytosol = ydot(19)+ydot(20)+ydot(22)+ydot(23)+ydot(25);

% Junctional and SL Ca Buffers
ydot(26) = kon_sll*y(36)*(Bmax_SLlowj-y(26))-koff_sll*y(26);       % SLLj      [mM/ms]
ydot(27) = kon_sll*y(37)*(Bmax_SLlowsl-y(27))-koff_sll*y(27);      % SLLsl     [mM/ms]
ydot(28) = kon_slh*y(36)*(Bmax_SLhighj-y(28))-koff_slh*y(28);      % SLHj      [mM/ms]
ydot(29) = kon_slh*y(37)*(Bmax_SLhighsl-y(29))-koff_slh*y(29);     % SLHsl     [mM/ms]
J_CaB_junction = ydot(26)+ydot(28);
J_CaB_sl = ydot(27)+ydot(29);

%% Ion concentrations
% SR Ca Concentrations
ydot(30) = kon_csqn*y(31)*(Bmax_Csqn-y(30))-koff_csqn*y(30); % Csqn [mM/ms]
ydot(31) = J_serca-(J_SRleak*Vmyo/Vsr+J_SRCarel)-ydot(30); % Ca_sr [mM/ms] %Ratio 3 leak current
% ydot(30) = 0;
% ydot(31) = 0;

% Sodium Concentrations
I_Na_tot_junc = I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc; % [uA/uF]
I_Na_tot_sl = I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl; % [uA/uF]
%I_Na_tot_sl2 = I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl*0; % [uA/uF]

ydot(32) = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(y(33)-y(32))-ydot(17);
ydot(33) = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(y(32)-y(33))...
   +J_na_slmyo/Vsl*(y(34)-y(33))-ydot(18); %FluxNaSL=ydot(33);
ydot(34) = J_na_slmyo/Vmyo*(y(33)-y(34)); % [mM/msec] 
% ydot(32) = 0;
% ydot(33) = 0;
% ydot(34) = 0;

% Potassium Concentration
I_K_tot = I_to+I_kr+I_ks+I_ki-2*I_nak+I_CaK+I_kp; % [uA/uF]
% ydot(35) = -I_K_tot*Cmem/(Vmyo*Frdy); % [mM/msec]
ydot(35) = 0; % -I_K_tot*Cmem/(Vmyo*Frdy); % [mM/msec]

% Calcium Concentrations
I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc; % [uA/uF]
I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl; % [uA/uF]
ydot(36) = -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(y(37)-y(36))...
    -J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc; % Ca_j
ydot(37) = -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(y(36)-y(37))...
    + J_ca_slmyo/Vsl*(y(38)-y(37))-J_CaB_sl;  % Ca_sl
ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol +J_ca_slmyo/Vmyo*(y(37)-y(38)); % [mM/msec]
% ydot(36) = 0;
% ydot(37) = 0;
% ydot(38) = 0;

if flag_Ba  == 1 || flag_EGTA == 1,
    ydot(38) = 0;
end
if flag_BAPTA == 1,
    ydot(36) = 0; ydot(37) = 0; ydot(38) = 0;
end

%% Simulation type
switch lower(protocol)
    
    case {'none',''},
        I_app = 0;

    case 'pace', % pace w/ current injection at rate 'rate'
		rate = (p_HR)*1e-3;
		if mod(t+0,1/rate) <= 5
            I_app = 9.5;
        else
            I_app = 0.0;
        end
        
    case 'step', % 200-ms voltage step (from -80 to 10 mV) at rate 'rate'
		rate = (p_HR)*1e-3;
        if mod(t,1/rate) < 10
            V_clamp = -80;
        elseif mod(t,1/rate) >= 10 && mod(t,1/rate) < 210
            V_clamp = 10;
        elseif mod(t,1/rate) >= 210
            V_clamp = -80;
        end
        R_clamp = 0.02;
		I_app = (V_clamp-y(39))/R_clamp;
        
end 

%% Membrane Potential
I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;          % [uA/uF]
I_Cl_tot = I_ClCa+I_Clbk;                        % [uA/uF]
I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;            % [uA/uF]
I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
ydot(39) = -(I_tot-I_app);
vmax = ydot(39);

%% Adjust output depending on the function call
if (nargin == 3)
    output = ydot;
elseif (nargin == 4) && strcmp(runType,'ydot')
    output = ydot;
elseif (nargin == 4) && strcmp(runType,'rates')
    output = r;
elseif (nargin == 4) && strcmp(runType,'currents')
    %currents = [I_Catot I_Natot I_nabk I_nak I_kr I_ks I_kp I_tos I_tof I_ki I_ClCa I_Clbk I_ncx I_pca I_cabk];
    currents = [I_Catot I_Natot];
    output = currents;
end
% ----- END E-C COUPLING MODEL --------------

%% Calculate timecourse for currents and other intermediates
function currents = calcCurrents(t,y,p)
% After running a simulation, feed the time vector and state variables into
% this function to compute ionic currents, etc.
% currents = [I_Catot I_Natot];
currents = [];
for i=1:size(t)
    if ceil(i/1000)==i/1000
        disp(['t = ',num2str(ceil(t(i)))]);
    end
    currents = [currents;f(t(i),y(i,:),p,'currents')];
end
% ----- END calcCurrents FUNC ---------------