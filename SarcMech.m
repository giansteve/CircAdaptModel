function ParSarc=SarcMech(ParSarc);
%function ParSarc=SarcMech(ParSarc);
%Theo Arts, University of Maastricht / Eindhoven University of Technology, Mar 2004.
%Sarcomere mechanics, G=Stress/sarcomere length= function of time and sarcomere length
%contractility and other parameters may be varied by varying ParSarc-values
%i=INPUT, o=OUTPUT
%i Lsi: State variable unloaded active sarcomere length
%i C  : State variable contractility, follows [Ca](tc)
%i tc : time after activation
%i Ls : actual sarcomere length
%i Contractility: limitation of contractility to a maximum
%i TimeAct : activation duration
%i GAct    : active stress scaling
%i vMax, LsStress0Pass, GPas, LenSeriesElement
%o LsiDot  : d/dt State variable Lsi
%o CDot    : d/dt State variable C
%o G       : stress/sarcomere length
%o DGDLs   : stiffness

%==== Input variables
t            =ParSarc.tc;
Ls           =ParSarc.Ls;
Lsi          =ParSarc.Lsi;
C            =ParSarc.C;
TR           =ParSarc.TR;
Contractility=ParSarc.Contractility;
Tau          =ParSarc.TimeAct; %time scale of contraction pulse
Ls0          =ParSarc.Lsi0Act; %zero active stress sarcomere length

% series elasticity and sarcomere shortening velocity
LNormSe   = (Ls-Lsi)/ParSarc.LenSeriesElement; % normalized Series Elastic Element length
ParSarc.LsiDot= ParSarc.vMax*(LNormSe-1)     ; % time derivative of sarcomere length

% constants related to timing are mainly based on experimental findings
L  = max(0,Lsi-Ls0)     ; % safety to avoid negative active stress
tA = Tau.*(0.8+0.30*L)  ; % activation time lengthens with sarcomere length
tR = TR*Tau             ; % rise time
tD = 0.07*Tau           ; % decay time
ft1= (1-min((1-t/tR).^2,1)).^2/tR; % rise regulating hump for d stress/ dt
ParSarc.CDot= Contractility.*L.*ft1 - C./((0.3+exp((tA-t)/tD))*tD); % 1st order control of [Ca++]
Giso= ParSarc.GAct * C.*(L.^0.5)*1.374; % Regulates wall thickness

%=== Passive stress ans stiffness characteristics
y=max(0,Ls-ParSarc.LsStress0Pas)/ParSarc.dLsPas;
a2=2.5; a2y2=a2.*y.^2; Ea2y2=exp(a2y2); % stiffness to bucking is low
GPas    =ParSarc.GPas*Ea2y2.* y       ; %close to y/(1-y) curve or zero-pole stress model
DGPasDLs=ParSarc.GPas*Ea2y2.*(1+2*a2y2); % stiffness dstress/dLs

%=== Stress/Ls and stiffness, collected for output
ParSarc.G     = GPas + Giso.*LNormSe; % adding active and passive stress
ParSarc.DGDLs = DGPasDLs+Giso/ParSarc.LenSeriesElement; % estimate of sarcomere stiffness
return
