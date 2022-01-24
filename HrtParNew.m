function Par=HrtParNew
% Theo Arts, Maastricht University, Eindhoven University of Technology,
% April 3, 2004, email: t.arts@bf.unimaas.nl
% Init: Initialization of simulation of hemodynamics and cardiac mechanics
% Left atrium, mitral valve, left ventricle, aortic valve, arteries, venes,
% Right atrium, tricuspid valve, right ventricle, pulmonary arteries, lung
% venes
% venous outflow is a circular orifice
% Par is a structure, containing parameters and variables
% Par contains sufficient information to start the problem
% The structure of Par can be viewed by executing <MapStruct(Par)>

%% === Initialization of Par
Par.tEnd  = 1.0  ; % minimum duration of simulation (s), will be rounded to a complete cycle

%==== used specific functions

%% ==== GENERAL AT REST
Par.q0        =1.00e-4     ; % Cardiac Output, mean systemic flow
Par.p0        =12200       ; % mean systemic blood pressure
Par.pDropPulm =1500        ; % pressure drop in pulmonary circulation
Par.rhob      =1050        ; % blood density


%==== Time scaling with body size
tCycleRef  = 0.85*(Par.q0/1e-4)^(0.28); % Scaling time for body size/systemic flow
Par.Dt     = 0.0005*2^round(log(tCycleRef/0.2)/log(2)); % time discretization ODD-solution (s)
Par.tCycle = round(tCycleRef/Par.Dt)*Par.Dt; % Cycle time= 1/heart rate (s) at rest

%=== Scaling of state variables
ExpV=round(log10(Par.tCycle*Par.q0)); V=10^ExpV; % scaling factors of volume,
Ls=1;                                            % scaling of sarcomere length
Expq=round(log10(Par.q0)); q=10^Expq;            % scaling of flow
Par.Scal.Flow=q; Par.Scal.Volume=V;
Par.Scale=[1,V,V,V,V,q,q,q,q,q,q,q,q,q,10*V,10*V,10*V,10*V,Ls,1,Ls,1,Ls,1,Ls,1]';
%                                           column vector for scaling of SVariables

%% ==== DEFAULTS

%==== Default Valve properties
Par.Valve.rhob=Par.rhob; % blood density
Par.Valve.FunctionNr=1; %initialization
Par.Valve.q=0; %initialial condition

%==== Default Tube properties
Par.Tube.k               = 6.0     ; % stiffness exponent vessel wall fibers, determines ZWave+++++
Par.Tube.Adapt.WallStress= 700e3   ; % Wall stress at peak (systolic) pressure (Pa)
Par.Tube.Adapt.WallShear = 60      ; % Mean wall shear rate for diameter remodeling, aorta at rest (s-1)
Par.Tube.Adapt.vMax      = 5.0     ; % Velocity of possible whole body impact, blood shockwave
Par.Tube.q0              = Par.q0  ; % estimated mean tube flow for adaptation
Par.Tube.rhob            = Par.rhob; % blood density (kgm-3)

%=== default systemic arterial Tube
Par.Tube.p0  = Par.p0;
Par.Tube.Len=1             ; %dummy value
Par.Tube.Len  = 18.0       ; % scaling factor effective tube: length/sqrt(A0)
Par.Tube=TubeInit(Par.Tube); % Wall and lumen size adapted

%==== Default (LV-) Sarcomere mechanics
Par.Sarc.ActivationDelay=0; % time of depolarization set to zero (=p-wave)

%==== Adaptation parameters for sarcomere behavior
Par.Sarc.Adapt.LsBe         = 2.30  ; % end-diastolic sarcomere length with excercise
Par.Sarc.Adapt.LsEe         = 1.75  ; % end-systolic sarcomere length with excercise
Par.Sarc.Adapt.GPasMax      = 2800  ; % maximum passive G=stress/Ls

Par.Sarc.LsStress0Pas    = 1.95 ; % [um] sarcomere length with zero passive stress
Par.Sarc.dLsPas          = 0.60 ; % [um] LsiStress0Pas+dLsPas equals stress pole
Par.Sarc.Lsi0Act         = 1.51 ; % Zero force sarcomere length
Par.Sarc.LenSeriesElement= 0.04 ; % [um] isometric dLs of series elastic element
Par.Sarc.GAct            = 60E3 ; % [Pa/um] max isometric stress*Par.Sarc.Adapt.LsMax
Par.Sarc.vMax            = 8.5 * 0.8/tCycleRef; % [um/s] shortening velocity scaled to HR(Body weight)
Par.Sarc.TimeAct         = 1e-3 ; % [s]dummy, force pulse duration at Ls=LsStress0Act (s), (Arts, AJP, 1982)
Par.Sarc.TR              = 0.25 ; % ratio rise time to TimeAct
Par.Sarc.C               = 0.0  ; % [Ca++] concentration factor
Par.Sarc.Contractility   = 1    ; % Contractility relative to normal

%==== Initialization of Sarcomere default factors for passive and active stress
%==== Initialization of passive GPas [Pa/um] level Par.Sarc.GPas
Par.Sarc.GPas= 1;
Par.Sarc.Lsi = Par.Sarc.Adapt.LsBe; Par.Sarc.Ls=Par.Sarc.Lsi; Par.Sarc.tc=0;
Par.Sarc     = SarcMech(Par.Sarc); % result ParSarc.G= maximal passive stress/Ls
Par.Sarc.GPas= Par.Sarc.GPas*Par.Sarc.Adapt.GPasMax/Par.Sarc.G; % setting GPas

%==== Default Cavity properties
Par.Cavity.Sarc = Par.Sarc;
Par.Cavity.rhob = Par.rhob; % blood density (kgm-3)


%% ======================= SPECIFIC PROPERTIES =============================

%==== Mean Valve Flows for initial adaptation
Par.ValveLVen= Par.Valve; % default
Par.ValveLAv = Par.Valve; % default
Par.ValveLArt= Par.Valve; % default
Par.ValveRVen= Par.Valve; % default
Par.ValveRAv = Par.Valve; % default
Par.ValveRArt= Par.Valve; % default
Par.ValveASD = Par.Valve; % default
Par.ValveVSD = Par.Valve; % default
Par.ValveDUCT= Par.Valve; % default


Par.ValveRVen.q = Par.q0;
Par.ValveLVen.q = Par.q0;
Par.ValveASD.q  = 0;
Par.ValveVSD.q  = 0;

Par.ValveLAv.q=Par.ValveLVen.q-Par.ValveASD.q;
Par.ValveRAv.q=Par.ValveRVen.q+Par.ValveASD.q;
Par.ValveLArt.q=Par.ValveLAv.q-Par.ValveVSD.q;
Par.ValveRArt.q=Par.ValveRAv.q+Par.ValveVSD.q;
Par.ValveDUCT.q=Par.ValveLVen.q-Par.ValveRArt.q;

%% ==== TUBES
%==== Left Arteries= Aorta
Par.TubeLArt      = Par.Tube           ; % Aorta as a tube
Par.TubeLArt.p0   = Par.p0             ; % estimated mean pressure
Par.TubeLArt.q0   = Par.ValveLArt.q    ; % flow to adapt to
Par.TubeLArt.Len  = 18.0               ; % scaling factor effective tube: length/sqrt(A0)
Par.TubeLArt=TubeInit(Par.TubeLArt)    ; % Wall and lumen size adapted

%==== Left Venes (pulmonary)
Par.TubeLVen      = Par.Tube           ; % Venes as a tube with open end to heart
Par.TubeLVen.p0   = 0.63*Par.Sarc.GPas ; % estimated mean venous pressure
Par.TubeLVen.q0   = Par.ValveLVen.q    ; % flow to adapt to
Par.TubeLVen.Len  = 10                 ;
Par.TubeLVen=TubeInit(Par.TubeLVen)    ; % Wall and lumen size adapted

%==== Right Venes (systemic)
Par.TubeRVen      = Par.TubeLVen       ; % Venes as a tube with open end to heart
Par.TubeRVen.p0   = 0.35*Par.Sarc.GPas ; % estimated mean venous pressure
Par.TubeRVen.q0   = Par.ValveRVen.q    ; % flow to adapt to
Par.TubeRVen.Len  = 18.0               ;
Par.TubeRVen = TubeInit(Par.TubeRVen)  ; % Wall and lumen size adapted

%==== Right Arteries (Pulmonary)
Par.TubeRArt      = Par.TubeLArt       ; % Pulmonary trunc as a tube
Par.TubeRArt.p0   = Par.pDropPulm+Par.TubeRVen.p0 ; % estimated mean pressure
Par.TubeRArt.q0   = Par.ValveRArt.q    ; % flow to adapt to
Par.TubeRArt.Len  = 10.0               ;
Par.TubeRArt = TubeInit(Par.TubeRArt)  ; % Wall and lumen size adapted

%==== Peripheral resistance
Par.LRp.R=(Par.TubeLArt.p0-Par.TubeLVen.p0)/Par.TubeLArt.q0; %set to normal
Par.RRp.R=Par.pDropPulm/Par.TubeRArt.q0; %pulmonary resistance

%% ==== VALVES
%==== Pulmonary venous outlet orifice
Par.ValveLVen.AOpen   = Par.TubeLVen.A0; % venous orifice area
Par.ValveLVen.ALeak   = Par.ValveLVen.AOpen; % orifice, slight valve function
Par.ValveLVen.Len     = 0.5*sqrt(Par.ValveLVen.AOpen);

%==== Left Atrio-ventricular Valve= Mitral valve
Par.ValveLAv.AOpen   = 2.5*Par.TubeLArt.A0; % mitral valve diameter
Par.ValveLAv.ALeak   = 1e-6*Par.ValveLAv.AOpen;
Par.ValveLAv.Len     = 0.5*sqrt(Par.ValveLAv.AOpen);

%==== Left Ventricular-arterial Valve= Aortic valve
Par.ValveLArt.AOpen   = Par.TubeLArt.A0; % aortic valve cross-section .. aortic cross-section
Par.ValveLArt.ALeak   = 1e-6*Par.ValveLArt.AOpen;
Par.ValveLArt.Len     = 0.5*sqrt(Par.ValveLArt.AOpen); % estimate of moving inertial mass

%==== Systemic venous outlet orifice
Par.ValveRVen.AOpen   = Par.TubeRVen.A0; % venous orifice diameter= venous diameter
Par.ValveRVen.ALeak   = Par.ValveRVen.AOpen; % orifice, slight valve function
Par.ValveRVen.Len     = 0.5*sqrt(Par.ValveRVen.AOpen);

%==== Right Atrio-ventricular Valve= Tricuspid valve
Par.ValveRAv.AOpen   = 2.5*Par.TubeRArt.A0; % tricuspid area
Par.ValveRAv.ALeak   = 1e-6*Par.ValveRAv.AOpen;
Par.ValveRAv.Len     = 0.5*sqrt(Par.ValveRAv.AOpen);

%==== Right Ventricular-arterial Valve= Pulmonary valve
Par.ValveRArt.AOpen   = Par.TubeRArt.A0; % pulmonary valve cross-section .. pulm. cross-section
Par.ValveRArt.ALeak   = 1e-6*Par.ValveRArt.AOpen;
Par.ValveRArt.Len     = 0.5*sqrt(Par.ValveRArt.AOpen); % estimate of moving inertial mass

%==== Ductus Arteriosis
Par.ValveDUCT.AOpen   = 1e-6*Par.Tube.A0; % Ductus cross-section
Par.ValveDUCT.ALeak   = Par.ValveDUCT.AOpen;
Par.ValveDUCT.Len     = 0.5*sqrt(Par.Tube.A0); % estimate of moving inertial mass

%==== VSD
Par.ValveVSD.AOpen   = 1e-6*Par.Tube.A0; % VSD cross-section
Par.ValveVSD.ALeak   = Par.ValveVSD.AOpen;
Par.ValveVSD.Len     = 0.5*sqrt(Par.Tube.A0); % estimate of moving inertial mass

%==== ASD
Par.ValveASD.AOpen   = 1e-6*Par.Tube.A0; % ASD cross-section
Par.ValveASD.ALeak   = Par.ValveASD.AOpen;
Par.ValveASD.Len     = 0.5*sqrt(Par.Tube.A0); % estimate of moving inertial mass


%% ==== CAVITIES
%==== Left Atrial cavity, estimate for excercise
Par.La        = Par.Cavity; % default
Par.La.VStroke= 0.50*Par.ValveLAv.q*Par.tCycle;
Par.La.pMax   = 0.28*Par.TubeLVen.pMax   ; % maximum La pressure for adaptation
Par.La        = CavityGeometry(Par.La)  ; % Adaptation of wall and cavity size
%==== Atrium specific sarcomere behavior
Par.La.Sarc.vMax         = 1.5*Par.Sarc.vMax ; % faster contraction possible
Par.La.Sarc.TR           = 0.5               ; % longer rise time relative to pulse duration
Par.La.Sarc.GAct=2.0*Par.Sarc.GAct;
Par.La.Sarc.GPas=6.0*Par.Sarc.GPas;
Par.La.Sarc.Adapt.LsEe=Par.Sarc.Adapt.LsEe-0.1;

%==== Right atrium
Par.Ra        = Par.La; % default
Par.Ra.VStroke= 0.40*Par.q0*Par.tCycle;
Par.Ra.pMax   = 0.18*Par.TubeRVen.pMax; % maximum Ra pressure, systolic pressure is quite high
Par.Ra        = CavityGeometry(Par.Ra); % Adaptation of wall and cavity size

%==== Left Ventricular cavity, inner shell
Par.Lv        = Par.Cavity;
Par.Lv.VStroke= Par.q0*Par.tCycle                    ; % stroke volume for adaptation
Par.Lv.pMax   = 1.1*(Par.TubeLArt.p0-Par.TubeRArt.p0); % maximum LV pressure for adaptation
Par.Lv        = CavityGeometry(Par.Lv)               ; % Adaptation of wall and cavity size

%==== Right Ventricular cavity, outer shell
Par.Rv        = Par.Lv;
Par.Rv.VStroke= 2.0*Par.Lv.VStroke    ;
Par.Rv.pMax   = 2.1*Par.TubeRArt.p0   ; % maximum RV pressure for adaptation
Par.Rv        = CavityGeometry(Par.Rv); % Adaptation of wall and cavity size

%% =========== SETTING INITIAL STATE
Par.t=0;

%==== TUBES
%Par.Tube.V volume are already initialized to Par.Tube.A0*Par.Tube.Len

%==== VALVES
%Par.Valve.q are already initialized to 0

%% ==== SARCOMERES
%==== Estimate of initial sarcomere length in wall of cavity
Par.La         = CavityMech(Par.La)   ; % Estimate of initial sarcomere length and cavity volume
Par.La.Sarc.Lsi= Par.La.Sarc.Ls-Par.La.Sarc.LenSeriesElement;

Par.Ra         = CavityMech(Par.Ra)   ; % Estimate of initial sarcomere length and cavity volume
Par.Ra.Sarc.Lsi= Par.Ra.Sarc.Ls-Par.Ra.Sarc.LenSeriesElement;

Par.Lv         = CavityMech(Par.Lv) ; % Estimate of sarcomere length
Par.Lv.Sarc.Lsi= Par.Lv.Sarc.Ls-Par.Lv.Sarc.LenSeriesElement;

Par.Rv         = CavityMech(Par.Rv) ; % Estimate of sarcomere length
Par.Rv.Sarc.Lsi= Par.Rv.Sarc.Ls-Par.Rv.Sarc.LenSeriesElement;

%=== Aditionals, maybe overwritten
Par.AdaptationType='t';

return


% ============== Auxiliary functions ================
function ParTube=TubeInit(ParTube)
ParTube.pMax = ParTube.p0+0.5*ParTube.rhob*ParTube.Adapt.vMax^2;
ParTube=TubeAdapt(ParTube);
ParTube.Len  = ParTube.Len*sqrt(ParTube.A0) ; % estimated effective tube windkessel length
ParTube.pMax= ParTube.p0 + abs(ParTube.Z.*ParTube.A0*ParTube.Adapt.vMax);
ParTube=TubeAdapt(ParTube)        ; % Wall and lumen size adapted
ParTube.V = ParTube.A0*ParTube.Len; %initial volume
return

function Par=CavityGeometry(Par)
%function Par=CavityGeometry(Par);
%Simulates result of initial adaptation of cavity geometry to pump load
%Par
%  i   .VStroke
%  i   .pMax
%  o   .VWall
%  o   .LsV0
%  o   .V   % estimate V for begin ejection
%      .Sarc.Adapt
%  i           .LsBe
%  i           .LsEe

wStroke=Par.Sarc.GAct*(Par.Sarc.Adapt.LsBe-Par.Sarc.Adapt.LsEe)*0.232; %estimate
%    estimate of stroke work
VWallDVStroke= Par.pMax/wStroke;
a=1.5; % excercise factor for stroke volume, rest->adaptation
VWall= a*Par.VStroke*VWallDVStroke; %multiplied by excercise stroke volume factor
Aux= 1-(Par.Sarc.Adapt.LsEe/Par.Sarc.Adapt.LsBe)^3;
%                       Aux increases with ejection fraction
Par.VWall= VWall;
Par.LsV0 = Par.Sarc.Adapt.LsBe * (Aux*VWallDVStroke/3)^(1/3);
Par.V    = Par.VStroke/Aux - VWall/3; % Cavity volume at Begin Ejection

return
