function Par=Adapt(Par);
%function Par=Adapt(Par);
% Theo Arts, Maastricht University, Eindhoven University of Technology,
% April 3, 2004, email: t.arts@bf.unimaas.nl
% Simulates adaptation of vessels and chambers to current hemodynamics
%Par=structue
%INPUT  i
%  SVar: State Variables, solution differential equations
%  AdaptationType: 'r'=rest, 'e'=exercise, 't'=contractility only
%  tCycle
%  Dt
%  Par.La, Ra, Lv, Rv
%  Par.TubeLArt, TubeRArt, TubeLVen, TubeRVen
%  Par.LRp, RRp
%  Par.PDropPulm, p0
%OUTPUT o
%  for chambers: VWall, LsV0, Sarc.Contractility
%  for vessels: AWall, A0

save ParTemp Par; %saves intermediate solution

SVar=Par.SVar;
type=Par.AdaptationType;

nt=size(SVar,1);
nC=ceil(Par.tCycle/Par.Dt); %number of samples per cycle
SVar=SVar([-nC+1:0]+end,:);
[SVarDot,Par]=HrtSVarDot(0,SVar',[],Par);

%===contractility, adaptation always performed
disp('sequence La, Ra, Lv, Rv')
Par.La=AdaptContractility(Par.La,'A');
Par.Ra=AdaptContractility(Par.Ra,'A');
Par.Lv=AdaptContractility(Par.Lv,'V');
Par.Rv=AdaptContractility(Par.Rv,'V');
% fine tuning of peripheral resistances
Par.LRp.R= Par.LRp.R * abs(Par.p0/mean(Par.TubeLArt.p)).^0.3; %slower adaptation of LRp.R
Par.RRp.R= Par.RRp.R * Par.pDropPulm/mean(Par.TubeRArt.p-Par.TubeLVen.p);
%Par.LRp.R= Par.p0/mean(Par.LRp.q);%++++++++++++
%Par.RRp.R= Par.pDropPulm/mean(Par.RRp.q);%+++++++++++++
disp('1e3*mean([Par.ValveLArt.q,Par.ValveRArt.q,Par.ValveLVen.q,Par.ValveRVen.q])');
disp( 1e3*mean([Par.ValveLArt.q,Par.ValveRArt.q,Par.ValveLVen.q,Par.ValveRVen.q]) );

%=== Adaptation of contractility, always
% filling in adaptation parameter values
%=== vessels
%disp('sequence adaptation: LArt,RArt,LVen,RVen');
Par.TubeLArt=AdaptVessel(Par.TubeLArt,type); 
Par.TubeRArt=AdaptVessel(Par.TubeRArt,type);
Par.TubeLVen=AdaptVessel(Par.TubeLVen,type);
Par.TubeRVen=AdaptVessel(Par.TubeRVen,type);

%==== Exercise
if lower(type)=='e';
    % ATRIA
    
    Parav=Par.La;
    disp('La');
    AuxDisp=[1e6*Parav.VWall,Parav.LsV0,1e-3*Parav.Sarc.GPas];
    Parav=AdaptVA(Parav,'A');
    disp('VWall(ml), LsV0(um), GPas(kPa/um)');
    disp([AuxDisp;[1e6*Parav.VWall,Parav.LsV0,1e-3*Parav.Sarc.GPas]]);
    Par.La=Parav;
    
    Parav=Par.Ra;
    disp('Ra');
    AuxDisp=[1e6*Parav.VWall,Parav.LsV0,1e-3*Parav.Sarc.GPas];
    Parav=AdaptVA(Parav,'A');
    disp('VWall(ml), LsV0(um), GPas(kPa/um)');
    disp([AuxDisp;[1e6*Parav.VWall,Parav.LsV0,1e-3*Parav.Sarc.GPas]]);
    Par.Ra=Parav;
    
    % VENTRICLES
    
    Parav=Par.Lv;
    disp('Lv');
    AuxDisp=[1e6*Parav.VWall,Parav.LsV0,1e-3*Parav.Sarc.GPas];
    Parav=AdaptVA(Parav,'V');
    disp(['VWall [l], LsV0 [um]']);
    disp([AuxDisp;[1e6*Parav.VWall,Parav.LsV0,1e-3*Parav.Sarc.GPas]]);
    Par.Lv=Parav;
    
    Parav=Par.Rv; 
    disp('Rv');
    AuxDisp=[1e6*Parav.VWall,Parav.LsV0,1e-3*Parav.Sarc.GPas];
    Parav=AdaptVA(Parav,'V');
    disp(['VWall [l], LsV0 [um]']);
    disp([AuxDisp;[1e6*Parav.VWall,Parav.LsV0,1e-3*Parav.Sarc.GPas]]);
    Par.Rv=Parav;
end

[SVarDot,Par]=HrtSVarDot(0,SVar(end,:)',[],Par);
Par.ErrContractility=sqrt(mean([...
        Par.La.Sarc.ContractilityIncrease,...
        Par.Ra.Sarc.ContractilityIncrease,...
        Par.Lv.Sarc.ContractilityIncrease,...
        Par.Rv.Sarc.ContractilityIncrease].^2));
%Par.AdaptationType='t'; %reset to standard adaptation
return

function ParCav=AdaptContractility(ParCav,VA);
ParSarc=ParCav.Sarc;

LsMin=min(ParSarc.Ls); % shortest sarcomere length in the cycle
LsMax=max(ParSarc.Ls); %  longest sarcomere length in the cycle

%==== contractility regulation for next beat
aux= LsMax-ParSarc.Adapt.LsBe; % increases with contractility, aux=0->contr=1
switch VA
    case 'V'; a=4.0*aux;
    case 'A'; a=30.0*aux;
end;
Contractility=1.2/(1+0.2*exp(-a));
ParCav.Sarc.ContractilityIncrease=Contractility-ParCav.Sarc.Contractility;

aux=abs(ParCav.Sarc.Contractility);
a=0.3*(0.5+Contractility); %slower feedback for low contractility
ParCav.Sarc.Contractility=aux*(Contractility/aux)^a; %a=factor of renewed fraction
disp(['Contractility: ',num2str([aux,ParCav.Sarc.Contractility])]);
return

%===== Auxiliary functions ======

function ParTube=AdaptVessel(ParTube,RestExc);
AuxDisplay=[ParTube.p0*1e-3,ParTube.q0*1e3,ParTube.pMax*1e-3,ParTube.A0*1e+4,ParTube.AWall*1e+4];

switch RestExc;
    case 'r'
        FracNew=0.5;%========================
        ParTube.q0  = ParTube.q0  * abs( (mean(ParTube.qRemod)/2) /ParTube.q0 )^FracNew;
        ParTube.p0  = ParTube.p0  * abs( max(0.2*ParTube.pMax, mean(abs(ParTube.p))) /ParTube.p0  )^FracNew;
        ParTube=TubeAdapt(ParTube)   ; % Wall and lumen size adapted
        disp('    p0(kPa), qRemod*(L/s), pMax(kPa), A0(cm2), AWall(cm2)');
        disp([AuxDisplay;[ParTube.p0*1e-3,ParTube.q0*1e3,ParTube.pMax*1e-3,...
                    ParTube.A0*1e+4,ParTube.AWall*1e+4]]);
    case 'e'
        % ParTube.p0  = ParTube.p0  *( max(0.2*ParTube.pMax, mean(abs(ParTube.p))) /ParTube.p0  )^FracNew;
        FracNew=0.3;%=========================
        aux=max(ParTube.pIn + ParTube.Z.*ParTube.A*ParTube.Adapt.vMax); aux=max(aux,0.1*abs(aux));
        ParTube.pMax= ParTube.pMax* abs(aux/ParTube.pMax)^FracNew;
        ParTube=TubeAdapt(ParTube)   ; % Wall and lumen size adapted
        disp('    p0(kPa), qRemod*(L/s), pMax(kPa), A0(cm2), AWall(cm2)');
        disp([AuxDisplay;[ParTube.p0*1e-3,ParTube.q0*1e3,ParTube.pMax*1e-3,...
                    ParTube.A0*1e+4,ParTube.AWall*1e+4]]);
    otherwise %no action
end

Aux=abs(max(ParTube.p(end),0.01)/ParTube.p0)^(1/(ParTube.k/3-1))/(1+3*ParTube.A0/ParTube.AWall);
ParTube.V(end)=((Aux-1)/3)*ParTube.A0*ParTube.Len; %new estimate of blood vessel volume

return

function ParCav=AdaptVA(ParCav,VA);
%Simulates adaptation of cavity geometry to pump load
%ParCav
%  i   .VStroke
%  i   .pMax
%  o   .VWall
%  o   .LsV0
%  o   .VBe % estimate V for begin ejection
%  (i) .Sarc.Adapt
%  i           .wStroke
%  i           .LsBe
%  i           .LsEe
%  i           .dLsPas

%=== purely active myofiber G as fu of time

ParSarc=ParCav.Sarc;
LsMax=max(ParSarc.Ls); % largest stretch, maybe in diastole
LsMin=min(ParSarc.Ls); % lowest stretch, maybe in diastole

FracNew=0.3; %renewal fraction of adaptation per beat

if VA=='V'; %=== VENTRICLE
    LStrain=(LsMax/LsMin)/(ParSarc.Adapt.LsBe/ParSarc.Adapt.LsEe);
end
if VA=='A'; %=== ATRIUM
    LStrain=(LsMax/LsMin)/((ParSarc.Adapt.LsBe)/ParSarc.Adapt.LsEe); %(max-->Ee)  
end;

disp(['LsMax,LsMin: ',num2str([LsMax,LsMin])]);
VWallFac= ParSarc.Contractility ^ FracNew;
ParCav.VWall=ParCav.VWall*VWallFac;
dFac=abs(LStrain)^FracNew;
ParCav.LsV0=ParCav.LsV0/dFac; % Cavity shrink for subnormal strain, strain softening

%==== estimate new cavity volume with given sarcomere length
ParCav.V(end)=ParCav.VWall*(abs(ParSarc.Ls(end)/ParCav.LsV0)^3-1)/3;

ParCav.Sarc=ParSarc;

return

%VWall= 178 ml for standard LV
