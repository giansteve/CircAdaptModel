function [SVarDot,Par]=SVarHrtDot(tDummy,SVar,flag,Par);
%function [SVarDot,Par]=SVarHrtDot(tDummy,SVar,flag,Par);
% Theo Arts, Maastricht University, Eindhoven University of Technology,
% April 3, 2004, email: t.arts@bf.unimaas.nl
% Differential equation in the format to be eaten by ode23-function
% DynACirc, model of the dynamic auto-adaptating circulation
%INPUT
% tDummy  = a dummy variable, irrelevant (was time)
% SVar= column vector of State variables,
% !!! if SVar is a matrix, different colums of SVar define
%     different states. This is very useful for display purposes
% flag    = irrelevant
% Par= Parameter structure
%OUTPUT
% SVarDot= column vector (or matrix) of derivatives of State Variables
% Par        = structure of variables

ColOnes=ones(size(SVar,2),1); % auxilary column vector
S=SVar.*(Par.Scale*ColOnes'); % unscaling to SI-units

%=== Convert State variables to understandable physiological quantities
i=1;

Par.t           =S(i,:)';i=i+1; %time
Par.TubeLArt.V  =S(i,:)';i=i+1; %systemic arterial volume
Par.TubeLVen.V  =S(i,:)';i=i+1; %pulmonary venous volume
Par.TubeRArt.V  =S(i,:)';i=i+1; %pulmonary arterial volume
Par.TubeRVen.V  =S(i,:)';i=i+1; %systemic venous volume
Par.ValveLArt.q =S(i,:)';i=i+1; %aortic valve flow
Par.ValveLAv.q  =S(i,:)';i=i+1; %mitral valve flow
Par.ValveLVen.q =S(i,:)';i=i+1; %pulmonary venous outflow
Par.ValveRArt.q =S(i,:)';i=i+1; %pulmonary valve flow
Par.ValveRAv.q  =S(i,:)';i=i+1; %tricuspid valve flow
Par.ValveRVen.q =S(i,:)';i=i+1; %systemic venous outflow
Par.ValveDUCT.q =S(i,:)';i=i+1; %ductus arteriosis flow
Par.ValveVSD.q  =S(i,:)';i=i+1; %VSD flow
Par.ValveASD.q  =S(i,:)';i=i+1; %ASD flow
Par.La.V        =S(i,:)';i=i+1; %La cavity volume
Par.Lv.V        =S(i,:)';i=i+1; %Lv cavity volume
Par.Ra.V        =S(i,:)';i=i+1; %Ra cavity volume
Par.Rv.V        =S(i,:)';i=i+1; %Rv cavity volume
Par.La.Sarc.Lsi =S(i,:)';i=i+1; %contractile unit length of sarcomere La
Par.La.Sarc.C   =S(i,:)';i=i+1; %contractility of sarcomere La
Par.Lv.Sarc.Lsi =S(i,:)';i=i+1; %contractile unit length of sarcomere Lv
Par.Lv.Sarc.C   =S(i,:)';i=i+1; %contractility of sarcomere La
Par.Ra.Sarc.Lsi =S(i,:)';i=i+1; %contractile unit length of sarcomere Ra
Par.Ra.Sarc.C   =S(i,:)';i=i+1; %contractility of sarcomere La
Par.Rv.Sarc.Lsi =S(i,:)';i=i+1; %contractile unit length of sarcomere Rv
Par.Rv.Sarc.C   =S(i,:)';i=i+1; %contractility of sarcomere La
%disp(num2str(Par.t));%FETUS+++++++++++++++++

%=== Electrical activation time
Par.La.Sarc.tc= mod(Par.t-Par.La.Sarc.ActivationDelay,Par.tCycle); % left atrium
Par.Lv.Sarc.tc= mod(Par.t-Par.Lv.Sarc.ActivationDelay,Par.tCycle); % left ventricle
Par.Ra.Sarc.tc= mod(Par.t-Par.Ra.Sarc.ActivationDelay,Par.tCycle); % left atrium
Par.Rv.Sarc.tc= mod(Par.t-Par.Rv.Sarc.ActivationDelay,Par.tCycle); % left ventricle

%=== Cavity pressures and cross-sections
Par.La        =CavityMech(Par.La)         ; % sarcomere length and pressures
Par.Lv        =CavityMech(Par.Lv)         ; % inner shell
Par.Ra        =CavityMech(Par.Ra)         ; % 
Par.Rv        =CavityMech(Par.Rv)         ; % outer shell
%=== Left to right cross-talk factor, trial, not optimal yet

%LR interaction Ventricles
y= Par.Lv.VWall/Par.Rv.VWall  ;
LAdapt=Par.Sarc.Adapt.LsBe;
%x=((LAdapt/Par.Rv.LsV0)^3-2/3)/((LAdapt/Par.Lv.LsV0)^3-2/3)*y; %normally>0, ~3.5
x=(Par.Lv.LsV0/Par.Rv.LsV0)^3*y; %normally>0, ~3.5
LrxV=0.5-1/(x^2+1); %now pressure dependent cross talk
%   Ratio, determining right ventricle exterior to left
%   based on expected pressure ratio
%LrxV=0.5;%Conventional
%LrxV=0.0;%Symmetric
MlrV=[1,(0.5-LrxV)^2;(0.5+LrxV)^2,1];

%Alternative LR-interaction
%xLmR=log(Par.Lv.VWall/Par.Rv.VWall)+log(Par.Lv.LsV0/Par.Rv.LsV0); %normally ~0.6
%x=tanh(3.0*xLmR); %clamping abs(x)<1
%a1=0*x.*(1-x.^2); a2=0.25*(1+x).^2; a3=0.25*(1-x).^2;
%MlrV=[1-a1,a3;a2,1+a1];

%LR interaction Atria (optional)
%y= Par.La.VWall/Par.Ra.VWall  ;
%LAdapt=Par.Sarc.Adapt.LsBe;
%x=((LAdapt/Par.Ra.LsV0)^3-2/3)/((LAdapt/Par.La.LsV0)^3-2/3)*y;
%LrxA=0.5-1/(x^2+1); %now pressure dependent cross talk
%MlrA=[1,(0.5-LrxA)^2;(0.5+LrxA)^2,1];
MlrA=[1,0;0,1]; %No cross-talk between atria

%=== Tube pressures, wave impedances and cross-sections
Par.TubeLArt  =Tube(Par.TubeLArt)        ; % pressures and wave impedance of Aorta
Par.TubeLVen  =Tube(Par.TubeLVen)        ; % pressures and wave impedance of Pulmonary Venes
Par.TubeRArt  =Tube(Par.TubeRArt)        ; % pressures and wave impedance of Pulmonary Artery
Par.TubeRVen  =Tube(Par.TubeRVen)        ; % pressures and wave impedance of Systemic Venes

%=== Systemic peripheral resistance with blood pressure control
Par.LRp.q= Par.TubeLArt.p/Par.LRp.R;
Par.TubeLArt.VDot=Par.TubeLArt.VDot-Par.LRp.q; %outflow should be equal to inflow, if stationary
Par.TubeRVen.VDot=Par.TubeRVen.VDot+Par.q0; %inflow of constant systemic flow
Par.TubeLArt.qRemod=Par.TubeLArt.qRemod+abs(Par.LRp.q);
Par.TubeRVen.qRemod=Par.TubeRVen.qRemod+abs(Par.q0);

%--- nearly constant pressure over pulmonary circulation
Par.RRp.q=(Par.TubeRArt.p-Par.TubeLVen.p)/Par.RRp.R;
Par.TubeRArt.VDot=Par.TubeRArt.VDot-Par.RRp.q;
Par.TubeLVen.VDot=Par.TubeLVen.VDot+Par.RRp.q;
Par.TubeRArt.qRemod=Par.TubeRArt.qRemod+abs(Par.RRp.q);
Par.TubeLVen.qRemod=Par.TubeLVen.qRemod+abs(Par.RRp.q);

%=== FLOWS AROUND VALVES
% i=1: Flow continuity
% i=2: Time derivatives of volumes and boundary pressures pIn and pOut
% i=3: Time derivatives of flows
for i=1:3,
    [Par.TubeLVen,Par.ValveLVen ,Par.La      ]=ValveFlow(Par.TubeLVen,Par.ValveLVen ,Par.La      );
    [Par.La      ,Par.ValveLAv  ,Par.Lv      ]=ValveFlow(Par.La      ,Par.ValveLAv  ,Par.Lv      );
    [Par.Lv      ,Par.ValveLArt ,Par.TubeLArt]=ValveFlow(Par.Lv      ,Par.ValveLArt ,Par.TubeLArt);
    [Par.TubeRVen,Par.ValveRVen ,Par.Ra      ]=ValveFlow(Par.TubeRVen,Par.ValveRVen ,Par.Ra      );
    [Par.Ra      ,Par.ValveRAv  ,Par.Rv      ]=ValveFlow(Par.Ra      ,Par.ValveRAv  ,Par.Rv      );
    [Par.Rv      ,Par.ValveRArt ,Par.TubeRArt]=ValveFlow(Par.Rv      ,Par.ValveRArt ,Par.TubeRArt);
    [Par.TubeLArt,Par.ValveDUCT,Par.TubeRArt ]=ValveFlow(Par.TubeLArt,Par.ValveDUCT ,Par.TubeRArt);
    [Par.Lv      ,Par.ValveVSD ,Par.Rv       ]=ValveFlow(Par.Lv      ,Par.ValveVSD  ,Par.Rv      );
    [Par.La      ,Par.ValveASD ,Par.Ra       ]=ValveFlow(Par.La      ,Par.ValveASD  ,Par.Ra      );
    
    if i==1,
        Aux=[Par.Lv.p,Par.Rv.p]*MlrV;
        Par.Lv.p=Aux(:,1); Par.Rv.p=Aux(:,2); % from myo-shells to hemodynamic pressures
        Aux=[Par.La.p,Par.Ra.p]*MlrA;
        Par.La.p=Aux(:,1); Par.Ra.p=Aux(:,2); % from my0-shells to hemodynamic pressures
    end
    if i==2,
        %--- from hemodynamic flows to shell volume derivatives, LR
        %ventricular interaction
        Aux=[Par.Lv.VDot,Par.Rv.VDot]*MlrV';
        Par.Lv.VDot=Aux(:,1); Par.Rv.VDot=Aux(:,2);
        Aux=[Par.Lv.V,Par.Rv.V]*inv(MlrV');
        Par.Lv.VHem=Aux(:,1); Par.Rv.VHem=Aux(:,2);
        
        %--- from hemodynamic flows to shell volume derivatives, LR
        %atrial interaction
        Aux=[Par.La.VDot,Par.Ra.VDot]*MlrA';
        Par.La.VDot=Aux(:,1); Par.Ra.VDot=Aux(:,2);
        Aux=[Par.La.V,Par.Ra.V]*inv(MlrA');
        Par.La.VHem=Aux(:,1); Par.Ra.VHem=Aux(:,2);
     end
end


%=== Collecting derivatives of State Variables
SDot=zeros(size(S));
i=1;

SDot(i,:)=ColOnes'           ;i=i+1; %time
SDot(i,:)=Par.TubeLArt.VDot' ;i=i+1; %systemic arterial volume
SDot(i,:)=Par.TubeLVen.VDot' ;i=i+1; %pulmonary venous volume
SDot(i,:)=Par.TubeRArt.VDot' ;i=i+1; %pulmonary arterial volume
SDot(i,:)=Par.TubeRVen.VDot' ;i=i+1; %systemic venous volume
SDot(i,:)=Par.ValveLArt.qDot';i=i+1; %aortic valve flow
SDot(i,:)=Par.ValveLAv.qDot' ;i=i+1; %mitral valve flow
SDot(i,:)=Par.ValveLVen.qDot';i=i+1; %pulmonary venous outflow
SDot(i,:)=Par.ValveRArt.qDot';i=i+1; %pulmonary valve flow
SDot(i,:)=Par.ValveRAv.qDot' ;i=i+1; %tricuspid valve flow
SDot(i,:)=Par.ValveRVen.qDot';i=i+1; %systemic venous outflow
SDot(i,:)=Par.ValveDUCT.qDot';i=i+1; %ductus arteriosis
SDot(i,:)=Par.ValveVSD.qDot' ;i=i+1; %VSD
SDot(i,:)=Par.ValveASD.qDot' ;i=i+1; %ASD
SDot(i,:)=Par.La.VDot'       ;i=i+1; %La cavity volume
SDot(i,:)=Par.Lv.VDot'       ;i=i+1; %Lv cavity volume
SDot(i,:)=Par.Ra.VDot'       ;i=i+1; %Ra cavity volume
SDot(i,:)=Par.Rv.VDot'       ;i=i+1; %Rv cavity volume
SDot(i,:)=Par.La.Sarc.LsiDot';i=i+1; %contractile unit length of sarcomere La
SDot(i,:)=Par.La.Sarc.CDot'  ;i=i+1; %contractile unit length of sarcomere La
SDot(i,:)=Par.Lv.Sarc.LsiDot';i=i+1; %contractile unit length of sarcomere Lv
SDot(i,:)=Par.Lv.Sarc.CDot'  ;i=i+1; %contractile unit length of sarcomere La
SDot(i,:)=Par.Ra.Sarc.LsiDot';i=i+1; %contractile unit length of sarcomere Ra
SDot(i,:)=Par.Ra.Sarc.CDot'  ;i=i+1; %contractile unit length of sarcomere La
SDot(i,:)=Par.Rv.Sarc.LsiDot';i=i+1; %contractile unit length of sarcomere Rv
SDot(i,:)=Par.Rv.Sarc.CDot'  ;i=i+1; %contractile unit length of sarcomere La

SVarDot=SDot./(Par.Scale*ColOnes'); %transpose to column vector

return




