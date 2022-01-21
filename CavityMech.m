function ParCavity=CavityMech(ParCavity);
%function ParCavity=CavityMech(ParCavity);
% Theo Arts, Maastricht University, Eindhoven University of Technology,
% April 3, 2004, email: t.arts@bf.unimaas.nl
% Calculates cavity pressure from volume, given sarcomere conditions
% Needs SarcMech
%INPUT/OUTPUT i/o
%ParCavity = , fields marked with 'i' of input
% i      V: 1.3000e-004
% i   Sarc: [1x1 struct]
%     i     Lsi: 2.2172
%     i      tc: 0
%     o      Ls: 2.2572
%     o  LsiDot: 6.2172e-015
%     o       G: 1.3019e+003
%     o   DGDLs: 1.0484e+004
% o      p: 862.4522
% i  VWall: 1.6200e-004
% i   LsV0: 1.5000
% o      Z:

VWall=ParCavity.VWall;

V=max(0.001,ParCavity.V/VWall); % normalized cavity volume, safety for negative volume

H= VWall^(1/3)      ; % scaling length for LV
A= VWall*V/H; % crude estimate of cross-section of a (left) ventricle
ParCavity.A   = A; % cross-sectional area for valve inflow and outflow pressure
ParCavity.VDot= 0   ; % initialization
SDp= 1+3*V          ; % stress/pressure
Ef = (1+3*V).^(1/3)  ; % stress/pressure

Ls= ParCavity.LsV0 * Ef; % sarcomere length
ParCavity.Sarc.Ls= Ls;
ParCavity.Sarc=SarcMech(ParCavity.Sarc); % Sarcomere mechanics, stress, stiffness calculation

G = ParCavity.Sarc.G    ; % stress/Ls
GD= ParCavity.Sarc.DGDLs; % derivative dG/dLs defines incremental stiffness
C = VWall*(Ef.^6)./(Ls.*(Ls.*GD-2*G)); %  dVdp derived with Mathematica
Cgrav= A / 10300; % compliance related to gravity effects: area/(rho g)
C = C+Cgrav;

% estimate of wave resistance
rho=1050; L= 0.5*H*rho./A; % estimate of inertia of cavity filling
ParCavity.Z= sqrt(L./C)     ; % characteristic wave impedance of cavity

ParCavity.p= ParCavity.Sarc.G .* Ls ./Ef.^3; ...
%   -0.1*ParCavity.Sarc.dLsPas*ParCavity.Sarc.GPas*exp(-2*V); % safety for low volumes by suction
%   -ParCavity.Sarc.dLsPas*ParCavity.Sarc.GPas*exp(-2*V); % safety for low volumes by suction

ParCavity.V=V*VWall;
ParCavity.VDot=0;
ParCavity.qRemod=0;