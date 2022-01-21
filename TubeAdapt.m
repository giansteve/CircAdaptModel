function ParTube=TubeAdapt(ParTube);
% function ParTube=TubeAdapt(ParTube);
% Theo Arts, Maastricht University, Eindhoven University of Technology,
% April 3, 2004, email: t.arts@bf.unimaas.nl
% Vessel diameter adapts to flow, wall thickness to maximum pressure
% ParTube. 
%INPUT  i
%OUTPUT o
% i  q0    mean tube flow
% i  pMax  maximum pressure
% i  p0    mean working pressure
% i  Len   effective length of vessel bed
% o  A0    reference cross-section
% o  AWall cross-section of wall
% o  V     volume of vessel for windkessel function
% i  Adapt
% i      WallStress
% i      WallShear

A0  =(pi*(4*ParTube.q0/ParTube.Adapt.WallShear)^2)^(1/3);
pMax=abs(ParTube.pMax);
SfMax=ParTube.Adapt.WallStress;
p0  =ParTube.p0;
k   =ParTube.k;

p0=max(p0,0.1*pMax); %safety for collapse
ParTube.AWall=3*A0/((SfMax/pMax)*(p0/pMax)^(3/(k-3))-1);
ParTube.A0=A0;

ParTube.V=A0*ParTube.Len;
ParTube=Tube(ParTube);

return

%Tested OK:SfMax=pMax^(k/(k-3))*(1+3*A0/ParTube.AWall)/p0^(3/(k-3))

%EQUATIONS likely to be about OK
%ex^3=(1+3AMax/AWall)/(1+3A0/AWall);
%s0=(1+3A0/AWall)*p0
%sMax=s0*ex^k
%sMax=(1+3AMax/AWall).pMax
%known: k,A0,pMax,p0,sMax
%to be solved: AWall
