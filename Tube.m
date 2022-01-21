function Par=Tube(Par);
%function Par=Tube(Par);
%Theo Arts, University of Maastricht / Technological University of Eindhoven, Feb 2004.
%Elastic Tube hemodynamics, related to the pressure/flow wave (e.g. aorta)
%
%INPUT/OUTPUT i/o
%Par.
% i  k     stiffness wall material
% i  rhob  density of blood
% i  p0    reference pressure
% i  Len   length of vessel for windkessel function
% i  A0    cross-section of lumen at reference
% i  AWall cross-section of wall
% i  V     volume of vessel for windkessel function
% o  A     cross-section of lumen
% o  p     pressure of windkessel compliance
% o  Z     wave impedance

AWalld3=Par.AWall/3;
Par.A=max(0.01*Par.A0,Par.V/Par.Len); %mean cross-section of tube, safe for A0<0
Par.p=Par.p0*((Par.A+AWalld3)./(Par.A0+AWalld3)).^(Par.k/3-1);
Par.Z=sqrt( Par.rhob*(Par.k-3) * abs(Par.p)./(Par.A.*(3*Par.A+Par.AWall)));

%=== Initializations
Par.VDot=0; %inflow
Par.qRemod=0; %Through flow
