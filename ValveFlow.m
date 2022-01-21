function [ParProx,ParValve,ParDist]=ValveFlow(ParProx,ParValve,ParDist);
%function [ParValve,ParProx,ParDist]=Valve(ParValve,ParProx,ParDist);

%Theo Arts, University of Maastricht / Technological University of Eindhoven, Feb 2004.
%Valve hemodynamics
%Calculates derivative of flow (qDot) as a function of flow q and pressure drop Dp
%This function can also be used for an orifice, and backward valve by making AOpen<ALeak
%INPUT/OUTPUT i/o
% i ParValve,ParProx,ParDist
% o ParValve.qDot

q=ParValve.q;

switch ParValve.FunctionNr
    
    case 1, %flow continuity to elements in environment
        ParProx.VDot=ParProx.VDot-q; % outflowflow to calculate through flow
        ParDist.VDot=ParDist.VDot+q;
        ParProx.qRemod=ParProx.qRemod+abs(q);
        ParDist.qRemod=ParDist.qRemod+abs(q);
        ParProx.FunctionNr=1;
        ParDist.FunctionNr=1;
        ParValve.FunctionNr=2;
        
    case 2, %Time derivatives of volumes + boundary pressures
        ParProx.pIn =ParProx.p+ParProx.VDot .*ParProx.Z;
        ParDist.pIn =ParDist.p+ParDist.VDot .*ParDist.Z;
        ParValve.FunctionNr=3;
        
    otherwise, %Time derivative of flow, VALVE function, bi-directionally symmetric
        Dp=ParProx.pIn-ParDist.pIn; %pressure drop
        Q=(q>0); P=(Dp>0);
        AFw=( Q| P)*ParValve.AOpen;
        ABw=(~Q|~P)*ParValve.ALeak;
        A=max(AFw,ABw); %opening area
        
        v1=q./ParProx.A; v2=q./A; rho=ParValve.rhob; v3=q./ParDist.A;
        v2max=max([v1.^2,v2.^2,v3.^2],[],2);
        DpB=(0.5*rho)*(Q.*(v2max-v1.^2)-(~Q).*(v2max-v3.^2));
        %                       Bernouilli pressure drop     
        ParValve.qDot=(Dp-DpB).*A/(rho*ParValve.Len); %flow derivative
             
        ParValve.FunctionNr=1;
end

return
