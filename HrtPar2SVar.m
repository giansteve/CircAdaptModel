function SVar=HrtPar2SVar(Par)
%function SVar=Par2SVar(Par);
% Theo Arts, Maastricht University, Eindhoven University of Technology,
% April 3, 2004, email: t.arts@bf.unimaas.nl

ColOnes=ones(size(Par.Lv.V)); % auxilary column vector
S= Par.t                ; %time
S=[S,Par.TubeLArt.V    ]; %left arterial volume
S=[S,Par.TubeLVen.V    ]; %left venous volume
S=[S,Par.TubeRArt.V    ]; %right arterial volume
S=[S,Par.TubeRVen.V    ]; %right venous volume
S=[S,Par.ValveLArt.q   ]; %aortic valve flow
S=[S,Par.ValveLAv.q    ]; %mitral valve flow
S=[S,Par.ValveLVen.q   ]; %venous outflow
S=[S,Par.ValveRArt.q   ]; %pulmonary valve flow
S=[S,Par.ValveRAv.q    ]; %tricuspid valve flow
S=[S,Par.ValveRVen.q   ]; %right venous outflow
S=[S,Par.ValveDUCT.q  ]; %ductus arteriosis flow
S=[S,Par.ValveVSD.q    ]; %ductus arteriosis flow
S=[S,Par.ValveASD.q   ]; %ductus arteriosis flow
S=[S,Par.La.V          ]; %LA cavity volume
S=[S,Par.Lv.V          ]; %LA cavity volume
S=[S,Par.Ra.V          ]; %LA cavity volume
S=[S,Par.Rv.V          ]; %LA cavity volume
S=[S,Par.La.Sarc.Lsi   ]; %contractile unit length of sarcomere
S=[S,Par.La.Sarc.C     ]; %contractility of sarcomere
S=[S,Par.Lv.Sarc.Lsi   ]; %contractile unit length of sarcomere
S=[S,Par.Lv.Sarc.C     ]; %contractility of sarcomere
S=[S,Par.Ra.Sarc.Lsi   ]; %contractile unit length of sarcomere
S=[S,Par.Ra.Sarc.C     ]; %contractility of sarcomere
S=[S,Par.Rv.Sarc.Lsi   ]; %contractile unit length of sarcomere
S=[S,Par.Rv.Sarc.C     ]; %contractility of sarcomere
SVar=S./(ColOnes*Par.Scale'); % Conversion to SVar matrix

return