function Par=CircAdapt(Par)
%function Par=CircAdapt(Par);
% Theo Arts, Maastricht University, Eindhoven University of Technology,
% April 3, 2004, email: t.arts@bf.unimaas.nl
% === Solution of differential equations, defined in parameter structure 'Par'
% Solution stored in Par.SVar

%===== timing of contraction
tCycleRef=0.85*(Par.Tube.q0/1e-4)^(0.28); %Scaling time for body size/systemic flow
Par.tCycle=round(Par.tCycle/Par.Dt)*Par.Dt; % to avoid starting in the middle of a Dt-interval
Par.t=round(Par.t/Par.tCycle)*Par.tCycle; % to avoid starting in the middle of the cycle
Par.La.Sarc.ActivationDelay = Par.Ra.Sarc.ActivationDelay + 0.01*tCycleRef/0.85; %Ra->La delay
Par.Rv.Sarc.ActivationDelay = Par.Ra.Sarc.ActivationDelay +0.13*tCycleRef/0.85; %AV-node delay
Par.Lv.Sarc.ActivationDelay = Par.Rv.Sarc.ActivationDelay;
Par.Lv.Sarc.TimeAct         = 0.12*tCycleRef +0.30*Par.tCycle ; % [s] moment of decay at Ls=LsStress0Act (s), (Arts, AJP, 1982)
Par.Rv.Sarc.TimeAct         = Par.Lv.Sarc.TimeAct;
Par.La.Sarc.TimeAct         = 0.16*tCycleRef; % [s] moment of decay at Ls=LsStress0Act (s)
Par.Ra.Sarc.TimeAct         = Par.La.Sarc.TimeAct;

%==== corrections

% solving differential equations
SVar=HrtPar2SVar(Par); % setting initial condition of differential equation
SVar=SVar(end,:);
t=SVar(end,1); tEnd=t+Par.tEnd;
while t<=tEnd % solving 1 cardiac cycle each time
    disp(['t=',num2str(t),';  Time to go= ',num2str(tEnd-t)]); pause(0.01);
    TimePoints=[0:Par.Dt:Par.tCycle]; % time interval to be added to the already existing solution
    %opt = odeset('RelTol',1e-4,'AbsTol',1e-5);
    %             somewhat higher relative accuracy and lower absolute accuracy than default
    opt = odeset('RelTol',1e-4,'AbsTol',1e-3);%low accuracy, faster
    [tc1,SVarAppend]=ode23('HrtSVarDot',...
        TimePoints,SVar(end,:),opt,Par); % solving of Differential Equations
    Par.SVar=SVarAppend;
    CircAdaptDisplay(Par);
    
    %===ADAPTATION after each cardiac cycle
    %Par.AdaptationType determines type of adaptation
    Par=Adapt(Par);
    
    SVar=[SVar;SVarAppend(2:end,:)]; % appends solution to existing solution
    SVarLast=SVar(end,:);
    [SVarDot,Par]=HrtSVarDot(0,transpose(SVarLast),[],Par);
    %     set Par to last SVar solution for continuation to next beat
    t=Par.t;
end
Par.SVar=SVar;

return
