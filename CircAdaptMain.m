% CircAdaptMain
% CircAdaptMain, a model of the dynamic, adapting circulation
% Theo Arts, Maastricht University, Eindhoven University of Technology,
% April 3, 2004, email: t.arts@bf.unimaas.nl
% All parameters are stored in structure 'Par'. The solution of the differential equation
% is stored in a data matrix Par.SVar with columns of time functions.
% Some data are displayed graphically
% Map of structure 'Par' obtained by execution of 'mapstruct(Par)'

%% GM settings
clearvars
close all
clc

% graphic settings
set(0,'DefaultFigureWindowStyle','default')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex')
set(0,'defaultAxesFontSize',11)

%% Code

if isfile('Par.mat') 
    load Par
    if ~isfield(Par,'AdaptationType')
        Par.AdaptationType='o'; %empty
    end
    
    ShowMenu=1;
    while ShowMenu
        c2=input('[N]ew, [I]ntervention, [R]eference, [U]ndo last, [C]ontinue <Ret>: ','s');
        a=[c2,'c']; c=a(1); % convert <Ret> to 'c'
        switch lower(c)
            case('u') % undo last
                load ParLast;
                c='i';
                ShowMenu=false;
            case('i') % intervention
                load Par;
                c='i';
                ShowMenu=false;
            case('r') % reference
                load ParRef;
                c='i';
                ShowMenu=false;
            case('n') % new
                Par=HrtParNew;
                c='i';
                ShowMenu=false;
            case('c') % continue
                load Par;
                save ParLast Par;
                ShowMenu=false;
            otherwise
                ShowMenu = true;
        end
    end
else
    if isfile('ParRef.mat')
        c2=input('[N]ew, [R]eference: ','s');
        a=[c2,'c']; c=a(1); % convert <Ret> to 'c'
        switch lower(c)
            case('r') % reference
                load ParRef;
                c='i';
            otherwise % new
                Par=HrtParNew;
                c='i';
        end
    else
        Par=HrtParNew; % Parameter initialization, some remodeling rules inclusive
        % Generates parameter structure Par and initial conditions of the variables Par
        c='i';
    end
end

%tCycle=850*(q/100)^(log(350/850)/log(4)); %estimate of tCycle


Par.AdaptationType='t'; %standard contractility adaptation
switch lower(c)
    case 'i'
        OK=1;
        AdaptStr='No';
        while OK
            disp(' '); 
            disp(['[P]ressure                  (kPa): ',num2str(Par.p0/1e3)]); 
            disp(['[F]low                     (ml/s): ',num2str(Par.q0*1e6)]);
            disp(['cycle [T]ime                 (ms): ',num2str(Par.tCycle*1e3)]);
            disp(['[D]uration simulation         (s): ',num2str(Par.tEnd)]);
            disp(['Adaptation [R]est ot [E]xercise  : ',AdaptStr]);
            disp('<Ret> Continue');
            c1=input('Choose PFTDRE<Ret>: ','s');
            switch lower(c1)
                case 'p'
                    Par.p0=input('Mean Arterial Pressure (kPa): ')*1e3;
                    Par.LRp.R = Par.p0/Par.q0; %initial estimates
                    Par.RRp.R = Par.pDropPulm/Par.q0; %initial estimates
                case 'f'
                    Par.q0=input('Systemic Flow (ml/s): ')/1e6;
                    Par.LRp.R = Par.p0/Par.q0; %initial estimates
                    Par.RRp.R = Par.pDropPulm/Par.q0; %initial estimates
                case 't'
                    Par.tCycle = input('Cycle Time (ms): ')*1e-3;
                case 'd'
                    Par.tEnd = input('Duration of simulation (s): ');
                case 'e'
                    Par.AdaptationType = 'e';
                    AdaptStr = 'Exercise';
                case 'r'
                    Par.AdaptationType = 'r';
                    AdaptStr = 'Rest';
                otherwise
                    OK = 0;
            end
        end
    otherwise
end


% === Solves SVar for problem defined in parameter structure 'Par'
Par = CircAdapt(Par);
Par.AdaptationType = 't';%back to normal non-adaptation, acute

%=== Saving State Variables and Model Parameters
save Par Par;
disp('Differential equation has been solved');

ParDisplay=CircAdaptDisplay(Par); % graphical display

