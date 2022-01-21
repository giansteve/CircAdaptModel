The CircAdapt model of heart and circulation dynamics with adaptation
(T. Arts, Maastricht University / Eindhoven University of Technology, Nov, 2004)

The setup of the CircAdapt model has been described in the paper "Adaptation to mechanical load determines shape and properties of heart and circulation, the CircAdapt model". The name of the main executable program script is 'CircAdaptMain'.

Data structure

The state of the model has been listed in a tree-like structure 'Par', in which parameters are sorted per module. After the program has been executed, this structure is in memory, visible by the MatLab instruction <whos>.The structure can be listed with instruction <MapStruct(Par)>. The names of the variables are quite logic, and can be reconstructed with the help of the article.

Normally the program gets its starting information from the file Par.mat, in which the last Par-structure has been stored. Starting conditions can be manipulated by loading a *(mat file with structure Par in it. Par may be changed. By storing this new Par structure with the instruction >> save Par Par  , the program starts with this new Par-structure


The program

CircAdaptMain decides whether to start from scratch (New), from the Reference settings stored in ParRef.mat,  to continue from the last simulation stored in Par.mat (Continue or just <Ret>), perform an Intervention, or Undo the last simulation result (Undo).

Before starting a simulation after 'Intervention' or 'Reference', interventions may be carried out like: Pressure, Flow, Cycle Time, Duration of simulation, or adaptation of at Rest or Exercise.

Fetus:
For the fetus, there was a slight modification of the pulmonary flow equation structure, leading to the files with 'Foet' in the name


Running the program
The simplest start of the program: put all files in one directory, add this directory to the MatLab path, press:

>> CircAdaptMain
[N]ew, [I]ntervention, [R]eference, [U]ndo last, [C]ontinue <Ret>: r
 
[P]ressure                  (kPa): 12.2
[F]low                     (ml/s): 100
cycle [T]ime                 (ms): 850
[D]uration simulation         (s): 0.1
Adaptation [R]est ot [E]xercise  : No
<Ret> Continue
Choose PFTDRE<Ret>: 
t=0.85;  Time to go= 0.1
sequence La, Ra, Lv, Rv
Contractility: 0.078189     0.07827
Contractility: 0.048948    0.048955
Contractility: 0.82572     0.82586
Contractility: 0.7581     0.75811
1e3*mean([Par.ValveLArt.q,Par.ValveRArt.q,Par.ValveLVen.q,Par.ValveRVen.q])
    0.1000    0.1001    0.1000    0.1000

sequence adaptation: LArt,RArt,LVen,RVen
Differential equation has been solved
>>

Pressing return starts the simulation. Graphic windows will pop up.

All hemodynamics are stored in the combination of the state variables StateVar and the parameter setting Par.


