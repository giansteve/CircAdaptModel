function strTot=MapStruct(Struct,varargin);
%function MapStruct(Struct,varargin);
%Commonly used as function MapStructure(Struct)
%Nested search in structure tree to map sizes of branch arrays within structure
%Single records (end branches) are not printed
%The varargin facility is needed for the specific nesting in this procedure

%INPUT
% Struct=stucture to map
% varargin{1}=nr, current branch number
% varargin{2}=depth, depth in structure tree
% varargin{3}=string array of output
%OUTPUT
% cell array of strings, containing output to print

%==== initialization of Root of structure
if nargin<4; % initialization of Root of structure
   nr=1; depth=0; strTot={};
   if isfield(Struct,'Name');
      strTot=[strTot;{Struct.Name}];
      disp(Struct.Name);
   end
else
   nr=varargin{1}; depth=varargin{2}; strTot=varargin{3};
end

if ~isstruct(Struct); return; end; % Escape if not a structure
Fn1=fieldnames(Struct);
n1=length(Fn1);

for i1=1:n1,
   Twig=getfield(Struct,{1},Fn1{i1});
   Len=length(Twig);
   str1=[]; if nr>1, str1=num2str(nr); end;
   str2=[]; if Len>1, str2=['(',num2str(Len),')']; end;
   str=[repmat(' ',[1,4*depth]), str1 ,'  ',Fn1{i1},str2];
   strTot=[strTot;{str}];
   %==== consider only first element of an array within a field
   i2=1; strTot=MapStruct(Twig(i2),i2,depth+1,strTot);
end

return
%==== part developed for Patient tree/ MRI tagging data
%for i1=1:n1,
%   Twig=getfield(Struct,{1},Fn1{i1});
%   Len=length(Twig);
%   if (~ischar(Twig))&(Len>1);
%      str=[repmat(' ',[1,4*depth]), num2str(nr),'  ',Fn1{i1},'(',num2str(Len),')'];
%      strTot=[strTot;{str}];
%      disp(str);
%      for i2=1:Len
%         strTot=MapStructure(Twig(i2),i2,depth+1,strTot);
%      end
%   end
%end


