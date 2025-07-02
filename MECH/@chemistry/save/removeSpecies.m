function   removeSpecies(obj,varargin)
%  removeSpecies(mechFileOut,deadSpeciesIndexList)
%  removeSpecies(baseMech,mechFileOut,deadSpeciesIndexList,deadSpeciesNamesList)
%   removes species given by deadSpeciesIndexList and deadSpeciesNamesList
%   from a mechanism file. If baseMech is a cti file, removeSpeciesCTI
%   will used, is mechFIleIn a inp-file, removeSpeciesCK will be used. The
%   type of the format of the input file is detected by the ending inp or
%   cti.
%   New Version.
%   (C) 2013 Uwe Prüfert 

switch nargin
    case 3
        baseMech = obj.mechanism;        
        mechFileOut = varargin{1};
        if strcmp(class(varargin{2}),'cell')
             deadSpeciesNamesList = varargin{2};
             deadSpeciesIndexList = obj.getSpeciesIndexByName(...
                 deadSpeciesNamesList);
        else
            deadSpeciesIndexList = varargin{2}
            deadSpeciesNamesList = obj.speciesNames(...
                deadSpeciesIndexList);
        end
    case 5
        baseMech =  varargin{1};       
        mechFileOut = varargin{2};
        deadSpeciesIndexList = varargin{3};
        deadSpeciesNamesList = varargin{4};
    otherwise
end


     
% identify the input file type
dotloc = max(strfind(baseMech,'.'));
if dotloc > 1
   filetype = baseMech(dotloc+1:end);
end
switch filetype
    case 'inp'
        % call CK version
        chemistry.removeSpeciesCK(baseMech,mechFileOut,...
            deadSpeciesNamesList)
    case 'cti'
        % call CTI version
        chemistry.removeSpeciesCTI(baseMech,mechFileOut,...
            deadSpeciesIndexList,deadSpeciesNamesList)
    otherwise
        error('removeSpecies:wrongFileType','Input-Output can be *.inp or *.cti')
end

end

