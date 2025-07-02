% function species = scanSpeciesList(elementsLine)
% % utility function for removeSpecies
% % it reads the species from the elementsLine and gives back a cell array
% % containing the species names
% 
% species = [];
% l = 0;
% candidate = [];
% spec = [];
% for k = 1:length(elementsLine)
%     
%     
%     if  ~strcmp(elementsLine(k),' ')
%         spec = [spec,elementsLine(k)];
%     end
%     if strcmp(elementsLine(k),' ')
%         l = l+1;
%         if ~isempty(spec)
%             candidate{l}= spec;
%             spec = [];
%         end
%     end 
% end
% l = 1;
% for k = 1:length(candidate)
%     if ~(isempty(candidate{k})||strcmp(candidate{k},'""",'))
%         species{l}=candidate{k};
%         l = l+1;
%     end
% end
% 
% 
% end



