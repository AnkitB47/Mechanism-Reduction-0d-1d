function blocks = indexCtiFile(ctiFile)
% blocks = indexCtiFile(ctiFile)
% indexes a cti file
% 
% blocks.nlines  : no of lines;
% blocks.species : startline of species block
% blocks.ideal_gas :  startline ideal_gas_block
% blocks.reactions :  startline of a reaction_block
% blocks.three_bodies :  startline of thethree_body_reaction_block
% blocks.falloffs     :  startlinefalloff_reaction_block
% blocks.units        :  startline of the unit_block

fid = fopen(ctiFile,'r');
nline = 0;
ideal_gas_block = [];
species_block = [];
reaction_block = [];
falloff_reaction_block = [];
three_body_reaction_block = [];
unit_block = [];
reaction_index = [];
while 1
    cline = fgetl(fid);     %read line
    nline = nline+1;            % count line
    
    % identify the keywords and the associate blocks
    if length(cline)>=10 && strcmp(cline(1:9),'ideal_gas')
        ideal_gas_block = [ideal_gas_block,nline];
    elseif length(cline)>=7 && strcmp(cline(1:7),'species')
        species_block = [species_block,nline];
    elseif length(cline)>=8 && strcmp(cline(1:8),'reaction')
        reaction_block = [reaction_block,nline];
        reaction_index = [reaction_index,nline];
    elseif length(cline)>=16 && strcmp(cline(1:16),'falloff_reaction')
        falloff_reaction_block = [falloff_reaction_block,nline];
        reaction_index = [reaction_index,nline];
    elseif length(cline)>=19 && strcmp(cline(1:19),'three_body_reaction')
        three_body_reaction_block = [three_body_reaction_block,nline];
        reaction_index = [reaction_index,nline];
    elseif length(cline)>=5 && strcmp(cline(1:5),'units')
        unit_block = nline;
    end
    
    if ~ischar(cline), 
        err = fclose(fid);
        break,
    end
end
blocks.nlines = nline-1;
blocks.species = species_block;
blocks.nspecies = length(blocks.species);
blocks.ideal_gas = ideal_gas_block;
blocks.reaction_index = reaction_index; 
blocks.reactions = reaction_block;
blocks.three_bodies = three_body_reaction_block;
blocks.falloffs = falloff_reaction_block;
blocks.units = unit_block;
blocks.nreactions = length(blocks.reactions)+length(blocks.three_bodies)+length(blocks.falloffs);

end