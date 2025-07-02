function  guessFixSpecies(obj,nsp)
            switch nargin
                case 1
                    nsp = 6;
                otherwise
            end
            [m] = max(obj.reference.moleFractions,[],2); 
            [~,indx] = sort(m,1,'descend');
            speciesList = obj.reference.speciesNames(indx(1:nsp));
            fprintf('ga:guessFixSpecies: Guess the following species for fixing:\n\t ')
            for k = 1:nsp-1
                fprintf([speciesList{k},', ']);
            end
            fprintf([speciesList{k+1},'\n']);
            obj.genotype = obj.genotype.fix(speciesList);
        end