function  solve(obj,options)
            % g.solve runs the GA to optimize the mechanism
            % g.solve(n) uses n generations
            % if n is not given, n is set internally to 5
            if nargin<2
                n = 5;
                tresh = 6;
                select_the = 'survivors';
                nind = [];
                prob = 0.3;
            else
                n = options.ngen;
                select_the = options.selectionStrategie;
                tresh = options.selectionThreshold;
                nind = options.nind;  
                prob = options.prob;
            end
            if isempty(obj.currentGeneration)
                 obj.makeFirstGeneration(nind);
            end
            try
                for k = 1:n
                    try
                        obj.evaluate; 
                    catch ME
                        ME.rethrow
                    end
                    fprintf(['Evaluation of the ',...
                        num2str(k),' generation successful.\n']);  
                    obj.select(select_the,tresh);
                    best  = obj.currentGeneration(1:min(length(obj.currentGeneration),3));
                    obj.xing;
                    purex = obj.currentGeneration;
                    obj.mutateOne; 
                    xandmutate1 = obj.currentGeneration;
                    obj.mutate(prob);
                    obj.currentGeneration = ...
                        [purex xandmutate1 obj.currentGeneration best ];                      
                end
                fprintf('Last generation...\n')
                try 
                    obj.evaluate;                     
                catch ME
                    ME.rethrow
                end
                obj.select(select_the,tresh); 
            catch ME
                fprintf(['ga.solve:evaluateFails: ',...
                    'Catch an error send from evaluate with message: ',ME.message,...
                    'Try to return.\n']);
                ME.rethrow
            end
            fprintf('\n\t\tsolved GA sucessfuly...\n\t\t       *******\n')
        end