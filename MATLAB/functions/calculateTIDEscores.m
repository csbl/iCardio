function [minRxnList, random_total] = calculateTIDEscores(model, minRxnList, data, num_iterations)
% calculateTIDEscores
% Calculates TIDE scores and associated significance for each task in the
% minRxnList. Details for the analysis can be found in the manuscript. 
%
% In brief, each TIDE score is calculated as the average reaction weight
% for reactions necessary for that metabolic task (excluding reactions
% without GPR rules). The significance is calculated by randomly shuffling
% the original gene expression data (data) and counting the number of
% random iterations that fall either greater or less than the task score,
% depending on how the task score fell relative to the mean task score for
% random samples. 
% 
% This function can be parallelized by running parpool() before calling the
% function. 
%
% INPUTS:
%    model:              a model structure
%    minRxnList:         a MATLAB structure generated using the 
%                        generateMinRxnList function where each entry 
%                        represents a metabolic task 
%    data:               log2 fold changes for genes in the original
%                        dataset, interpretation of results depends on
%                        whether data from all genes is included or only
%                        metabolic genes in the model (due to the
%                        randomized shuffling step)
%                        Data should be formated as a structure contaning a
%                        .gene and .value entry.
%    num_iterations:     the number of random shuffles that should be
%                        performed
%
% OUTPUTS:
%    minRxnList:         structure with the results:
%                        * id: corresponding to the task provided in the
%                        taskStructure
%                        * description: corresponding to the task provided
%                        in the taskStructure
%                        * rxns: minimum set of reactions necessary
%                        for that particular task to pass
%                        * taskScore: the calculated taskScore
%                        * rxnID: location of each reaction in the model
%                        * genesUsed: the gene used to weight each reaction
%                        in the model, 0 indicates no data or no change
%                        * significance: the calculated significance for
%                        the metabolic task
%    random_total:       the randomly generated taskScores that were used
%                        to calculate the final significance
% 
% .. Authors:
%       - Originally written by Bonnie Dougherty, 2020-04-08


% Replace mapExpressionToReactions function with individual function calls
% This step takes the longest, only need to run once
minSum = false;
parsedGPR = GPRparser(model,minSum);% Extracting GPR data from model

% Find wich genes in expression data are used in the model
[gene_id, gene_expr] = findUsedGenesLevels(model,data.gene, data.value);

% Change -1 values (no data) to 0
gene_expr(gene_expr == -1) = 0;

% Link the gene to the model reactions
[expressionRxns, genes_used] = selectGeneFromGPR(model, gene_id, gene_expr, parsedGPR, minSum);
for rxn = 1:length(expressionRxns)
    if expressionRxns(rxn) == -1
        expressionRxns(rxn) = 0;
    end
end

% Mark all genes_used that have a reaction score of 0 as 0
for k = 1:length(genes_used)
    if expressionRxns(k) == 0
        genes_used{1,k} = '0';
    end
end

% Calculate the actual TIDE score for the real data
for task = 1:length(minRxnList)
    % Translate reaction names into location IDs in the model
    current_rxns = [];
    for rxn = 1:length(minRxnList(task).rxns)
        [ans model_location] = max(strcmp(minRxnList(task).rxns{rxn}, model.rxns));
        current_rxns(end+1,1) = model_location;
    end
    current_task_reactions = expressionRxns(current_rxns);
    minRxnList(task).taskScore = sum(current_task_reactions)/length(minRxnList(task).rxns);
    minRxnList(task).rxnID = current_rxns;
    minRxnList(task).genesUsed = genes_used(current_rxns)';
end


% Randomize data and re-calculate a distribution of task scores

random_total = zeros(length(minRxnList),num_iterations);
genes = data.gene;
orig_data = data.value;
parfor x = 1:num_iterations
    % Randomize data within all genes that have data
    random_data = orig_data(randperm(length(orig_data)));
    
    % Find wich genes in expression data are used in the model
    [gene_id, gene_expr] = findUsedGenesLevels(model, genes, random_data);

    % Change -1 values (no data) to 0
    gene_expr(gene_expr == -1) = 0;

    % Link the gene to the model reactions
    [expressionRxns, gene_used] = selectGeneFromGPR(model, gene_id, gene_expr, parsedGPR, minSum);

    expressionRxns(:,x) = expressionRxns;
    expressionRxns(expressionRxns == -1) = 0;

    random_total_temp = zeros(length(minRxnList),1);
    for task = 1:length(minRxnList)
        current_task_reactions = expressionRxns(minRxnList(task).rxnID);
        random_total_temp(task) = sum(current_task_reactions)/length(minRxnList(task).rxnID);
    end
    random_total(:,x) = random_total_temp;
end

% Calculate significance for each task based off the random task scores
for task = 1:length(minRxnList)
    % Calculate p-value from ranked method, dependent on sign of the data
    % Do one-sided test based on value of expressionRxns
    
    % One-sided p-value from mean of the distribution
    if minRxnList(task).taskScore >= mean(random_total(task,:))
        minRxnList(task).significance = sum(random_total(task,:) >= minRxnList(task).taskScore)/num_iterations;
    else 
        minRxnList(task).significance = -sum(random_total(task,:) <= minRxnList(task).taskScore)/num_iterations;
    end
    
    % Check to see if the significance value is == 0
    if minRxnList(task).significance == 0
        if minRxnList(task).taskScore >= mean(random_total(task,:))
            minRxnList(task).signifiance = 0.0001;
        else
            minRxnList(task).significance = -0.0001;
        end
    end
end
end

function [gene_id, gene_expr] = findUsedGenesLevels(model, gene, value, printLevel)
% Returns vectors of gene identifiers and corresponding gene expression
% levels for each gene present in the model ('model.genes').
%
% USAGE:
%    [gene_id, gene_expr] = findUsedGenesLevels(model, exprData)
%
% INPUTS:
%
%   model:               input model (COBRA model structure)
%
%   exprData:            mRNA expression data structure
%       .gene                cell array containing GeneIDs in the same
%                            format as model.genes
%       .value               Vector containing corresponding expression value (FPKM)
%
% OPTIONAL INPUTS:
%    printLevel:         Printlevel for output (default 0);
%
% OUTPUTS:
%
%   gene_id:             vector of gene identifiers present in the model
%                        that are associated with expression data
%
%   gene_expr:           vector of expression values associated to each
%                        'gened_id'
%
%   
% Authors: - S. Opdam & A. Richelle May 2017

exprData.gene = gene;
exprData.value = value;

if ~exist('printLevel','var')
    printLevel = 0;
end

gene_expr=[];
gene_id = model.genes;

for i = 1:numel(gene_id)
        
    cur_ID = gene_id{i};
	dataID=find(ismember(exprData.gene,cur_ID));
	if isempty (dataID)
    	gene_expr(i)=-1;        
    elseif length(dataID)==1
    	gene_expr(i)=exprData.value(dataID);
    elseif length(dataID)>1    	
        if printLevel > 0
            disp(['Double for ',num2str(cur_ID)])
        end
    	gene_expr(i)=mean(exprData.value(dataID));
    end
end
           
end


function [expressionCol, gene_used] = selectGeneFromGPR(model, gene_names, gene_exp, parsedGPR, minSum)
% Map gene expression to reaction expression using the GPR rules. An AND
% will be replaced by MIN and an OR will be replaced by MAX.
%
% USAGE:
%   expressionCol = selectGeneFromGPR(model, gene_names, gene_exp, parsedGPR, minMax)
%
% INPUTS:
%   model:          COBRA model struct
%   gene_names:     gene identifiers corresponding to gene_exp. Names must
%                   be in the same format as model.genes (column vector)
%                   (as returned by "findUsedGeneLevels.m")
%   gene_exp:       gene FPKM/expression values, corresponding to names (column vector)
%                   (as returned by "findUsedGeneLevels.m")
%   parsedGPR:      GPR matrix as returned by "GPRparser.m"
%
% OPTIONAL INPUTS:
%   minSum:         instead of using min and max, use min for AND and Sum
%                   for OR
%
% OUTPUTS:
%   expressionCol:  reaction expression, corresponding to model.rxns.
%                   No gene-expression data and orphan reactions will
%                   be given a value of -1.
%
% AUTHOR: Anne Richelle, May 2017


if ~exist('minSum','var')
    minSum = false;
end
gene_used={};
for i=1:length(model.rxns)
	gene_used{i}='';
end
    
% -1 means unknown/no data
expressionCol = -1*ones(length(model.rxns),1); 
for i = 1:length(model.rxns)
    curExprArr=parsedGPR{i};
    curExpr= [];
    gene_potential=[];
    for j=1:length(curExprArr)
        if length(curExprArr{j})>=1
            geneID = find(ismember(gene_names,curExprArr{j}));
            %geneID = find(ismember(gene_names,str2num(curExprArr{j}{1})));
            % if the gene is measured
            if ~isempty(geneID) 
                if minSum
                    % This is an or rule, so we sum up all options.
                    curExpr= [curExpr, sum(gene_exp(geneID))]; 
                    gene_potential=[gene_potential, gene_names(geneID)'];
                else
                    % If there is data for any gene in 'AND' rule, take the minimum value
                    [minGenevalue, minID]=min(gene_exp(geneID));
                    curExpr= [curExpr, minGenevalue]; %If there is data for any gene in 'AND' rule, take the minimum value
                    gene_potential=[gene_potential, gene_names(geneID(minID))];
                end
            end
        end
    end

    if ~isempty(curExpr)
        if minSum
            % in case of min sum these are and clauses that are combined, so its the minimum.
            [expressionCol(i), ID_min]=min(curExpr);
            gene_used{i}=gene_potential(ID_min);
        else
            % Instead of taking the max, take the max of the absolute value
            if(max(curExpr) > 0 && min(curExpr) < 0)
                [maximum location] = max(curExpr);
                %[maximum location] = max(abs(curExpr));
                expressionCol(i) = curExpr(location);
                ID_max = location;
            else
                [maximum location] = max(abs(curExpr));
                expressionCol(i) = curExpr(location);
                ID_max = location;
            end

            % if there is data for any gene in the 'OR' rule, take the maximum value
            gene_used{i}=gene_potential(ID_max);
        end
    end
end

end