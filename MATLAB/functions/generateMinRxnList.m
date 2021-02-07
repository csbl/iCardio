function [minRxnList] = generateMinRxnList(model, taskStructure, removeNoGPR)
% generateMinRxnList
% Iterates through each task in the task list, setting bounds as are done
% in the checkMetabolicTasks functions (Richelle et al, 2017). Code for
% setting these bounds was directly taken from the checkMetabolicTasks
% function. 
% 
% After setting bounds for a particular task, pFBA is run to identify the
% minimum set of reactions that are necessary in order for that task to
% proceed. These reactions are saved in a structure with a list of
% reactions necessary for each task. 
%
% INPUTS:
%    model:              a model structure
%    taskStructure:      MATLAB representation of the taskStructure used to
%                        check tasks. Run the generateTaskStructure
%                        function in the COBRA toolbox. 
%    removeNoGPR:        'true' or 'false' for removing reactions that do  
%                        not have an associated GPR 
%
% OUTPUTS:
%    minRxnList:         structure with the results:
%                        * ID: corresponding to the task provided in the
%                        taskStructure
%                        * description: corresponding to the task provided
%                        in the taskStructure
%                        * reactions: minimum set of reactions necessary
%                        for that particular task to pass
% .. Authors:
%       - Originally written for RAVEN toolbox by Rasmus Agren, 2013-11-17
%       (checkMetabolicTasks)
%       - Adapted for cobratoolbox and modified to rely only on flux 
%       constraints by Richelle Anne, 2017-05-18 (checkMetabolicTasks)
%       - Adapted to include pFBA constraints and identify reactions
%       necessary for each task by Bonnie Dougherty, 2020-04-08

minRxnList = struct;

printOutput=true;
printOnlyFailed=false;
getEssential=false;

% Create new irreversible model to compare for pFBA solution
compare_model = convertToIrreversible(model);

% CHECK the format of the model
if size(model.rxns,2)>size(model.rxns,1)
    model.rxns=model.rxns';
end
if size(model.rxnNames,2)>size(model.rxnNames,1)
    model.rxnNames=model.rxnNames';
end
if size(model.rules,2)>size(model.rules,1)
    model.rules=model.rules';
end

if isfield(model,'grRules') && size(model.grRules,2)>size(model.grRules,1)
    model.grRules=model.grRules';
end

%Find all exchange/demand/sink reactions
Exchange = {};
for k=1:length(model.rxns)
    if sum(abs(model.S(:,k))) == 1
        Exchange(end+1) = model.rxns(k);
    end
end
Exchange=unique(Exchange);

%Close all exchange reactions
model.lb(findRxnIDs(model,Exchange))=0;
model.ub(findRxnIDs(model,Exchange))=0;

for i= 1:length(taskStructure)
     
    clear tModel
    tModel=model;
    modelMets=upper(tModel.mets);

    %%SETUP of the input model
    %suppress objective function if any
    tModel.c(tModel.c==1)=0;
    tModel.csense(1:length(model.b),1) = 'E';

    taskReport{i,1}=taskStructure(i).id;
    taskReport{i,2}=taskStructure(i).description;

    %Set the inputs
    if ~isempty(taskStructure(i).inputs)

        rxn_Subs={};
        for n=1:length(taskStructure(i).inputs)
            INPUT=taskStructure(i).inputs(n);
            INPUT=INPUT{1};
            match_INPUTS = strncmpi(INPUT,modelMets,length(INPUT(1:end-3)));
            match_INPUTS = modelMets(match_INPUTS==1);

            % Identify all compartments where the metabolite is located
            compSymbol={};
            for k=1:length(match_INPUTS)
                [tokens] = regexp(match_INPUTS{k},'(.+)\[(.+)\]','tokens');
                Symb = tokens{1}{2};
                compSymbol{end+1} = Symb;
            end

            % Set the exchange reactions for the inputs
            % If the metabolites already exist extracellularly
            AddExchange=0;
            if ismember('S',compSymbol)==1
                Tsp_ID=findRxnIDs(tModel,findRxnsFromMets(tModel,INPUT));
                Tsp_rxn = full(tModel.S(:,Tsp_ID));
                Nb_React=sum(abs(Tsp_rxn),1);

                % If an exchange reaction already exist
                if isempty(Nb_React==1)==0
                    ID_exc=find(Nb_React==1);
                    % Remove the existing exchange reaction
                    K=ismember(Exchange,tModel.rxns(Tsp_ID(ID_exc)));
                    tModel=removeRxns(tModel,tModel.rxns(Tsp_ID(ID_exc)));
                    if (sum(K) > 0)
                        % Remove the exchange reaction that was removed
                        % from the model from the list of exchange rxns
                        [temp location] = max(K);
                        Exchange(location) = [];
                    end
                    AddExchange=1;
                else
                    AddExchange=1;
                end
            else
            	AddExchange=1;
            end

            % Add a temporary exchange reaction that allows the import of
            % the metabolite
            if AddExchange==1
                % If the input is also member of the outputs, let the exchange reversible
                warning off
                if ismember(INPUT,taskStructure(i).outputs)==1
                	[tModel]=addReaction(tModel,['temporary_exchange_',INPUT(1:end-3),'_INPUT'],[' <=> ',INPUT],[],[],-1000,1000,[], [], [], [], [], [],0);
                else
                	[tModel]=addReaction(tModel,['temporary_exchange_',INPUT(1:end-3),'_INPUT'],[' => ',INPUT],[],[],taskStructure(i).LBin(n),taskStructure(i).UBin(n),[], [], [], [], [], [],0);
                end
                warning on
                rxn_Subs(end+1) = {['temporary_exchange_',INPUT(1:end-3),'_INPUT']};
            end

            % Definition of the compartment for the transport reaction
            if ischar(taskStructure(i).COMP)==1
                comp_used=taskStructure(i).COMP;
                if strcmpi(comp_used,'[s]')==1
                    continue
                end
            elseif ismember('C',compSymbol)==1
                comp_used='[c]';
            elseif ismember('M',compSymbol)==1
                comp_used='[m]';
            elseif ismember('N',compSymbol)==1
                comp_used='[n]';
            elseif ismember('X',compSymbol)==1
                comp_used='[x]';
            elseif ismember('L',compSymbol)==1
                comp_used='[l]';
            elseif ismember('R',compSymbol)==1
                comp_used='[R]';
            end

            % Set the transport reactions for the inputs
            % Find existing transporters associated with the input
            AddTransport=0;
            Tsp_ID=findRxnIDs(tModel,findRxnsFromMets(tModel,INPUT));
            Tsp_rxn = full(tModel.S(:,Tsp_ID));
            Nb_React=sum(abs(Tsp_rxn),1);
        end
    end

    modelMets=upper(tModel.mets);
    [I J]=ismember(upper(taskStructure(i).inputs),modelMets);
    J=J(I);
    %Check that all metabolites are either real metabolites
    if ~all(I)
        fprintf(['ERROR: Could not find all inputs in "[' taskStructure(i).id '] ' taskStructure(i).description '"\n']);
    	taskReport{i,3}='Could not find all inputs';
    	notPresent=notPresent+1;
    end
    if numel(J)~=numel(unique(J))
    	dispEM(['The constraints on some input(s) in "[' taskStructure(i).id '] ' taskStructure(i).description '" are defined more than one time']);
    end

    %Set the outputs
    if ~isempty(taskStructure(i).outputs)

        rxn_Prod={};
        for n=1:length(taskStructure(i).outputs)
            OUTPUT=taskStructure(i).outputs(n);
        	OUTPUT=OUTPUT{1};
            
            if strcmp(OUTPUT,'ALLMETS')
                % exchange = findRxnIDs(tModel, Exchange);
                for p = 1:length(Exchange)
                    tModel.ub(findRxnIDs(tModel,Exchange(p)))=1000;
                end
            end

        	%skip the setup if output is also input as it has already been
        	%setup
            if ismember(upper(OUTPUT),upper(taskStructure(i).inputs))==1
            	continue
            end

            match_OUTPUTS = strncmpi(OUTPUT,modelMets,length(OUTPUT(1:end-3)));
            match_OUTPUTS = modelMets(match_OUTPUTS==1);
            compSymbol={};
            for k=1:length(match_OUTPUTS)
                [tokens] = regexp(match_OUTPUTS{k},'(.+)\[(.+)\]','tokens');
            	Symb = tokens{1}{2};
            	compSymbol{end+1} = Symb;
            end

            % Set the exchange reactions for the outputs
            % If the metabolites already exist extracellularly
            AddExchange=0;
            if ismember('S',compSymbol)==1
            	Tsp_ID=findRxnIDs(tModel,findRxnsFromMets(tModel,OUTPUT));
                Tsp_rxn = full(tModel.S(:,Tsp_ID));
                Nb_React=sum(abs(Tsp_rxn),1);
                % If an exchange reaction already exist
                if isempty(Nb_React==1)==0
                	ID_exc=find(Nb_React==1);
                    % Remove the existing exchange reaction
                    tModel=removeRxns(tModel,tModel.rxns(Tsp_ID(ID_exc)));
                    AddExchange=1;
                else
                    AddExchange=1;
                end
            else
                AddExchange=1;
            end

            % Add a temporary exchange reaction that allows the export of
            % the metabolite
            if AddExchange==1
                warning off
            	[tModel]=addReaction(tModel,['temporary_exchange_',OUTPUT(1:end-3)],[OUTPUT,' => '],[],[],taskStructure(i).LBout(n),taskStructure(i).UBout(n),[], [], [], [], [], [],0);
                warning on
                rxn_Prod(end+1) = {['temporary_exchange_',OUTPUT(1:end-3)]};
            end

            % Definition of the compartment for the transport reaction
            if ischar(taskStructure(i).COMP)==1
            	comp_used=taskStructure(i).COMP;
                if strcmpi(comp_used,'[e]')==1
                	continue
                end
            elseif ismember('C',compSymbol)==1
            	comp_used='[c]';
            elseif ismember('M',compSymbol)==1
            	comp_used='[m]';
            elseif ismember('N',compSymbol)==1
            	comp_used='[n]';
           	elseif ismember('X',compSymbol)==1
                comp_used='[x]';
            elseif ismember('L',compSymbol)==1
             	comp_used='[l]';
            elseif ismember('R',compSymbol)==1
                comp_used='[R]';
            end

            % Set the transport reactions for the outputs
            % Find existing transporters associated with the ouput
            AddTransport=0;
            Tsp_ID=findRxnIDs(tModel,findRxnsFromMets(tModel,OUTPUT));
            Tsp_rxn = full(tModel.S(:,Tsp_ID));
            Nb_React=sum(abs(Tsp_rxn),1);
                    end
    end

    if ~isempty(taskStructure(i).EQU)
        currentTask = taskStructure(i).EQU;
        for p = 1:length(taskStructure(i).EQU)
            equation = currentTask{p,1};
            UB = taskStructure(i).EQULB(p,1);
            LB = taskStructure(i).EQULB(p,1);
            if UB == 0 
                UB = 1000;
            end
            if LB == 0
                LB = -1000;
            end
            % add reaction to the model
            [tModel, rxnIDexists] = addReaction(tModel, strcat('TEMPORARY_',taskStructure(i).id,'_',num2str(p)), 'reactionFormula', equation, 'lowerBound', LB, 'upperBound', UB);
            
            % Check to see if the reaction is already in the model
            % temp = 
        end
    end

    modelMets=upper(tModel.mets);
    [I J]=ismember(upper(taskStructure(i).outputs),modelMets);
    J=J(I);

    %Check that all metabolites are either real metabolites
    if ~all(I)
        fprintf(['ERROR: Could not find all outputs in "[' taskStructure(i).id '] ' taskStructure(i).description '"\n']);
        taskReport{i,3}='Could not find all outputs';
        notPresent=notPresent+1;
    end
    if numel(J)~=numel(unique(J))
        dispEM(['The constraints on some output(s) in "[' taskStructure(i).id '] ' taskStructure(i).description '" are defined more than one time']);
    end


	%Solve the constrained problem
	tModel.csense(1:length(tModel.mets),1) = 'E';
	tModel.osense = -1;
	tModel.A=tModel.S;
	sol=solveCobraLP(tModel);
    
    % Check to see if a solution exists:
    if ~isempty(sol.full) && sum(abs(sol.full))~=0
         % Identify the reactions that are necessary for the model
        [MinimizedFlux modelIrrev]= pFBA_edit(tModel,0);
        carry_temp = MinimizedFlux.full ~= 0;
        carry_temp = modelIrrev.rxns(carry_temp);
        
        carrying_rxns = {};
        for a = 1:length(carry_temp)
            % remove the b/f tag to the reaction
            if length(carry_temp{a}) > 7
                [ans model_location] = max(strcmp(carry_temp{a}(1:8), model.rxns));
                if ans ~= 0
                    carrying_rxns{end+1,1} = carry_temp{a}(1:8);
                end
            end
        end
        
        % Create a new entry in minRxnList
        minRxnList(end+1).id = taskStructure(i).id;
        minRxnList(end).description = taskStructure(i).description;
        minRxnList(end).rxns = carrying_rxns;
    end
    
end

% Remove the first entry which initialized the data structure
minRxnList(1) = [];

if strcmp(removeNoGPR,'true')
    % Remove reactions in the task list that don't have a GPR
    % Identify indexes of reactions without GPR rules
    for k = 1:length(model.rxns)
        noGPR(k) = isempty(model.grRules{k});
    end
    noGPR = model.rxns(noGPR);

    % Remove these reactions from the list of total_carrying_rxns
    for task = 1:length(minRxnList)
        common = intersect(noGPR, minRxnList(task).rxns);
        minRxnList(task).rxns = setdiff(minRxnList(task).rxns, common);
    end

    % Remove tasks that now have no reactions (i.e. all carrying reactions had
    % no GPR rule)
    task = 1;
    while task ~= length(minRxnList)
        if isempty(minRxnList(task).rxns)
            minRxnList(task) = [];
        else
            task = task+1;
        end
    end
end

end

