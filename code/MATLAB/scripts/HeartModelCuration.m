% Exploring different integration algorithm for building a draft heart
% metabolic network model
% Curation of chosen draft cardiomyocyte-specific metabolic model
% Analysis of curated cardiomyocyte-specific model
% 5/21/2019 BVD
clear all
close all
clc

initCobraToolbox
changeCobraSolver('ibm_cplex', 'all');

% Load in the draft models
load('data/draft_heart_models.mat')

% Assign the draft model as the CORDA model
heart_model_draft = tissueModel_CORDA;

%% Load draft model and test with task list
% Load tasks from original hepatocyte model
% Convert from original format to COBRA format
inputFile = ['data/AllTasks_CardiomyocyteSpecific.xlsx'];
[FINAL] = generateCobraTaskList(inputFile, hsa_cobra_load);
xlswrite('data/AllTasks_CardiomyocyteSpecific_COBRA.xlsx', FINAL, 'TASKS')
inputFile = ['data/AllTasks_CardiomyocyteSpecific_COBRA.xlsx'];

draft_tasks = checkMetabolicTasks_BVD(heart_model_draft,inputFile);

%% Add reactions back to ensure that cardiomyocyte tasks pass
add_rxns = {
    % For checking ATP hydryolsis
    'RCR11017'
    % For dATP synthesis
    'RCR10480'
    'RCR10087'
    'RCR10100'
    % For dGTP synthesis
	'RCR10481'
	'RCR10482'
    'RCR10483'
	'RCR14272'
	'RCR10099'
    % For dCTP synthesis
    'RCR10484'
	'RCR10127'
    % For dUTP synthesis
	'RCR14273'
    'RCR10096'
	% ATP salvage from hypoxanthine and adenosine
    'RCR40226'
	% fructose degradation
    'RCR40148'
	% glucuronate de novo synthesis
    'RCR10458'
    'RCR10459'
	% cardiolipin synthesis
    'RCR10495'
    'RCR10488'
    'RCR10854'
    'RCR13155'
    'RCR13156'
    'RCR13157'
    'RCR13158'
    'RCR13159'
    'RCR14157'
    'RCR20542'
    'RCR20543'
    'RCR20616'
	% PS, SM, ceramine, lactosylceramide de novo synthesis
    'RCR11289'
    'RCR20666'
    'RCR41633'
	% NAD and NADP de novo synthesis
    'RCR40217'
	% FAD de novo synthesis
    'RCR13931'
    'RCR40301'
 	% Ubiquinol-to-ATP
    'RCR30165'
	% arginine to nitric oxide
    'RCR11265' 
    'RCR11263' 
    'RCR40988' 
	% glyogen storage from glucose
    'RCR90154'
	% asparagine breadown to OAA
    'RCR10090'
	% proline degradation to Krebs cycle intermediate
    'RCR10316'
    % glycine degradation to Krebs cycle intermediate
    'RCR20075'
	% cysteine to Krebs cycle intermediate
    'RCR11646'
    'RCR11645'
	% beta-alanine uptake
    'RCR40108'
	% pyroxidal-P uptake
    'RCR40430'
	% cholesterol synthesis
    'RCR12925'
    'RCR13838'
	% ATP from histidine
    'RCR11325'
    'RCR11326'
    'RCR13803'
    'RCR14328'
    'RCR14329'
    % methionine ATP yield
    'RCR13862'
    'RCR11648'
    % phenylalanine ATP yield
    'RCR11388'
    % tryptophan ATP yield
    'RCR10534'
    'RCR10535'
    'RCR10929'
    'RCR11389'
    'RCR11390'
    'RCR11393'
    'RCR11396'
    'RCR11397'
    'RCR20259'
	% Complete oxidation of a number of fatty acids
    'RCR20304'
	% thioredoxin metabolism
    'RCR11505'
    'RCR11504'
    % NAD+ de novo synthesis from tryptophan
    'RCR11353'
    'RCR14666'
    'RCR41027'
    'RCR41476'
 
    % Added for general model functionality
	% fix carbon source utilization/ATP yields - H2O transporters, H2O2
    % transporter
    'RCR40013'
    'RCR20015'
    'RCR20257'
	% protein synthesis and demand
    'RCR90143'
    'RCR90134'
    % add glutamate exchange
    'RCR30005'
    };

% deleted contains rxn_ids for reactions removed from the model
deleted = [];
for i = 1:length(hsa_cobra_load.rxns)
    temp = strcmp(hsa_cobra_load.rxns{i}, heart_model_draft.rxns());
    if max(temp) == 0
        deleted(end+1,1) = i;
    end
end

deleted_update = deleted;
for i = 1:length(add_rxns)
    [trash rxn_ID] = max(strcmp(add_rxns{i}, hsa_cobra_load.rxns()));
    location = find(deleted_update == rxn_ID);
    deleted_update(location) = [];
end

% Re-create heart model based off new deleted reaction list
deleted_rxns = {};
for k = 1:length(deleted_update)
    deleted_rxns{end+1} = hsa_cobra_load.rxns{deleted_update(k)};
end

heart_model_curation = removeRxns(hsa_cobra_load, deleted_rxns);

% prevent gluconeogenesis by restricting release of glucose
temp = findRxnIDs(heart_model_curation,'RCR40464');
heart_model_curation.ub(temp) = 0;

% Remove reactions for catalase activity
remove = {'RCR10165' 'RCR10607' 'RCR11007'};
heart_model_curation = removeRxns(heart_model_curation, remove);

% Removing reactions to correct false positive tasks
remove = {
    % acetoacetate production
    'RCR20138'
    % B-hydroxybutanoate production
    'RCR20267'  
	% bile acid synthesis (taurocholate and taucochendeoxycholate)
    'RCR12926'
	% glutamine to citrulline conversion
    'RCR10504'
    % CDCA synthesis
    'RCR40170'
    % CA synthesis
    'RCR10217'
    'RCR13611'
	% gluconeogenesis
    'RCR14184'
    'RCR14208'
    };

heart_model_curation = removeRxns(heart_model_curation, remove);

%% Re-check metabolic tasks after updates
% Convert from original format to COBRA format
inputFile = ['data/AllTasks_CardiomyocyteSpecific.xlsx'];
[FINAL] = generateCobraTaskList(inputFile, hsa_cobra_load);
xlswrite('data/AllTasks_CardiomyocyteSpecific_COBRA.xlsx', FINAL, 'TASKS')
inputFile = ['data/AllTasks_CardiomyocyteSpecific_COBRA.xlsx'];

curated_tasks = checkMetabolicTasks_BVD(heart_model_curation,inputFile);

%% Write model to SBML format
% Change the objective function to protein demand
heart_model_curation = changeObjective(heart_model_curation, 'RCR90143');
outmodel = writeCbModel(heart_model_curation, 'format','sbml');

save('data/HeartModel.mat', 'heart_model_curation')

%% Identify reactions carrying flux for each of the metabolic tasks
% Generate a taskStructure for determining the min set of reactions
inputFile = ['data/AllTasks_CardiomyocyteSpecific_COBRA.xlsx'];
taskStructure=generateTaskStructure_BVD(inputFile);

% calculate the reactions necessary for each task
removeNoGPR = 'false';
minRxnList = generateMinRxnList(heart_model_curation, taskStructure, removeNoGPR);

% Step through the minRxnList to identify reactions covered under different
% tasks
cardio = [1:4 221:308];
iHsa = [5:222];

cardio_rxns = [];
for k = 1:length(cardio)
    rxns = findRxnIDs(heart_model_curation, minRxnList(cardio(k)).rxns);
    cardio_rxns = [cardio_rxns; rxns];
end

% Remove repeats in the cardio_rxns list
k = 1;
while k < length(cardio_rxns)
    to_remove = cardio_rxns(k) == cardio_rxns;
    to_remove(k) = 0;
    cardio_rxns(to_remove) = [];
    k = k+1;
end

% Identify carrying reactions in the iHsa task list
iHsa_rxns = [];
for k = 1:length(iHsa)
    rxns = findRxnIDs(heart_model_curation, minRxnList(iHsa(k)).rxns);
    iHsa_rxns = [iHsa_rxns; rxns];
end

% Remove repeats in the iHsa_rxns list
k = 1;
while k < length(iHsa_rxns)
    to_remove = iHsa_rxns(k) == iHsa_rxns;
    to_remove(k) = 0;
    iHsa_rxns(to_remove) = [];
    k = k+1;
end

% Identify shared reactions between the two
shared_rxns = cardio_rxns;
k = 1;
while k < length(shared_rxns)
    if sum(shared_rxns(k) == iHsa_rxns) == 0
        shared_rxns(k) = [];
    else 
        k = k+1;
    end
end

% For heart model curation,
% shared reactions: 753
% cardio reactions: 874
% iHsa reactions: 1593

%% Identify the average number of reactions covered in each task
num_rxns = [];
for task = 1:length(minRxnList)
    num_rxns(task) = length(minRxnList(task).rxns);
end
%% ATP yields for common carbon sources
clc

% Load in exchange reactions from iHsa model
% Since exchange reactions were changed to [-1000, 0] for building context
% specific models, need to change to [0,0] to make ATP predictions
exchange_model = ncomm_blais_xls2model('data/ncomms14250-s10, iHsa COBRA.xlsx');
lbx = find(exchange_model.lb ~= -1000 & exchange_model.lb ~= 0);
lbx_rxns = findRxnIDs(heart_model_curation, exchange_model.rxns(lbx));
lbx_rxns(lbx_rxns == 0) = [];
heart_model_curation.lb(lbx_rxns) = 0;

model_base = heart_model_curation;

% Change to correct carbon substrate uptake constraints
reactions_in_model = {
    'RCR30315' %string('glucose exchange')
    'RCR90123' %string('palmitate exchange')
    'RCR30024' %string('L-lactate exchange')
    'RCR30046' %string('(R)-3-hydroxybutanoate exchange')
    'RCR30045' %string('acetoacetate exchange')
    'RCR30032' %string('acetate exchange')
    'RCR30013' %string('alanine exchange')
    'RCR30015' %string('arginine exchange')
    'RCR30014' %string('asparagine exchange')
    'RCR30145' %string('aspartate exchange')
    'RCR30095' %string('cysteine exchange')
    'RCR30005' %string('glutamate exchange')
    'RCR30321' %string('glutamine exchange')
    'RCR30116' %string('glycine exchange')
    'RCR30527' %string('histidine exchange')
    'RCR30528' %string('isoleucine exchange')
    'RCR30529' %string('leucine exchange')
    'RCR30530' %string('lysine exchange')
    'RCR30531' %string('methionine exchange')
    'RCR90202' %string('oletate exchange')
    'RCR30184' %string('phenylalanine exchange')
    'RCR30117' %string('proline exchange')
    'RCR30042' %string('serine exchange')
    'RCR30317' %string('threonine exchange')
    'RCR30318' %string('tryptophan exchange')
    'RCR30031' %string('tyrosine exchange')
    'RCR30185' %string('valine exchange')
    };

f = findRxnIDs(model_base, reactions_in_model);
ATP_production = zeros(length(reactions_in_model),1);
for i = 1:length(reactions_in_model)
    model = model_base;
    % Set constraint inputs
    model.ub(f(i)) = 0;
    model.lb(f(i)) = -1;
    % Add ATP demand as a reaction
    model = changeObjective(model, 'RCR11017');
    sol = optimizeCbModel(model);
    ATP_production(i,1) = sol.f;
    RxnID = findRxnIDs(model, 'RCR30011');
    ATP_production(i,2) = abs(sol.full(RxnID));
    
    % Run pFBA to find the number of reactions used
    model.lb(findRxnIDs(model, 'RCR11017')) = 0.99*sol.f;
    [MinimizedFlux modelIrrev]= pFBA_edit(model,0);
    ATP_production(i,3) = sum(MinimizedFlux.full ~= 0);
end
