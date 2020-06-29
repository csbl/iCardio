%% Comparing different integration algorithms using the same set of HPA data
clear all
close all
clc

% Load the original model and make the updates to the model
initCobraToolbox(false)
changeCobraSolver('gurobi6', 'all');

% Loading a COBRA file
hsa_cobra_load = xls2model('data\ncomms14250-s10, iHsa COBRA.xlsx');
% Change objective function for ATP use
hsa_cobra_load = changeObjective(hsa_cobra_load, 'RCR11017');

% Make changes to the original reconstruction based off new metabolic tasks
% 1: Change the DHAP reaction to be reveresed
f = findRxnIDs(hsa_cobra_load,'RCR21050');
hsa_cobra_load.lb(f) = -1000;
hsa_cobra_load.ub(f) = 0;

% 2: Remove superoxide movement from the model (currently just exported out)
f = findRxnIDs(hsa_cobra_load,'RCR40428');
hsa_cobra_load.lb(f) = 0;
hsa_cobra_load.ub(f) = 0;

% 3: All catalase reactions have bounds of [0,0] rchange to [0,1000]
rxns = findRxnIDs(hsa_cobra_load, {'RCR10165' 'RCR10607' 'RCR11007' 'RCR11029' 'RCR14124'});
hsa_cobra_load.lb(rxns) = 0;
hsa_cobra_load.ub(rxns) = 1000;

% 4: Change the Complex I reaction to include superoxide production
hsa_cobra_load = addReaction(hsa_cobra_load, 'RCR21048', {'m03103[m]','m02553[m]','m02039[m]', 'm02630[m]', 'm03102[m]','m02552[m]','m02039[c]','m02631[m]'}, [-1 -1 -5 -0.001 1 1 4 0.001], false);

% 5: Change the number of protons moved to generate one ATP
hsa_cobra_load = addReaction(hsa_cobra_load, 'RCR20085', {'m02751[m]','m01285[m]','m02039[c]', 'm01371[m]', 'm02039[m]','m02040[m]'}, [-1 -1 -2.7 1 2.7 1], true);

% 6: Change all exchange reaction bounds that are not 1000, 0, or -1000
lbx = find(hsa_cobra_load.lb ~= -1000 & hsa_cobra_load.lb ~= 0);
hsa_cobra_load.lb(lbx) = -1000;
ubx = find(hsa_cobra_load.ub ~= 1000 & hsa_cobra_load.ub ~= 0);
hsa_cobra_load.ub(ubx) = 1000;

%% Check tasks with the updated models
% Load tasks from original hepatocyte model
% Convert from original format to COBRA format
inputFile = ['data\ncomms14250-s5, metabolic tasks.xls'];
[FINAL] = generateCobraTaskList(inputFile, hsa_cobra_load);
xlswrite('data\RATCONTasks_COBRA.xlsx', FINAL, 'TASKS')
inputFile_RATCON = ['RATCONTasks_COBRA.xlsx'];

% Load cardiomyocyte tasks
inputFile = ['data\CardiomyocyteTasks.xlsx'];
[FINAL] = generateCobraTaskList(inputFile, hsa_cobra_load);
xlswrite('data\CardiomyocyteTasks_COBRA.xlsx', FINAL, 'TASKS')
inputFile_cardio = ['data\CardiomyocyteTasks_COBRA.xlsx'];

% Check that tasks from the original model are passing correctly 
model = hsa_cobra_load;
RATCON = checkMetabolicTasks_BVD(model,inputFile_RATCON);
cardio = checkMetabolicTasks_BVD(model,inputFile_cardio);

% These two variable interfere with CORDA implementation
hsa_cobra_load = rmfield(hsa_cobra_load, 'osenseStr');
hsa_cobra_load = rmfield(hsa_cobra_load, 'csense');

%% Load cardiomyocyte-specific HPA data
heart_data = readtable('data\20190625 -- MATLAB_integrate_HPA.csv');
              
% table contains conservative and liberal calls for proteins
% conservative: high: 2, medium: 1, low/not detected/NA: -1
% liberal: high: 3, medium: 2, low: 1, not detected/NA: -1

express_data_heart.gene = heart_data.GENE_NAME;
express_data_heart.value = heart_data.conservative;

% Map to reactions in the model
[expressionRxns parsedGPR] = mapExpressionToReactions(hsa_cobra_load, express_data_heart);

%% Running algorithms using the createTissueSpecificModel function
% GIMME
clear options
model = hsa_cobra_load;
options.solver = 'GIMME';
%   for GIMME
%       options.expressionRxns       reaction expression, expression data corresponding to model.rxns.
%                                    Note : If no gene-expression data are
%                                    available for the reactions, set the
%                                    value to -1
%       options.threshold            expression threshold, reactions below this are minimized
%       options.obj_frac*            minimum fraction of the model objective function
%                                    (default - 0.9)

% Since GIMME wants -1 for no presence reactions, change conservative
% mappings to 2, 1, and 0
express_data.gene = heart_data.GENE_NAME;
express_data.value = heart_data.conservative;
express_data.value(express_data.value == -1) = 0;

% Map to reactions in the model
[expressionRxns_GIMME parsedGPR] = mapExpressionToReactions(hsa_cobra_load, express_data);

options.expressionRxns = expressionRxns_GIMME;
options.threshold = 0.9;
options.obj_frac = 0.1;

[tissueModel_GIMME] = createTissueSpecificModel(model, options);

%% iMAT or Shlomi et al algorithm
clear options
%    for iMAT
%       options.expressionRxns       reaction expression, expression data corresponding to model.rxns.
%                                    Note : If no gene-expression data are
%                                    available for the reactions, set the value to -1
%       options.threshold_lb         lower bound of expression threshold, reactions with
%                                    expression below this value are "non-expressed"
%       options.threshold_ub         upper bound of expression threshold, reactions with
%                                    expression above this value are
%                                    "expressed"
%       options.tol*                 minimum flux threshold for "expressed" reactions
%                                    (default 1e-8)
%       options.core*                cell with reaction names (strings) that are manually put in
%                                    the high confidence set (default - no core reactions)
%       options.logfile*             name of the file to save the MILP log (defaut - 'MILPlog')
%       options.runtime*             maximum solve time for the MILP (default - 7200s)
%       options.epsilon*             small value to consider when modeling
%                                    flux (default 1)

% Since iMAT wants -1 for no presence reactions, change conservative
% mappings to 2, 1, and 0
express_data.gene = heart_data.GENE_NAME;
express_data.value = heart_data.conservative;
express_data.value(express_data.value == -1) = 0;

% Map to reactions in the model
[expressionRxns_iMAT parsedGPR] = mapExpressionToReactions(hsa_cobra_load, express_data);

model = hsa_cobra_load;
options.solver = 'iMAT';
options.expressionRxns = expressionRxns_iMAT;
options.threshold_lb = 0;
options.threshold_ub = 0.9;

[tissueModel_iMAT] = createTissueSpecificModel(model, options);

%% Create a flux consistent base model
clear options

% Create globally consistent network:
% checking for consistency of the network
% A - n x 1 boolean vector constaining consistent reactions
epsilon = 1e-6; % smallest flux considered non-zero
printLevel = 2; % summary print level
[A,modelFlipped,V] = fastcc(hsa_cobra_load, epsilon, printLevel);

% A - 5837 reactions/8336 reactions - number of reactions doesn't change
% based on epsilon
remove = setdiff(1:numel(hsa_cobra_load.rxns), A);
rxnRemoveList = hsa_cobra_load.rxns(remove);
base_model_FluxConsistent = removeRxns(hsa_cobra_load, rxnRemoveList);

% epsilon = 1e-8: 5864 reactions
% epsilon = 1e-6: 5850 reactions
% epsilon = 1e-4: 5850 reactions

%% fastCore - using createTissueSpecificModel

% Including only high and medium data, the algorithm doesn't run
% Include high, medium, and low proteins
express_data_heart_FastCore.gene = heart_data.GENE_NAME;
express_data_heart_FastCore.value = heart_data.liberal;

% Map to reactions in the model
[expressionRxns_FluxConsistent parsedGPR] = mapExpressionToReactions(base_model_FluxConsistent, express_data_heart_FastCore);
core = 1:length(base_model_FluxConsistent.rxns);
core = core(expressionRxns_FluxConsistent >= 0);

% FastCore
%       options.core                 indices of reactions in cobra model that are part of the
%                                    core set of reactions
%       options.epsilon*             smallest flux value that is considered
%                                    nonzero (default 1e-4)
%       options.printLevel*          0 = silent, 1 = summary, 2 = debug (default 0)
model = base_model_FluxConsistent;
options.solver = 'fastCore';
options.core = core;
% lower epsilon gives more core reactions without throwing an error
options.epsilon = 1e-8;
options.printLevel = 2;

[tissueModel_fastCore_liberal] = createTissueSpecificModel(model, options);

% Include only high reactions
core = 1:length(base_model_FluxConsistent.rxns);
core = core(expressionRxns_FluxConsistent > 1);

options.core = core;
[tissueModel_fastCore_conservative] = createTissueSpecificModel(model, options);

%% MBA - using createTissueSpecificModel
% Use the FastCore model from above since MBA calls fastcc during algorithm
clear options
%   for MBA
%       options.medium_set           list of reaction names with medium confidence
%       options.high_set             list of reaction names with high confidence
%       options.tol*                 minimum flux threshold for "expressed" reactions
%                                    (default - 1e-8)
model = base_model_FluxConsistent;
options.solver = 'MBA';
options.medium_set = []; % base_model_FastCore.rxns(expressionRxns_FastCore >= 0);
options.high_set = base_model_FluxConsistent.rxns(expressionRxns_FluxConsistent == 2);
options.tol = 1e-8;

[tissueModel_MBA_conservative] = createTissueSpecificModel(model, options);

options.medium_set = base_model_FluxConsistent.rxns(expressionRxns_FluxConsistent == 1); % | expressionRxns_FluxConsistent == 0);
[tissueModel_MBA_liberal] = createTissueSpecificModel(model,options);

%% CORDA 
% Function not in the new COBRA toolbox
%   metTests - metabolic tests to be included in the reconstruction. This
%       argument should be a cell array of strings of size nx2, where n is 
%       the number of metabolic tests to be performed. Column 1 should be 
%       the name of the reaction to be included and the corresponding column 
%       2 should be the reaction to be included. For example, to test for
%       the production of pep and pyruvate, metTests should be equal to
%       {'DM_pep[c]' 'pep[c] -> ';'DM_pyr[c]' 'pyr[c] -> '}. May be left
%       empty.
%   ES - High confidence reactions. Reactions to be included in the model. Cell
%       array of strings.
%   PR - Medium confidence reactions to be included in the model if they do 
%       not depend on too many NP reactions to carry a flux. Cell array of 
%       strings.
%   NP - Negatice confidence reactions not to be included in the model. 
%       These reactions will be included in the tissue model only if they 
%       are necessary for the flux of ES reactions or for the flux of PRtoNP 
%       or more PR reactions. Cell array of strings.

express_data_heart.gene = heart_data.GENE_NAME;
express_data_heart.value = heart_data.conservative;

% Map to reactions in the model
[expressionRxns parsedGPR] = mapExpressionToReactions(hsa_cobra_load, express_data_heart);

% Define medium reactions as those associated with medium-level proteins
% Classifications for reactions: -1, 1, 2
% ES: 2
% PR: 1
% NP: -1
ES = {};
PR = {};
NP = {};
for i = 1:length(expressionRxns)
    temp = hsa_cobra_load.rxns(i);
    if expressionRxns(i) == 2
        ES{end+1,1} = temp{1,1};
    elseif expressionRxns(i) == 1
        PR{end+1,1} = temp{1,1};
    elseif expressionRxns(i) == -1
        NP{end+1,1} = temp{1,1};
    end
end

constraint = 1e-4;
[tissueModel_CORDA, rescue, HCtoMC, HCtoNC, MCtoNC] = CORDA(hsa_cobra_load, {}, ES,PR,NP, 2, constraint);

%% Save models
% save('draft_heart_models.mat', tissueModel_CORDA, tissueModel_fastCore_conservative, tissueModel_fastCore_liberal, tissueModel_GIMME, tissueModel_iMAT, tissueModel_MBA_conservative, tissueModel_MBA_liberal)

%% Compare all the developed models based on task completion
% Load previously run models
load('draft_heart_models.mat')

% Load tasks from original hepatocyte model
% Convert from original format to COBRA format
inputFile = ['data\AllTasks_CardiomyocyteSpecific.xlsx'];
[FINAL] = generateCobraTaskList(inputFile, hsa_cobra_load);
xlswrite('data\AllTasks_CardiomyocyteSpecific_COBRA.xlsx', FINAL, 'TASKS')
inputFile = ['data\AllTasks_CardiomyocyteSpecific_COBRA.xlsx'];

% Check tasks for individual models
model = tissueModel_GIMME;
GIMME_tasks = checkMetabolicTasks_BVD(model,inputFile);

model = tissueModel_iMAT;
iMAT_tasks = checkMetabolicTasks_BVD(model,inputFile);

model = tissueModel_MBA_conservative;
MBA_conservative_tasks = checkMetabolicTasks_BVD(model,inputFile);

model = tissueModel_MBA_liberal;
MBA_liberal_tasks = checkMetabolicTasks_BVD(model,inputFile);

model = tissueModel_fastCore_conservative;
FastCore_conservative_tasks = checkMetabolicTasks_BVD(model,inputFile);

model = tissueModel_fastCore_liberal;
FastCore_liberal_tasks = checkMetabolicTasks_BVD(model,inputFile);

model = tissueModel_CORDA;
CORDA_tasks = checkMetabolicTasks_BVD(model,inputFile);
