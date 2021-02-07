% Testing additional heart failure datasets with TIDEs
% Additional data sets: GSE76701, GSE26887, GSE71613, GSE133054
% Consider running additional GSE141910 (over 360 samples)
clear all
clc

% initialize the toolbox
initCobraToolbox
changeCobraSolver('ibm_cplex', 'all');

% Load in the heart model
load('HeartModel.mat')

% Generate a taskStructure for determining the min set of reactions
% Convert from original format to COBRA format
inputFile = ['AllTasks_CardiomyocyteSpecific_COBRA.xlsx'];
taskStructure=generateTaskStructure_BVD(inputFile);

% calculate the reactions necessary for each task
removeNoGPR = 'true';
minRxnList = generateMinRxnList(heart_model_curation, taskStructure, removeNoGPR);

% Add on analysis by reaction subsystem
% Find a unique list of subsystems in this model
reactions = {};
for k = 1:length(heart_model_curation.subSystems)
    reactions{k,1} = char(heart_model_curation.subSystems{k});
end
subsystem = unique(reactions);
subsystem(1,:) = [];

for k = 1:length(subsystem)
    minRxnList(end+1).id = strcat('S',num2str(k));
    minRxnList(end).description = subsystem{k};
    minRxnList(end).rxns = heart_model_curation.rxns(strcmp(subsystem{k}, reactions));
end

% Remove reactions from the subsystem list that don't have GPR rules
for k = 1:length(heart_model_curation.rxns)
    noGPR(k) = isempty(heart_model_curation.grRules{k});
end
noGPR = heart_model_curation.rxns(noGPR);

% Remove these reactions from the list of total_carrying_rxns
for task = 1:length(minRxnList)
    common = intersect(noGPR, minRxnList(task).rxns);
    minRxnList(task).rxns = setdiff(minRxnList(task).rxns, common);
end

% Remove tasks that no longer have associated reactions
task = 1;
while task ~= length(minRxnList)
    if isempty(minRxnList(task).rxns)
        minRxnList(task) = [];
    else
        task = task+1;
    end
end

% Specify the number of random iterations of data
model = heart_model_curation;
parpool(4);
num_iterations = 1000;

%% Debugging reading in table problems
opts = detectImportOptions('C:\Users\bvd5nq\Documents\R scripts\Cardiomyocyte Model\GSEdata\GSE5406_DEGs.csv');
opts.Delimiter = {','};
opts.VariableNames = {'EntrezID','failure','failure_logFC'};

%% GSE26887
% Load example data set
humanHF = readtable('C:\Users\bvd5nq\Documents\R scripts\Cardiomyocyte Model\GSEdata\GSE26887_DEGs.csv', opts);
data.gene = humanHF.EntrezID;
data.value = humanHF.failure_logFC;

% Run test expression mapping to make sure EntrezIDs are formatted correctly
[expressionRxns parsedGPR] = mapExpressionToReactions(model, data);

% calculate TIDE scores based on a dataset
model = heart_model_curation;
[GSE26887_failure, GSE26887_failure_random] = calculateTIDEscores(model, minRxnList, data, num_iterations);

%% GSE71613
% Load example data set
humanHF = readtable('C:\Users\bvd5nq\Documents\R scripts\Cardiomyocyte Model\GSEdata\GSE71613_failure_DEGs.csv', opts);
data.gene = humanHF.EntrezID;
data.value = humanHF.failure_logFC;

% Run test expression mapping to make sure EntrezIDs are formatted correctly
[expressionRxns parsedGPR] = mapExpressionToReactions(model, data);

% calculate TIDE scores based on a dataset
model = heart_model_curation;
[GSE71613_failure, GSE71613_failure_random] = calculateTIDEscores(model, minRxnList, data, num_iterations);
 
%% GSE141910 - failure only
% Also read in only the failure data
humanHF = readtable('C:\Users\bvd5nq\Documents\R scripts\Cardiomyocyte Model\GSEdata\GSE141910_failure_DEGs.csv', opts);
data.gene = humanHF.EntrezID;
data.value = humanHF.failure_logFC;

% Run test expression mapping to make sure EntrezIDs are formatted correctly
[expressionRxns parsedGPR] = mapExpressionToReactions(model, data);

% calculate TIDE scores based on a dataset
model = heart_model_curation;
[GSE141910_failure, GSE141910_failure_random] = calculateTIDEscores(model, minRxnList, data, num_iterations);

%% GSE133054 - failure only
% Change options to reflect failure and hypertrophic data
humanHF = readtable('C:\Users\bvd5nq\Documents\R scripts\Cardiomyocyte Model\GSEdata\GSE133054_failure_DEGs.csv', opts);
data.gene = humanHF.EntrezID;
data.value = humanHF.failure_logFC;

% Run test expression mapping to make sure EntrezIDs are formatted correctly
[expressionRxns parsedGPR] = mapExpressionToReactions(model, data);

% calculate TIDE scores based on a dataset
model = heart_model_curation;
[GSE133054_failure, GSE133054_failure_random] = calculateTIDEscores(model, minRxnList, data, num_iterations);

%% Save all the datasets to Excel spreadsheets
% GSE26887_failure
data = GSE26887_failure;
% Save individual variable names to xlsx files in R/data folder
data_save = {};
for k = 1:length(data)
    data_save{k,2} = data(k).description;
    data_save{k,1} = data(k).id;
    data_save{k,4} = data(k).significance;
    data_save{k,3} = data(k).taskScore;
end
xlswrite('C:\Users\bvd5nq\Documents\R scripts\Cardiomyocyte Model\data\TIDEs\GSE26887_allGenes_Failure.xlsx', data_save)

% GSE71613_failure
data = GSE71613_failure;
% Save individual variable names to xlsx files in R/data folder
data_save = {};
for k = 1:length(data)
    data_save{k,2} = data(k).description;
    data_save{k,1} = data(k).id;
    data_save{k,4} = data(k).significance;
    data_save{k,3} = data(k).taskScore;
end
xlswrite('C:\Users\bvd5nq\Documents\R scripts\Cardiomyocyte Model\data\TIDEs\GSE71613_allGenes_Failure.xlsx', data_save)

% GSE141910_failure
data = GSE141910_failure;
% Save individual variable names to xlsx files in R/data folder
data_save = {};
for k = 1:length(data)
    data_save{k,2} = data(k).description;
    data_save{k,1} = data(k).id;
    data_save{k,4} = data(k).significance;
    data_save{k,3} = data(k).taskScore;
end
xlswrite('C:\Users\bvd5nq\Documents\R scripts\Cardiomyocyte Model\data\TIDEs\GSE141910_allGenes_Failure.xlsx', data_save)

% GSE133054_failure
data = GSE133054_failure;
% Save individual variable names to xlsx files in R/data folder
data_save = {};
for k = 1:length(data)
    data_save{k,2} = data(k).description;
    data_save{k,1} = data(k).id;
    data_save{k,4} = data(k).significance;
    data_save{k,3} = data(k).taskScore;
end
xlswrite('C:\Users\bvd5nq\Documents\R scripts\Cardiomyocyte Model\data\TIDEs\GSE133054_allGenes_Failure.xlsx', data_save)

%% Save relevant variables into a file
save('TIDEs_additional.mat','GSE26887_failure','GSE71613_failure','GSE141910_failure','GSE133054_failure')