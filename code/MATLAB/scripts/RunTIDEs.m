% Testing various data sets with TIDEs
clear all
clc

% initialize the toolbox
initCobraToolbox
changeCobraSolver('ibm_cplex', 'all');

% Load in the heart model
load('data/HeartModel.mat')

% Generate a taskStructure for determining the min set of reactions
% Convert from original format to COBRA format
inputFile = ['data/AllTasks_CardiomyocyteSpecific_COBRA.xlsx'];
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

%% Run TIDEs analysis for datasets presented in Figure 4 and 5
% GSE5406
% Load example data set
humanHF = readtable('data/GSE5406_DEGs.csv');
data.gene = humanHF.GENE_NAME;
data.value = humanHF.Ischemic_logFC;

% calculate TIDE scores based on a dataset
model = heart_model_curation;
[GSE5406_Ischemic, GSE5406_Ischemic_random] = calculateTIDEscores(model, minRxnList, data, num_iterations);

% Map idiopathic data to the model
data.value = humanHF.Idiopathic_logFC;
[GSE5406_Idiopathic, GSE5406_Idiopathic_random] = calculateTIDEscores(model, minRxnList, data, num_iterations);

%% GSE57345
humanHF = readtable('data/GSE57345_DEGs.csv');
data.gene = humanHF.GENE_NAME;
data.value = humanHF.Ischemic_logFC;

model = heart_model_curation;
[GSE57345_Ischemic, GSE57345_Ischemic_random] = calculateTIDEscores(model, minRxnList, data, num_iterations);

% Map idiopathic data to the model
data.value = humanHF.Idiopathic_logFC;
[GSE57345_Idiopathic, GSE57345_Idiopathic_random] = calculateTIDEscores(model, minRxnList, data, num_iterations);

%% GSE1869
humanHF = readtable('data/GSE1869_DEGs.csv');
data.gene = humanHF.ENTREZ_GENE_ID;
data.value = humanHF.Ischemic_logFC;

model = heart_model_curation;
[GSE1869_Ischemic, GSE1869_Ischemic_random] = calculateTIDEscores(model, minRxnList, data, num_iterations);

% Map idiopathic data to the model
data.value = humanHF.Idiopathic_logFC;
[GSE1869_Dilated, GSE1869_Dilated_random] = calculateTIDEscores(model, minRxnList, data, num_iterations);


%% Save all the datasets to Excel spreadsheets
% GSE5406_Idiopathic
data = GSE5406_Ischemic;
% Save individual variable names to xlsx files in R/data folder
data_save = {};
for k = 1:length(data)
    data_save{k,2} = data(k).description;
    data_save{k,1} = data(k).id;
    data_save{k,4} = data(k).significance;
    data_save{k,3} = data(k).taskScore;
end
xlswrite('data/GSE5406_allGenes_Ischemic.xlsx', data_save)

% GSE5406_Idiopathic
data = GSE5406_Idiopathic;
% Save individual variable names to xlsx files in R/data folder
data_save = {};
for k = 1:length(data)
    data_save{k,2} = data(k).description;
    data_save{k,1} = data(k).id;
    data_save{k,4} = data(k).significance;
    data_save{k,3} = data(k).taskScore;
end
xlswrite('data/GSE5406_allGenes_Idiopathic.xlsx', data_save)

%% GSE57345
% GSE57345_Idiopathic
data = GSE57345_Ischemic;
% Save individual variable names to xlsx files in R/data folder
data_save = {};
for k = 1:length(data)
    data_save{k,2} = data(k).description;
    data_save{k,1} = data(k).id;
    data_save{k,4} = data(k).significance;
    data_save{k,3} = data(k).taskScore;
end
xlswrite('data/GSE57345_allGenes_Ischemic.xlsx', data_save)

% GSE57345_Idiopathic
data = GSE57345_Idiopathic;
% Save individual variable names to xlsx files in R/data folder
data_save = {};
for k = 1:length(data)
    data_save{k,2} = data(k).description;
    data_save{k,1} = data(k).id;
    data_save{k,4} = data(k).significance;
    data_save{k,3} = data(k).taskScore;
end
xlswrite('data/GSE57345_allGenes_Idiopathic.xlsx', data_save)

%% GSE1869
% GSE1869
data = GSE1869_Dilated;
% Save individual variable names to xlsx files in R/data folder
data_save = {};
for k = 1:length(data)
    data_save{k,2} = data(k).description;
    data_save{k,1} = data(k).id;
    data_save{k,4} = data(k).significance;
    data_save{k,3} = data(k).taskScore;
end
xlswrite('data/GSE1869_allGenes_Dilated.xlsx', data_save)

% GSE1869_Idiopathic
data = GSE1869_Ischemic;
% Save individual variable names to xlsx files in R/data folder
data_save = {};
for k = 1:length(data)
    data_save{k,2} = data(k).description;
    data_save{k,1} = data(k).id;
    data_save{k,4} = data(k).significance;
    data_save{k,3} = data(k).taskScore;
end
xlswrite('data/GSE1869_allGenes_Ischemic.xlsx', data_save)

%% Save some of the random data sets
data = GSE5406_Ischemic_random; 
xlswrite('data/GSE5406_Ischemic_random.xlsx', data);

%% Save the variables to a MATLAB file
save('data/TIDEs_original.mat', 'GSE5406_Ischemic', 'GSE5406_Idiopathic', 'GSE57345_Idiopathic', 'GSE57345_Ischemic', 'GSE1869_Dilated', 'GSE1869_Ischemic')