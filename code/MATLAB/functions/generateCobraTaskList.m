function [FINAL] = generateCobraTaskList(inputFile, model)
% generateCobraTaskList
% Generate a new Excel file that translates between the RAVEN and COBRA
% formats for tasks. Uses a metabolic model (model) to translate between
% metabolite names and IDs. 
% 
% Adapted from the generateTaskStructure function in the COBRA toolbox. 
%
% INPUTS:
%    inputFile:          Original task list in the RAVEN format
%    model:              a COBRA model
%
% OUTPUTS:
%    FINAL:              a structure for the translated task list, can be
%                        saved to a new Excel file
%
% .. Authors:
%    - Originally written for RAVEN toolbox by Rasmus Agren, 2013-08-01
%    - Adapted for cobratoolbox and modified to rely only on flux constraints by Richelle Anne, 2017-04-18
%    - Adapted for translating between RAVEN and COBRA namespace by Bonnie
%    Dougherty, 2020-04-08


% Read in the task list: 
[crap,crap,raw]=xlsread(inputFile,'TASKS');

% Remove the first column and first row
raw(:,1) = [];
%raw(1,:) = [];

HEADERS = raw(1,:);
columns_COBRA={'ID';'DESCRIPTION';'IN';'IN LB';'IN UB';'OUT';'OUT LB';'OUT UB';'SHOULD FAIL';'COMP';'SYSTEM';'SUBSYSTEM';'EQU';'EQU LB';'EQU UB'};

FINAL = {};

for i = 1:length(columns_COBRA)
    current_header = columns_COBRA(i);

    [I colI]=ismember(current_header, HEADERS);
    
    if sum(colI) == 0
        % the header isn't in the original file
        header_present(i) = 0;
    else
        header_present(i) = 1;
        FINAL(:,i) = raw(:,colI);
    end
end

% Convert the IN catgeory to metabolite IDs rather than metabolite names
for i=2:length(FINAL(:,3))
    % Parse IN
    % Check to see if NaN, if NaN skip
    IN_string = [];
    if ~isnan(FINAL{i,3})
        inputs = regexp(FINAL(i,3),';','split');
        inputs = inputs{1,1};
        for k = 1:length(inputs)
            location = find(strcmp(inputs{1,k}(1:end-3), model.metNames));
            metID = model.mets(location(1));
            metID = metID{1,1}(1:end-3);
            metID = strcat(metID, inputs{1,k}(end-2:end));
            if length(IN_string) == 0
                IN_string = metID;
            else
                IN_string = strcat(IN_string,';', metID);
            end
        end
        FINAL{i,3} = IN_string;
        IN_string = [];
    end
    
    % Parse OUT
    % Check to see if NaN, if NaN skip
    OUT_string = [];
    if ~isnan(FINAL{i,6})
        if strcmp(FINAL{i,6} ,'ALLMETS')
            OUT_string = 'ALLMETS';
        else
            inputs = regexp(FINAL(i,6),';','split');
            inputs = inputs{1,1};
            for k = 1:length(inputs)
                location = find(strcmp(inputs{1,k}(1:end-3), model.metNames));
                metID = model.mets(location(1));
                metID = metID{1,1}(1:end-3);
                metID = strcat(metID, inputs{1,k}(end-2:end));
                if length(OUT_string) == 0
                    OUT_string = metID;
                else
                    OUT_string = strcat(OUT_string,';', metID);
                end
            end
        end
        FINAL{i,6} = OUT_string;
        OUT_string = [];
    end
    
    % Parse the equations
    if ~isnan(FINAL{i,13})
        equation = FINAL{i,13};
        new_equation{1,1} = '';
        % Figure out the number of reactants in the reaction
        num_reactants = strsplit(equation, ' <=> ');
        if (length(num_reactants) ~= 2) 
            num_reactants = strsplit(num_reactants{1,1}, ' => ');
        end
        num_reactants = strsplit(num_reactants{1,1}, ' + ');
        num_reactants = length(num_reactants);
        equation = strsplit(equation, {' + ', ' <=> ', ' => '});

        for k = 1:length(equation)
            if k == num_reactants+1
                % insert reaction sign
                new_equation{1,1}{1,1} = new_equation{1,1}{1,1}(1:end-3);
                new_equation{1,1} = (strcat(new_equation{1,1}, {' -> '}));
            end

            % Find the corresponding metabolite
            location = find(strcmp(equation{1,k}(1:end-3), model.metNames));
            metID = model.mets(location(1));
            metID = metID{1,1}(1:end-3);
            
            metID = strcat(metID, equation{1,k}(end-2:end));
            new_equation{1,1} = (strcat(new_equation{1,1}, metID, {' + '}));
        end
        new_equation{1,1}{1,1} = new_equation{1,1}{1,1}(1:end-3);
        
        % Have the final equation
        % Write equation
        FINAL{i,13} = new_equation{1,1}{1,1};
    end
end

% Create the SHOULD FAIL columns
for i = 2:length(FINAL(:,9))
    if isnan(FINAL{i,9})
        FINAL{i,9} = 0;
    else
        FINAL{i,9} = 1;
    end
end

% Create the COMP column
FINAL{1,10} = 'COMP';
% for i = 2:length(FINAL(:,10))
%     if isnan(FINAL{i,1})
%     else
%         FINAL{i,10} = 'c';
%     end
% end

FINAL{1,11} = 'SYSTEM';
FINAL{1,12} = 'SUBSYSTEM';

