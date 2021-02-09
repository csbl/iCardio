function [MinimizedFlux modelIrrev]= pFBA_edit(model,GeneOption)
% This function finds the minimum flux through the network and returns the
% minimized flux and an irreversible model convert model to irrev

    modelIrrev = convertToIrreversible(model);
    % add pseudo-metabolite to measure flux through network

    if nargin==1
        GeneOption=0;
    end
    
    %Add the metabolite
    modelIrrev = addMetabolite(modelIrrev,'fluxMeasure');
    %Then set all the Stoichiometric entries.
    
    if GeneOption==0 % signal that you want to minimize the sum of all gene and non-gene associated fluxes        
        modelIrrev.S(end,:) = ones(size(modelIrrev.S(1,:)));
    elseif GeneOption==1 % signal that you want to minimize the sum of only gene-associated fluxes
        %find all reactions which are gene associated
        Ind=find(sum(modelIrrev.rxnGeneMat,2)>0);
        modelIrrev.S(end,:) = zeros(size(modelIrrev.S(1,:)));
        modelIrrev.S(end,Ind) = 1;
    elseif GeneOption==2 % signal that you want to minimize the sum of only NON gene-associated fluxes
        %find all reactions which are gene associated
        Ind=find(sum(modelIrrev.rxnGeneMat,2)==0);
        modelIrrev.S(end,:) = zeros(size(modelIrrev.S(1,:)));
        modelIrrev.S(end,Ind) = 1;
    end   

    % add a pseudo reaction that measures the flux through the network
    modelIrrev = addReaction(modelIrrev,'netFlux',{'fluxMeasure'},[-1],false,0,inf,0,'','');

    % set the flux measuring demand as the objective
    modelIrrev.c = zeros(length(modelIrrev.rxns),1);
    modelIrrev = changeObjective(modelIrrev, 'netFlux');

    % minimize the flux measuring demand (netFlux)
    MinimizedFlux = optimizeCbModel(modelIrrev,'min');
end