%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
                                        % BULK %
% % %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%  
% Model autorship: Andreea Stan, TU Delft
%% Clear workspace and memory
% Press Ctrl+enter to evaluate it
clear; clc;
%% Define system parameters
%pb=5e-2; %probability of parasite emergence per cycle
%dil=10;  %dilution
vol=20e-8; %volume in liters (=20 ul)
NA=6.023e23; %avogadro
%tim= 16; %integration time(h)

probabilities = logspace(-4,1,1000.);
probabilities(probabilities>1)= 0.9;

dilutions_r=[];
times_r=[];
prob_r=[];
orip2p3_r=[];

dilutions_g=[];
times_g=[];
prob_g=[];
orip2p3_g=[];

figure

%%Monte Carlo
for s=1:1000 %number of simulations
    disp(s)

    save=[]; %matrix to save concentrations after each cycle

    dna=0.05;
    dnapdna= 0;
    dnap= 0;
    tp= 0;
    para=0;
    dnappara= 0;
    res=1000; 
    initialConditions = [dna dnapdna dnap tp para dnappara res];

    dil = randi([2 10],1);
    tim = randi([2 10],1);
    pb = probabilities(s);

    %% Create mex model Define initial concentrations (nM)

    model=IQMmodel('ODE_220616_2.txt');
    time=[0:0.1:tim];
    paramvec = [];
    options=[];
    options.reltol=10^-25;
    IQMmakeMEXmodel(model,'mexmodel');

    %% Simulate experiment
    for i=1:10 %number of rounds

        %Generate parasites
        chance=rand(1,1000); %in liposomes we have 1000 (=comp) chances of getting a parasite.
            for j=1:length(chance)
            if chance(j)>pb*(initialConditions(1)+initialConditions(2))
                initialConditions(5)= initialConditions(5)+0;
            else
                initialConditions(5)=initialConditions(5)+(1/(NA*vol*1e9)); %N --> C
            end
            end
            %% 
            %Run ODE
            initialConditions(initialConditions<1e-25)=0; 
            z=IQMPsimulate('mexmodel',[0:0.1:tim],initialConditions,paramvec,options);
            initialConditions = z.statevalues(end,:);
            save=[save,z.variablevalues(1,:)',z.variablevalues(end,:)'];
          
            %dilution + feed
            initialConditions =initialConditions/dil;

            initialConditions(3)=0; %Unless bound to DNA, DNAP inactivates after 1 cycle
            initialConditions(4)=0; %Unless bound to DNA, the TP inactivates after 1 cycle
            initialConditions(7)=initialConditions(7)+(res-res/dil);
    end
if save(1, 18) <= save(1,20) & save(1, 16) <= save(1,18)
    plot(save(1,:), 'g');
    orip2p3_g = [orip2p3_g; save(1,:)];
    dilutions_g=[dilutions_g, dil];
    times_g=[times_g, tim];
    prob_g=[prob_g, pb];
else 
    plot(save(1,:), 'r');
    orip2p3_r = [orip2p3_r; save(1,:)];
    dilutions_r=[dilutions_r, dil];
    times_r=[times_r, tim];
    prob_r=[prob_r, pb];
hold on;
end
end
%% Visualization

%% 3D SCATTER
bulks = xlsread("Bulk_data2.xlsx",'Parameter space suc');
prob_g = bulks(:,1);
times_g = bulks(:,2);
dilutions_g = bulks(:,3);
bulkf = xlsread("Bulk_data2.xlsx",'Parameter space fai')
prob_r = bulkf(:,1);
times_r = bulkf(:,2);
dilutions_r = bulkf(:,3);

figure
scatter3(prob_r, times_r, dilutions_r, 'r', 'filled')
hold on;
scatter3(prob_g, times_g, dilutions_g, 'g', 'filled')
xlabel('probability')
ylabel('time')
zlabel('dilution fold')
title('Parameter variation space')
hold off;

%% SIMULATIONS

orip2p3_g = readtable("Bulk_data2.xlsx",'Sheet','Sheet2');
orip2p3_g = table2array(orip2p3_g);
figure
for r=1:height(orip2p3_g)
    plot(orip2p3_g(r,:), 'g')
    hold on;
end
orip2p3_r = readtable("Bulk_data2.xlsx",'Sheet','Sheet4')
orip2p3_r = table2array(orip2p3_r);
for r=1:height(orip2p3_r)
    plot(orip2p3_r(r,:), 'r')
    hold on;
end
xlabel('round') 
ylabel('average concetration (nM)') 
title('Bulk model Monte Carlo simulations')
hold off;
