%Author: Kanadech Jirapongtanavech, Z5176970

%Program: Solution for AAS, T1 2019, Project4.part A....

close all; clear; clc;
global N_G N_C P_m P_c locNumber
N_G = 2000; %Number of generations
N_C = 300; %Population // Number of Chromosomes
P_m = 0.80; %Probablity of mutation
P_c = 0.50; %Probablity of crossover
locNumber = 20; %number of locations
map = zeros(3,locNumber);
history = zeros(2,N_G);
for i = 1:length(map)
   %Generate locNumber landmarks with coordinates randomly distributed
   %in the interval of (-40,40)
   map(:,i) = [i;-40 + (40+40)*rand([2,1])];
end
startPoint = -40 + (40+40)*rand([2,1]);
endPoint = -40 + (40+40)*rand([2,1]);
%disp(map);
handle = plot(0,0);
zoom on;
text(startPoint(1),startPoint(2),'Start');
text(endPoint(1),endPoint(2),'End');
for i = 1:locNumber
    text(map(2,i),map(3,i),int2str(map(1,i)));
end
population = zeros(N_C,locNumber);
for i = 1:N_C
   noDup = false;
   while(~noDup) 
       population(i,:) = randperm(locNumber);
       c = unique(population(1:i,:),'rows');
       [m,~] = size(c);
       if (m == i)
           noDup = true;
       else
           disp('Found duplicate sequence');
       end
   end
end

for i = 1:N_G
    % Obtain Fitness score of every chromosome
    disp(i)
    fitness = calcFitness(startPoint,endPoint,map,population);
    history(:,i) = [i;min(fitness(:,1))];
    [~,eliteIndex] = min(fitness(:,1));
    elite = population(eliteIndex,:);
    tempMap = map(:,elite);
    x = [startPoint(1),tempMap(2,:),endPoint(1)];
    y = [startPoint(2),tempMap(3,:),endPoint(2)];
    set(handle,'xdata',x,'ydata',y,'marker','*');
    pause(0.001);
    % disp(fitness);
    meanWeight = mean(fitness(:,2));
    fitnessScore = fitness(:,2)/meanWeight;
    probArray = fitnessScore/sum(fitnessScore);
    %sum(probArray)
    %% Cross Over
    nextGen = [];
    for j = 1:N_C
        if rand()>P_c
           crossOverPoint = randi([2,10],1,1);
           parent = randsample(length(probArray),2,true,probArray);
           parentA = population(parent(1),:);
           parentB = population(parent(2),:);
           child1 = [parentB(1:crossOverPoint),parentA(crossOverPoint+1:locNumber)];
           child2 = [parentA(1:crossOverPoint),parentB(crossOverPoint+1:locNumber)];
           nextGen = [nextGen;child1];
           nextGen = [nextGen;child2];
        end
    end
    population = nextGen;
    %% mutation
    populationSize = size(population,1);
    for j = 1:populationSize
        if rand()>(P_m)
            index1 = 0;
            index2 = 0;
            while index1 == index2
               index1 = randi([1,locNumber],1,1);
               index2 = randi([1,locNumber],1,1);
            end
            temp = population(j,index1);
            population(j,index1) = population(j,index2);
            population(j,index2) = temp;
        end
        % remove duplicate
        diff = setdiff(map(1,:),population(j,:));
        if diff
            [~,ia,~] = unique(population(j,:),'first');     %unique indices
            dupInd = setdiff(map(1,:),ia);                  %dup indices
            for k=1:length(dupInd)
                population(j,dupInd(k)) = diff(k);
            end
        end
        diff = setdiff(map(1,:),population(j,:));
        assert(isempty(diff));
    end
    %% Elitism
    population(randsample(populationSize,1),:) = elite;
end
figure()
plot(history(1,:),history(2,:));
xlabel('Generation');
ylabel('Distance');

function fitness = calcFitness(startP,endP,map,population)
    global locNumber
    fitness = zeros(size(population,1),2);
    for i=1:size(population,1)
        individual = population(i,:);
        distance = 0;
        distance = distance + norm(startP - map(2:3,individual(1))); %distance fromt start to first element
        for j=1:locNumber - 1
            distance = distance + norm(map(2:3,individual(j)) - map(2:3,individual(j+1)));
        end
        distance = distance + norm(endP - map(2:3,individual(locNumber))); %distance fromt last element to end
        fitness(i,:) = [distance,1/distance];
    end
end