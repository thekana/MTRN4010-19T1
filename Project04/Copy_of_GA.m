global N_G N_C P_m P_c
N_G = 1000; %Number of generations
N_C = 300; %Population // Number of Chromosomes
P_m = 0.99; %Probablity of mutation
P_c = 0.50; %Probablity of crossover
map = zeros(3,11);
history = zeros(2,N_G);
for i = 1:length(map)
   %Generate 11 landmarks with coordinates randomly distributed
   %in the interval of (-40,40)
   map(:,i) = [i;-40 + (40+40)*rand([2,1])];
end
startPoint = -40 + (40+40)*rand([2,1]);
endPoint = -40 + (40+40)*rand([2,1]);
%disp(map);
population = zeros(N_C,11);
for i = 1:N_C
   noDup = false;
   while(~noDup) 
       population(i,:) = randperm(11);
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
    % disp(fitness);
    meanWeight = mean(fitness(:,2));
    fitnessScore = fitness(:,2)/meanWeight;
    % remove fitnessScore < 1
    parents = population(fitnessScore>=1.0,:);
    % selection completed
    % Shuffle
    parents = parents(randperm(size(parents, 1)), :);
    
    %% Cross Over
    populationSize = size(population,1);
    nextGen = [];
    for j = 1:N_C
        %disp(j);
        index1 = 0;
        index2 = 0;
        if rand()>P_c
           crossOverPoint = randi([2,10],1,1);
           while index1 == index2
               index1 = randi([1,size(parents,1)],1,1);
               index2 = randi([1,size(parents,1)],1,1);
           end
           parent1 = parents(index1,:);
           parent2 = parents(index2,:);
           child1 = [parent2(1:crossOverPoint),parent1(crossOverPoint+1:11)];
           child2 = [parent1(1:crossOverPoint),parent2(crossOverPoint+1:11)];
           nextGen = [nextGen;child1;child2];
        end
%        if size(population,1) >= N_C
%            break;
%        end
    end
    %% mutation
    populationSize = size(nextGen,1);
    for j = 1:populationSize
        if rand()>(P_m)
            index1 = 0;
            index2 = 0;
            while index1 == index2
               index1 = randi([1,11],1,1);
               index2 = randi([1,11],1,1);
            end
            temp = nextGen(j,index1);
            nextGen(j,index1) = nextGen(j,index2);
            nextGen(j,index2) = temp;
        end
        % remove duplicate
        diff = setdiff(map(1,:),nextGen(j,:));
        if diff
            [~,ia,~] = unique(nextGen(j,:),'first');     %unique indices
            dupInd = setdiff(map(1,:),ia);                  %dup indices
            for k=1:length(dupInd)
                nextGen(j,dupInd(k)) = diff(k);
            end
        end
        diff = setdiff(map(1,:),nextGen(j,:));
        assert(isempty(diff));
    end
    %% Elitism
    if size(parents,1) < 10
        nextGen = [nextGen;parents];
    else
        nextGen = [nextGen;parents(1:10,:)];
    end
    population = nextGen;
end

plot(history(1,:),history(2,:));

function fitness = calcFitness(startP,endP,map,population)

    fitness = zeros(size(population,1),2);
    for i=1:size(population,1)
        individual = population(i,:);
        distance = 0;
        distance = distance + norm(startP - map(2:3,individual(1))); %distance fromt start to first element
        for j=1:10
            distance = distance + norm(map(2:3,individual(j)) - map(2:3,individual(j+1)));
        end
        distance = distance + norm(endP - map(2:3,individual(11))); %distance fromt last element to end
        fitness(i,:) = [distance,1/distance];
    end
end