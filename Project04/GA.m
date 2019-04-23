global N_G N_C P_m P_c
N_G = 5000; %Number of generations
N_C = 300; %Population // Number of Chromosomes
P_m = 0.5; %Probablity of mutation
P_c = 0.5; %Probablity of crossover
map = zeros(3,11);
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
    
    fitness = calcFitness(startPoint,endPoint,map,population);
    % disp(fitness);
    meanWeight = mean(fitness(:,2));
    fitnessScore = fitness(:,2)/meanWeight;
    % find fitnessScore > 1
end

function fitness = calcFitness(startP,endP,map,population)
    global N_C
    fitness = zeros(N_C,2);
    for i=1:N_C
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