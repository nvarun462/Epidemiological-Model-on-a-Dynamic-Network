%This code creates a configuration network in which nodes(humans) with
%varying levels of contacts(degree) are put under SIR dynamics.

%CLEARING SPECIFIC VARIABLES FROM PREVIOUS RUNS
clear D_avg E2E E2I I2I E_connections Rand c connections e edge index indexedge;
clear infected infectedE infectedI infecE infecI;
clear recovered recoveredE recoveredI recE recI;
clear susceptible susceptibleE susceptibleI susE susI;
clear adjList distr distrE distrI;
clear histogramE histogramI;

max_days=200; %days to run epidemic

N = 10000; % Number of Nodes
A = 5000; x = 10; %A Extrovert Nodes each initially with degre x
B = 5000; y = 5; %B Introvert Nodes each initially with degree y

if ((A+B)~= N) %Confirms Correct Node Count
    disp 'Error: All types of Nodes should add to N';
    return;
end

S = A*x + B*y; %Total number of Initial Stubs(two stubs connect to make a contact pair)
if (mod(S,2)==1)
    disp 'Error: Total Stubs must be even';
    return;
end

%Creates an Matrix to keep track of which edge belongs to which node.
%The designation of edge E can be found with e(E,:)
k = 1;
for i=1:A
    for j=1:x
        q = strlength(['N',num2str(i),'E',num2str(j)]);
        e(k,1:q)= ['N',num2str(i),'E',num2str(j)];
        k= k+1;
    end
end
for i=(A+1):(B+A)
    for j=1:y
        q = strlength(['N',num2str(i),'E',num2str(j)]);
        e(k,1:q)= ['N',num2str(i),'E',num2str(j)];
        k= k+1;
    end
end

%Randomly Pairing Edges Together
Rand = randperm(S);  %Generates a random permutation of all Stubs, with every 2 representing a pair
c(1:2,1:(S/2)) = 0;  %Matrix of all Connections between 2 Stubs made by segmenting Rand into pieces of length 2
for i = 1:(S/2)  %Populates c matrix
    c(1,i) = Rand((2*i)-1);
    c(2,i) = Rand(2*i);
end

edge(1:2, 1:(S/2)) = 0; %matrix of all connections
for i = 1:S
    edge(i) = getNode(e(c(i),:)); %see getNode function below
end

%START OF EPIDEMICS CODE - - - - - - - - - - - - - - - - - - - - - - - - -

%Definitions
dmax=25; %upper limit to one node's degree
dmin=5;  %lower limit to one node's degree
transmission=0.01;
alpha = 0.00875;
duration=15;
beta=0.2;
d = x; %initial D value = (A*x)/A
E2E = 0; % number of extrovert to extrovert connections only -- will be determined below
E2I = 0; % number of extrovert to introvert connections only -- will be determined below
I2I = 0; % number of introvert to introvert connections only -- will be determined below
ULmt = round(dmax*2.5);

%Creating an adjacency list to monitor all connections
adjList(1:N,1:(ULmt))=0;
for i=1:size(edge,2)
    n1 = edge(1,i);
    n2 = edge(2,i);
    temp1 = 1;
    while (adjList(n1,temp1) ~= 0)
        temp1 = temp1+1;
    end
    adjList(n1,temp1) = n2;
    temp1 = 1;
    while (adjList(n2,temp1) ~= 0)
        temp1 = temp1+1;
    end
    adjList(n2,temp1) = n1;
end

for i=1:size(edge,2) % Find Initial Number of extrovert only connections
    if((edge(1,i)<= A )&& (edge(2,i)<= A))
        E2E = E2E+1;
    end
end

for i=1:size(edge,2) % Find Initial Number of extrovert to introvert connections
    if((edge(1,i)<= A )&& (edge(2,i) > A))
        E2I = E2I+1;
    end
end

for i=1:size(edge,2) % Find Initial Number of extrovert to introvert connections
    if((edge(1,i) > A )&& (edge(2,i) <= A))
        E2I = E2I+1;
    end
end

for i=1:size(edge,2) % Find Initial Number of introvert only connections
    if((edge(1,i)> A )&& (edge(2,i) > A))
        I2I = I2I+1;
    end
end

%Pre-allocate up all tracking matrices
susceptible(1:max_days) = 0;
infected(1:max_days) = 0;
recovered(1:max_days) = 0;
susceptibleI(1:max_days) = 0;
infectedI(1:max_days) = 0;
recoveredI(1:max_days) = 0;
susceptibleE(1:max_days) = 0;
infectedE(1:max_days) = 0;
recoveredE(1:max_days) = 0;
E_connections(1:max_days) = 0; %E_to_E_connectionCount
D_avg(1:max_days) = 0; %D_avg as a function of days




%1 - susceptible, 2 - infected, 3 - recovered
%initial conditions
nodes(1:N) = 1;

% NO LONGER USED - will infect one each of first 3 days
% for i=1:1
%     index = floor(N*rand)+1; % randomly infect a single node
%     nodes(index) = 2;
% end


% start the dynamics
internalclock(1:N) = 0; % represents the internal clock of infected

%DAILY LOOP
inQuarantine = 0; quarTimer = 0;
for days=1:max_days
    
    % Infect one node each of first 3 days
    if (days <=3)
        index = floor(N*rand)+1; % randomly infect a single node
        nodes(index) = 2;
    end
    
    for i=1:size(edge,2) %size(edge,2) refers to current total no. of connections (this will change every day)
        indexedge = floor((size(edge,2))*rand)+1; % randomly choose an edge
        n1 = edge(1,indexedge);
        n2 = edge(2,indexedge);
        %determine if the virus is passed along according to transmission
        if ((nodes(n1)==1)&&(nodes(n2)==2)&&(rand < transmission))
            nodes(n1) = 2;
            internalclock(n1) = 1;
        end
        if ((nodes(n2)==1)&&(nodes(n1)==2)&&(rand < transmission))
            nodes(n2) = 2;
            internalclock(n2) = 1;
        end
    end
    
    % updating the recovery
    for i=1:N
        if (nodes(i)==2)
            internalclock(i) = internalclock(i) + 1;
            if (internalclock(i)>duration)
                nodes(i) = 3;
            end
        end
    end
    
    %Updating all tracking matrices
    susE = 0; infecE = 0; recE = 0; %SIR for Extroverts
    susI = 0; infecI = 0; recI = 0; %SIR for Introverts
    for i=1:A
        if (nodes(i)==1) susE = susE + 1;
        elseif (nodes(i)==2) infecE = infecE + 1;
        else             recE = recE +1;
        end
    end
    for i=(A+1):N
        if (nodes(i)==1) susI = susI + 1;
        elseif (nodes(i)==2) infecI = infecI + 1;
        else             recI = recI +1;
        end
    end
    %Totals
    susceptible(days) = susE+susI;
    infected(days) = infecE+infecI;
    recovered(days) = recE+recI;
    susceptibleI(days) = susI;
    infectedI(days) = infecI;
    recoveredI(days) = recI;
    susceptibleE(days) = susE;
    infectedE(days) = infecE;
    recoveredE(days) = recE;
    
    %MODIFYING THE NETWORK---------------------------------
    %Compute current average degree of extroverts
    d = ((E2E * 2)+E2I)/A;
    D_avg(days) = d;
    
    if (days == 50) inQuarantine = 1; end %Triggers Quarantine on Day 50 (similar to UMich)
    
    if (inQuarantine == 1) 
        disp 'in Quarantine'
        quarTimer = quarTimer + 1;
        if (quarTimer == 1) %first day of quarantine removes every extroverts connection to other extroverts
            E2E = 0;
            disp 'start of Quarantine'
            remover(1)=0;
            for iii = 1:size(edge,2)
                rem1 = edge(1,iii);
                rem2 = edge(2,iii);
                if ((rem1 <= A) && (rem2 <= A))
                    if (remover(1)==0)
                        remover(1) = iii;
                    else
                        remover = [remover iii];
                    end
                    temp1 = 1;
                    while (adjList(rem1,temp1) ~= rem2)
                        temp1 = temp1+1;
                    end
                    adjList(rem1,temp1) = 0;
                    for kk=temp1:ULmt-1
                        adjList(rem1,kk) = adjList(rem1,kk+1);
                    end
                    adjList(rem1,ULmt) = 0;
                    temp1=1;
                    while (adjList(rem2,temp1) ~= rem1)
                        temp1 = temp1+1;
                    end
                    adjList(rem2,temp1) = 0;
                    for kk=temp1:ULmt-1
                        adjList(rem2,kk) = adjList(rem2,kk+1);
                    end
                    adjList(rem2,ULmt) = 0;
                end
            end
            edge(:,remover) = [];
            clear remover;
        end
        
        if (quarTimer > 14) %quarantine ends in two weeks
            inQuarantine = 0;
            disp 'done with Quarantine'
        end
        connections = 0;
    else
        quarTimer = 0;
        
        %Determine dC/dT
        %Later can implement a Time Delay (exposed but not yet infected)
        connections = round(alpha*A*(dmax-d) - beta*infectedE(days)*(d-dmin)); %No Time Delay


        %adding connections
        if (connections>0)
            for i = (1:connections)
                add1 = randi(A);
                add2 = randi(A);
                edge = [edge [add1, add2]']; %randomly selects two Extrovert Nodes only WITH REPETITION ALLOWED
                temp1 = 1;
                while (adjList(add1,temp1) ~= 0)
                    temp1 = temp1+1;
                end
                adjList(add1,temp1) = add2;
                temp1 = 1;
                while (adjList(add2,temp1) ~= 0)
                    temp1 = temp1+1;
                end
                adjList(add2,temp1) = add1;
            end
        end
        
        %makes sure if we are removing connections, there are enough to remove
        if(connections <= -1*E2E)
            disp 'Error: Not enough E2E connections to remove';
            return;
        end
        
        %removing connections
        if (connections<0)
            disp 'I am removing';
            remover(1)=0;
            for iii = 1:size(edge,2)
                if ((edge(1,iii) <=A) && (edge(2,iii) <=A))
                    if (remover(1)==0)
                        remover(1) = iii;
                    else
                        remover = [remover iii];
                    end
                end
            end
            toBeRemoved = randsample(remover,abs(connections));
            for iii = 1:abs(connections)
                rem1 = edge(1,toBeRemoved(iii));
                rem2 = edge(2,toBeRemoved(iii));
                temp1 = 1;
                while (adjList(rem1,temp1) ~= rem2)
                    temp1 = temp1+1;
                end
                adjList(rem1,temp1) = 0;
                for kk=temp1:ULmt-1
                    adjList(rem1,kk) = adjList(rem1,kk+1);
                end
                adjList(rem1,ULmt) = 0;
                temp1=1;
                while (adjList(rem2,temp1) ~= rem1)
                    temp1 = temp1+1;
                end
                adjList(rem2,temp1) = 0;
                for kk=temp1:ULmt-1
                    adjList(rem2,kk) = adjList(rem2,kk+1);
                end
                adjList(rem2,ULmt) = 0;
            end
            edge(:,toBeRemoved) = [];
            clear remover toBeRemoved;
        end
        
    end
    
    %Tracking extrovert connections only with time
    E2E = E2E+connections;
    
    days

end

%Collecting histigrams of Extroverts/Introverts' remaining susceptible neighbors
distrE(1:A)=0;
for i=1:A
    for j=1:(ULmt)
        if ((adjList(i,j) ~= 0) && (nodes(adjList(i,j)) < 1.5))
            distrE(i) = distrE(i)+1;
        end
    end
end

histogramE(1:(ULmt))=0;
for kk = 1:A
    histogramE(distrE(kk)+1) = histogramE(distrE(kk)+1) + 1;
end

distrI(1:B)=0;
for i=A+1:N
    for j=1:(ULmt)
        if ((adjList(i,j) ~= 0) && (nodes(adjList(i,j)) < 1.5))
            distrI(i-A) = distrI(i-A)+1;
        end
    end
end

histogramI(1:(ULmt))=0;
for kk = 1:B
    histogramI(distrI(kk)+1) = histogramI(distrI(kk)+1) + 1;
end






function Nreturn = getNode(Esel) %Parses the Designation to get just the Node number
P = find(Esel == 'E') - 1;
Nreturn = str2double(Esel(:,2:P));
end