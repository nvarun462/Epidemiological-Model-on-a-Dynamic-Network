%This code will run College_Network_Final multiple times and automatically
%average only the successful runs and plot the data


clear all; clc;
iterations = 3; %How many good runs are wanted
max_days = 200; %used to pre-allocate below matrices. Make sure is equal to max_days parameter in College_Network2
dmax = 25;
ULmt = round(dmax*2.5);

%pre-allocate matrices
susceptibleOVR(1:max_days) = 0;
infectedOVR(1:max_days) = 0;
recoveredOVR(1:max_days) = 0;
susceptibleIOVR(1:max_days) = 0;
infectedIOVR(1:max_days) = 0;
recoveredIOVR(1:max_days) = 0;
susceptibleEOVR(1:max_days) = 0;
infectedEOVR(1:max_days) = 0;
recoveredEOVR(1:max_days) = 0;
E_connectionsOVR(1:max_days) = 0; %E_to_E_connectionCount
D_avgOVR(1:max_days) = 0; %D_avg as a function of days
HistogramEALL(iterations,1:ULmt) = 0;
HistogramEAVG(1:ULmt) = 0;
HistogramIALL(iterations,1:ULmt) = 0;
HistogramIAVG(1:ULmt) = 0;
kkk = 0;
while (kkk<iterations) %keep going until we have (iterations) good runs
College_Network_Final;
if (susceptible(max_days)/N < .9) %a run is considered 'good' if the virus reached at least 10% of the population before dying off
    susceptibleOVR(1:max_days) = susceptibleOVR(1:max_days) + susceptible(1:max_days);
    infectedOVR(1:max_days) = infectedOVR(1:max_days) + infected(1:max_days);
    recoveredOVR(1:max_days) = recoveredOVR(1:max_days) + recovered(1:max_days);
    susceptibleIOVR(1:max_days) = susceptibleIOVR(1:max_days) + susceptibleI(1:max_days);
    infectedIOVR(1:max_days) = infectedIOVR(1:max_days) + infectedI(1:max_days);
    recoveredIOVR(1:max_days) = recoveredIOVR(1:max_days) + recoveredI(1:max_days);
    susceptibleEOVR(1:max_days) = susceptibleEOVR(1:max_days) + susceptibleE(1:max_days);
    infectedEOVR(1:max_days) = infectedEOVR(1:max_days) + infectedE(1:max_days);
    recoveredEOVR(1:max_days) = recoveredEOVR(1:max_days) + recoveredE(1:max_days);
    E_connectionsOVR(1:max_days) = E_connectionsOVR(1:max_days) + E_connections(1:max_days);
    D_avgOVR(1:max_days) = D_avgOVR(1:max_days) + D_avg(1:max_days);
    kkk = kkk+1 %kkk will only increment if it is a good run
    HistogramEALL(kkk,1:ULmt) = histogramE;
    HistogramIALL(kkk,1:ULmt) = histogramI;
else
    disp 'Failed Run';
end
end

%PLOTTING DATA INTO GRAPHS----------------------------------------

% for iii = 1:iterations
% HistogramEAVG(1:ULmt) = HistogramEAVG(1:ULmt) + HistogramEALL(iii,(1:ULmt));
% HistogramIAVG(1:ULmt) = HistogramIAVG(1:ULmt) + HistogramIALL(iii,(1:ULmt));
% HistogramALL(iii,1:ULmt) = HistogramEALL(iii,1:ULmt) + HistogramIALL(iii,1:ULmt);
% end
% 
% 
% 
% HistogramAVG(1:ULmt) = HistogramEAVG(1:ULmt) + HistogramIAVG(1:ULmt);


% h=figure;
% subplot(1,3,1);
% plot(HistogramAVG/iterations/N);grid;hold;
% xlabel('Number of Susceptible Neighbors');
% ylabel('Proportions');
% title('Degree Distribution for Entire Network');
% subplot(1,3,2);
% plot(HistogramEAVG/iterations/A);grid;hold;
% xlabel('Number of Susceptible Neighbors');
% ylabel('Proportions');
% title('Degree Distribution for Extroverts');
% subplot(1,3,3);
% plot(HistogramIAVG/iterations/B);grid;hold;
% xlabel('Number of Susceptible Neighbors');
% ylabel('Proportions');
% title('Degree Distribution for Introverts');


% save histograms.mat HistogramIALL HistogramEALL HistogramALL

% %Plot SIR Entire Network
% h=figure;
% plot(susceptibleOVR/iterations/N,'b');grid;hold;
% plot(infectedOVR/iterations/N,'r');
% plot(recoveredOVR/iterations/N,'g');
% title('SIR for Entire Network');
% xlabel('Days');
% ylabel('Proportion');
% saveas(h,'SIR_Entire_Network_Fig3.fig');

% %Plot SIR for Introverts
% h1=figure;
% plot(susceptibleIOVR/iterations/B,'b');grid;hold;
% plot(infectedIOVR/iterations/B,'r');
% plot(recoveredIOVR/iterations/B,'g');
% title('SIR for Introverts Only');
% xlabel('Days');
% ylabel('Proportion');
% saveas(h1,'SIR_Introverts_newForm.fig');
% 
% %Plot SIR for Extroverts
% h2=figure;
% plot(susceptibleEOVR/iterations/A,'b');grid;hold;
% plot(infectedEOVR/iterations/A,'r');
% plot(recoveredEOVR/iterations/A,'g');
% title('SIR for Extroverts Only');
% xlabel('Days');
% ylabel('Proportion');
% saveas(h2,'SIR_Extroverts_newForm.fig');

% %Plot Average Degree of Extroverts
% h3=figure;
% plot(D_avgOVR/iterations,'k');grid;
% title('Average degree of Extroverts');
% xlabel('Days');
% ylabel('Degree');
% saveas(h3,'Avg_degree_Extroverts_Fig3.fig');
% 
% save Fig3Data.mat susceptibleOVR infectedOVR recoveredOVR iterations N D_avgOVR;
