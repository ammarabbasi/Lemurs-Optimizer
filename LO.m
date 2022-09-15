%=======================================================================
%            Lemurs Optimizer: A New Metaheuristic Algorithm 
%               for Global Optimization (LO)

% This work is published in Journal of "Applied Sciences"
% URL:
% DOI: 

% Copyright (c) 2022, Ammar Kamal Abasi and Sharif Naser Makhadmeh and
% Mohammed Azmi Al-Betar and Osama Ahmad Alomari and
%Mohammed A. Awadallah and Zaid Abdi Alkareem Alyasseri and
% Iyad Abu Doush and Ashraf Elnagar and Eman H.Alkhammash12 and
% Myriam Hadjoun.
% All rights reserved.
%=======================================================================

clear all
close all
clc

PopSize = 30; %/* The number of Solutions*/
jumping_rate_min=0.1;   % JR0.1  0.1, 0.2, 0.3, 0.4, 0.5
jumping_rate_max=0.5;   %     0.9, 0.8, 0.7, 0.6, 0.5

Max_iter = 100000; %/*The number of cycles for foraging {a stopping criteria}*/

runs = 2; %/*Algorithm can be run many times in order to see its robustness*/

ObjVal = zeros(1, PopSize);

BestResults = zeros(runs, 1); % saving the best solution at each run

for funNum = 1:23 % fun#1 to fun#23
    if (funNum == 1)
        Function_name = 'F1';
    elseif (funNum == 2)
        Function_name = 'F2';
    elseif (funNum == 3)
        Function_name = 'F3';
    elseif (funNum == 4)
        Function_name = 'F4';
    elseif (funNum == 5)
        Function_name = 'F5';
    elseif (funNum == 6)
        Function_name = 'F6';
    elseif (funNum == 7)
        Function_name = 'F7';
    elseif (funNum == 8)
        Function_name = 'F8';
    elseif (funNum == 9)
        Function_name = 'F9';
    elseif (funNum == 10)
        Function_name = 'F10';
    elseif (funNum == 11)
        Function_name = 'F11';
    elseif (funNum == 12)
        Function_name = 'F12';
    elseif (funNum == 13)
        Function_name = 'F13';
    elseif (funNum == 14)
        Function_name = 'F14';
    elseif (funNum == 15)
        Function_name = 'F15';
    elseif (funNum == 16)
        Function_name = 'F16';
    elseif (funNum == 17)
        Function_name = 'F17';
    elseif (funNum == 18)
        Function_name = 'F18';
    elseif (funNum == 19)
        Function_name = 'F19';
    elseif (funNum == 20)
        Function_name = 'F20';
    elseif (funNum == 21)
        Function_name = 'F21';
    elseif (funNum == 22)
        Function_name = 'F22';
    elseif (funNum == 23)
        Function_name = 'F23';
    end
 
    % Load details of the selected benchmark function
    [lb, ub, dim, fobj] = Get_Functions_details(Function_name);
 
    for run = 1:runs
        % Initializing arrays
        swarm = zeros(PopSize, dim);
     
        % Initialize the population/solutions
        swarm = initialization(PopSize, dim, ub, lb);
     
        for i = 1:PopSize,
            ObjVal(i) = fobj(swarm(i, :));
        end
     
        Fitness = calculateFitness(ObjVal);
     
        %===================== loop ===================================
        tic
     
        itr = 0; % Loop counter
     
        while itr < Max_iter
         
            jumping_rate = jumping_rate_max - itr * ((jumping_rate_max - jumping_rate_min) / Max_iter);
           % swarm looking to go away from killer "free risk started too high " 
            [sorted_objctive, sorted_indexes] = sort(Fitness);
         
            current_solution = find(sorted_indexes == i);  % possition of the current solution 
         
            near_solution_postion = current_solution - 1;  % 
            if near_solution_postion == 0
                near_solution_postion = 1;
            end
            near_solution = sorted_indexes(near_solution_postion); 
         
            [cost, best_solution_Index] = min(ObjVal);
         
            for i = 1:PopSize,
             
                NewSol = swarm(i, :);
             
                for j = 1: dim,
                 
                    r = rand(); % select a number within range 0 to 1.
                 
                    if (r < jumping_rate)   
                     
                        NewSol(j) = swarm(i, j) + abs(swarm(i, j) - swarm(near_solution, j)) * (rand - 0.5) * 2;
                     
                        % manipulate range between lb and ub
                      if (size(lb,2)~=1)
                          NewSol(j) = min(max(NewSol(j), lb(j)), ub(j));
                      else
                          
                        NewSol(j) = min(max(NewSol(j), lb), ub);
                      end
                    else
                        % for long jumbing will take from best solution
                        NewSol(j) = swarm(i, j) + abs(swarm(i, j) - swarm(best_solution_Index, j)) * (rand - 0.5) * 2;
                        % manipulate range between lb and ub
                        
                       if (size(lb,2)~=1)
                          NewSol(j) = min(max(NewSol(j), lb(j)), ub(j));
                      else
                          
                        NewSol(j) = min(max(NewSol(j), lb), ub);
                      end
                     
                    end
                 
                end
             
                %evaluate new solution
                ObjValSol = fobj(NewSol);
                FitnessSol = calculateFitness(ObjValSol);
             
                % Update the curent solution  & Age of the current solution
                if (ObjVal(i) > ObjValSol)
                    swarm(i, :) = NewSol;
                    Fitness(i) = FitnessSol;
                    ObjVal(i) = ObjValSol;
                 
                end
             
            end
         
            if (mod(itr, 100) == 0)
             
                display(['Fun#', num2str(funNum), ' Run#', num2str(run), ', Itr ', num2str(itr), ' Results ', num2str(min(ObjVal))]);
            end
            
            itr = itr + 1;
            conv(:,itr)=min(ObjVal);
        end
     
        toc;
     
        % Save the best results at each iteration
        BestResults(run) = min(ObjVal);
		cg_curve(run,:)=conv;
    end % run
    
   Best=min(BestResults);
   worst=max(BestResults);
   avg=mean(BestResults);
   mstd=std(BestResults);
   Overall=[Best,worst,avg,mstd];

 ave_cg_curve(1,:)=mean(cg_curve);
 min_cg_curve(1,:)=min(cg_curve);
 
figure('Position',[290   206   648   287])

%Draw the search space
subplot(1,2,1);
func_plot(Function_name);
title('Test function')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])
grid off
shading interp;
light;
lighting phong;
shading interp;

%Draw the convergence curve
subplot(1,2,2);
semilogy( ave_cg_curve(1:200),'Color','r')
hold on 
semilogy( min_cg_curve(1:200),'Color','b')
title('Mean and best convergence curve')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid off
box on
legend('Mean-LO','Best-LO' )
     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Write the best results for all runs in excel file
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fileName = strcat('LO',',Fun_',num2str(funNum),',Sel_', datestr(now, ', dd-mmmm-yyyy,HH-MM-SS_FFF_AM'),'.xlsx');
  
  fprintf(1,fileName);
  
  sheetName = strcat('fun_',num2str(funNum),'_Best_Results');
  sheetName2 = strcat('fun_',num2str(funNum),'_Overall_Results');
  xlswrite(fileName,BestResults,sheetName);% write results to excel file
  xlswrite(fileName,Overall,sheetName2);% write results to excel file
%  sheetName = strcat('fun_',num2str(func_num),'_Error_Rate');
  
    %  xlswrite(fileName,ErrorRate,sheetName);
 saveas(gcf,[num2str(funNum),'.fig']);

    fprintf(1, '\n\n Done \n\n');
 
end
