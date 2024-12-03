% clearing the command Window
clear; home
clc;
% path to the location of the input files 
addpath('\\MSU-COSM-01.labs.msu\data$\Jarette Greene\GOlgatha cave Github\Golgotha-cave-drip-analysis-master\Golgotha-cave-drip-analysis-master\Data')
% asking user to enter acronymn of the cave data they would like to be
% anaylzed. NOTE: repeats input prompt until valid response is recieved 
prompt3 = "\nWhich cave data would you like to use Natural Bridge Cavern(NBC) or Harrie Wood Cave(HWC)? Enter NBC or HWC?: ";
cave_data = input(prompt3,'s');
cave_data = lower(cave_data);
while (cave_data ~= "nbc") && (cave_data~= "hwc")
    prompt = "\nInvalid response please enter valid response. NBC for Natural bridge cavern or HWC for Harrie wood cave.: ";
    cave_data = input(prompt,'s');
    cave_data = lower(cave_data);
end
%collecting input data for Natural Bridge cavern and Harrie Wood Cave based
%on the string previously entered
%% Collecting input data 
if cave_data == "nbc" 
    [D24,text24] = xlsread('Logger_daily_sum_full_ year.xlsx'); %#ok<XLSRD> daily driplogger data 
    [D24_15,text24_15] = xlsread('COMBINATIONJANMAY2024_clean.xlsx'); %#ok<XLSRD> 15 min drip logger data
    start_time= datetime(text24{2,1},'Locale','en_US','InputFormat','MM/dd/uuuu'); %start date of data collection
    end_time = datetime(text24{size(text24,1),1},'Locale','en_US','InputFormat','MM/dd/uuuu'); % end date of the data collection
    entire_time = start_time:end_time; %entire daily time
    entire_time = entire_time(:); %flipping to rows

    %getting the entire time in 15 minute intervals
    start_time_15 = datetime(text24_15{2,1});
    end_time_15 = datetime(text24_15{size(text24_15,1),1});
    end_time_15 = end_time_15 +minutes(15);
    entire_time_15 = start_time_15 : +minutes(15) : end_time_15;
    entire_time_15 = entire_time_15(:);
  
elseif cave_data == "hwc" 
    [D24_15,text24_15] = xlsread("HWC_All_Data_clean.xlsx") ; %#ok<XLSRD> 15 minutes data  
    [D24,text24] = xlsread("dripdata_Harrie_Wood_Cave_2014-2017.xlsx") ; %#ok<XLSRD> daily data 
    start_time = datetime(text24{3,1},'Locale','en_US','InputFormat','dd/MM/uuuu'); %daily start time
    end_time = datetime(text24{size(text24,1),1},'Locale','en_US','InputFormat','dd/MM/uuuu'); % daily end time
    entire_time = start_time:end_time; % entire daily time of study 
    entire_time = entire_time(:);
    start_time_15 = datetime(text24_15{2,1},'InputFormat','dd/MM/yyyy HH:mm'); %start time of 15min data
    end_time_15 = datetime(text24_15{size(text24_15,1),1},"InputFormat", 'dd/MM/yyyy HH:mm'); %end time of 15min data
    entire_time_15 = start_time_15 : +minutes(15) : end_time_15;
    entire_time_15 = entire_time_15(:);
    % removing the blank spaces from the daily data set 
    j = 1;
    %locating blank spaces
    for i = 2:2:size(D24,1)
       r(j) = i;
       j = j+1;
    end
    %removing the locations of the blank spaces.
    D24(r,:) = [];
end
%% interpolating missing data 
%finding locations of missing data points
[row,col] = find(isnan(D24));
[row2,col2] = find(isnan(D24_15));

%getting the names of the drip loggers
drip_log_name = text24(1,:);
drip_log_name(:,1) = [];

%compiling the missing data in the 15 minute data
temp_count = 1;
temp_count_15 = 1;
for i = 1:size(col2,1)
    temp_15 = col2(i,1);
    if i == 1
        missing_col2(1,temp_count_15) = temp_15;
        temp_count_15 = temp_count_15 + 1;
    else 
        if ~ismember(temp_15,missing_col2)
            missing_col2(1,temp_count_15) = temp_15;
            temp_count_15 = temp_count_15 + 1;
        end
    end
end
%converting dates into numbers
date_numbers_15 = convertTo(entire_time_15,'excel');

%putting 15 minute data together (dates with drip data)
all_data_15(:,1) = date_numbers_15;
all_data_15(:,2:size(D24_15,2)+1) = D24_15;

%performing interpolation
for i = 1:size(missing_col2,2) 
    temp2_15 = D24_15(:,missing_col2(:,i));
    missing_15 = isnan(D24_15(:,missing_col2(:,i)));
    %peforming interpolation based on the data being analyzed (changing
    %interpolation method based on what works best for the data)
    if cave_data == "hwc" 
        temp2_15(missing_15)= interp1(all_data_15(~missing_15),D24_15(~missing_15),all_data_15(missing_15),'nearest','extrap');
    elseif cave_data == "nbc" 
        temp2_15(missing_15)= interp1(all_data_15(~missing_15),D24_15(~missing_15),all_data_15(missing_15),'makima');
    end
    
    D24_15(:,missing_col2(:,i)) = temp2_15;
end

%performing the same interpolation process but this time with the daily
%data
for i = 1:size(col,1)
    temp = col(i,1);
    if i == 1
        missing_col(1,temp_count) = temp;
        temp_count = temp_count + 1;
    else 
        if ~ismember(temp,missing_col)
            missing_col(1,temp_count) = temp;
            temp_count = temp_count + 1;
        end
    end
end
date_numbers= convertTo(entire_time,'excel');
all_data(:,1) = date_numbers;
all_data(:,2:size(D24,2)+1) = D24;
for i = 1:size(missing_col,2)
    temp2 = D24(:,missing_col(:,i));
    missing = isnan(D24(:,missing_col(:,i)));
    temp2(missing) = interp1(all_data(~missing),D24(~missing),all_data(missing),'pchip');
    D24(:,missing_col(:,i)) = temp2;
end

%% Converting data into further intervals (Hourly, Weekly, Monthly) 
%converting 15 minutes data to hourly intervals
for i =  1:size(D24_15,2) 
   Data_temp = zeros(size(D24_15,1),1);
   Data_temp(:,1) = filter(ones(5,1),1,D24_15(:,i));
   Data_temp = Data_temp(5:5:end,:);
   D24_Hr(:,i) = Data_temp(:,1);
end


%sconverting daily data to weekly time intervals
D24_week = zeros([fix((size(D24,1)/7)),size(D24,2)]);
tempsum = zeros([1,size(D24,2)]);
count = 1;
count2 = 1;
for i = 1:size(D24,1)
    for j = 1:size(D24,2)
        tempsum(1,j) = tempsum(1,j) + D24(i,j);
    end
    count = count + 1;
    if count == 8
        D24_week(count2,:) = tempsum(1,:);
        count2 = count2 + 1;
        tempsum(:,:) = 0;
        count = 1;
    end
end
j = 1;
for i = 1:size(D24_week) 
    week_dates(i,1) = entire_time(j,1);
    j = j+6;
end
%converting daily data into monthly intervals 
previous_time = start_time;
tempsum2 = zeros([1,size(text24,2)-1]);
count1 = 1;
for i = 1:size(entire_time,1)
   current_time = datetime(entire_time(i,1));
    if(month(current_time) == month(previous_time))
        for j = 1:size(text24,2)-1
            tempsum2(1,j) = tempsum2(1,j) + D24(i,j);
        end
        previous_time = current_time;
    else
        D24_month(count1,:) = tempsum2(1,:);
        count1 = count1 + 1;
        tempsum2(:,:) = D24(i,:);
        previous_time = current_time;
    end
end
D24_month(count1,:) = tempsum2(1,:);
temp_date = start_time;
j =1;
for i = 1:size(entire_time,1)
  current = datetime(entire_time(i,1));
  if (month(current) ~= month(temp)) 
    month_date(j,1) = current;
    j= j+1;
    temp = current;
  else 
    temp = current;
  end
end

%% Calculating Distance Matrix

%calculating Distance Matrix 
for i = 1:size(D24,2)
    for j = 1:i
        ind1 = isfinite(D24(:,j));
        ind2 = isfinite(D24(:,i));
        ind1_2 = ind1 & ind2;
        v1 = D24(ind1_2,j); v2 = D24(ind1_2,i);
        
        %Correlation Coefficient
        c1 = corrcoef([D24(ind1_2,j),D24(ind1_2,i)]);
        corr1(i,j)= c1(1,2);
        corr1(j,i)= corr1(i,j);

        %Euclidean Distance
        edist1(i,j) = sqrt(mean(v1-v2).^2);
        edist1(j,i) = edist1(i,j);

        %Offset Distance
        [xcf1,lags1] = xcorr(v1,v2);
        m1 = find(xcf1 == max(xcf1));
        offsetdist1(i,j) = abs(lags1(m1(1)));
        offsetdist1(j,i) = offsetdist1(i,j);

        % New Distance
        NewDist1(i,j)= offsetdist1(i,j)*(1 - corr1(i,j));
        NewDist1(j,i)= NewDist1(i,j);
    end
    corr1(i,i) = 1;
end
  cdist1 = -corr1+1;
  % vector of colors to be used on graphs in the future
  colors24 = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 0 0 0];
  color = 1;
    

% asking the user how many cluster groups they would like and also which
% sampling frequency thry would like to see
prompt = "/nHow many cluster groups would you like to make (2/3/4 cluster groups)?: ";
no_of_cgroups = input(prompt);
while (no_of_cgroups ~= 2) && (no_of_cgroups ~= 3) && (no_of_cgroups ~=4 ) 
    prompt = "\nInvalid cluster grouping. Please enter valid cluster grouping (2,3,4).: ";
    no_of_cgroups = input(prompt);
end
prompt2 = "\nPick sampling frequncy by typing daily, weekly or monthly: ";
frequency = input(prompt2,'s');
frequency = lower(frequency);
while (frequency ~= "daily") && (frequency ~= "weekly") && (frequency ~= "monthly") 
    prompt2 = "\nInvalid sampling frequency selected. Please enter valid sampling frequency (daily, weekly, monthly): ";
    frequency = input(prompt2,'s');
    frequency = lower(frequency);
end
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));

%% Performing clustering
%beginning clustering using daily sampling frequency
if frequency == "daily" 

    %Performing Multidimensional Scaling (MDS) using Correlation
    %Coefficient distance matrix
    [Y1cc,E1cc] = cmdscale(cdist1);
    % Performing Kmeans cluster using results of MDS using CC
    [IDX1cc,C1cc]= kmeans(Y1cc,no_of_cgroups,'OnlinePhase','on','replicates',10);
    
    %Performing Multidimensional Scaling (MDS) using Euclidean distance matrix
    [Y1ed,E1ed] = cmdscale(edist1);
     % Performing Kmeans cluster using results of MDS using ED
    [IDX1ed,C1ed]= kmeans(Y1ed,no_of_cgroups,'OnlinePhase','on','replicates',10);
    
     %Performing Multidimensional Scaling (MDS) using Offset Distance matrix
    [Y1osd,E1osd] = cmdscale(offsetdist1);
     % Performing Kmeans cluster using results of MDS using OSD
    [IDX1osd,C1osd]= kmeans(Y1osd,no_of_cgroups,'OnlinePhase','on','replicates',10);
    
    %Performing Multidimensional Scaling (MDS) using New Distance matrix
    [Y1ND,E1ND] = cmdscale(NewDist1);
    % Performing Kmeans cluster using results of MDS using ND
    [IDX1ND,C1ND]= kmeans(Y1ND,no_of_cgroups,'OnlinePhase','on','replicates',10);
 

    j = 1;
    n = 1;
    m = 1;
    p = 1; 
    % getting the loggers that are in cluster 1 for each distance matrix
    for i = 1:size(D24,2)
        %cluster 1 of Correlation Coefficient
        if IDX1cc(i,1) == 1
            Cluster_1cc(1,j) = i;
            j = j+1;
        end
        %cluster 1 of Euclidean Distance
        if IDX1ed(i,1) == 1
            Cluster_1ed(1,n) = i;
            n = n+1;
        end
        %cluster 1 of Offset Distance
        if IDX1osd(i,1) == 1
            Cluster_1osd(1,m) = i;
            m = m+1;
        end
        %cluster 1 of New Distance
        if IDX1ND(i,1) == 1
            Cluster_1nd(1,p) = i;
            p = p+1;
        end
    end
    
    j = 1;
    n = 1;
    m = 1;
    p = 1;
%% Plotting Cluster Groupings 
    plot1cc = round(size(Cluster_1cc,2)/2);
    plot1ed = round(size(Cluster_1ed,2)/2);
    plot1osd = round(size(Cluster_1osd,2)/2);
    plot1nd = round(size(Cluster_1nd,2)/2);
    % graphing the loggers in cluster 1 for Correlation Coefficient
    for i = 1:size(Cluster_1cc,2)
        if color == 7
            color = 1;
        end
        figure(1)
        subplot (plot1cc,2,i)
        plot(entire_time(:,1),D24(:,Cluster_1cc(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_1cc(1,i)))
        sgtitle('Cluster1: Correlation Coeffiecient (Daily)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        % changing the axis based on the cave data being analyzed 
        if cave_data == "nbc"
            xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
            set(gca,'XTick',linspace(start_time,end_time,13))
            set(gca,'XTickLabel',xtick)
        elseif cave_data == "hwc"
             xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
             set(gca,'XTick',linspace(start_time,end_time,11))
             set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    color = 1;
    %graphing cluster 1 of Euclidean Distance
    for i = 1:size(Cluster_1ed,2)
        if color == 7
            color = 1;
        end
        figure(2)
        subplot (plot1ed,2,i)
        plot(entire_time(:,1),D24(:,Cluster_1ed(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_1ed(1,i)))
        sgtitle('Cluster1: Ecludiean Distance (Daily)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
         if cave_data == "nbc"
            xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
            set(gca,'XTick',linspace(start_time,end_time,13))
            set(gca,'XTickLabel',xtick)
        elseif cave_data == "hwc"
             xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
             set(gca,'XTick',linspace(start_time,end_time,11))
             set(gca,'XTickLabel',xtick)
         end
        color = color+1;
    end
    color = 1;
    % graphing cluster 1 of Offset Distance
    for i = 1:size(Cluster_1osd,2)
        if color == 7
            color = 1;
        end
        figure(3)
        subplot (plot1osd,2,i)
        plot(entire_time(:,1),D24(:,Cluster_1osd(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_1osd(1,i)))
        sgtitle('Cluster1: Offset Distance (Daily)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
         if cave_data == "nbc"
            xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
            set(gca,'XTick',linspace(start_time,end_time,13))
            set(gca,'XTickLabel',xtick)
        elseif cave_data == "hwc"
             xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
             set(gca,'XTick',linspace(start_time,end_time,11))
             set(gca,'XTickLabel',xtick)
         end
        color = color+1;
    end
    color = 1;
    % graphing cluster 1 of New Distance
    for i = 1:size(Cluster_1nd,2)
        if color == 7
            color = 1;
        end
        figure(4)
        subplot (plot1nd,2,i)
        plot(entire_time(:,1),D24(:,Cluster_1nd(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_1nd(1,i)))
        sgtitle('Cluster1: New Distance (Daily)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
         if cave_data == "nbc"
            xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
            set(gca,'XTick',linspace(start_time,end_time,13))
            set(gca,'XTickLabel',xtick)
        elseif cave_data == "hwc"
             xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
             set(gca,'XTick',linspace(start_time,end_time,11))
             set(gca,'XTickLabel',xtick)
         end
        color = color+1;
    end
    color = 1;
    
    j = 1;
    n = 1;
    m = 1;
    p = 1;
    % gathering loggers in cluster 2 from the different distance matrix
    for i = 1:size(D24,2)
        % gathering cluster 2 for Correlation Coefficient 
        if IDX1cc(i,1) == 2
            Cluster_2cc(1,j) = i;
            j = j+1;
        end
        % gathering cluster 2 for Euclidiean Distance
        if IDX1ed(i,1) == 2
            Cluster_2ed(1,n) = i;
            n = n+1;
        end
        % gathering cluster 2 for Offset Distance
        if IDX1osd(i,1) == 2
            Cluster_2osd(1,m) = i;
            m = m+1;
        end
        % gathering cluster 2 for New Distance
        if IDX1ND(i,1) == 2
            Cluster_2nd(1,p) = i;
            p = p+1;
        end
    end
    
    j = 1;
    n = 1;
    m = 1;
    p = 1;
    
    plot2cc = round(size(Cluster_2cc,2)/2);
    plot2ed = round(size(Cluster_2ed,2)/2);
    plot2osd = round(size(Cluster_2osd,2)/2);
    plot2nd = round(size(Cluster_2nd,2)/2);
    color = 1;
    %Plotting Cluster 2 of Correlation Coefficients
    for i = 1:size(Cluster_2cc,2)
        if color == 7
            color = 1;
        end
        figure(5)
        subplot (plot2cc,2,i)
        plot(entire_time(:,1),D24(:,Cluster_2cc(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_2cc(1,i)))
        sgtitle('Cluster2: Correlation Coeffiecient (Daily)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
         if cave_data == "nbc"
            xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
            set(gca,'XTick',linspace(start_time,end_time,13))
            set(gca,'XTickLabel',xtick)
        elseif cave_data == "hwc"
             xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
             set(gca,'XTick',linspace(start_time,end_time,11))
             set(gca,'XTickLabel',xtick)
         end
        color = color+1;
    end
    color = 1;
    % plotting Cluster 2 of Euclidean Distance
    for i = 1:size(Cluster_2ed,2)
        if color == 7
            color = 1;
        end
        figure(6)
        subplot (plot2ed,2,i)
        plot(entire_time(:,1),D24(:,Cluster_2ed(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_2ed(1,i)))
        sgtitle('Cluster2: Ecludiean Distance (Daily)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
         if cave_data == "nbc"
            xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
            set(gca,'XTick',linspace(start_time,end_time,13))
            set(gca,'XTickLabel',xtick)
        elseif cave_data == "hwc"
             xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
             set(gca,'XTick',linspace(start_time,end_time,11))
             set(gca,'XTickLabel',xtick)
         end
        color = color+1;
    end
    color = 1;
    %plotting Cluster 2 of Offset Distance 
    for i = 1:size(Cluster_2osd,2)
        if color == 7
            color = 1;
        end
        figure(7)
        subplot (plot2osd,2,i)
        plot(entire_time(:,1),D24(:,Cluster_2osd(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_2osd(1,i)))
        sgtitle('Cluster2: Offset Distance (Daily)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
            xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
            set(gca,'XTick',linspace(start_time,end_time,13))
            set(gca,'XTickLabel',xtick)
        elseif cave_data == "hwc"
             xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
             set(gca,'XTick',linspace(start_time,end_time,11))
             set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    color = 1;
    %Plotting Cluster 2 of New Distance
    for i = 1:size(Cluster_2nd,2)
        if color == 7
            color = 1;
        end
        figure(8)
        subplot (plot2nd,2,i)
        plot(entire_time(:,1),D24(:,Cluster_2nd(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_2nd(1,i)))
        sgtitle('Cluster2: New Distance (Daily)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
         if cave_data == "nbc"
            xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
            set(gca,'XTick',linspace(start_time,end_time,13))
            set(gca,'XTickLabel',xtick)
        elseif cave_data == "hwc"
             xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
             set(gca,'XTick',linspace(start_time,end_time,11))
             set(gca,'XTickLabel',xtick)
         end
        color = color+1;
    end
    color = 1;
    % if the no_of_cgroups is 3 or greater run this next session if not
    % skip
    % Gathering the loggers in Cluster 3 for all distance matrix 
    if no_of_cgroups >= 3 
        j = 1;
        n = 1;
        m = 1;
        p = 1;
     
        for i = 1:size(D24,2)
            % gathering cluster 3 for Correlation Coefficient 
            if IDX1cc(i,1) == 3
                Cluster_3cc(1,j) = i;
                j = j+1;
            end
            % gathering cluster 3 for Euclidean Distance 
            if IDX1ed(i,1) == 3
                Cluster_3ed(1,n) = i;
                n = n+1;
            end
            % gathering cluster 3 for Offset Distance 
            if IDX1osd(i,1) == 3
                Cluster_3osd(1,m) = i;
                m = m+1;
            end
            % gathering cluster 3 for New Distance
            if IDX1ND(i,1) == 3
                Cluster_3nd(1,p) = i;
                p = p+1;
            end
        end
        
        j = 1;
        n = 1;
        m = 1;
        p = 1;
        
        plot3cc = round(size(Cluster_3cc,2)/2);
        plot3ed = round(size(Cluster_3ed,2)/2);
        plot3osd = round(size(Cluster_3osd,2)/2);
        plot3nd = round(size(Cluster_3nd,2)/2);
        color = 1;
        % plotting the cluster 3 of Correlation Coefficient 
        for i = 1:size(Cluster_3cc,2)
            if color == 7
                color = 1;
            end
            figure(9)
            subplot (plot3cc,2,i)
            plot(entire_time(:,1),D24(:,Cluster_3cc(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_3cc(1,i)))
            sgtitle('Cluster3: Correlation Coeffiecient (Daily)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
             if cave_data == "nbc"
                xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
                set(gca,'XTick',linspace(start_time,end_time,13))
                set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
                xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
                set(gca,'XTick',linspace(start_time,end_time,11))
                set(gca,'XTickLabel',xtick)
             end
            color = color+1;
        end
        color = 1;
        % Plotting cluster 3 of the Euclidean Distance 
        for i = 1:size(Cluster_3ed,2)
            if color == 7
                color = 1;
            end
            figure(10)
            subplot (plot3ed,2,i)
            plot(entire_time(:,1),D24(:,Cluster_3ed(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_3ed(1,i)))
            sgtitle('Cluster3: Ecludiean Distance (Daily)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
             if cave_data == "nbc"
                xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
                set(gca,'XTick',linspace(start_time,end_time,13))
                set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
                xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
                set(gca,'XTick',linspace(start_time,end_time,11))
                set(gca,'XTickLabel',xtick)
             end
            color = color+1;
        end
        color = 1;
        % Plotting cluster 3 of the Offset Distance
        for i = 1:size(Cluster_3osd,2)
            if color == 7
                color = 1;
            end
            figure(11)
            subplot (plot3osd,2,i)
            plot(entire_time(:,1),D24(:,Cluster_3osd(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_3osd(1,i)))
            sgtitle('Cluster3: Offset Distance (Daily)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
                xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
                set(gca,'XTick',linspace(start_time,end_time,13))
                set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
                xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
                set(gca,'XTick',linspace(start_time,end_time,11))
                set(gca,'XTickLabel',xtick)
             end
            color = color+1;
        end
        color = 1;
        % Plotting cluster 3 of the New Distance 
        for i = 1:size(Cluster_3nd,2)
            if color == 7
                color = 1;
            end
            figure(12)
            subplot (plot3nd,2,i)
            plot(entire_time(:,1),D24(:,Cluster_3nd(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_3nd(1,i)))
            sgtitle('Cluster3: New Distance (Daily)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
                xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
                set(gca,'XTick',linspace(start_time,end_time,13))
                set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
                xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
                set(gca,'XTick',linspace(start_time,end_time,11))
                set(gca,'XTickLabel',xtick)
             end
            color = color+1;
        end
    end
    %if the no_of_cgroup is 4 or greater run this section or else skip it 
    % Gathering the loggers that are in cluster 4 for each Distance matrix
    if no_of_cgroups >=4 
        j = 1;
        n = 1;
        m = 1;
        p = 1;
        for i = 1:size(D24,2)
            % gathering cluster 4 for Correlation Coefficient
            if IDX1cc(i,1) == 4
                Cluster_4cc(1,j) = i;
                j = j+1;
            end
            % gathering cluster 4 for Euclidean Distance
            if IDX1ed(i,1) == 4
                Cluster_4ed(1,n) = i;
                n = n+1;
            end
            % gathering cluster 4 for Offset Distance
            if IDX1osd(i,1) == 4
                Cluster_4osd(1,m) = i;
                m = m+1;
            end
            % gathering cluster 4 for New Distance
            if IDX1ND(i,1) == 4
                Cluster_4nd(1,p) = i;
                p = p+1;
            end
        end
        
        j = 1;
        n = 1;
        m = 1;
        p = 1;
        
        plot4cc = round(size(Cluster_4cc,2)/2);
        plot4ed = round(size(Cluster_4ed,2)/2);
        plot4osd = round(size(Cluster_4osd,2)/2);
        plot4nd = round(size(Cluster_4nd,2)/2);
        color = 1;
        %Plotting Cluster 4 for Correlation Coefficient 
         for i = 1:size(Cluster_4cc,2)
            if color == 7
                color = 1;
            end
            figure(13)
            subplot (plot4cc,2,i)
            plot(entire_time(:,1),D24(:,Cluster_4cc(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_4cc(1,i)))
            sgtitle('Cluster4: Correlation Coeffiecient (Daily)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
                xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
                set(gca,'XTick',linspace(start_time,end_time,13))
                set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
                xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
                set(gca,'XTick',linspace(start_time,end_time,11))
                set(gca,'XTickLabel',xtick)
             end
            color = color+1;
        end
        color = 1;
        %Plotting Cluster 4 for Euclidean Distance
        for i = 1:size(Cluster_4ed,2)
            if color == 7
                color = 1;
            end
            figure(14)
            subplot (plot4ed,2,i)
            plot(entire_time(:,1),D24(:,Cluster_4ed(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_4ed(1,i)))
            sgtitle('Cluster4: Ecludiean Distance (Daily)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
                xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
                set(gca,'XTick',linspace(start_time,end_time,13))
                set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
                xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
                set(gca,'XTick',linspace(start_time,end_time,11))
                set(gca,'XTickLabel',xtick)
             end
            color = color+1;
        end
        color = 1;
        %Plotting Cluster 4 for Offset Distance 
        for i = 1:size(Cluster_4osd,2)
            if color == 7
                color = 1;
            end
            figure(15)
            subplot (plot4osd,2,i)
            plot(entire_time(:,1),D24(:,Cluster_4osd(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_4osd(1,i)))
            sgtitle('Cluster4: Offset Distance (Daily)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
                xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
                set(gca,'XTick',linspace(start_time,end_time,13))
                set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
                xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
                set(gca,'XTick',linspace(start_time,end_time,11))
                set(gca,'XTickLabel',xtick)
             end
            color = color+1;
        end
        color = 1;
        %Plotting Cluster 4 for New Distance 
        for i = 1:size(Cluster_4nd,2)
            if color == 7
                color = 1;
            end
            figure(16)
            subplot (plot4nd,2,i)
            plot(entire_time(:,1),D24(:,Cluster_4nd(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_4nd(1,i)))
            sgtitle('Cluster4: New Distance (Daily)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
                xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
                set(gca,'XTick',linspace(start_time,end_time,13))
                set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
                xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
                set(gca,'XTick',linspace(start_time,end_time,11))
                set(gca,'XTickLabel',xtick)
             end
            color = color+1;
        end
    end
%repeating the previous section but instead weekly data 
%only accesses this section if the user enters weekly 
elseif frequency == "weekly" 
    for i = 1:size(D24_week,2) 
        for j= 1:i
            ind11 = isfinite(D24_week(:,j));
            ind22 = isfinite(D24_week(:,i));
            ind11_22 = ind11 & ind22;
            v11 = D24_week(ind11_22,j); v22 = D24_week(ind11_22,i);
    
            c11 = corrcoef([D24_week(ind11_22,j),D24_week(ind11_22,i)]);
            corr11(i,j) = c11(1,2);
            corr11(j,i) = corr11(i,j);
    
            edist11(i,j) = sqrt(mean(v11-v22).^2);
            edist11(j,i) = edist11(i,j);
    
            [xcf11,lags11] = xcorr(v11,v22);
            m11 = find(xcf11 == max(xcf11));
            offsetdist11(i,j) = abs(lags11(m11(1)));
            offsetdist11(j,i) = offsetdist11(i,j);
    
            NewDist11(i,j) = offsetdist11(i,j)*(1 - corr11(i,j));
            NewDist11(j,i) = NewDist11(i,j);
        end
        corr11(i,i) = 1;
    end
    cdist11 = -corr11+1;
    CC11 = cdist11;
    ED11 = edist11;
    OSD11 = offsetdist11;
    NEWD11 = NewDist11;
    
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));
    
    [Y11cc,E11cc] = cmdscale(CC11);
    [IDX11cc,C11cc]= kmeans(Y11cc,no_of_cgroups,'OnlinePhase','on','replicates',10);
    
    [Y11ed,E11ed] = cmdscale(ED11);
    [IDX11ed,C11ed]= kmeans(Y11ed,no_of_cgroups,'OnlinePhase','on','replicates',10);
    
    [Y11osd,E11osd] = cmdscale(OSD11);
    [IDX11osd,C11osd]= kmeans(Y11osd,no_of_cgroups,'OnlinePhase','on','replicates',10);
    
    [Y11ND,E11ND] = cmdscale(NEWD11);
    [IDX11ND,C11ND]= kmeans(Y11ND,no_of_cgroups,'OnlinePhase','on','replicates',10);
    
    j = 1;
    n = 1;
    m = 1;
    p = 1;
    
    for i = 1:size(D24_week,2)
        if IDX11cc(i,1) == 1
            Cluster_11cc(1,j) = i;
            j = j+1;
        end
        if IDX11ed(i,1) == 1
            Cluster_11ed(1,n) = i;
            n = n+1;
        end
        if IDX11osd(i,1) == 1
            Cluster_11osd(1,m) = i;
            m = m+1;
        end
        if IDX11ND(i,1) == 1
            Cluster_11nd(1,p) = i;
            p = p+1;
        end
    end
    
    j = 1;
    n = 1;
    m = 1;
    p = 1;
    
    plot11cc = round(size(Cluster_11cc,2)/2);
    plot11ed = round(size(Cluster_11ed,2)/2);
    plot11osd = round(size(Cluster_11osd,2)/2);
    plot11nd = round(size(Cluster_11nd,2)/2);
    
    for i = 1:size(Cluster_11cc,2)
        if color == 7
            color = 1;
        end
        figure(100)
        subplot (plot11cc,2,i)
        plot(week_dates(:,1),D24_week(:,Cluster_11cc(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_11cc(1,i)))
        sgtitle('Cluster 1: Correlation Coeffiecient (Weekly)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
            xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
            set(gca,'XTick',linspace(start_time,end_time,16))
            set(gca,'XTickLabel',xtick)
       elseif cave_data == "hwc"
            xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
            set(gca,'XTick',linspace(start_time,end_time,13))
            set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    color = 1;
    
    for i = 1:size(Cluster_11ed,2)
        if color == 7
            color = 1;
        end
        figure(101)
        subplot (plot11ed,2,i)
        plot(week_dates(:,1),D24_week(:,Cluster_11ed(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_11ed(1,i)))
        sgtitle('Cluster 1: Ecludiean Distance (Weekly)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
            xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
            set(gca,'XTick',linspace(start_time,end_time,16))
            set(gca,'XTickLabel',xtick)
       elseif cave_data == "hwc"
            xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
            set(gca,'XTick',linspace(start_time,end_time,13))
            set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    color = 1;
    
    for i = 1:size(Cluster_11osd,2)
        if color == 7
            color = 1;
        end
        figure(102)
        subplot (plot11osd,2,i)
        plot(week_dates(:,1),D24_week(:,Cluster_11osd(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_11osd(1,i)))
        sgtitle('Cluster 1: Offset Distance (Weekly)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
            xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
            set(gca,'XTick',linspace(start_time,end_time,16))
            set(gca,'XTickLabel',xtick)
       elseif cave_data == "hwc"
            xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
            set(gca,'XTick',linspace(start_time,end_time,13))
            set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    color = 1;
    
    for i = 1:size(Cluster_11nd,2)
        if color == 7
            color = 1;
        end
        figure(103)
        subplot (plot11nd,2,i)
        plot(week_dates(:,1),D24_week(:,Cluster_11nd(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_11nd(1,i)))
        sgtitle('Cluster 1: New Distance (Weekly)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
            xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
            set(gca,'XTick',linspace(start_time,end_time,16))
            set(gca,'XTickLabel',xtick)
       elseif cave_data == "hwc"
            xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
            set(gca,'XTick',linspace(start_time,end_time,13))
            set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    color = 1;
    j = 1;
    n = 1;
    m = 1;
    p = 1;
    
    for i = 1:size(D24_week,2)
        if IDX11cc(i,1) == 2
            Cluster_22cc(1,j) = i;
            j = j+1;
        end
        if IDX11ed(i,1) == 2
            Cluster_22ed(1,n) = i;
            n = n+1;
        end
        if IDX11osd(i,1) == 2
            Cluster_22osd(1,m) = i;
            m = m+1;
        end
        if IDX11ND(i,1) == 2
            Cluster_22nd(1,p) = i;
            p = p+1;
        end
    end
    
    j = 1;
    n = 1;
    m = 1;
    p = 1;
    
    plot22cc = round(size(Cluster_22cc,2)/2);
    plot22ed = round(size(Cluster_22ed,2)/2);
    plot22osd = round(size(Cluster_22osd,2)/2);
    plot22nd = round(size(Cluster_22nd,2)/2);
    color = 1;
    for i = 1:size(Cluster_22cc,2)
        if color == 7
            color = 1;
        end
        figure(104)
        subplot (plot22cc,2,i)
        plot(week_dates(:,1),D24_week(:,Cluster_22cc(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_22cc(1,i)))
        sgtitle('Cluster 2: Correlation Coeffiecient (Weekly)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
            xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
            set(gca,'XTick',linspace(start_time,end_time,16))
            set(gca,'XTickLabel',xtick)
       elseif cave_data == "hwc"
            xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
            set(gca,'XTick',linspace(start_time,end_time,13))
            set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    color = 1;
    
    for i = 1:size(Cluster_22ed,2)
        if color == 7
            color = 1;
        end
        figure(105)
        subplot (plot22ed,2,i)
        plot(week_dates(:,1),D24_week(:,Cluster_22ed(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_22ed(1,i)))
        sgtitle('Cluster 2: Ecludiean Distance (Weekly)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
            xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
            set(gca,'XTick',linspace(start_time,end_time,16))
            set(gca,'XTickLabel',xtick)
       elseif cave_data == "hwc"
            xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
            set(gca,'XTick',linspace(start_time,end_time,13))
            set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    color = 1;
    
    for i = 1:size(Cluster_22osd,2)
        if color == 7
            color = 1;
        end
        figure(106)
        subplot (plot22osd,2,i)
        plot(week_dates(:,1),D24_week(:,Cluster_22osd(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_22osd(1,i)))
        sgtitle('Cluster 2: Offset Distance (Weekly)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
            xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
            set(gca,'XTick',linspace(start_time,end_time,16))
            set(gca,'XTickLabel',xtick)
       elseif cave_data == "hwc"
            xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
            set(gca,'XTick',linspace(start_time,end_time,13))
            set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    color = 1;
    
    for i = 1:size(Cluster_22nd,2)
        if color == 7
            color = 1;
        end
        figure(107)
        subplot (plot22nd,2,i)
        plot(week_dates(:,1),D24_week(:,Cluster_22nd(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_22nd(1,i)))
        sgtitle('Cluster 2: New Distance (Weekly)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
            xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
            set(gca,'XTick',linspace(start_time,end_time,16))
            set(gca,'XTickLabel',xtick)
       elseif cave_data == "hwc"
            xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
            set(gca,'XTick',linspace(start_time,end_time,13))
            set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    color = 1;
    if no_of_cgroups >= 3 
        j = 1;
        n = 1;
        m = 1;
        p = 1;
        
        for i = 1:size(D24_week,2)
            if IDX11cc(i,1) == 3
                Cluster_33cc(1,j) = i;
                j = j+1;
            end
            if IDX11ed(i,1) == 3
                Cluster_33ed(1,n) = i;
                n = n+1;
            end
            if IDX11osd(i,1) == 3
                Cluster_33osd(1,m) = i;
                m = m+1;
            end
            if IDX11ND(i,1) == 3
                Cluster_33nd(1,p) = i;
                p = p+1;
            end
        end
        
        j = 1;
        n = 1;
        m = 1;
        p = 1;
        
        plot33cc = round(size(Cluster_33cc,2)/2);
        plot33ed = round(size(Cluster_33ed,2)/2);
        plot33osd = round(size(Cluster_33osd,2)/2);
        plot33nd = round(size(Cluster_33nd,2)/2);
        color = 1;
        for i = 1:size(Cluster_33cc,2)
            if color == 7
                color = 1;
            end
            figure(108)
            subplot (plot33cc,2,i)
            plot(week_dates(:,1),D24_week(:,Cluster_33cc(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_33cc(1,i)))
            sgtitle('Cluster 3: Correlation Coeffiecient (Weekly)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
                xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
                set(gca,'XTick',linspace(start_time,end_time,16))
                set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
                xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
                set(gca,'XTick',linspace(start_time,end_time,13))
                set(gca,'XTickLabel',xtick)
            end
            color = color+1;
        end
        color = 1;
        
        for i = 1:size(Cluster_33ed,2)
            if color == 7
                color = 1;
            end
            figure(109)
            subplot (plot33ed,2,i)
            plot(week_dates(:,1),D24_week(:,Cluster_33ed(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_33ed(1,i)))
            sgtitle('Cluster 3: Ecludiean Distance (Weekly)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
                xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
                set(gca,'XTick',linspace(start_time,end_time,16))
                set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
                xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
                set(gca,'XTick',linspace(start_time,end_time,13))
                set(gca,'XTickLabel',xtick)
            end
            color = color+1;
        end
        color = 1;
        
        for i = 1:size(Cluster_33osd,2)
            if color == 7
                color = 1;
            end
            figure(110)
            subplot (plot33osd,2,i)
            plot(week_dates(:,1),D24_week(:,Cluster_33osd(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_33osd(1,i)))
            sgtitle('Cluster 3: Offset Distance (Weekly)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
                xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
                set(gca,'XTick',linspace(start_time,end_time,16))
                set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
                xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
                set(gca,'XTick',linspace(start_time,end_time,13))
                set(gca,'XTickLabel',xtick)
            end
            color = color+1;
        end
        color = 1;
        
        for i = 1:size(Cluster_33nd,2)
            if color == 7
                color = 1;
            end
            figure(111)
            subplot (plot33nd,2,i)
            plot(week_dates(:,1),D24_week(:,Cluster_33nd(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_33nd(1,i)))
            sgtitle('Cluster 3: New Distance (Weekly)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
                xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
                set(gca,'XTick',linspace(start_time,end_time,16))
                set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
                xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
                set(gca,'XTick',linspace(start_time,end_time,13))
                set(gca,'XTickLabel',xtick)
            end
            color = color+1;
        end
    end
    if no_of_cgroups >=4 
        j = 1;
        n = 1;
        m = 1;
        p = 1;
        
        for i = 1:size(D24_week,2)
            if IDX11cc(i,1) == 4
                Cluster_44cc(1,j) = i;
                j = j+1;
            end
            if IDX11ed(i,1) == 4
                Cluster_44ed(1,n) = i;
                n = n+1;
            end
            if IDX11osd(i,1) == 4
                Cluster_44osd(1,m) = i;
                m = m+1;
            end
            if IDX11ND(i,1) == 4
                Cluster_44nd(1,p) = i;
                p = p+1;
            end
        end
        
        j = 1;
        n = 1;
        m = 1;
        p = 1;
        
        plot44cc = round(size(Cluster_44cc,2)/2);
        plot44ed = round(size(Cluster_44ed,2)/2);
        plot44osd = round(size(Cluster_44osd,2)/2);
        plot44nd = round(size(Cluster_44nd,2)/2);
        color = 1;
         for i = 1:size(Cluster_44cc,2)
            if color == 7
                color = 1;
            end
            figure(112)
            subplot (plot44cc,2,i)
            plot(week_dates(:,1),D24_week(:,Cluster_44cc(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_44cc(1,i)))
            sgtitle('Cluster 4: Correlation Coeffiecient (Weekly)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
                xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
                set(gca,'XTick',linspace(start_time,end_time,16))
                set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
                xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
                set(gca,'XTick',linspace(start_time,end_time,13))
                set(gca,'XTickLabel',xtick)
            end
            color = color+1;
        end
        color = 1;
        
        for i = 1:size(Cluster_44ed,2)
            if color == 7
                color = 1;
            end
            figure(113)
            subplot (plot44ed,2,i)
            plot(week_dates(:,1),D24_week(:,Cluster_44ed(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_44ed(1,i)))
            sgtitle('Cluster 4: Ecludiean Distance (Weekly)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
                xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
                set(gca,'XTick',linspace(start_time,end_time,16))
                set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
                xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
                set(gca,'XTick',linspace(start_time,end_time,13))
                set(gca,'XTickLabel',xtick)
            end
            color = color+1;
        end
        color = 1;
        
        for i = 1:size(Cluster_44osd,2)
            if color == 7
                color = 1;
            end
            figure(114)
            subplot (plot44osd,2,i)
            plot(week_dates(:,1),D24_week(:,Cluster_44osd(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_44osd(1,i)))
            sgtitle('Cluster 4: Offset Distance (Weekly)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
                xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
                set(gca,'XTick',linspace(start_time,end_time,16))
                set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
                xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
                set(gca,'XTick',linspace(start_time,end_time,13))
                set(gca,'XTickLabel',xtick)
            end
            color = color+1;
        end
        color = 1;
        
        for i = 1:size(Cluster_44nd,2)
            if color == 7
                color = 1;
            end
            figure(115)
            subplot (plot44nd,2,i)
            plot(week_dates(:,1),D24_week(:,Cluster_44nd(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_44nd(1,i)))
            sgtitle('Cluster 4: New Distance (Weekly)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
                xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
                set(gca,'XTick',linspace(start_time,end_time,16))
                set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
                xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
                set(gca,'XTick',linspace(start_time,end_time,13))
                set(gca,'XTickLabel',xtick)
            end
            color = color+1;
        end
    end
% performing the same process as clustering as the last sections but
% instead using monthly data 
elseif frequency == "monthly" 
    for i = 1:size(D24_month,2) 
        for j= 1:i
            ind111 = isfinite(D24_month(:,j));
            ind222 = isfinite(D24_month(:,i));
            ind111_222 = ind111 & ind222;
            v111 = D24_month(ind111_222,j); v222 = D24_month(ind111_222,i);
    
            c111 = corrcoef([D24_month(ind111_222,j),D24_month(ind111_222,i)]);
            corr111(i,j) = c111(1,2);
            corr111(j,i) = corr111(i,j);
    
            edist111(i,j) = sqrt(mean(v111-v222).^2);
            edist111(j,i) = edist111(i,j);
    
            [xcf111,lags111] = xcorr(v111,v222);
            m111 = find(xcf111 == max(xcf111));
            offsetdist111(i,j) = abs(lags111(m111(1)));
            offsetdist111(j,i) = offsetdist111(i,j);
    
            NewDist111(i,j) = offsetdist111(i,j)*(1 - corr111(i,j));
            NewDist111(j,i) = NewDist111(i,j);
        end
        corr111(i,i) = 1;
    end
    cdist111 = -corr111+1;
    CC111 = cdist111;
    ED111 = edist111;
    OSD111 = offsetdist111;
    NEWD111 = NewDist111;
    
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));
    
    [Y111cc,E111cc] = cmdscale(CC111);
    [IDX111cc,C111cc]= kmeans(Y111cc,no_of_cgroups,'OnlinePhase','on','replicates',10);
    
    [Y111ed,E111ed] = cmdscale(ED111);
    [IDX111ed,C111ed]= kmeans(Y111ed,no_of_cgroups,'OnlinePhase','on','replicates',10);
    
    [Y111osd,E111osd] = cmdscale(OSD111);
    [IDX111osd,C111osd]= kmeans(Y111osd,no_of_cgroups,'OnlinePhase','on','replicates',10);
    
    [Y111ND,E111ND] = cmdscale(NEWD111);
    [IDX111ND,C111ND]= kmeans(Y111ND,no_of_cgroups,'OnlinePhase','on','replicates',10);
    
    j = 1;
    n = 1;
    m = 1;
    p = 1;
    
    for i = 1:size(D24_month,2)
        if IDX111cc(i,1) == 1
            Cluster_111cc(1,j) = i;
            j = j+1;
        end
        if IDX111ed(i,1) == 1
            Cluster_111ed(1,n) = i;
            n = n+1;
        end
        if IDX111osd(i,1) == 1
            Cluster_111osd(1,m) = i;
            m = m+1;
        end
        if IDX111ND(i,1) == 1
            Cluster_111nd(1,p) = i;
            p = p+1;
        end
    end
    
    j = 1;
    n = 1;
    m = 1;
    p = 1;
    
    plot111cc = round(size(Cluster_111cc,2)/2);
    plot111ed = round(size(Cluster_111ed,2)/2);
    plot111osd = round(size(Cluster_111osd,2)/2);
    plot111nd = round(size(Cluster_111nd,2)/2);
    
    for i = 1:size(Cluster_111cc,2)
        if color == 7
            color = 1;
        end
        figure(200)
        subplot (plot111cc,2,i)
        plot(month_date(:,1),D24_month(:,Cluster_111cc(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_111cc(1,i)))
        sgtitle('Cluster 1: Correlation Coeffiecient (Monthly)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
           xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
           set(gca,'XTick',linspace(start_time,end_time,14))
           set(gca,'XTickLabel',xtick)
        elseif cave_data == "hwc"
           xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
           set(gca,'XTick',linspace(start_time,end_time,12))
           set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    color = 1;
    
    for i = 1:size(Cluster_111ed,2)
        if color == 7
            color = 1;
        end
        figure(201)
        subplot (plot111ed,2,i)
        plot(month_date(:,1),D24_month(:,Cluster_111ed(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_111ed(1,i)))
        sgtitle('Cluster 1: Ecludiean Distance (Monthly)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
           xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
           set(gca,'XTick',linspace(start_time,end_time,14))
           set(gca,'XTickLabel',xtick)
        elseif cave_data == "hwc"
           xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
           set(gca,'XTick',linspace(start_time,end_time,12))
           set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    color = 1;
    
    for i = 1:size(Cluster_111osd,2)
        if color == 7
            color = 1;
        end
        figure(202)
        subplot (plot111osd,2,i)
        plot(month_date(:,1),D24_month(:,Cluster_111osd(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_111osd(1,i)))
        sgtitle('Cluster 1: Offset Distance (Monthly)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
           xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
           set(gca,'XTick',linspace(start_time,end_time,14))
           set(gca,'XTickLabel',xtick)
        elseif cave_data == "hwc"
           xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
           set(gca,'XTick',linspace(start_time,end_time,12))
           set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    color = 1;
    
    for i = 1:size(Cluster_111nd,2)
        if color == 7
            color = 1;
        end
        figure(203)
        subplot (plot111nd,2,i)
        plot(month_date(:,1),D24_month(:,Cluster_111nd(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_111nd(1,i)))
        sgtitle('Cluster 1: New Distance (Monthly)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
           xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
           set(gca,'XTick',linspace(start_time,end_time,14))
           set(gca,'XTickLabel',xtick)
        elseif cave_data == "hwc"
           xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
           set(gca,'XTick',linspace(start_time,end_time,12))
           set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    
    color = 1;
    j = 1;
    n = 1;
    m = 1;
    p = 1;
    
    for i = 1:size(D24_month,2)
        if IDX111cc(i,1) == 2
            Cluster_222cc(1,j) = i;
            j = j+1;
        end
        if IDX111ed(i,1) == 2
            Cluster_222ed(1,n) = i;
            n = n+1;
        end
        if IDX111osd(i,1) == 2
            Cluster_222osd(1,m) = i;
            m = m+1;
        end
        if IDX111ND(i,1) == 2
            Cluster_222nd(1,p) = i;
            p = p+1;
        end
    end
    
    j = 1;
    n = 1;
    m = 1;
    p = 1;
    
    plot222cc = round(size(Cluster_222cc,2)/2);
    plot222ed = round(size(Cluster_222ed,2)/2);
    plot222osd = round(size(Cluster_222osd,2)/2);
    plot222nd = round(size(Cluster_222nd,2)/2);
    color = 1;
    for i = 1:size(Cluster_222cc,2)
        if color == 7
            color = 1;
        end
        figure(204)
        subplot (plot222cc,2,i)
        plot(month_date(:,1),D24_month(:,Cluster_222cc(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_222cc(1,i)))
        sgtitle('Cluster 2: Correlation Coeffiecient (Monthly)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
           xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
           set(gca,'XTick',linspace(start_time,end_time,14))
           set(gca,'XTickLabel',xtick)
        elseif cave_data == "hwc"
           xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
           set(gca,'XTick',linspace(start_time,end_time,12))
           set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    color = 1;
    
    for i = 1:size(Cluster_222ed,2)
        if color == 7
            color = 1;
        end
        figure(205)
        subplot (plot222ed,2,i)
        plot(month_date(:,1),D24_month(:,Cluster_222ed(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_222ed(1,i)))
        sgtitle('Cluster 2: Ecludiean Distance (Monthly)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
           xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
           set(gca,'XTick',linspace(start_time,end_time,14))
           set(gca,'XTickLabel',xtick)
        elseif cave_data == "hwc"
           xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
           set(gca,'XTick',linspace(start_time,end_time,12))
           set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    color = 1;
    
    for i = 1:size(Cluster_222osd,2)
        if color == 7
            color = 1;
        end
        figure(206)
        subplot (plot222osd,2,i)
        plot(month_date(:,1),D24_month(:,Cluster_222osd(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_222osd(1,i)))
        sgtitle('Cluster 2: Offset Distance (Monthly)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
           xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
           set(gca,'XTick',linspace(start_time,end_time,14))
           set(gca,'XTickLabel',xtick)
        elseif cave_data == "hwc"
           xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
           set(gca,'XTick',linspace(start_time,end_time,12))
           set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    color = 1;
    
    for i = 1:size(Cluster_222nd,2)
        if color == 7
            color = 1;
        end
        figure(207)
        subplot (plot222nd,2,i)
        plot(month_date(:,1),D24_month(:,Cluster_222nd(1,i)),'Color',colors24(color,:),'LineStyle','--')
        title(drip_log_name(1,Cluster_222nd(1,i)))
        sgtitle('Cluster 2: New Distance (Monthly)')
        axis tight
        xlabel('Time')
        ylabel('Drips')
        if cave_data == "nbc"
           xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
           set(gca,'XTick',linspace(start_time,end_time,14))
           set(gca,'XTickLabel',xtick)
        elseif cave_data == "hwc"
           xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
           set(gca,'XTick',linspace(start_time,end_time,12))
           set(gca,'XTickLabel',xtick)
        end
        color = color+1;
    end
    
    color = 1;
    if no_of_cgroups >= 3 
        j = 1;
        n = 1;
        m = 1;
        p = 1;
        
        for i = 1:size(D24_month,2)
            if IDX111cc(i,1) == 3
                Cluster_333cc(1,j) = i;
                j = j+1;
            end
            if IDX111ed(i,1) == 3
                Cluster_333ed(1,n) = i;
                n = n+1;
            end
            if IDX111osd(i,1) == 3
                Cluster_333osd(1,m) = i;
                m = m+1;
            end
            if IDX111ND(i,1) == 3
                Cluster_333nd(1,p) = i;
                p = p+1;
            end
        end
        
        j = 1;
        n = 1;
        m = 1;
        p = 1;
        
        plot333cc = round(size(Cluster_333cc,2)/2);
        plot333ed = round(size(Cluster_333ed,2)/2);
        plot333osd = round(size(Cluster_333osd,2)/2);
        plot333nd = round(size(Cluster_333nd,2)/2);
        color = 1;
        for i = 1:size(Cluster_333cc,2)
            if color == 7
                color = 1;
            end
            figure(208)
            subplot (plot333cc,2,i)
            plot(month_date(:,1),D24_month(:,Cluster_333cc(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_333cc(1,i)))
            sgtitle('Cluster 3: Correlation Coeffiecient (Monthly)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
               xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
               set(gca,'XTick',linspace(start_time,end_time,14))
               set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
               xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
               set(gca,'XTick',linspace(start_time,end_time,12))
               set(gca,'XTickLabel',xtick)
            end
            color = color+1;
        end
        color = 1;
        
        for i = 1:size(Cluster_333ed,2)
            if color == 7
                color = 1;
            end
            figure(209)
            subplot (plot333ed,2,i)
            plot(month_date(:,1),D24_month(:,Cluster_333ed(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_333ed(1,i)))
             sgtitle('Cluster 3: Ecludiean Distance (Monthly)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
               xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
               set(gca,'XTick',linspace(start_time,end_time,14))
               set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
               xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
               set(gca,'XTick',linspace(start_time,end_time,12))
               set(gca,'XTickLabel',xtick)
            end
            color = color+1;
        end
        color = 1;
        
        for i = 1:size(Cluster_333osd,2)
            if color == 7
                color = 1;
            end
            figure(210)
            subplot (plot333osd,2,i)
            plot(month_date(:,1),D24_month(:,Cluster_333osd(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_333osd(1,i)))
            sgtitle('Cluster 3: Offset Distance (Monthly)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
               xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
               set(gca,'XTick',linspace(start_time,end_time,14))
               set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
               xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
               set(gca,'XTick',linspace(start_time,end_time,12))
               set(gca,'XTickLabel',xtick)
            end
            color = color+1;
        end
        color = 1;
        
        for i = 1:size(Cluster_333nd,2)
            if color == 7
                color = 1;
            end
            figure(211)
            subplot (plot333nd,2,i)
            plot(month_date(:,1),D24_month(:,Cluster_333nd(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_333nd(1,i)))
            sgtitle('Cluster 3: New Distance (Monthly)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
               xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
               set(gca,'XTick',linspace(start_time,end_time,14))
               set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
               xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
               set(gca,'XTick',linspace(start_time,end_time,12))
               set(gca,'XTickLabel',xtick)
            end
            color = color+1;
        end
    end
    
    if no_of_cgroups >=4 
        j = 1;
        n = 1;
        m = 1;
        p = 1;
        
        for i = 1:size(D24_month,2)
            if IDX111cc(i,1) == 4
                Cluster_444cc(1,j) = i;
                j = j+1;
            end
            if IDX111ed(i,1) == 4
                Cluster_444ed(1,n) = i;
                n = n+1;
            end
            if IDX111osd(i,1) == 4
                Cluster_444osd(1,m) = i;
                m = m+1;
            end
            if IDX111ND(i,1) == 4
                Cluster_444nd(1,p) = i;
                p = p+1;
            end
        end
        
        j = 1;
        n = 1;
        m = 1;
        p = 1;
        
        plot444cc = round(size(Cluster_444cc,2)/2);
        plot444ed = round(size(Cluster_444ed,2)/2);
        plot444osd = round(size(Cluster_444osd,2)/2);
        plot444nd = round(size(Cluster_444nd,2)/2);
        color = 1;
         for i = 1:size(Cluster_444cc,2)
            if color == 7
                color = 1;
            end
            figure(212)
            subplot (plot444cc,2,i)
            plot(month_date(:,1),D24_month(:,Cluster_444cc(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_444cc(1,i)))
            sgtitle('Cluster 4: Correlation Coeffiecient (Monthly)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
               xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
               set(gca,'XTick',linspace(start_time,end_time,14))
               set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
               xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
               set(gca,'XTick',linspace(start_time,end_time,12))
               set(gca,'XTickLabel',xtick)
            end
            color = color+1;
        end
        color = 1;
        
        for i = 1:size(Cluster_444ed,2)
            if color == 7
                color = 1;
            end
            figure(213)
            subplot (plot444ed,2,i)
            plot(month_date(:,1),D24_month(:,Cluster_444ed(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_444ed(1,i)))
            sgtitle('Cluster 4: Ecludiean Distance (Monthly)') 
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
               xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
               set(gca,'XTick',linspace(start_time,end_time,14))
               set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
               xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
               set(gca,'XTick',linspace(start_time,end_time,12))
               set(gca,'XTickLabel',xtick)
            end
            color = color+1;
        end
        color = 1;
        
        for i = 1:size(Cluster_444osd,2)
            if color == 7
                color = 1;
            end
            figure(214)
            subplot (plot444osd,2,i)
            plot(month_date(:,1),D24_month(:,Cluster_444osd(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_444osd(1,i)))
            sgtitle('Cluster 4: Offset Distance (Monthly)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
               xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
               set(gca,'XTick',linspace(start_time,end_time,14))
               set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
               xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
               set(gca,'XTick',linspace(start_time,end_time,12))
               set(gca,'XTickLabel',xtick)
            end
            color = color+1;
        end
        color = 1;
        
        for i = 1:size(Cluster_444nd,2)
            if color == 7
                color = 1;
            end
            figure(215)
            subplot (plot444nd,2,i)
            plot(month_date(:,1),D24_month(:,Cluster_444nd(1,i)),'Color',colors24(color,:),'LineStyle','--')
            title(drip_log_name(1,Cluster_444nd(1,i)))
            sgtitle('Cluster 4: New Distance (Monthly)')
            axis tight
            xlabel('Time')
            ylabel('Drips')
            if cave_data == "nbc"
               xtick = {'May23';'Jun23';'Jul23';'Aug23';'Sep23';'Oct23';'Nov23';'Dec23';'Jan24';'Feb24';'Mar24';'Apr24';'May24'};
               set(gca,'XTick',linspace(start_time,end_time,14))
               set(gca,'XTickLabel',xtick)
            elseif cave_data == "hwc"
               xtick = {'Jul14';'Oct14';'Jan15';'Apr15';'Jul15';'Oct15';'Jan16';'Apr16';'Jul16';'Oct16';'Jan17'};
               set(gca,'XTick',linspace(start_time,end_time,12))
               set(gca,'XTickLabel',xtick)
            end
            color = color+1;
        end
    end
end


%% Coefficient of Variance vs Sampling freqency

%size of different sampling frequencies 
k24 = [1 round(size(D24_15,1)/size(D24_Hr,1)) round(size(D24_15,1)/size(D24,1)) round(size(D24_15,1)/size(D24_week,1)) round(size(D24_15,1)/size(D24_month,1))];
flow_class24 = [1 0 0];
% Calculating the COV for each sampling frequency 
for i = 1:size(D24_15,2) 
    %Calculating the average drips for each sampling frequency
    % 15min, hourly, daily, weekly, monthly respectively
    avg_drip_15 = mean(D24_15(:,i));
    avg_drip_hr = mean(D24_Hr(:,i));
    avg_drip_daily = mean(D24(:,i));
    avg_drip_weekly = mean(D24_week(:,i));
    avg_drip_month = mean(D24_month(:,i));

    %Calculating the standard deviation for each sampling frequency
    %15min, hourly, daily, weekly, monthly respectively
    std_drip_15 = std(D24_15(:,i));
    std_drip_hr = std(D24_Hr(:,i));
    std_drip_daily = std(D24(:,i));
    std_drip_weekly = std(D24_week(:,i));
    std_drip_month = std(D24_month(:,i));

    %Calculating the Coefficient of Variance for each sampling frequency
    % 15min, hourly, daily, weekly, monthly respectively
    cov_24(i,1) = std_drip_15/avg_drip_15*100;
    cov_24(i,2) = std_drip_hr/avg_drip_hr*100;
    cov_24(i,3) = std_drip_daily/avg_drip_daily*100;
    cov_24(i,4) = std_drip_weekly/avg_drip_weekly*100;
    cov_24(i,5) = std_drip_month/avg_drip_month*100;
end
log_cov_24 = log(cov_24);

%Plotting the graph of COV but using daily cluster group the categorize
%each logger.
if frequency == "daily" 
    figure(20); clf;
    hold on;
    %plotting COV of Cluster 1 (Correlation Coefficient) represented by O on the graph
    for i = 1:size(Cluster_1cc,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_1cc(1,i),j),50,flow_class24,'o','LineWidth',1.5,'MarkerEdgeColor','r');
        end
    end
    %plotting COV of Cluster 2 (Correlation Coefficient) represented by a diamond shape
    for i = 1:size(Cluster_2cc,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_2cc(1,i),j),50,flow_class24,'diamond','LineWidth',1.5,'MarkerEdgeColor','b');
        end
    end
    %plotting COV of Cluster 3 (Correlation Coefficient) represented by a star 
    if no_of_cgroups >= 3 
        for i = 1:size(Cluster_3cc,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_3cc(1,i),j),50,flow_class24,'pentagram','LineWidth',1.5,'MarkerEdgeColor','m');
            end
        end
    end
    %plotting COV of Cluster 4 (Correlation Coefficient)represented by square
    if no_of_cgroups >= 4 
        for i = 1:size(Cluster_4cc,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_4cc(1,i),j),50,flow_class24,'square','LineWidth',1.5,'MarkerEdgeColor','k');
            end
        end
    end
    %performing line smoothing 
    for j = 1:size(D24_15,2) 
        % this function gives a warning, this warning can be ignored 
        [xx24,yy24] = smoothLine(log(k24(:)*15),log_cov_24(j,:));
        plot(xx24,yy24, '--g','LineWidth',1.2);
    end
    xlab = zeros(size(log_cov_24));
    for i = 1:size(xlab,2) 
        xlab(:,i) = log(k24(i)*15)*ones(size(D24_15,2),1);
    end
    xlab2 = zeros(size(1,size(log_cov_24,2)));
    for i = 1:size(log_cov_24,2) 
        xlab2(1,i) = xlab(1,i);
    end
    %Axis Labels
    axis equal square
    axis tight; box on;
    xlabel('ln [Sampling frequency(min)]')
    ylabel('ln [Coefficient of Variation (COV)]')
    title('Sampling frequency vs Coefficient (Correlation Coefficient) ');
    xtick = {'15mins';'1hr';'1day';'1week';'1month'};
    set(gca,'XTick',xlab2)
    set(gca,'XTickLabel',xtick)
    
    %Plotting the graph of COV but using daily cluster group the categorize
    %each logger.
    figure(21); clf;
    hold on;
    %plotting COV of Cluster 1 (Euclidean Distance) represented by O on the graph
    for i = 1:size(Cluster_1ed,2) 
        for j = 1:5 
      scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_1ed(1,i),j),50,flow_class24,'o','LineWidth',1.5,'MarkerEdgeColor','r');
           
        end
    end
    %plotting COV of Cluster 2 (Euclidean Distance) represented by a diamond shape
    for i = 1:size(Cluster_2ed,2) 
        for j = 1:5 
          scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_2ed(1,i),j),50,flow_class24,'diamond','LineWidth',1.5,'MarkerEdgeColor','b');
        
        end
    end
    %plotting COV of Cluster 3 (Euclidean Distance) represented by a star
    if no_of_cgroups >= 3 
        for i = 1:size(Cluster_3ed,2) 
            for j = 1:5 
         scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_3ed(1,i),j),50,flow_class24,'pentagram','LineWidth',1.5,'MarkerEdgeColor','m');
           
            end
        end
    end
    %plotting COV of Cluster 4 (Eucldiean distance)represented by square
    if no_of_cgroups >= 4 
        for i = 1:size(Cluster_4ed,2) 
            for j = 1:5 
           scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_4ed(1,i),j),50,flow_class24,'square','LineWidth',1.5,'MarkerEdgeColor','k');
         
            end
        end
    end
    for j = 1:size(D24_15,2) 
        [xx24,yy24] = smoothLine(log(k24(:)*15),log_cov_24(j,:));
        plot(xx24,yy24, '--g','LineWidth',1.2);
    end
    xlab = zeros(size(log_cov_24));
    for i = 1:size(xlab,2) 
        xlab(:,i) = log(k24(i)*15)*ones(size(D24_15,2),1);
    end
    xlab2 = zeros(size(1,size(log_cov_24,2)));
    for i = 1:size(log_cov_24,2) 
        xlab2(1,i) = xlab(1,i);
    end
    axis equal square
    axis tight; box on;
    xlabel('ln [Sampling frequency(min)]')
    ylabel('ln [Coefficient of Variation (COV)]')
    title('Sampling frequency vs Coefficient (Euclidiean Distance) ');
    xtick = {'15mins';'1hr';'1day';'1week';'1month'};
    set(gca,'XTick',xlab2)
    set(gca,'XTickLabel',xtick)
    
    figure(22); clf;
    hold on;
    %Plotting the graph of COV but using daily cluster group the categorize
    %each logger.
    %plotting COV of Cluster 1 (Offset Distance) represented by O on the graph
    for i = 1:size(Cluster_1osd,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_1osd(1,i),j),50,flow_class24,'o','LineWidth',1.5,'MarkerEdgeColor','r');
        end
    end
    %plotting COV of Cluster 2 (Offset Distance) represented by a diamond shape
    for i = 1:size(Cluster_2osd,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_2osd(1,i),j),50,flow_class24,'diamond','LineWidth',1.5,'MarkerEdgeColor','b');
        end
    end
    %plotting COV of Cluster 3 (Offset Distance) represented by a star
    if no_of_cgroups >= 3 
        for i = 1:size(Cluster_3osd,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_3osd(1,i),j),50,flow_class24,'pentagram','LineWidth',1.5,'MarkerEdgeColor','m');
            end
        end
    end
    %plotting COV of Cluster 4 (Offset Distance)represented by square
    if no_of_cgroups >= 4 
        for i = 1:size(Cluster_4osd,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_4osd(1,i),j),50,flow_class24,'square','LineWidth',1.5,'MarkerEdgeColor','k');
            end
        end
    end
    for j = 1:size(D24_15,2) 
        [xx24,yy24] = smoothLine(log(k24(:)*15),log_cov_24(j,:));
        plot(xx24,yy24, '--g','LineWidth',1.2);
    end
    xlab = zeros(size(log_cov_24));
    for i = 1:size(xlab,2) 
        xlab(:,i) = log(k24(i)*15)*ones(size(D24_15,2),1);
    end
    xlab2 = zeros(size(1,size(log_cov_24,2)));
    for i = 1:size(log_cov_24,2) 
        xlab2(1,i) = xlab(1,i);
    end
    axis equal square
    axis tight; box on;
    xlabel('ln [Sampling frequency(min)]')
    ylabel('ln [Coefficient of Variation (COV)]')
    title('Sampling frequency vs Coefficient (Offset Distance) ');
    xtick = {'15mins';'1hr';'1day';'1week';'1month'};
    set(gca,'XTick',xlab2)
    set(gca,'XTickLabel',xtick)
    %Plotting the graph of COV but using daily cluster group the categorize
    %each logger.
    figure(23); clf;
    hold on;
    %plotting COV of Cluster 1 (New Distance) represented by O on the graph
    for i = 1:size(Cluster_1nd,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_1nd(1,i),j),50,flow_class24,'o','LineWidth',1.5,'MarkerEdgeColor','r');
        end
    end
    %plotting COV of Cluster 2 (New Distance) represented by a diamond shape
    for i = 1:size(Cluster_2nd,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_2nd(1,i),j),50,flow_class24,'diamond','LineWidth',1.5,'MarkerEdgeColor','b');
        end
    end
    %plotting COV of Cluster 3 (New Distance) represented by a star
    if no_of_cgroups >= 3 
        for i = 1:size(Cluster_3nd,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_3nd(1,i),j),50,flow_class24,'pentagram','LineWidth',1.5,'MarkerEdgeColor','m');
            end
        end
    end
    %plotting COV of Cluster 4 (New Distance)represented by square
    if no_of_cgroups >= 4 
        for i = 1:size(Cluster_4nd,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_4nd(1,i),j),50,flow_class24,'square','LineWidth',1.5,'MarkerEdgeColor','k');
            end
        end
    end
    for j = 1:size(D24_15,2) 
        [xx24,yy24] = smoothLine(log(k24(:)*15),log_cov_24(j,:));
        plot(xx24,yy24, '--g','LineWidth',1.2);
    end
    xlab = zeros(size(log_cov_24));
    for i = 1:size(xlab,2) 
        xlab(:,i) = log(k24(i)*15)*ones(size(D24_15,2),1);
    end
    xlab2 = zeros(size(1,size(log_cov_24,2)));
    for i = 1:size(log_cov_24,2) 
        xlab2(1,i) = xlab(1,i);
    end
    axis equal square
    axis tight; box on;
    xlabel('ln [Sampling frequency(min)]')
    ylabel('ln [Coefficient of Variation (COV)]')
    title('Sampling frequency vs Coefficient (New Distance) ');
    xtick = {'15mins';'1hr';'1day';'1week';'1month'};
    set(gca,'XTick',xlab2)
    set(gca,'XTickLabel',xtick)

% Plotting the same graphs but using weekly data clustering
elseif frequency == "weekly" 
    figure(20); clf;
    hold on;
    for i = 1:size(Cluster_11cc,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_11cc(1,i),j),50,flow_class24,'o','LineWidth',1.5,'MarkerEdgeColor','r');
        end
    end
    for i = 1:size(Cluster_22cc,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_22cc(1,i),j),50,flow_class24,'diamond','LineWidth',1.5,'MarkerEdgeColor','b');
        end
    end
    if no_of_cgroups >= 3 
        for i = 1:size(Cluster_33cc,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_33cc(1,i),j),50,flow_class24,'pentagram','LineWidth',1.5,'MarkerEdgeColor','m');
            end
        end
    end
    if no_of_cgroups >= 4 
        for i = 1:size(Cluster_44cc,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_44cc(1,i),j),50,flow_class24,'square','LineWidth',1.5,'MarkerEdgeColor','k');
            end
        end
    end
    for j = 1:size(D24_15,2) 
        [xx24,yy24] = smoothLine(log(k24(:)*15),log_cov_24(j,:));
        plot(xx24,yy24, '--g','LineWidth',1.2);
    end
    xlab = zeros(size(log_cov_24));
    for i = 1:size(xlab,2) 
        xlab(:,i) = log(k24(i)*15)*ones(size(D24_15,2),1);
    end
    xlab2 = zeros(size(1,size(log_cov_24,2)));
    for i = 1:size(log_cov_24,2) 
        xlab2(1,i) = xlab(1,i);
    end
    axis equal square
    axis tight; box on;
    xlabel('ln [Sampling frequency(min)]')
    ylabel('ln [Coefficient of Variation (COV)]')
    title('Sampling frequency vs Coefficient (Correlation Coefficient) ');
    xtick = {'15mins';'1hr';'1day';'1week';'1month'};
    set(gca,'XTick',xlab2)
    set(gca,'XTickLabel',xtick)
    
    
    figure(21); clf;
    hold on;
    for i = 1:size(Cluster_11ed,2) 
        for j = 1:5 
      scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_11ed(1,i),j),50,flow_class24,'o','LineWidth',1.5,'MarkerEdgeColor','r');
           
        end
    end
    for i = 1:size(Cluster_22ed,2) 
        for j = 1:5 
          scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_22ed(1,i),j),50,flow_class24,'diamond','LineWidth',1.5,'MarkerEdgeColor','b');
        
        end
    end
    if no_of_cgroups >= 3 
        for i = 1:size(Cluster_33ed,2) 
            for j = 1:5 
         scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_33ed(1,i),j),50,flow_class24,'pentagram','LineWidth',1.5,'MarkerEdgeColor','m');
           
            end
        end
    end
    if no_of_cgroups >= 4 
        for i = 1:size(Cluster_44ed,2) 
            for j = 1:5 
           scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_44ed(1,i),j),50,flow_class24,'square','LineWidth',1.5,'MarkerEdgeColor','k');
         
            end
        end
    end
    for j = 1:size(D24_15,2) 
        [xx24,yy24] = smoothLine(log(k24(:)*15),log_cov_24(j,:));
        plot(xx24,yy24, '--g','LineWidth',1.2);
    end
    xlab = zeros(size(log_cov_24));
    for i = 1:size(xlab,2) 
        xlab(:,i) = log(k24(i)*15)*ones(size(D24_15,2),1);
    end
    xlab2 = zeros(size(1,size(log_cov_24,2)));
    for i = 1:size(log_cov_24,2) 
        xlab2(1,i) = xlab(1,i);
    end
    axis equal square
    axis tight; box on;
    xlabel('ln [Sampling frequency(min)]')
    ylabel('ln [Coefficient of Variation (COV)]')
    title('Sampling frequency vs Coefficient (Euclidiean Distance) ');
    xtick = {'15mins';'1hr';'1day';'1week';'1month'};
    set(gca,'XTick',xlab2)
    set(gca,'XTickLabel',xtick)
    
    figure(22); clf;
    hold on;
    for i = 1:size(Cluster_11osd,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_11osd(1,i),j),50,flow_class24,'o','LineWidth',1.5,'MarkerEdgeColor','r');
        end
    end
    for i = 1:size(Cluster_22osd,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_22osd(1,i),j),50,flow_class24,'diamond','LineWidth',1.5,'MarkerEdgeColor','b');
        end
    end
    if no_of_cgroups >= 3 
        for i = 1:size(Cluster_33osd,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_33osd(1,i),j),50,flow_class24,'pentagram','LineWidth',1.5,'MarkerEdgeColor','m');
            end
        end
    end
    if no_of_cgroups >= 4 
        for i = 1:size(Cluster_44osd,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_44osd(1,i),j),50,flow_class24,'square','LineWidth',1.5,'MarkerEdgeColor','k');
            end
        end
    end
    for j = 1:size(D24_15,2) 
        [xx24,yy24] = smoothLine(log(k24(:)*15),log_cov_24(j,:));
        plot(xx24,yy24, '--g','LineWidth',1.2);
    end
    xlab = zeros(size(log_cov_24));
    for i = 1:size(xlab,2) 
        xlab(:,i) = log(k24(i)*15)*ones(size(D24_15,2),1);
    end
    xlab2 = zeros(size(1,size(log_cov_24,2)));
    for i = 1:size(log_cov_24,2) 
        xlab2(1,i) = xlab(1,i);
    end
    axis equal square
    axis tight; box on;
    xlabel('ln [Sampling frequency(min)]')
    ylabel('ln [Coefficient of Variation (COV)]')
    title('Sampling frequency vs Coefficient (Offset Distance) ');
    xtick = {'15mins';'1hr';'1day';'1week';'1month'};
    set(gca,'XTick',xlab2)
    set(gca,'XTickLabel',xtick)
    
    figure(23); clf;
    hold on;
    for i = 1:size(Cluster_11nd,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_11nd(1,i),j),50,flow_class24,'o','LineWidth',1.5,'MarkerEdgeColor','r');
        end
    end
    for i = 1:size(Cluster_22nd,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_22nd(1,i),j),50,flow_class24,'diamond','LineWidth',1.5,'MarkerEdgeColor','b');
        end
    end
    if no_of_cgroups >= 3 
        for i = 1:size(Cluster_33nd,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_33nd(1,i),j),50,flow_class24,'pentagram','LineWidth',1.5,'MarkerEdgeColor','m');
            end
        end
    end
    if no_of_cgroups >= 4 
        for i = 1:size(Cluster_44nd,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_44nd(1,i),j),50,flow_class24,'square','LineWidth',1.5,'MarkerEdgeColor','k');
            end
        end
    end
    for j = 1:size(D24_15,2) 
        [xx24,yy24] = smoothLine(log(k24(:)*15),log_cov_24(j,:));
        plot(xx24,yy24, '--g','LineWidth',1.2);
    end
    xlab = zeros(size(log_cov_24));
    for i = 1:size(xlab,2) 
        xlab(:,i) = log(k24(i)*15)*ones(size(D24_15,2),1);
    end
    xlab2 = zeros(size(1,size(log_cov_24,2)));
    for i = 1:size(log_cov_24,2) 
        xlab2(1,i) = xlab(1,i);
    end
    axis equal square
    axis tight; box on;
    xlabel('ln [Sampling frequency(min)]')
    ylabel('ln [Coefficient of Variation (COV)]')
    title('Sampling frequency vs Coefficient (New Distance) ');
    xtick = {'15mins';'1hr';'1day';'1week';'1month'};
    set(gca,'XTick',xlab2)
    set(gca,'XTickLabel',xtick)
 
%plotting the same graphs but using Monthly clustering 
elseif frequency == "monthly" 
    figure(20); clf;
    hold on;
    for i = 1:size(Cluster_111cc,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_111cc(1,i),j),50,flow_class24,'o','LineWidth',1.5,'MarkerEdgeColor','r');
        end
    end
    for i = 1:size(Cluster_222cc,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_222cc(1,i),j),50,flow_class24,'diamond','LineWidth',1.5,'MarkerEdgeColor','b');
        end
    end
    if no_of_cgroups >= 3 
        for i = 1:size(Cluster_333cc,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_333cc(1,i),j),50,flow_class24,'pentagram','LineWidth',1.5,'MarkerEdgeColor','m');
            end
        end
    end
    if no_of_cgroups >= 4 
        for i = 1:size(Cluster_444cc,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_444cc(1,i),j),50,flow_class24,'square','LineWidth',1.5,'MarkerEdgeColor','k');
            end
        end
    end
    for j = 1:size(D24_15,2) 
        [xx24,yy24] = smoothLine(log(k24(:)*15),log_cov_24(j,:));
        plot(xx24,yy24, '--g','LineWidth',1.2);
    end
    xlab = zeros(size(log_cov_24));
    for i = 1:size(xlab,2) 
        xlab(:,i) = log(k24(i)*15)*ones(size(D24_15,2),1);
    end
    xlab2 = zeros(size(1,size(log_cov_24,2)));
    for i = 1:size(log_cov_24,2) 
        xlab2(1,i) = xlab(1,i);
    end
    axis equal square
    axis tight; box on;
    xlabel('ln [Sampling frequency(min)]')
    ylabel('ln [Coefficient of Variation (COV)]')
    title('Sampling frequency vs Coefficient (Correlation Coefficient) ');
    xtick = {'15mins';'1hr';'1day';'1week';'1month'};
    set(gca,'XTick',xlab2)
    set(gca,'XTickLabel',xtick)
    
    
    figure(21); clf;
    hold on;
    for i = 1:size(Cluster_111ed,2) 
        for j = 1:5 
      scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_111ed(1,i),j),50,flow_class24,'o','LineWidth',1.5,'MarkerEdgeColor','r');
           
        end
    end
    for i = 1:size(Cluster_222ed,2) 
        for j = 1:5 
          scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_222ed(1,i),j),50,flow_class24,'diamond','LineWidth',1.5,'MarkerEdgeColor','b');
        
        end
    end
    if no_of_cgroups >= 3 
        for i = 1:size(Cluster_333ed,2) 
            for j = 1:5 
         scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_333ed(1,i),j),50,flow_class24,'pentagram','LineWidth',1.5,'MarkerEdgeColor','m');
           
            end
        end
    end
    if no_of_cgroups >= 4 
        for i = 1:size(Cluster_444ed,2) 
            for j = 1:5 
           scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_444ed(1,i),j),50,flow_class24,'square','LineWidth',1.5,'MarkerEdgeColor','k');
         
            end
        end
    end
    for j = 1:size(D24_15,2) 
        [xx24,yy24] = smoothLine(log(k24(:)*15),log_cov_24(j,:));
        plot(xx24,yy24, '--g','LineWidth',1.2);
    end
    xlab = zeros(size(log_cov_24));
    for i = 1:size(xlab,2) 
        xlab(:,i) = log(k24(i)*15)*ones(size(D24_15,2),1);
    end
    xlab2 = zeros(size(1,size(log_cov_24,2)));
    for i = 1:size(log_cov_24,2) 
        xlab2(1,i) = xlab(1,i);
    end
    axis equal square
    axis tight; box on;
    xlabel('ln [Sampling frequency(min)]')
    ylabel('ln [Coefficient of Variation (COV)]')
    title('Sampling frequency vs Coefficient (Euclidiean Distance) ');
    xtick = {'15mins';'1hr';'1day';'1week';'1month'};
    set(gca,'XTick',xlab2)
    set(gca,'XTickLabel',xtick)
    
    figure(22); clf;
    hold on;
    for i = 1:size(Cluster_111osd,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_111osd(1,i),j),50,flow_class24,'o','LineWidth',1.5,'MarkerEdgeColor','r');
        end
    end
    for i = 1:size(Cluster_222osd,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_222osd(1,i),j),50,flow_class24,'diamond','LineWidth',1.5,'MarkerEdgeColor','b');
        end
    end
    if no_of_cgroups >= 3 
        for i = 1:size(Cluster_333osd,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_333osd(1,i),j),50,flow_class24,'pentagram','LineWidth',1.5,'MarkerEdgeColor','m');
            end
        end
    end
    if no_of_cgroups >= 4 
        for i = 1:size(Cluster_444osd,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_444osd(1,i),j),50,flow_class24,'square','LineWidth',1.5,'MarkerEdgeColor','k');
            end
        end
    end
    for j = 1:size(D24_15,2) 
        [xx24,yy24] = smoothLine(log(k24(:)*15),log_cov_24(j,:));
        plot(xx24,yy24, '--g','LineWidth',1.2);
    end
    xlab = zeros(size(log_cov_24));
    for i = 1:size(xlab,2) 
        xlab(:,i) = log(k24(i)*15)*ones(size(D24_15,2),1);
    end
    xlab2 = zeros(size(1,size(log_cov_24,2)));
    for i = 1:size(log_cov_24,2) 
        xlab2(1,i) = xlab(1,i);
    end
    axis equal square
    axis tight; box on;
    xlabel('ln [Sampling frequency(min)]')
    ylabel('ln [Coefficient of Variation (COV)]')
    title('Sampling frequency vs Coefficient (Offset Distance) ');
    xtick = {'15mins';'1hr';'1day';'1week';'1month'};
    set(gca,'XTick',xlab2)
    set(gca,'XTickLabel',xtick)
    
    figure(23); clf;
    hold on;
    for i = 1:size(Cluster_111nd,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_111nd(1,i),j),50,flow_class24,'o','LineWidth',1.5,'MarkerEdgeColor','r');
        end
    end
    for i = 1:size(Cluster_222nd,2) 
        for j = 1:5 
            scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_222nd(1,i),j),50,flow_class24,'diamond','LineWidth',1.5,'MarkerEdgeColor','b');
        end
    end
    if no_of_cgroups >= 3 
        for i = 1:size(Cluster_333nd,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_333nd(1,i),j),50,flow_class24,'pentagram','LineWidth',1.5,'MarkerEdgeColor','m');
            end
        end
    end
    if no_of_cgroups >= 4 
        for i = 1:size(Cluster_444nd,2) 
            for j = 1:5 
                scatter(log(k24(j)*15)*ones(size(D24_15,2),1),log_cov_24(Cluster_444nd(1,i),j),50,flow_class24,'square','LineWidth',1.5,'MarkerEdgeColor','k');
            end
        end
    end
    for j = 1:size(D24_15,2) 
        [xx24,yy24] = smoothLine(log(k24(:)*15),log_cov_24(j,:));
        plot(xx24,yy24, '--g','LineWidth',1.2);
    end
    xlab = zeros(size(log_cov_24));
    for i = 1:size(xlab,2) 
        xlab(:,i) = log(k24(i)*15)*ones(size(D24_15,2),1);
    end
    xlab2 = zeros(size(1,size(log_cov_24,2)));
    for i = 1:size(log_cov_24,2) 
        xlab2(1,i) = xlab(1,i);
    end
    axis equal square
    axis tight; box on;
    xlabel('ln [Sampling frequency(min)]')
    ylabel('ln [Coefficient of Variation (COV)]')
    title('Sampling frequency vs Coefficient (New Distance) ');
    xtick = {'15mins';'1hr';'1day';'1week';'1month'};
    set(gca,'XTick',xlab2)
    set(gca,'XTickLabel',xtick)
end
%% Collecting Statistical Data
Analysis_15min = cell(7,size(D24_15,2)+1);
Analysis_15min = cell2table(Analysis_15min);
% Adding driplogger labels to the approipriate columns (starting on the
% second column
for n = 2:size(Analysis_15min,2)
    Analysis_15min(1,n) = drip_log_name(1,n-1);
end
% Adding statistical  labels to the table 
Analysis_15min{2,1}{1} = 'Average Drip';
Analysis_15min{3,1}{1} = 'Drip Variance';
Analysis_15min{4,1}{1} = 'Standard Deviation';
Analysis_15min{5,1}{1} = 'Coefficient of Variance';
Analysis_15min{6,1}{1} = 'Total Drips';
Analysis_15min{7,1}{1} = 'Mode';
% Performing necessary statistical calculations
for m = 2:size(Analysis_15min,2) 
   %summation of all drips for a single logger
   total_drip24 = sum(D24_15(:,m-1));
   %total number of intervals of drip collects 
    total_time24 = size(D24_15,1)/4;
   % calculating avergae drip for specfic logger
    Analysis_15min{2,m} = num2cell(total_drip24/total_time24/4);
    % calculating drip variance
    Analysis_15min{3,m} = num2cell(var(D24_15(:,m-1)));
    %calculating Standard Deviation
    Analysis_15min{4,m} = num2cell(std(D24_15(:,m-1)));
    % Calculating Coeffiecient of Variance (COV)
    Analysis_15min{5,m} = num2cell(cell2mat(Analysis_15min{4,m})/cell2mat(Analysis_15min{2,m})*100);
    % placing total drips into the table
    Analysis_15min{6,m} = num2cell(total_drip24);
    % Calculating the Mode of the of the driplogger
    Analysis_15min{7,m} = num2cell(mode(D24_15(:,m-1)));
   % log_cov24(m-1,1) = log(cell2mat(Analysis_15min{4,m})/cell2mat(Analysis_15min{2,m})*100);
end
% converting Analysis_15min to an excel file
if cave_data == "nbc" 
    writetable(Analysis_15min,'NBC_Statistical_Data_15mins.xlsx')
else
    writetable(Analysis_15min,'HWC_Statistical_Data_15mins.xlsx')
end


Analysis_daily = cell(7,size(D24,2)+1);
Analysis_daily = cell2table(Analysis_daily);
% Adding driplogger labels to the approipriate columns (starting on the
% second column
for n = 2:size(Analysis_daily,2)
    Analysis_daily(1,n) = drip_log_name(1,n-1);
end
% Adding statistical  labels to the table 
Analysis_daily{2,1}{1} = 'Average Drip';
Analysis_daily{3,1}{1} = 'Drip Variance';
Analysis_daily{4,1}{1} = 'Standard Deviation';
Analysis_daily{5,1}{1} = 'Coefficient of Variance';
Analysis_daily{6,1}{1} = 'Total Drips';
Analysis_daily{7,1}{1} = 'Mode';
% Performing necessary statistical calculations
for m = 2:size(Analysis_daily,2) 
   %summation of all drips for a single logger
   total_drip24 = sum(D24(:,m-1));
   %total number of intervals of drip collects 
    total_time24 = size(D24,1)/4;
   % calculating avergae drip for specfic logger
    Analysis_daily{2,m} = num2cell(total_drip24/total_time24/4);
    % calculating drip variance
    Analysis_daily{3,m} = num2cell(var(D24(:,m-1)));
    %calculating Standard Deviation
    Analysis_daily{4,m} = num2cell(std(D24(:,m-1)));
    % Calculating Coeffiecient of Variance (COV)
    Analysis_daily{5,m} = num2cell(cell2mat(Analysis_daily{4,m})/cell2mat(Analysis_daily{2,m})*100);
    % placing total drips into the table
    Analysis_daily{6,m} = num2cell(total_drip24);
    % Calculating the Mode of the of the driplogger
    Analysis_daily{7,m} = num2cell(mode(D24(:,m-1)));
   % log_cov24(m-1,1) = log(cell2mat(Analysis_daily{4,m})/cell2mat(Analysis_daily{2,m})*100);
end
% converting Analysis_daily to an excel file
if cave_data == "nbc" 
    writetable(Analysis_15min,'NBC_Statistical_Data_daily.xlsx')
else
    writetable(Analysis_15min,'HWC_Statistical_Data_daily.xlsx')
end

Analysis_weekly = cell(7,size(D24_week,2)+1);
Analysis_weekly = cell2table(Analysis_weekly);
% Adding driplogger labels to the approipriate columns (starting on the
% second column
for n = 2:size(Analysis_weekly,2)
    Analysis_weekly(1,n) = drip_log_name(1,n-1);
end
% Adding statistical  labels to the table 
Analysis_weekly{2,1}{1} = 'Average Drip';
Analysis_weekly{3,1}{1} = 'Drip Variance';
Analysis_weekly{4,1}{1} = 'Standard Deviation';
Analysis_weekly{5,1}{1} = 'Coefficient of Variance';
Analysis_weekly{6,1}{1} = 'Total Drips';
Analysis_weekly{7,1}{1} = 'Mode';
% Performing necessary statistical calculations
for m = 2:size(Analysis_weekly,2) 
   %summation of all drips for a single logger
   total_drip24 = sum(D24_week(:,m-1));
   %total number of intervals of drip collects 
    total_time24 = size(D24_week,1)/4;
   % calculating avergae drip for specfic logger
    Analysis_weekly{2,m} = num2cell(total_drip24/total_time24/4);
    % calculating drip variance
    Analysis_weekly{3,m} = num2cell(var(D24_week(:,m-1)));
    %calculating Standard Deviation
    Analysis_weekly{4,m} = num2cell(std(D24_week(:,m-1)));
    % Calculating Coeffiecient of Variance (COV)
    Analysis_weekly{5,m} = num2cell(cell2mat(Analysis_weekly{4,m})/cell2mat(Analysis_weekly{2,m})*100);
    % placing total drips into the table
    Analysis_weekly{6,m} = num2cell(total_drip24);
    % Calculating the Mode of the of the driplogger
    Analysis_weekly{7,m} = num2cell(mode(D24_week(:,m-1)));
   % log_cov24(m-1,1) = log(cell2mat(Analysis_weekly{4,m})/cell2mat(Analysis_weekly{2,m})*100);
end
% converting Analysis_weekly to an excel file
if cave_data == "nbc" 
    writetable(Analysis_15min,'NBC_Statistical_Data_weekly.xlsx')
else
    writetable(Analysis_15min,'HWC_Statistical_Data_weekly.xlsx')
end

Analysis_monthly = cell(7,size(D24_month,2)+1);
Analysis_monthly = cell2table(Analysis_monthly);
% Adding driplogger labels to the approipriate columns (starting on the
% second column
for n = 2:size(Analysis_monthly,2)
    Analysis_monthly(1,n) = drip_log_name(1,n-1);
end
% Adding statistical  labels to the table 
Analysis_monthly{2,1}{1} = 'Average Drip';
Analysis_monthly{3,1}{1} = 'Drip Variance';
Analysis_monthly{4,1}{1} = 'Standard Deviation';
Analysis_monthly{5,1}{1} = 'Coefficient of Variance';
Analysis_monthly{6,1}{1} = 'Total Drips';
Analysis_monthly{7,1}{1} = 'Mode';
% Performing necessary statistical calculations
for m = 2:size(Analysis_monthly,2) 
   %summation of all drips for a single logger
   total_drip24 = sum(D24_month(:,m-1));
   %total number of intervals of drip collects 
    total_time24 = size(D24_month,1)/4;
   % calculating avergae drip for specfic logger
    Analysis_monthly{2,m} = num2cell(total_drip24/total_time24/4);
    % calculating drip variance
    Analysis_monthly{3,m} = num2cell(var(D24_month(:,m-1)));
    %calculating Standard Deviation
    Analysis_monthly{4,m} = num2cell(std(D24_month(:,m-1)));
    % Calculating Coeffiecient of Variance (COV)
    Analysis_monthly{5,m} = num2cell(cell2mat(Analysis_monthly{4,m})/cell2mat(Analysis_monthly{2,m})*100);
    % placing total drips into the table
    Analysis_monthly{6,m} = num2cell(total_drip24);
    % Calculating the Mode of the of the driplogger
    Analysis_monthly{7,m} = num2cell(mode(D24_month(:,m-1)));
   % log_cov24(m-1,1) = log(cell2mat(Analysis_monthly{4,m})/cell2mat(Analysis_monthly{2,m})*100);
end
% converting Analysis_monthly to an excel file
if cave_data == "nbc" 
    writetable(Analysis_15min,'NBC_Statistical_Data_monthly.xlsx')
else
    writetable(Analysis_15min,'HWC_Statistical_Data_monthly.xlsx')
end
