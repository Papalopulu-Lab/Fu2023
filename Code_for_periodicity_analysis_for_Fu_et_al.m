% The following is the code used for the data analysis used for results in
% Fu et. al.

clear all; close all; clc;

%% Preallocation of labels 

Condition = {'BMP4','PROLIF'};
Protein = {'GFP','SOX2mKate2'};

% names for the labels on figures

Condition_name = {'Quiescent','Proliferative'};
Protein_name = {'UbC-GFP','SOX2-mKate2'};

% number of experiments 

Expts = 3;


%% Importing the data

%Here is where the data must be imported, it needs to have each column as
%time series of gene expression with time running down.

%The data can then be separated into different conditions, proteins,
%experiments etc. as desired. 

% Here we separate the data into MATLAB arrays called "Data" in the following order
% {condition}{protien}{experiment}

% E.G. 
Data{1}{1}{1} = % matrix of all cells in quiescent condition, UbC-GFP expression in experiment 1.

%% Loop through conditions and proteins, experiments and cells, get statistics (mean, variance, etc.)

% Here we separate each cell further into its own array, first preallocate
% each statistic we want to derive, this is done at the beginning of each loop. We also provide cases where we combine
% all experiments.

% Within this section, we also develop arrays of detrended data via
% the subtraction of a Savitzky-Golay filter which shall then go into the
% oscillation analysis

Array_of_cells = cell(1,2);
Smoothened_array_of_cells = cell(1,2);
Variance = cell(1,2);
COV = cell(1,2);
COV_all_experiments = cell(1,2);
Fano_factor = cell(1,2);
Fano_factor_all_experiments = cell(1,2);
Variance_all_experiments = cell(1,2);
Mean_levels = cell(1,2);
Mean_levels_all_experiments = cell(1,2);
Max_fold_change = cell(1,2);
Max_fold_change_all_experiments = cell(1,2);
Max_amplitude = cell(1,2);
Max_amplitude_of_detrended_trace = cell(1,2);
Max_amplitude_all_experiments = cell(1,2);
Detrended_array_of_cells = cell(1,2);
Linear_trend = cell(1,2);
Linearly_detrended_array_of_cells = cell(1,2);

for condition_index = 1:2
    
    Array_of_cells{condition_index} = cell(1,2);
    Smoothened_array_of_cells{condition_index} = cell(1,2);
    Variance{condition_index} = cell(1,2);
    COV{condition_index} = cell(1,2);
    COV_all_experiments{condition_index} = cell(1,2);
    Fano_factor{condition_index} = cell(1,2);
    Fano_factor_all_experiments{condition_index} = cell(1,2);
    Variance_all_experiments{condition_index} = cell(1,2);
    Mean_levels{condition_index} = cell(1,2);
    Mean_levels_all_experiments{condition_index} = cell(1,2);
    Max_amplitude{condition_index} = cell(1,2);
    Max_amplitude_of_detrended_trace{condition_index} = cell(1,2);
    Max_amplitude_all_experiments{condition_index} = cell(1,2);
    Max_fold_change{condition_index} = cell(1,2);
    Max_fold_change_all_experiments{condition_index} = cell(1,2);
    Detrended_array_of_cells{condition_index} = cell(1,2);
    Linear_trend{condition_index} = cell(1,2);
    Linearly_detrended_array_of_cells{condition_index} = cell(1,2);
    
    
    for protein_index = 1:2
        
        Array_of_cells{condition_index}{protein_index} = cell(1,numel(expts));
        Smoothened_array_of_cells{condition_index}{protein_index} = cell(1,numel(expts));
        Variance{condition_index}{protein_index} = cell(1,numel(expts));
        COV{condition_index}{protein_index} = cell(1,numel(expts));
        COV_all_experiments{condition_index}{protein_index} = [];
        Fano_factor{condition_index}{protein_index} = cell(1,numel(expts));
        Fano_factor_all_experiments{condition_index}{protein_index} = [];
        Variance_all_experiments{condition_index}{protein_index} = [];
        Mean_levels{condition_index}{protein_index} = cell(1,numel(expts));
        Mean_levels_all_experiments{condition_index}{protein_index} = [];
        Max_fold_change{condition_index}{protein_index} = cell(1,numel(expts));
        Max_amplitude{condition_index}{protein_index} = cell(1,numel(expts));
        Max_amplitude_of_detrended_trace{condition_index}{protein_index} = cell(1,numel(expts));
        Max_amplitude_all_experiments{condition_index}{protein_index} = [];
        Max_fold_change_all_experiments{condition_index}{protein_index} = [];
        Detrended_array_of_cells{condition_index}{protein_index} = cell(1,numel(expts));
        Linear_trend{condition_index}{protein_index} = cell(1,numel(expts));
        Linearly_detrended_array_of_cells{condition_index}{protein_index} = cell(1,numel(expts));
        
        
        for expt_index = 1:numel(expts)
            
            Array_of_cells{condition_index}{protein_index}{expt_index} = cell(1,size(Data{condition_index}{protein_index}{expt_index},2));
            Smoothened_array_of_cells{condition_index}{protein_index}{expt_index} = cell(1,size(Data{condition_index}{protein_index}{expt_index},2));
            Variance{condition_index}{protein_index}{expt_index} = zeros(1,size(Data{condition_index}{protein_index}{expt_index},2));
            COV{condition_index}{protein_index}{expt_index} = zeros(1,size(Data{condition_index}{protein_index}{expt_index},2));
            Fano_factor{condition_index}{protein_index}{expt_index} = zeros(1,size(Data{condition_index}{protein_index}{expt_index},2));
            Mean_levels{condition_index}{protein_index}{expt_index} = zeros(1,size(Data{condition_index}{protein_index}{expt_index},2));
            Max_amplitude{condition_index}{protein_index}{expt_index} = zeros(1,size(Data{condition_index}{protein_index}{expt_index},2));
            Max_amplitude_of_detrended_trace{condition_index}{protein_index}{expt_index} = zeros(1,size(Data{condition_index}{protein_index}{expt_index},2));
            Max_fold_change{condition_index}{protein_index}{expt_index} = zeros(1,size(Data{condition_index}{protein_index}{expt_index},2));
            Detrended_array_of_cells{condition_index}{protein_index}{expt_index} = cell(1,size(Data{condition_index}{protein_index}{expt_index},2));
            Linear_trend{condition_index}{protein_index}{expt_index} = cell(1,size(Data{condition_index}{protein_index}{expt_index},2));
            Linearly_detrended_array_of_cells{condition_index}{protein_index}{expt_index} = cell(1,size(Data{condition_index}{protein_index}{expt_index},2));
            
            for cell_index = 1:size(Data{condition_index}{protein_index}{expt_index},2)
                
                Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index} = ...
                    Data{condition_index}{protein_index}{expt_index}(~isnan(Data{condition_index}{protein_index}{expt_index}(:,cell_index)),cell_index);
                
                Mean_levels{condition_index}{protein_index}{expt_index}(cell_index) = mean(Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index});
                
                Variance{condition_index}{protein_index}{expt_index}(cell_index) = var(Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index});
                
                COV{condition_index}{protein_index}{expt_index}(cell_index) = ...
                    std(Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index})/...
                    mean(Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index});
                
                Fano_factor{condition_index}{protein_index}{expt_index}(cell_index) = ...
                    var(Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index})/...
                    mean(Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index});
                
                
                Max_fold_change{condition_index}{protein_index}{expt_index}(cell_index) = ...
                    max(Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index})/min(Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index});
            
                Max_amplitude{condition_index}{protein_index}{expt_index}(cell_index) = ...
                    max(Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index})-min(Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index});
            
                if mod(length(Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index}),2) == 0
                    FL = length(Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index})-1;
                else 
                    FL = length(Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index});
                end
                
            Smoothened_array_of_cells{condition_index}{protein_index}{expt_index}{cell_index} =...
                smoothdata(Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index},'Gaussian',smoothing_window);
            
            Detrended_array_of_cells{condition_index}{protein_index}{expt_index}{cell_index} = ...
                Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index}-...
                Smoothened_array_of_cells{condition_index}{protein_index}{expt_index}{cell_index};
            
            Max_amplitude_of_detrended_trace{condition_index}{protein_index}{expt_index}(cell_index) = ...
                    max(Detrended_array_of_cells{condition_index}{protein_index}{expt_index}{cell_index})-min(Detrended_array_of_cells{condition_index}{protein_index}{expt_index}{cell_index});
            
            
            Linearly_detrended_array_of_cells{condition_index}{protein_index}{expt_index}{cell_index} = ...
                Array_of_cells{condition_index}{protein_index}{expt_index}{cell_index}-...
                Linear_trend{condition_index}{protein_index}{expt_index}{cell_index};
            
            
            end
                        
            Variance_all_experiments{condition_index}{protein_index} = [Variance_all_experiments{condition_index}{protein_index},Variance{condition_index}{protein_index}{expt_index}];
            Mean_levels_all_experiments{condition_index}{protein_index} = [Mean_levels_all_experiments{condition_index}{protein_index},Mean_levels{condition_index}{protein_index}{expt_index}];
            Max_fold_change_all_experiments{condition_index}{protein_index} = [Max_fold_change_all_experiments{condition_index}{protein_index},Max_fold_change{condition_index}{protein_index}{expt_index}];
            Max_amplitude_all_experiments{condition_index}{protein_index} = [Max_amplitude_all_experiments{condition_index}{protein_index},Max_amplitude{condition_index}{protein_index}{expt_index}];
            COV_all_experiments{condition_index}{protein_index} = [COV_all_experiments{condition_index}{protein_index},COV{condition_index}{protein_index}{expt_index}];
            Fano_factor_all_experiments{condition_index}{protein_index} = [Fano_factor_all_experiments{condition_index}{protein_index},Fano_factor{condition_index}{protein_index}{expt_index}];

        end  
    end
end


%% Loop through conditions, proteins, experiments and cells and perform oscillatory analysis

% Within this section, we use Lomb-Scargle periodogram to estimate the
% dominant periodicity of each time series trace

LSP_power = cell(1,2);
LSP_frequency = cell(1,2);
Dominant_period = cell(1,2);
Dominant_frequency = cell(1,2);
Dominant_power = cell(1,2);
VoS_freq = cell(1,2);
VOS_matrix = cell(1,2);
LSP_power_interpolated = cell(1,2);
LSP_power_interpolated_all_experiments = cell(1,2);
VOS_matrix_all_experiments = cell(1,2);
Dominant_period_all_experiments = cell(1,2);
Dominant_power_all_experiments = cell(1,2);
Normalised_data = cell(1,2);
Data_all_experiments = cell(1,2);
Normalised_data_all_experiments = cell(1,2);

max_length = 216;

for condition_index = 1:2
    
    LSP_power{condition_index} = cell(1,2);
    LSP_frequency{condition_index} = cell(1,2);
    Dominant_period{condition_index} = cell(1,2);
    Dominant_frequency{condition_index} = cell(1,2);
    Dominant_power{condition_index} = cell(1,2);
    VoS_freq{condition_index} = cell(1,2);
    VOS_matrix{condition_index} = cell(1,2);
    LSP_power_interpolated{condition_index} = cell(1,2);
    LSP_power_interpolated_all_experiments{condition_index} = cell(1,2);
    VOS_matrix_all_experiments{condition_index} = cell(1,2);
    Dominant_period_all_experiments{condition_index} = cell(1,2);
    Dominant_power_all_experiments{condition_index} = cell(1,2);   
    Normalised_data{condition_index} = cell(1,2); 
    Data_all_experiments{condition_index} = cell(1,2); 
    Normalised_data_all_experiments{condition_index} = cell(1,2); 
     
    for protein_index = 1:2
        
        LSP_power{condition_index}{protein_index} = cell(1,numel(expts));
        LSP_frequency{condition_index}{protein_index} = cell(1,numel(expts));
        Dominant_period{condition_index}{protein_index} = cell(1,numel(expts));
        Dominant_frequency{condition_index}{protein_index} = cell(1,numel(expts));
        Dominant_power{condition_index}{protein_index} = cell(1,numel(expts));
        VoS_freq{condition_index}{protein_index} = cell(1,numel(expts));
        VOS_matrix{condition_index}{protein_index} = cell(1,numel(expts));
        LSP_power_interpolated{condition_index}{protein_index} = cell(1,numel(expts));
        LSP_power_interpolated_all_experiments{condition_index}{protein_index} = [];
        VOS_matrix_all_experiments{condition_index}{protein_index} = [];
        Dominant_period_all_experiments{condition_index}{protein_index} = [];
        Dominant_power_all_experiments{condition_index}{protein_index} = [];
        Normalised_data{condition_index}{protein_index}= cell(1,numel(expts));
        Data_all_experiments{condition_index}{protein_index} = [];
        Normalised_data_all_experiments{condition_index}{protein_index} = [];

        for expt_index = 1:numel(expts)
            
            LSP_power{condition_index}{protein_index}{expt_index} = cell(1,size(Data{condition_index}{protein_index}{expt_index},2));
            LSP_frequency{condition_index}{protein_index}{expt_index} = cell(1,size(Data{condition_index}{protein_index}{expt_index},2));
            Dominant_period{condition_index}{protein_index}{expt_index} = zeros(1,size(Data{condition_index}{protein_index}{expt_index},2));
            Dominant_power{condition_index}{protein_index}{expt_index} = zeros(1,size(Data{condition_index}{protein_index}{expt_index},2));
            VoS_freq{condition_index}{protein_index}{expt_index} = cell(1,size(Data{condition_index}{protein_index}{expt_index},2));
            VOS_matrix{condition_index}{protein_index}{expt_index} = zeros(max_length,size(Data{condition_index}{protein_index}{expt_index},2));
            LSP_power_interpolated{condition_index}{protein_index}{expt_index}  = zeros(max_length,size(Data{condition_index}{protein_index}{expt_index},2));
            Normalised_data{condition_index}{protein_index}{expt_index}= nan(max_length,size(Data{condition_index}{protein_index}{expt_index},2));
            
            for cell_index = 1:size(Data{condition_index}{protein_index}{expt_index},2)
                
                 y = Detrended_array_of_cells{condition_index}{protein_index}{expt_index}{cell_index};

                time_vector = (1:length(y));
                
                [LSP_power{condition_index}{protein_index}{expt_index}{cell_index},...
                    LSP_frequency{condition_index}{protein_index}{expt_index}{cell_index}] = ...
                    plomb(y,time_vector,'normalized');
                
                Dominant_frequency{condition_index}{protein_index}{expt_index}(cell_index) =...
                    LSP_frequency{condition_index}{protein_index}{expt_index}{cell_index}...
                    (LSP_power{condition_index}{protein_index}{expt_index}{cell_index}...
                    ==max(LSP_power{condition_index}{protein_index}{expt_index}{cell_index}));
                
                Dominant_period{condition_index}{protein_index}{expt_index}(cell_index) =...
                    (1/Dominant_frequency{condition_index}{protein_index}{expt_index}(cell_index))/3;
                
                Dominant_power{condition_index}{protein_index}{expt_index}(cell_index) =...
                    max(LSP_power{condition_index}{protein_index}{expt_index}{cell_index});
                
                % interpolate power and freq in order to average correctly
                
                VoS_freq{condition_index}{protein_index}{expt_index}{cell_index}=...
                    linspace(LSP_frequency{condition_index}{protein_index}{expt_index}{cell_index}(1),...
                    LSP_frequency{condition_index}{protein_index}{expt_index}{cell_index}(end),max_length);
                
                VOS_matrix{condition_index}{protein_index}{expt_index}(:,cell_index) =...
                    VoS_freq{condition_index}{protein_index}{expt_index}{cell_index};
                
                LSP_power_interpolated{condition_index}{protein_index}{expt_index}(:,cell_index) =...
                    interp1(LSP_frequency{condition_index}{protein_index}{expt_index}{cell_index},...
                    LSP_power{condition_index}{protein_index}{expt_index}{cell_index},...
                    VoS_freq{condition_index}{protein_index}{expt_index}{cell_index});
                
                Normalised_data{condition_index}{protein_index}{expt_index}(1:length(y),cell_index) = zscore(y);
                
            end
            
            VOS_matrix_all_experiments{condition_index}{protein_index} = ...
                [VOS_matrix_all_experiments{condition_index}{protein_index},...
                VOS_matrix{condition_index}{protein_index}{expt_index}];
            
            LSP_power_interpolated_all_experiments{condition_index}{protein_index} = ...
                [LSP_power_interpolated_all_experiments{condition_index}{protein_index},...
                LSP_power_interpolated{condition_index}{protein_index}{expt_index}];
            
            Dominant_period_all_experiments{condition_index}{protein_index} = ...
                [Dominant_period_all_experiments{condition_index}{protein_index},...
                Dominant_period{condition_index}{protein_index}{expt_index}];
            
             Dominant_power_all_experiments{condition_index}{protein_index} = ...
                [Dominant_power_all_experiments{condition_index}{protein_index},...
                Dominant_power{condition_index}{protein_index}{expt_index}];
            
%              Data_all_experiments{condition_index}{protein_index}=...
%                  [Data_all_experiments{condition_index}{protein_index},...
%                 Data{condition_index}{protein_index}{expt_index}];
%             
%              Normalised_data_all_experiments{condition_index}{protein_index}=...
%                  [Normalised_data_all_experiments{condition_index}{protein_index},...
%                 Normalised_data{condition_index}{protein_index}{expt_index}];
            
            
        end
    end
end

%% Plot figures


% Figures can now be plotted using any statistics for example... 

% Mean leavels of all experiments


figure
boxplot([Mean_levels_all_experiments{1}{1},Mean_levels_all_experiments{1}{2},...
         Mean_levels_all_experiments{2}{1},Mean_levels_all_experiments{2}{2}],...
    [ones(1,length(Mean_levels_all_experiments{1}{1})),2*ones(1,length(Mean_levels_all_experiments{1}{2})),...
    3*ones(1,length(Mean_levels_all_experiments{2}{1})),4*ones(1,length(Mean_levels_all_experiments{2}{2}))])
xticks(1:4)
ylabel('Mean levels')


% Or LSP output, power/frequency for single cells or averages for each
% conditon/protein

figure
plot(VOS_matrix{1}{1}{1}{1}'.*3,mean(LSP_power_interpolated{1}{1}{1}{1}'),'LineWidth',2)
xlabel('Frequency')
ylabel('Power')
xlim([0,0.3])
legend([strcat(Condition_name{1},{' '},Protein_name{1}),...
    strcat(Condition_name{1},{' '},Protein_name{2}),...
    strcat(Condition_name{2},{' '},Protein_name{1}),...
    strcat(Condition_name{2},{' '},Protein_name{2})])

figure
plot(mean(VOS_matrix_all_experiments{1}{1}').*3,mean(LSP_power_interpolated_all_experiments{1}{1}'),'LineWidth',2)
hold on
plot(mean(VOS_matrix_all_experiments{1}{2}').*3,mean(LSP_power_interpolated_all_experiments{1}{2}'),'LineWidth',2)
plot(mean(VOS_matrix_all_experiments{2}{1}').*3,mean(LSP_power_interpolated_all_experiments{2}{1}'),'LineWidth',2)
plot(mean(VOS_matrix_all_experiments{2}{2}').*3,mean(LSP_power_interpolated_all_experiments{2}{2}'),'LineWidth',2)
xlabel('Frequency')
ylabel('Power')
xlim([0,0.3])
legend([strcat(Condition_name{1},{' '},Protein_name{1}),...
    strcat(Condition_name{1},{' '},Protein_name{2}),...
    strcat(Condition_name{2},{' '},Protein_name{1}),...
    strcat(Condition_name{2},{' '},Protein_name{2})])







