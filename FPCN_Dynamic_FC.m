function [Dynamic_corr temporal_data]=FPCN_Dynamic_FC(ROI_Timecourses, Ci)

%computes dynamic changes across 60-sec windows in mean between-network FC and clustering and their relationship
%Input: 
    %ROI_Timecourses; a file containing preprocessed ROI timecourse data for a given subject and condition. Data should be filtered for 60 sec window, and organized as timepoints (rows) x ROIs (columns)
    %Ci: vector of community assignments for all nodes, where Ci=1 for DAN nodes, Ci=2 for DN nodes, Ci=3 for FPCNa nodes, and Ci=4 for FPCNb nodes 
%Output: 
    %Dynamic_corr - contains correlation between temporal variation in between network FC and clustering; 
    %temporal data - contains the mean clustering and mean b/w network FC values for each timepoint as a reference 

            
%% Create correlation matrix within each window
Windows=147; %number of windows (based on number of timepoints and window size)

for t = 1:Windows
    window=t:t+29; %30 sliding window timepoints (60 second window)
    Corr_Matrix=corrcoef(ROI_Timecourses(window,:)); %Compute correlation matrix
    Corr_Matrix(1:length(Corr_Matrix)+1:end)=0; %remove self-connections
    Z(:,:,t)=.5.*log((1+Corr_Matrix)./(1-Corr_Matrix)); %Z transform correlations
    Z_pos=Z(:,:,t); Z_pos(Z_pos(:)<0)=0; %retain only pos weights for computing clustering
    
    %Specify subgraphs
    FPCNa_DN=Z(Ci==3,Ci==2,t); % subgraph of between network (FPCNa and DN) FC values
    FPCNb_DN=Z(Ci==4,Ci==2,t);
    FPCNa_DAN=Z(Ci==3,Ci==1,t);
    FPCNb_DAN=Z(Ci==4,Ci==1,t);
    DN_within(:,:,t)=Z_pos(Ci==2,Ci==2); %subgraph of within network FC values that will be used for clustering
    DAN_within(:,:,t)=Z_pos(Ci==1,Ci==1);
    
    % compute between network FC
    Mean_FPCNa_DN(t)=mean(reshape(FPCNa_DN,1,[])); %Vectorize subgraph and compute mean FC
    Mean_FPCNb_DN(t)=mean(reshape(FPCNb_DN,1,[]));
    Mean_FPCNa_DAN(t)=mean(reshape(FPCNa_DAN,1,[]));
    Mean_FPCNb_DAN(t)=mean(reshape(FPCNb_DAN,1,[]));
end % t (Timepoint counter)

%% compute clustering within each window

maxvalue_DAN=max(DAN_within(:)); %determine max z(r) value from set of all matrices for normalization
maxvalue_DN=max(DN_within(:));

for t=1:Windows
    DAN_nrm=DAN_within(:,:,t)./maxvalue_DAN; %normalize within-network values to [0 1]
    DN_nrm=DN_within(:,:,t)./maxvalue_DN;
    DAN_C=clustering_coef_wu(DAN_nrm); %Run BCT clustering script
    DN_C=clustering_coef_wu(DN_nrm);
    Mean_DAN_C(t)=mean(DAN_C); %Compute mean clustering across DAN nodes
    Mean_DN_C(t)=mean(DN_C);
end


%% Compute correlation between temporal variation in between-network FC and clustering

corr=corrcoef(Mean_FPCNa_DN,Mean_DN_C);Dynamic_corr(1,1)=corr(2);
corr=corrcoef(Mean_FPCNb_DN,Mean_DN_C); Dynamic_corr(1,2)=corr(2);
corr=corrcoef(Mean_FPCNa_DAN,Mean_DAN_C); Dynamic_corr(1,3)=corr(2);
corr=corrcoef(Mean_FPCNb_DAN,Mean_DAN_C); Dynamic_corr(1,4)=corr(2);

temporal_data(:,1)=Mean_FPCNa_DN;
temporal_data(:,2)=Mean_FPCNb_DN;
temporal_data(:,3)=Mean_FPCNa_DAN;
temporal_data(:,4)=Mean_FPCNb_DAN;
end
              




                   