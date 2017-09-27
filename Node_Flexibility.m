function Mean_SD=Node_Flexibility(W)
 
%Computes the standard deviation of FC across contexts (6 conditions in our data set) between each FPCN node and DN and DAN nodes
%Input: data file called W, containing correlation matrices for a given subject, where 3rd dimension is condition 
%output: Mean_SD, which is the mean variability of FC for each FPCN node

%Specify community assignments for each node in the data file
DN_DAN_index=1:17; %all DN and DAN nodes
FPCN_index=18:37;%all FPCN nodes
            
for FPCN_node=1:length(FPCN_index)
    for DN_DAN_node=1:length(DN_DAN_index)
        
        for i = 1:6; %condition
            FC_Value(i)=W(FPCN_index(FPCN_node),DN_DAN_index(DN_DAN_node),i); %extract FC value for given FPCN - DN/DAN node pair for each condition
        end 
        
        FC_SD(1,DN_DAN_node)=std(FC_Value); % compute SD of FC values across the 6 contexts, and store for every DN/DAN node
    end
    
    Mean_SD(1,FPCN_node)=mean(FC_SD); %Compute mean variability and store for each FPCN node  
end 

end


             
                  