function [Node_SI mean_SI_FPCNa mean_SI_FPCNb mean_FC_with_DAN mean_FC_with_DN] = Selectivity_Index(W,Ci)

% Computes the mean selectivity (differential coupling with the DN and DAN)
% for FPCNa and FPCNb nodes. 

%The SI is computed as mean FC with DN nodes - mean FC with DAN nodes for each FPCN node, and then averaged across
%relevant nodes. 

%Input: Ci; a vector consisisting of community assignments for each node
%       W; Weighted graph for a given participant and condition

%Node index based on Ci
DAN=[find(Ci==1) find(Ci==2)];
DN=[find(Ci==3) find(Ci==4)];
FPCNa=[find(Ci==5) find(Ci==6)];
FPCNb=[find(Ci==7) find(Ci==8)];
FPCN=[FPCNa FPCNb];

Node_SI=[];

for i = 1:length(FPCN)
    current_FPCN_ROI=FPCN(i);
    mean_FC_with_DAN(1,i)=mean(reshape(W(current_FPCN_ROI,DAN),1,[])); %Compute mean FC with DAN for each FPCN node
    mean_FC_with_DN(1,i)=mean(reshape(W(current_FPCN_ROI,DN),1,[])); %Compute mean FC with DN for each FPCN node
    Node_SI(1,i)=mean_FC_with_DN(1,i)-mean_FC_with_DAN(1,i); % Compute difference score (selectivity index) for node
end

mean_SI_FPCNa=mean(Node_SI(1,1:length(FPCNa))); %compute mean selectivity index across FPCNa nodes
mean_SI_FPCNb=mean(Node_SI(1,length(FPCNa)+1:end)); %compute mean selectivity index across FPCNb nodes

end