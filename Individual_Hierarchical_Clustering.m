function Mean_between_FC=Individual_Hierarchical_Clustering(W,subj)

%Function clusters FPCN nodes for each participant based on FC patterns
%with the DN and DAN, and then creates ROI masks for the two empirically derived FPCN
%subsystems

%Input:
%   W: weighted graph for a given individual
%   subj: subject number

%Index of relevant nodes within the FC matrix
DAN=[60:85 259:284];
DN=[149:187 358:384];
FPCN=[101:103 121:148 325:357 306:308]; 
DN_DAN=[DAN DN];

X=W(FPCN,DN_DAN); %subgraph of between network connections  

%cluster FPCN nodes into 2 families (or three for the subjects that only
%had 1 region separate from all others in the first division)

switch subj 
    case {'S10' 'S15' 'S16' 'S17'}
        cluster_results = clusterdata(X,'distance', 'spearman', 'linkage','average','maxclust',3);
    otherwise
        cluster_results = clusterdata(X,'distance', 'spearman', 'linkage','average','maxclust',2);
end

%get node index values for the two subsystems, arbitrarily called FPCN 1
%and FPCN 2 for now (determining whether FPCN 1 corresponded to FPCNa or
%FPCNb was done by manual inspection of each individual's data)
FPCN1=FPCN(cluster_results==1); 
FPCN2=FPCN(cluster_results==2);

%Vectorize and compute mean between network FC. 
FPCN1_DN=mean(reshape(W(FPCN1,DN),1,[]));
FPCN1_DAN=mean(reshape(W(FPCN1,DAN),1,[]));
FPCN2_DN=mean(reshape(W(FPCN2,DN),1,[]));
FPCN2_DAN=mean(reshape(W(FPCN2,DAN),1,[]));

%% Create FPCN ROI masks based on clustering results

roi_dir='\Yeo400\'; %directory holding ROIs
roi_output_dir=pwd;
cd(roi_dir);

for k = 1:2
    current_index=FPCN(cluster_results==k);
    ROI_Name=['Network' int2str(k)];
    clear roilist ROI_holder
    for n=1:length(current_index)
        ROI_holder{n}=['Yeo_400_' int2str(current_index(n)) '_roi.mat']; %Grab Yeo ROIs for current network
    end
    roilist=ROI_holder;
    
    % Combine ROIs with Marsbar
    clear func m b a c
    roilist=maroi('load_cell',roilist); % roilist contains the ROIs to be combined
    [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Combine ROIs');
    
    for c=1:length(roilist)
        a{c,:}=['r' int2str(c) '|'];
    end
    
    b=[];
    for m=1:size(a,1)
        b=[b a{m,:}];
    end
    func=b(1:end-1); 
    
    for v = 1:length(roilist)
        eval(sprintf('r%d = roilist{%d};', v, v));
    end
    
    o=eval(func);   
    o = label(o, [ROI_Name]);
    saveroi(o, fullfile(roi_output_dir, strcat(ROI_Name,'_roi.mat')));  % save as .mat
    save_as_image(o, fullfile(roi_output_dir, strcat(ROI_Name, '_roi.nii'))); % Save as image (.nii)  
end
