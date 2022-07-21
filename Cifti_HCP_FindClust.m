%This code is meant to automagically extract data from the tFMRI HCP data
%First, it reads in the full sample task data to find clusters, peaks of
%activation, and automagic labels for each peak
workbench = 'C:\Users\rawls017\Documents\workbench\bin_windows64\wb_command';
infile = 'HCP_S1200_997_tfMRI_ALLTASKS_level2_cohensd_hp200_s4_MSMAll.dscalar.nii';
clus_thresh = '.8 1 .8 1';
outfile = 'CLUSTER_HCP_S1200_997_tfMRI_ALLTASKS_level2_cohensd_hp200_s4_MSMAll.dscalar.nii';
L_surface = 'S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii';
R_surface = 'S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii';
labelout = 'LABEL_CLUSTER_HCP_S1200_997_tfMRI_ALLTASKS_level2_cohensd_hp200_s4_MSMAll.dlabel.nii';
ROIout = 'ROI_LABEL_CLUSTER_HCP_S1200_997_tfMRI_ALLTASKS_level2_cohensd_hp200_s4_MSMAll.dscalar.nii';

system([ workbench ' -cifti-find-clusters ' infile  ' ' clus_thresh ' ' 'COLUMN ' outfile ...
    ' -left-surface ' L_surface ' -right-surface ' R_surface]);

[cii_orig, cii_data_orig] = ciftiopen(infile,workbench);
[cii_clus, cii_data_clus] = ciftiopen(outfile,workbench);
tmp = size(cii_data_orig);
numcond = tmp(2);
d = [];
for i = 1:numcond
    tmp = cii_data_orig(:,i);
    tmp2 = cii_data_clus(:,i);
    tmp3 = unique(tmp2(tmp2>0));
    d = [d length(tmp3)];
end

numclustmax = max(d);

ClusNumOut = zeros(numclustmax,numcond);
ClusPeakVal = zeros(numclustmax,numcond);
ClusPeakVox = zeros(numclustmax,numcond);
ClusSize = zeros(numclustmax,numcond);

for i = 1:numcond
    tmp = cii_data_orig(:,i);
    tmp2 = cii_data_clus(:,i);
    tmp3 = unique(tmp2(tmp2>0));
    for j = 1:length(tmp3)
        %tmp4 = tmp(find(tmp2==tmp3(j)));
        C = find(tmp2==tmp3(j));
        [a,b] = max(tmp(C));
        
        ClusNumOut(j,i) = tmp3(j);
        ClusPeakVal(j,i) = a;
        ClusPeakVox(j,i) = C(b);
        ClusSize(j,i) = length(C);
    end
end 
%Now find coordinates of each peak voxel
ClusLoxX = zeros(size(ClusNumOut));
ClusLoxY = zeros(size(ClusNumOut));
ClusLoxZ = zeros(size(ClusNumOut));
PeakSpotLabel = cell(size(ClusNumOut));

system([ workbench ' -cifti-export-dense-mapping ' infile ' COLUMN' ' -volume-all' ' -structure vol_map_wLabel_S1200.txt']);
system([ workbench ' -cifti-export-dense-mapping ' infile ' COLUMN' ' -volume-all' ' vol_map_S1200.txt']);
system([ workbench ' -cifti-export-dense-mapping ' infile ' COLUMN' ' -surface' ' CORTEX_LEFT' ' left_cortex_map_S1200.txt']);
system([ workbench ' -cifti-export-dense-mapping ' infile ' COLUMN' ' -surface' ' CORTEX_RIGHT' ' right_cortex_map_S1200.txt']);

subcort_lookup = importdata('vol_map_S1200.txt',' ');
subcort_label = importdata('vol_map_wLabel_S1200.txt',' ');
leftcort_lookup = importdata('left_cortex_map_S1200.txt',' ');
rightcort_lookup = importdata('right_cortex_map_S1200.txt',' ');

L_g = gifti(L_surface);
R_g = gifti(R_surface);

system([ workbench ' -cifti-parcellate ' infile ' Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii COLUMN' ' parcellated_test.pscalar.nii']);
cortcomp = ft_read_cifti('parcellated_test.pscalar.nii');
cortcomp_vox2parc = cortcomp.brainordinate.parcellation(find(cortcomp.brainordinate.parcellation>0));

subcomp = ft_read_cifti(infile);
subcomp = subcomp.pos(64985:96854,:);

addpath('C:\Users\rawls017\Documents\MATLAB\spm12\spm12');
for i = 1:numcond
    dat = ClusPeakVox(:,i);
    dat  =dat(dat>0);
    for j = 1:length(dat)
        if dat(j)>59411
            index=find(subcort_lookup(:,1)==dat(j));
            ClusLoxX(j,i)=subcomp(index,1);
            ClusLoxY(j,i)=subcomp(index,2);
            ClusLoxZ(j,i)=subcomp(index,3);
            PeakSpotLabel{j,i}=subcort_label.textdata(index,2);
        elseif dat(j)>29695
            index = find(rightcort_lookup(:,1)==dat(j));
            truindex = rightcort_lookup(index,2)+1;
            ClusLoxX(j,i)=R_g.vertices(truindex,1);
            ClusLoxY(j,i)=R_g.vertices(truindex,2);
            ClusLoxZ(j,i)=R_g.vertices(truindex,3);
            vox_parc = cortcomp_vox2parc(dat(j));
            PeakSpotLabel{j,i}=cortcomp.brainordinate.parcellationlabel(vox_parc); %spm_atlas('query','neuromorphometrics',[ClusLoxX(j,i),ClusLoxY(j,i),ClusLoxZ(j,i)]');   
        else
            index = find(leftcort_lookup(:,1)==dat(j));
            truindex = leftcort_lookup(index,2)+1;
            ClusLoxX(j,i)=L_g.vertices(truindex,1);
            ClusLoxY(j,i)=L_g.vertices(truindex,2);
            ClusLoxZ(j,i)=L_g.vertices(truindex,3);
            vox_parc = cortcomp_vox2parc(dat(j));
            PeakSpotLabel{j,i}=cortcomp.brainordinate.parcellationlabel(vox_parc); %spm_atlas('query','neuromorphometrics',[ClusLoxX(j,i),ClusLoxY(j,i),ClusLoxZ(j,i)]');   
        end
    end
end
rmpath('C:\Users\rawls017\Documents\MATLAB\spm12\spm12');

GroupClustInfo.ClustLabel = PeakSpotLabel;
GroupClustInfo.Xlocs = ClusLoxX;
GroupClustInfo.Ylocs = ClusLoxY;
GroupClustInfo.Zlocs = ClusLoxZ;
GroupClustInfo.PeakVoxel = ClusPeakVox;
save('Group_Clust_taskdata.m','GroupClustInfo');

%Meaningful task contrasts - 
%11 - WM 2bk-0bk
%%%%14 - WM 0bk-2bk
%31 - gamble-punish
%32 - gamble-reward
%%%%33 - Gamble punish-reward
%%%%36 - Gamble reward-punish
%%%%65 - language math-story
%66 - language story-math
%%%%71 - social, random-TOM
%74 - social, TOM-random
%%%%77 - relational, match-rel
%%%%78 - relational, rel-match
%83 - emotion, face-shape
%%%%86 - emotion, shape-face
GroupRedClust.ClustLabel = GroupClustInfo.ClustLabel(:,[11,31,32,66,74,83]);
GroupRedClust.Xlocs = GroupClustInfo.Xlocs(:,[11,31,32,66,74,83]);
GroupRedClust.Ylocs = GroupClustInfo.Ylocs(:,[11,31,32,66,74,83]);
GroupRedClust.Zlocs = GroupClustInfo.Zlocs(:,[11,31,32,66,74,83]);
GroupRedClust.PeakVoxel = GroupClustInfo.PeakVoxel(:,[11,31,32,66,74,83]);

count=0;
[tmp1,tmp2]=size(GroupRedClust.ClustLabel);
for i=1:tmp1
    for j=1:tmp2
        if strcmp(class(GroupRedClust.ClustLabel{i,j}),'cell')
        count=count+1;
        end
    end
end
disp(['Data contains ' num2str(count) ' clusters of activation!']);
%Cool! We found and labeled our clusters! Now we need to extract the
%relevant info from each cluster for our subjects!
SubjectIDs = [105014];
for sub = 1:length(SubjectIDs)
    %Contrast list:
    %WM - 1:30
    %Gambling - 31:36
    %Motor - 37:62
    %Language - 63:68
    %Social - 69:74
    %Relational - 75:80
    %Emotion - 81:86
    
    SubjMeanVal=zeros(numclustmax,numcond);
    SubjPeakVal=zeros(numclustmax,numcond);
    SubjPeakVox=zeros(numclustmax,numcond);

    SingleSubData = [];
    
    for i = 1:7
        %load subject
        if i==1
            TaskName = 'WM';
        elseif i==2
            TaskName = 'GAMBLING';
        elseif i==3
            TaskName = 'MOTOR';
        elseif i==4
            TaskName = 'LANGUAGE';
        elseif i==5
            TaskName = 'SOCIAL';
        elseif i==6
            TaskName = 'RELATIONAL';
        elseif i==7
            TaskName = 'EMOTION';
        end
        subfile = ['/home/annaz/rawls017/Documents/MATLAB/TaskData/ZIP/' Taskname '/' num2str(SubjectIDs(sub)) '_3T_tfMRI_' Taskname '_analysis_s4'];
            
        C:/Users/rawls017/Documents/MATLAB/TaskData/' TaskName '/' num2str(SubjectIDs(sub)) '_3T_tfMRI_' TaskName ...
            '_analysis_s4/' num2str(SubjectIDs(sub)) '/MNINonLinear/Results/tfMRI_' ...
            TaskName '/tfMRI_' TaskName '_hp200_s4_level2_MSMAll.feat/' num2str(SubjectIDs(sub)) '_tfMRI_' TaskName '_level2_hp200_s4_MSMAll.dscalar.nii'];
        
        
        [cii_subj, cii_data_subj] = ciftiopen(subfile,workbench);
        SingleSubData = [SingleSubData cii_data_subj];
    end
    for i = 1:numcond 
        tmp = SingleSubData(:,i);
        tmp2 = cii_data_clus(:,i);
        tmp3 = unique(tmp2(tmp2>0));
        for j = 1:length(tmp3)
            %tmp4 = tmp(find(tmp2==tmp3(j)));
            C = find(tmp2==tmp3(j));
            %extract peak value within cluster
            [a,b] = max(tmp(C));
            %extract mean value
            c = mean(tmp(C));
            
            SubjPeakVox(j,i) = b;
            SubjPeakVal(j,i) = a;
            SubjMeanVal(j,i) = c;
        end
    end
    Sub_Data_Out.PeakVal = SubjPeakVal;
    Sub_Data_Out.MeanVal = SubjMeanVal;
    Sub_Data_Out.PeakVox = SubjPeakVox;
    save([num2str(subjectID) '_taskdata_extracted.m'],'Sub_Data_Out');
end
