% edit those for path handling
lead path
niicopies='/Users/NetStim/Documents/HCP1000_2/copied_nii';
root='/Users/NetStim/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/HCP1200/hcp-openaccess';
storagepath='/Users/NetStim/Documents/HCP1000_2/'; % folder to which the matrix will be written.
% storagepath='/Volumes/Backup Plus/HCP1000/result/';
backuppath='/Volumes/Backup Plus/HCP1000/';
TR=0.72;
%
try
    delete('bucket');
end
system(['ln -s ',ea_path_helper(root),' bucket']);

root='';
load('subIDs.txt')
load('runinfo.mat')
load('idx.mat');
d=load('dataset_info.mat');
runNames={'rfMRI_REST1_LR','rfMRI_REST1_RL','rfMRI_REST2_LR','rfMRI_REST2_RL'};

cd bucket
% cd hcp-openaccess/HCP_1200
cd HCP_1200

% runlength=1200;
chunk=7000;
chunk_y=7000;
dim=285903;


if ~exist([storagepath,'HCP1200.mat'],'file')
    X=zeros(1000,dim,'single'); % smaller dim is fine.
    save([storagepath,'HCP1200.mat'],'X','-v7.3');
    db=matfile([storagepath,'HCP1200.mat'],'Writable',true);
else
    db=matfile([storagepath,'HCP1200.mat'],'Writable',true);
end

if ~exist([storagepath,'nanmask.mat'],'file')
    M=zeros(1000,dim,'uint16'); % smaller dim is fine.
    save([storagepath,'nanmask.mat'],'M','-v7.3');
    nanmask=matfile([storagepath,'nanmask.mat'],'Writable',true);
else
    nanmask=matfile([storagepath,'nanmask.mat'],'Writable',true);
end

if ~exist([storagepath,'nextsub.mat'],'file') % always start from the next subject after the last one that was done. This is done in case the script breaks.
    nextsub=1;
    nextvox=1;
    nextvox_y=1;
else
    load([storagepath,'nextsub']);
    if ~exist('nextvox','var')
        nextvox=1;
    end
    if ~exist('nextvox_y','var')
        nextvox_y=1;
    end
end

backupcnt=1; % backup every 100 subjects.

for sub=nextsub:length(subIDs)
    nextsub=sub;
    save([storagepath,'nextsub.mat'],'nextsub','nextvox','nextvox_y');
    if ismember(nextsub,idx)
        continue
    end
    ea_dispercent(0,['Gathering data for sub ',num2str(sub)]);
    runs = runNames(logical(runinfo(sub,:)));
    
    
    if sum(runinfo(sub,:))==4
        
     tc = cell(1,length(runs));   
     for run=1:length(runs)
       
         nii=ea_load_untouch_nii([root,num2str(subIDs(sub)),filesep,'MNINonLinear',filesep,'Results',filesep,runs{run},filesep,runs{run},'_hp2000_clean.nii.gz']);
        
        tc{run} = zeros(length(d.dataset.vol.outidx), size(nii.img,4),'single');
        sampleLength=size(nii.img,4);
        for vol=1:sampleLength
            thisimg=nii.img(:,:,:,vol);
            tc{run}(:,vol)=thisimg(d.dataset.vol.outidx);
        end
        tc{run} = tc{run} - repmat(ea_nanmean(tc{run},2),1,size(nii.img,4));
        tc=ea_bpfilter(tc{run},TR,sampleLength);
        ea_dispercent(run/length(runs));

     end
    elseif sum(runinfo(sub,1:2))==2 && sum(runinfo(sub,:))~=4
        tc = cell(1,2);   
        for run=1:2

            nii=ea_load_untouch_nii([root,num2str(subIDs(sub)),filesep,'MNINonLinear',filesep,'Results',filesep,runNames{run},filesep,runNames{run},'_hp2000_clean.nii.gz']);
        
            tc{run} = zeros(length(d.dataset.vol.outidx), size(nii.img,4),'single');
            sampleLength=size(nii.img,4);
            for vol=1:size(nii.img,4)
                thisimg=nii.img(:,:,:,vol);
                tc{run}(:,vol)=thisimg(d.dataset.vol.outidx);
            end
            tc{run} = tc{run} - repmat(ea_nanmean(tc{run},2),1,size(nii.img,4));
            tc=ea_bpfilter(tc{run},TR,sampleLength);
            ea_dispercent(run/length(runs));
        end
    elseif sum(runinfo(sub,3:4))==2 && sum(runinfo(sub,:))~=4
        tc = cell(1,2);  
        for run=3:4
            nii=ea_load_untouch_nii([root,num2str(subIDs(sub)),filesep,'MNINonLinear',filesep,'Results',filesep,runNames{run},filesep,runNames{run},'_hp2000_clean.nii.gz']);
        
            tc{run} = zeros(length(d.dataset.vol.outidx), size(nii.img,4),'single');
            sampleLength=size(nii.img,4);
            for vol=1:size(nii.img,4)
                thisimg=nii.img(:,:,:,vol);
                tc{run}(:,vol)=thisimg(d.dataset.vol.outidx);
            end
            tc{run} = tc{run} - repmat(ea_nanmean(tc{run},2),1,size(nii.img,4));
            tc=ea_bpfilter(tc{run},TR,sampleLength);
            ea_dispercent(run/length(runs));
        end
    end
    ea_dispercent(1,'end');
    
    tc = horzcat(tc{:});
    tc=tc';
    
    % now write out matrix in chunks.
    ea_dispercent(0,['Writing out sub ',num2str(sub)]);
    from=nextvox;
    from_y=nextvox_y;
    for vox=nextvox:chunk:dim
        nextvox=vox;
        if from>dim
            from = vox;
        end
        to=from+(chunk-1);
        if to>dim
            to=dim;
        end
        for vox_y=nextvox_y:chunk_y:dim

            nextvox_y=vox_y;

            save([storagepath,'nextsub.mat'],'nextsub','nextvox','nextvox_y');
            if from_y>dim
                from_y=vox_y;
            end
            to_y=from_y+(chunk_y-1);

            if to_y>dim
                to_y=dim;
            end
            mat=corr(tc(:,from:to),tc(:,from_y:to_y),'rows','pairwise');

            nanchunk = uint16(isnan(mat));
            if sub==1
                nanmask.M(from:to,from_y:to_y)=nanchunk; % initialize matrix with first subject
            else
                nanmask.M(from:to,from_y:to_y)=uint16(sum(cat(3,nanmask.M(from:to,from_y:to_y),nanchunk),3)); % add entry.
            end

            if sub==1
                db.X(from:to,from_y:to_y)=mat; % initialize matrix with first subject
            else
                db.X(from:to,from_y:to_y)=ea_nansum(cat(3,db.X(from:to,from_y:to_y),mat),3); % add entry.
            end
            from_y=from_y+chunk;
        end
        nextvox_y=1;
        from=from+chunk;
        ea_dispercent(from/dim);

    end
    ea_dispercent(1,'end');
    nextvox=1;
    backupcnt=backupcnt+1;
%     if backupcnt==20
%         disp('Backing up data');
%         backupcnt=0;
%         copyfile([storagepath,'HCP1200.mat'],[backuppath,'HCP1200.mat']);
%         copyfile([storagepath,'nanmask.mat'],[backuppath,'nanmask.mat']);
%         copyfile([storagepath,'nextsub.mat'],[backuppath,'nextsub.mat']);
%         clear mountainduck cache
%         try
%             rmdir('/Users/NetStim/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Cache','s');
%         end
%     end


end