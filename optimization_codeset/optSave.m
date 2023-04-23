function [] = optSave(prepath,name,batchtype,batchscen, ...
    batchloc,batchc)

%clearvars -except name prepath batchtype scen loc c

%unpack ARRAY TASK ID into batchloc
if isequal(batchtype,'ssm') %ssm
    loc_id = batchloc;
    switch loc_id
        case 1
            batchloc = 'argBasin';
        case 2
            batchloc = 'cosEndurance_wa';
        case 3
            batchloc = 'irmSea';
    end
end
optScript

if isequal(batchtype,'ssm') %ssm- THE SENS FUNCTION WILL CHANGE ONCE THE HYBRID VERSION IS CREATED
        save([prepath name '_' num2str(loc_id) '.mat'], ...
        'tiv','tcm','twf','cis','rsp','cos','szo', 'pmm',...
        'pvd','psr','pcm','pwf','pve','rai', ...
        'wiv','wcm','whl','ect','rhs', ...
        'giv','fco','fca','fsl','oci','gcm', ...
        'lft','dtc','osv','spv','tmt','eol','dep','bcc', ...
        'bhc','utp','ild','sdr','s0','-v7.3')

else %save single structure
    if exist('multStruct','var')
        stru.(name) = multStruct;
    elseif exist('allLocUses','var')
        stru.(name) = allLocUses;
    elseif exist('allScenUses','var')
        stru.(name) = allScenUses;
    else
        stru.(name) = optStruct;
    end
    save([prepath name '.mat'], '-struct','stru','-v7.3')
end

end

