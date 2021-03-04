chunk=7000;
chunk_y=7000;
dim=285903;
nextvox=1;
nextvox_y=1;
from=nextvox;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
from_y=nextvox_y;
% load using -matfile- two matrices from two folders as a try and name one of them as "db" and take X
% of the other one to name it "mat"
    for vox=nextvox:chunk:dim
        nextvox=vox;
        if from>dim
%             from = vox;
            break
        end
        to=from+(chunk-1);
        if to>dim
            to=dim;
        end
        for vox_y=nextvox_y:chunk_y:dim

            nextvox_y=vox_y;

%             save([storagepath,'nextsub.y.X'],'nextsub','nextvox','nextvox_y');
            if from_y>dim
                from_y=vox_y;
            end
            to_y=from_y+(chunk_y-1);

            if to_y>dim
                to_y=dim;
            end
%             y.X=corr(tc(:,from:to),tc(:,from_y:to_y),'rows','pairwise');

%             nanchunk = uint16(isnan(y.X));
%             if sub==501
%                 nanmask.M(from:to,from_y:to_y)=nanchunk; % initialize y.Xrix with first subject
%             else
%                 nanmask.M(from:to,from_y:to_y)=uint16(sum(cat(3,nanmask.M(from:to,from_y:to_y),nanchunk),3)); % add entry.
%             end

%             if sub==501
%                 db.X(from:to,from_y:to_y)=y.X; % initialize y.Xrix with first subject
%             else
                db.X(from:to,from_y:to_y)=ea_nansum(cat(3,db.X(from:to,from_y:to_y),y.X(from:to,from_y:to_y)),3); % add entry.
%             end
            from_y=from_y+chunk;
        end
        nextvox_y=1;
        from=from+chunk;
        ea_dispercent(from/dim);

    end
