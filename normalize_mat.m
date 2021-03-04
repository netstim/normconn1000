% Small script to normalize summed up HCP matrix to make diagonal entries
% == 1
lead path
out=matfile('HCP1200_R.mat','Writable',true);
chunk=5000;


numel=285903;

ea_dispercent(0,'iterating');

for vox=1:chunk:numel
    if (vox+chunk-1)>numel % last iter
        out.X(:,vox:numel)=(out.X(:,vox:numel)./repmat(nanmax(out.X(:,vox:numel)),numel,1));
        ea_dispercent(1,'end');
    else
        out.X(:,vox:(vox+chunk-1))=(out.X(:,vox:(vox+chunk-1))./repmat(nanmax(out.X(:,vox:(vox+chunk-1))),numel,1));
        ea_dispercent((vox+chunk-1)/numel);
    end
end
