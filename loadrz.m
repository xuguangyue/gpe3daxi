function [r,z,var]=loadrz(name)

data = dlmread([name,'.txt']);
r  = unique(data(:,1));
Nr = length(r);
z  = unique(data(:,2));
Nz = length(z);
[~,ss]=size(data);
if ss==4
    var=reshape(data(:,ss-1)+1i*data(:,ss),Nz,Nr);
else
    if ss==3 
        var=reshape(data(:,ss),Nz,Nr);
    end
end

end