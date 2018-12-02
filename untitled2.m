
clear
clc

ff = './';

for j=50
    
    if j==0        
        file = [ff,'realaxi-den-init'];
    else
        file =[ff,'realaxi-tmp_',num2str(j)];
    end
    
    [r,z,den]=loadrz(file);
    figure
    mesh(r,z,den)

end