
clear
clc

ff = './realaxi-tmp_';

for j=1

    [r,z,den]=loadrz([ff,num2str(j)]);
    imagesc(r,z,den)
    
end