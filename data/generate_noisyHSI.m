load('simu_indian.mat')
Ohsi=simu_indian;
load('indian_nlevel.mat')
load('band_s_d.mat')
%%

if max(Ohsi(:))>1
    Ohsi=my_normalized(Ohsi);
end

Nway = size(Ohsi);


%% noise_case1: G:0.1-0.2
Nhsi = zeros(Nway);
for j=1:Nway(3)
    Nhsi(:,:,j) = Ohsi(:,:,j)+sigma_n3(j)*randn(Nway(1),Nway(2));
end
save('indian_case1','Ohsi','Nhsi'); % save results  

%% noise_case2 G:0.1-0.2 S;0.1-0.2
for j=1:Nway(3)
    Nhsi(:,:,j) = imnoise(Ohsi(:,:,j),'salt & pepper',p_n3(j))+sigma_n3(j)*randn(Nway(1),Nway(2));
end
save('indian_case2','Ohsi','Nhsi'); % save results  


%% noise_case3 G:0.1-0.2 S;0.1-0.2 stripes: 40%的band, 每个band为6-15条 
for j=1:Nway(3)
    Nhsi(:,:,j) = imnoise(Ohsi(:,:,j),'salt & pepper',p_n3(j))+sigma_n3(j)*randn(Nway(1),Nway(2));
end

for i = 1:length(add_band_s)
    stripenum = randperm(10,1)+5;
    locolumn    = randperm(Nway(2),stripenum);
    Nhsi(:,locolumn,add_band_s(i))=0.2*rand(1)+0.6;
end
save('indian_case3','Ohsi','Nhsi'); % save results  



%% noise_case4 G:0.1-0.2 S;0.1-0.2 deadline: 20%的band, 每个band为6-10条， 宽度为1-3
for j=1:Nway(3)
    Nhsi(:,:,j) = imnoise(Ohsi(:,:,j),'salt & pepper',p_n3(j))+sigma_n3(j)*randn(Nway(1),Nway(2));
end

for i = 1:length(add_band_d)
    deadlinenum = randperm(5,1)+5;
    locolumn    = randperm(Nway(2)-2,deadlinenum);
    an          = funrand(3,deadlinenum);
    loc1=find(an==1);loc2=find(an==2);loc3=find(an==3);
    Nhsi(:,locolumn(loc1),add_band_d(i))=0; 
    Nhsi(:,locolumn(loc2),add_band_d(i))=0;
    Nhsi(:,locolumn(loc2)+1,add_band_d(i))=0;
    Nhsi(:,locolumn(loc3),add_band_d(i))=0;
    Nhsi(:,locolumn(loc3)+1,add_band_d(i))=0; 
    Nhsi(:,locolumn(loc3)+2,add_band_d(i))=0; 
end
save('indian_case4','Ohsi','Nhsi'); % save results  



%% noise_case5 G:0.1-0.2 S;0.1-0.2  stripes: 40%的band, 每个band为6-15条  deadline: 20%的band, 每个band为5-10条， 宽度为1-3
for j=1:Nway(3)
    Nhsi(:,:,j) = imnoise(Ohsi(:,:,j),'salt & pepper',p_n3(j))+sigma_n3(j)*randn(Nway(1),Nway(2));
end

for i = 1:length(add_band_s)
    stripenum = randperm(10,1)+5;
    locolumn    = randperm(Nway(2),stripenum);
    Nhsi(:,locolumn,add_band_s(i))=0.2*rand(1)+0.6;
end

for i = 1:length(add_band_d)
    deadlinenum = randperm(5,1)+5;
    locolumn    = randperm(Nway(2)-2,deadlinenum);
    an          = funrand(3,deadlinenum);
    loc1=find(an==1);loc2=find(an==2);loc3=find(an==3);
    Nhsi(:,locolumn(loc1),add_band_d(i))=0; 
    Nhsi(:,locolumn(loc2),add_band_d(i))=0;
    Nhsi(:,locolumn(loc2)+1,add_band_d(i))=0;
    Nhsi(:,locolumn(loc3),add_band_d(i))=0;
    Nhsi(:,locolumn(loc3)+1,add_band_d(i))=0; 
    Nhsi(:,locolumn(loc3)+2,add_band_d(i))=0; 
end
save('indian_case5','Ohsi','Nhsi'); % save results  
    