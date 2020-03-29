if 0
    %%
    K_kf2 = zeros(size(K_kf,1)*size(K_kf,2),size(K_kf,3));
    for p = 1 : size(K_kf,3)
        K_ = K_kf(:,:,p);
        K_kf2(:,p) = K_(:);
    end
    
    %%
    for p = 1 : size(K_kf,3)
        eval(strcat(['Kstruct.K',num2str(p) ,'= K_kf(:,:,p);']))
    end
end
% 
% for i = 1 : 10
%     K_ = K_kf(:,:,i);
% end

dim1 = size(K_kf,1);
dim2 = size(K_kf,2);
%%
tic
for i = 1 : 30
    K_ = reshape(K_kf2(:,i), [dim1, dim2]);
end
toc

%%
tic
for i = 1 : 30
    eval(strcat(['K_ = Kstruct.K', num2str(p),';']))
end
toc

%%
tic

for i = 1 : 30
    K_ = Kstruct.(fns{i});
    if (i == 15)
        fprintf('best one?\n')
    end
end
toc
%% 
tic
for i = 1 : 30
    switch randi([1 11])
        case 1 
            K_ = K1;
        case 2 
            K_ = K2;
        case 3 
            K_ = K3;
        case 4
            K_ = K4;
        case 5
            K_ = K5;
        case 6
            K_ = K6;
        case 7
            K  = K7;
        case 8
            K_ = K8;
        case 9 
            K_ = K9;
        case 10 
            K_ = K10;
        otherwise
            K_ = K20;
    end
end
toc


%% 
tic
for i = 1 : 30
    switch randi([1 11])
        case 1 
            K_ = Kstruct.K1;
        case 2 
            K_ = Kstruct.K2;
        case 3 
            K_ = Kstruct.K3;
        case 4
            K_ = Kstruct.K4;
        case 5
            K_ = Kstruct.K5;
        case 6
            K_ = Kstruct.K6;
        case 7
            K  = Kstruct.K7;
        case 8
            K_ = Kstruct.K8;
        case 9 
            K_ = Kstruct.K9;
        case 10 
            K_ = Kstruct.K10;
        otherwise
            K_ = Kstruct.K20;
    end
end
toc