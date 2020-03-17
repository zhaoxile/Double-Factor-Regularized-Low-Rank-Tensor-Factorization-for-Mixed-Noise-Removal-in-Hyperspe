function [X, A, B, S, Out] = Mixed_LRTF_DFR_GPU(Y, opts)

max_it   = gpuArray(opts.max_it);  
Bmax_it  = gpuArray(opts.Bmax_it); 
tol      = gpuArray(opts.tol);              
R        = gpuArray(opts.R);
rho      = gpuArray(opts.rho);
tau      = gpuArray(opts.tau);
lambda   = gpuArray(opts.lambda);
beta     = gpuArray(opts.beta);
mu       = gpuArray(opts.mu);

Out.Res  = gpuArray([]); Out.PSNR=gpuArray([]);
Y        = gpuArray(Y);

%% Initiation
Nway     = gpuArray(size(Y));
NwayB    = gpuArray([Nway(1),Nway(2),R]);
A        = gpuArray.rand(Nway(3), R);
B        = gpuArray.rand(NwayB);
Y3       = Unfold(Y,Nway,3); Y3=gpuArray(Y3);
X        = Y;
S        = gpuArray.zeros(Nway);
%% Difference operator 
D= cell(1,3);     
for i = 1:3
diaga = gpuArray.ones(Nway(i),1);  diagb = gpuArray.ones(Nway(i)-1,1);
D{i}  = diag(-diaga)+diag(diagb,1);
D{i}(end,1) = 1;
end
%% D3 
d3         = D{3}(:,1);
deig3      = fft(d3);
SigD3tD3   = 2*lambda*(abs(deig3).^2);
%% D1 and D2
d1         = gpuArray.zeros(Nway(1),Nway(2));
d1(end,1)  = 1;  d1(1,1) = -1;
d2         = gpuArray.zeros(Nway(1),Nway(2));
d2(1,end)  = 1;  d2(1,1) = -1;
d1eig      = fft2(d1);
d2eig      = fft2(d2);
SigD12tD12 = beta*((abs(d1eig)).^2+(abs(d2eig)).^2);
SigD12tD12 = SigD12tD12(:);


for k=1: max_it
     Ak = A;   Bk = B;  Xk = X;  Sk=S;
     
    %% update A
    Bk3       = Unfold(Bk,NwayB,3);
    Sk3       = Unfold(Sk,Nway,3);  
    [U1,S1,~] = svd(Bk3*Bk3');
    Sig1      = diag(S1);
    Sig       = repmat(Sig1',Nway(3),1)+repmat(SigD3tD3,1,R)+rho;
    Sig       = 1./Sig;
    G         = (Y3-Sk3)*Bk3'+ rho*Ak;
    temp      = Sig.*(fft(G)*U1);
    A         = real(ifft(temp))*U1';
    
    %% update B
    B = GroupTV_B(Y, A, Bk, Sk, D, NwayB, SigD12tD12, rho, tau, beta, Bmax_it);
    X = my_ttm3(B,A);

    %% update S
    W = 1/(abs((Y-X+rho*Sk)/(1+rho))+eps);
    S = prox_l1((Y-X+rho*Sk)/(1+rho),W*mu/(1+rho));
    
    %%
    Res = norm(X(:)-Xk(:))/norm(Xk(:));
    Out.Res = [Out.Res,Res];
    
    if isfield(opts, 'Xtrue')
        XT=opts.Xtrue;
        psnr = PSNR3D(XT * 255, X * 255);
        Out.PSNR = [Out.PSNR,psnr];
    end
    
    if mod(k, 5) == 0
         if isfield(opts, 'Xtrue')
            fprintf('LRTF-DFR: iter = %d   PSNR = %f   res = %f \n', k, psnr, Res);
         else
            fprintf('LRTF-DFR: iter = %d   res = %f \n', k, Res);
         end 
    end
    
    %% check stopping criterion
    if Res<tol
        break
    end
  
end
end



function B = GroupTV_B(Y, A, Bk, Sk, D, NwayB, SigD12tD12, rho, tau, beta, Bmax_it)
% auxiliary variable
Z1 = gpuArray.zeros(NwayB);
Z2 = Z1;
% multiplier
P1 = Z1;
P2 = Z1;

for i=1:Bmax_it
    %% B subproblem
    [U2,S2,~]  = svd(A'*A);
    Sig2       = diag(S2);
    Sigg       = repmat(Sig2',NwayB(1)*NwayB(2),1)+repmat(SigD12tD12,1,NwayB(3))+rho;    
    Sigg       = 1./Sigg;
    K          = my_ttm3(Y-Sk,A')+beta*(my_ttm1(Z1-P1/beta,D{1}')+my_ttm2(Z2-P2/beta,D{2}'))+rho*Bk;
    K3         = Unfold(K,NwayB,3);
    temp       = Sigg.*(calF(K3',NwayB(1),NwayB(2))*U2);
    B3t        = real(calFt(temp,NwayB(1),NwayB(2)))*U2';
    B          = Fold(B3t',NwayB,3);
    
    %% Z subproblem
    Z1         = Thres_21(my_ttm1(B,D{1})+P1/beta, tau/beta);   
    Z2         = Thres_21(my_ttm2(B,D{2})+P2/beta, tau/beta);  
    %% updating P
    P1 =P1+beta*(my_ttm1(B,D{1})-Z1); 
    P2 =P2+beta*(my_ttm2(B,D{2})-Z2);  
end
end


function X=my_ttm1(B,A) %%%calculate B \times_1 A, matrix A, tensor B.
         sizeB = size(B);
         tempX = A*Unfold(B,sizeB,1);
         X     = Fold(tempX,[size(tempX,1),sizeB(2),sizeB(3)],1);
end

function X=my_ttm2(B,A) %%%calculate B \times_2 A, matrix A, tensor B.
         sizeB = size(B);
         tempX = A*Unfold(B,sizeB,2);
         X     = Fold(tempX,[sizeB(1),size(tempX,1),sizeB(3)],2);
end

function X=my_ttm3(B,A) %%%calculate B \times_3 A, matrix A, tensor B.
         sizeB = size(B);
         tempX = A*Unfold(B,sizeB,3);
         X     = Fold(tempX,[sizeB(1),sizeB(2),size(tempX,1)],3);
end



