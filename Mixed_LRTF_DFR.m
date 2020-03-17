function [X, A, B, S, Out] = Mixed_LRTF_DFR(Y, opts)

max_it   = opts.max_it;  
Bmax_it  = opts.Bmax_it; 
tol      = opts.tol;              
R        = opts.R;
rho      = opts.rho;
tau      = opts.tau;
lambda   = opts.lambda;
beta     = opts.beta;
mu       = opts.mu;

Out.Res=[]; Out.PSNR=[];

%% Initiation
Nway     = size(Y);
NwayB    = [Nway(1),Nway(2),R];
A        = rand(Nway(3), R);
B        = rand(NwayB);
Y3       = Unfold(Y,Nway,3);
X        = Y;
S        = zeros(Nway);
%% Difference operator 
D= cell(1,3);     
for i = 1:3
diaga = ones(Nway(i),1);  diagb = ones(Nway(i)-1,1);
D{i}  = diag(-diaga)+diag(diagb,1);
D{i}(end,1) = 1;
end
%% D3 
d3         = D{3}(:,1);
deig3      = fft(d3);
SigD3tD3   = 2*lambda*(abs(deig3).^2);
%% D1 and D2
d1         = zeros(Nway(1),Nway(2));
d1(end,1)  = 1;  d1(1,1) = -1;
d2         = zeros(Nway(1),Nway(2));
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
    X = double(ttm(tensor(B),A,3));

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
            fprintf('LRTF-DFR: iter = %d   res = %f\n', k, Res);
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
Z1 = zeros(NwayB);
Z2 = zeros(NwayB);
% multiplier
P1=zeros(NwayB);
P2=zeros(NwayB);

for i=1:Bmax_it
    %% B subproblem
    [U2,S2,~]  = svd(A'*A);
    Sig2       = diag(S2);
    Sigg       = repmat(Sig2',NwayB(1)*NwayB(2),1)+repmat(SigD12tD12,1,NwayB(3))+rho;
    Sigg       = 1./Sigg;
    K          = double(ttm(tensor(Y-Sk),A',3))+beta*(double(ttm(tensor(Z1-P1/beta),D{1}',1))+double(ttm(tensor(Z2-P2/beta),D{2}',2)))+rho*Bk;
    K3         = Unfold(K,NwayB,3);
    temp       = Sigg.*(calF(K3',NwayB(1),NwayB(2))*U2);
    B3t        = real(calFt(temp,NwayB(1),NwayB(2)))*U2';
    B          = Fold(B3t',NwayB,3);
    
    %% Z subproblem
    Z1         = Thres_21(double(ttm(tensor(B),D{1},1))+P1/beta, tau/beta);
    Z2         = Thres_21(double(ttm(tensor(B),D{2},2))+P2/beta, tau/beta);
    %% updating P
    P1 =P1+beta*(double(ttm(tensor(B),D{1},1))-Z1);
    P2 =P2+beta*(double(ttm(tensor(B),D{2},2))-Z2);
end
end


