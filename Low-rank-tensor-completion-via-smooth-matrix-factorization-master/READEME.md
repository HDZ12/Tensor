# Framelet_X3
函数首先初始化一些参数和变量，然后进入一个最大迭代次数为`maxit`的循环。在每次迭代中，函数首先解决`X`，然后解决`Z`，最后更新`Theta`。如果误差小于设定的阈值`tol`，则跳出循环。
## [D,R]=GenerateFrameletFilter(frame);
`GenerateFrameletFilter`函数用于生成帧滤波器的分解和重构系数。参数`frame`决定了所使用的滤波器类型。例如，当`frame`等于1时，使用的是Haar小波；当`frame`等于2时，使用的是分段线性帧滤波器；当`frame`等于3时，使用的是分段立方帧滤波器。

具体来说，`D`和`R`分别代表分解和重构系数。这些系数用于帧滤波器的分解和重构过程，这是信号处理和图像处理中的常见步骤。在这些过程中，原始信号或图像被转换为一种可以更容易分析和处理的形式。
```Matlab
function [D,R]=GenerateFrameletFilter(frame)

% function [D,R]=GenerateFrameletFilter(frame)
%
% This function generate the Decomposition and Reconstruction
% (D and R respectively) coefficients of the framelet filters
% The available filters are:
% frame=1 : Haar wavelet
% frame=2 : Piecewise Linear Framelet
% frame=3 : Piecewise Cubic Framelet
%
% Written by Jian-Feng Cai.
% email: tslcaij@nus.edu.sg
%
% See functions FraDec FraRec FraDecMultilevel FraRecMultiLevel

if frame==1          %Haar Wavelet
    D{1}=[1 1 1]/2;
    D{2}=[0 1 -1]/2;
    D{3}='cc';
    R{1}=[1 1 0]/2;
    R{2}=[-1 1 0]/2;
    R{3}='cc';
elseif frame==2      %Piecewise Linear Framelet
    D{1}=[1 2 1]/4;
    D{2}=[1 0 -1]/4*sqrt(2);
    D{3}=[-1 2 -1]/4;
    D{4}='ccc';
    R{1}=[1 2 1]/4;
    R{2}=[-1 0 1]/4*sqrt(2);
    R{3}=[-1 2 -1]/4;
    R{4}='ccc';
elseif frame==3      %Piecewise Cubic Framelet
    D{1}=[1 4 6 4 1]/16;
    D{2}=[1 2 0 -2 -1]/8;
    D{3}=[-1 0 2 0 -1]/16*sqrt(6);
    D{4}=[-1 2 0 -2 1]/8;
    D{5}=[1 -4 6 -4 1]/16;
    D{6}='ccccc';
    R{1}=[1 4 6 4 1]/16;
    R{2}=[-1 -2 0 2 1]/8;
    R{3}=[-1 0 2 0 -1]/16*sqrt(6);
    R{4}=[1 -2 0 2 -1]/8;
    R{5}=[1 -4 6 -4 1]/16;
    R{6}='ccccc';
end
```
## Z=FraDecMultiLevel(X_transition,D,Level)
`FraDecMultiLevel`函数用于实现多级帧分解。参数`X_transition`是需要被分解的数据，`D`是分解滤波器，`Level`是分解的级别。函数的返回值`Z`是帧分解的系数。

具体来说，`FraDecMultiLevel`函数首先获取`D`的长度，然后对`A`进行分解，得到`kDec`。接着，函数进入一个循环，在每次循环中，函数调用`FraDec`函数对`kDec`进行分解，并将结果存储在`Dec`中，然后更新`kDec`为`Dec`的第一个元素。这个过程会重复`Level`次
```matlab
unction DecVal=FraDecMultiLevel(A,D,L,idx) %#codegen

nD=length(D);

kDec=A;

coder.varsize('kDec');

Dec = cell(1,L);

for k = 1:L

temp = FraDec(kDec,D,k);

Dec{k} = temp;

kDec = temp{1,1};

end

assert(idx<=L);

DecVal = Dec{idx};
```
## summary
函数的主要步骤如下：

1. **初始化参数**：从`opts2`中获取`mu`、`beta`和`maxit`等参数，并设置一些初始值。
    
2. **生成Framelet滤波器**：调用`GenerateFrameletFilter`函数生成Framelet滤波器。
    
3. **主循环**：在`maxit`次迭代中，执行以下步骤：
    
    - **求解X**：计算临时变量`temp`和`C`，然后使用这些变量更新`X_out`。
    - **求解Z**：更新`X_transition`，然后使用`FraDecMultiLevel`函数和`wthresh`函数更新`Z`。
    - **更新Theta**：计算误差`error`，并根据`C`和`Z`的差值更新`Theta`。

如果误差`error`小于阈值`tol`，则跳出循环。最后，函数返回`X_out`。
```matlab
function X_out=Framelet_X3(A,Y,X_k,rho,opts2,x_size)

mu = opts2.mu;

beta = opts2.beta;

maxit = opts2.F_it;

frame = 1;

Level = 1;

wLevel = 1/2;

tol = 1e-6;

[D,R]=GenerateFrameletFilter(frame);

nD=length(D);

m=x_size(1);

n=x_size(2);

normg=norm(Y,'fro');

e=size(A,2);

X_way=[m,n,e];

X_transition=zeros(m,n*e);

X_out=zeros(e,m*n);

Z=FraDecMultiLevel(X_transition,D,Level);

Theta=Z;

for nstep=1:maxit

%% solve X

temp=mu*A'*Y;

C=CoeffOper('-',Z,Theta);

tempC=Unfold(Fold(FraRecMultiLevel(C,R,Level),X_way,1),X_way,3);

X_out=pinv(mu*A'*A+(rho+beta)*eye(size(A'*A)))*(rho*X_k+temp+beta*tempC); %#ok<MHERM>

%% solve Z

X_transition=Unfold(Fold(X_out,X_way,3),X_way,1);

C=FraDecMultiLevel(X_transition,D,Level);

Thresh=1/beta;

for ki=1:Level

for ji=1:nD-1

for jj=1:nD-1

Z{ki}{ji,jj}=wthresh(C{ki}{ji,jj}+Theta{ki}{ji,jj},'s',Thresh);

end

end

if wLevel<=0

Thresh=Thresh*norm(D{1});

else

Thresh=Thresh*wLevel;

end

end

%% update Theta

error=0;

for ki=1:Level

for ji=1:nD-1

for jj=1:nD-1

if ((ji~=1)||(jj~=1))||(ki==Level)

deltab=C{ki}{ji,jj}-Z{ki}{ji,jj};

error=error+norm(deltab,'fro')^2;

Theta{ki}{ji,jj}=Theta{ki}{ji,jj}+deltab;

end

end

end

end

error=sqrt(error)/normg;

if error<tol

break;

end

end
```
# TV_A3
函数的主要步骤如下：

1. **初始化参数**：从`opts5`中获取`mu`、`beta`和`maxit`等参数，并设置一些初始值。
    
2. **计算SVD**：对`X*X'`进行奇异值分解，并计算`Sig1`。
    
3. **生成D矩阵**：生成一个特殊的矩阵`D`，并计算`Sig2`。
    
4. **计算Sig**：根据`Sig1`、`Sig2`和`rho`计算`Sig`。
    
5. **初始化Theta和W**：初始化两个零矩阵`Theta`和`W`。
    
6. **主循环**：在`maxit`次迭代中，执行以下步骤：
    
    - **求解A**：计算临时变量`M`和`temp`，然后使用这些变量更新`A_out`。
    - **求解W**：使用`wthresh`函数更新`W`。
    - **更新Theta**：根据`D*A_out`和`W`的差值更新`Theta`。
```matlab
function A_out=TV_A3(A_p,Y,X,rho,opts5)

if isfield(opts5,'mu'); mu = opts5.mu; else mu = 10; end

if isfield(opts5,'beta'); beta= opts5.beta; else beta = 1000; end

if isfield(opts5,'F_it'); maxit = opts5.F_it; else maxit = 15; end

m=size(A_p,1);

n=size(A_p,2);

[U1,S1,~]=svd(X*X');

Sig1=mu*diag(S1);

diaga=ones(m,1);diagb=ones(m-1,1);

D=diag(-diaga)+diag(diagb,1);

D(end,1)=1;

d=D(:,1);

deig=fft(d);

Sig2=beta*(abs(deig).^2);

Sig=repmat(Sig1',m,1)+repmat(Sig2,1,n)+rho;

Sig=1./Sig;

Theta=zeros(m,n);

W=zeros(m,n);

for i=1:maxit

%% A subproblem

M=mu*Y*X'+beta*D'*(W-Theta)+rho*A_p;

temp=Sig.*(fft(M)*U1);

A_out=real(ifft(temp))*U1';

%% W subproblem

Thresh=1/beta;

W=wthresh(D*A_out+Theta,'s',Thresh);

%% updating Theta

Theta=Theta+D*A_out-W;

end

end
```
# inc_SMF_LRTC
函数的主要步骤如下：

1. **初始化参数**：从`opts`、`opts2`和`opts5`中获取各种参数，并设置一些初始值。
    
2. **初始化张量和矩阵**：对输入张量进行展开，生成初始的`Y`、`A`和`X`。
    
3. **初始化秩增加方案**：计算`YY1`、`YY2`和`YY3`，并计算初始的秩残差。
    
4. **主循环**：在`maxit`次迭代中，执行以下步骤：
    
    - **更新X**：根据`A_p`、`Y_p`和`X_p`更新`X`。
    - **更新A**：根据`X`、`Y_p`和`A_p`更新`A`。
    - **更新Y**：根据`X`、`A`和`Y_p`更新`Y`，然后将`Y_tensor`更新为`Y1`、`Y2`和`Y3`的加权和。
    - **更新秩**：计算新的秩残差，并根据条件判断是否需要增加秩。
    - **检查停止准则**：如果残差小于阈值，则跳出循环
```matlab
function [Y_tensor, A, X, Out]= inc_SMF_LRTC(Y_tensorP, known, opts, opts2, opts5)

maxit = opts.maxit;

tol = opts.tol;

alpha = opts.alpha;

rho1 = opts.rho1;

rho2 = opts.rho2;

rho3 = opts.rho3;

rhoX = opts2.rho;

rhoA = opts5.rho;

max_rank = opts.max_rank;

R = opts.R;

rank_inc = 5;

Out.Res=[];

%% Initiation

Y_tensor0 = Y_tensorP;

Nway = size(Y_tensorP);

coNway = zeros(1,3);

Y0 = cell(1,3);

A0 = cell(1,3);

X0 = cell(1,3);

for n = 1:3

coNway(n) = prod(Nway)/Nway(n);

end

for i = 1:3

Y0{i} = Unfold(Y_tensor0,Nway,i);

Y0{i} = Y0{i}';

X0{i} = rand(coNway(i), R(i));

A0{i} = rand(R(i),Nway(i));

end

Y_p=Y0;X_p=X0;A_p=A0;

%% Initiation for rank increasing scheme

YY1=Fold(A0{1}'*X0{1}',Nway,1); YY1(known)=Y_tensor0(known);

YY2=Fold(A0{2}'*X0{2}',Nway,2); YY2(known)=Y_tensor0(known);

YY3=Fold(A0{3}'*X0{3}',Nway,3); YY3(known)=Y_tensor0(known);

resrank1=norm(Y_tensor0(:)-YY1(:));

resrank2=norm(Y_tensor0(:)-YY2(:));

resrank3=norm(Y_tensor0(:)-YY3(:));

%%

for k=1: maxit

%% update X

temp=A_p{1}*A_p{1}'+rho1*(eye(size(A_p{1}*A_p{1}')));

X{1}=(Y_p{1}*A_p{1}'+rho1*X_p{1})*pinv(temp);

temp=A_p{2}*A_p{2}'+rho2*(eye(size(A_p{2}*A_p{2}')));

X{2}=(Y_p{2}*A_p{2}'+rho2*X_p{2})*pinv(temp);

x_size=Nway(1:2);

XX=Framelet_X3(A_p{3}',Y_p{3}',X_p{3}',rhoX,opts2,x_size);

X{3}=XX';

%% update A

temp=X{1}'*X{1}+rho1*eye(size(X{1}'*X{1}));

A{1}= pinv(temp)*(X{1}'*Y_p{1}+rho1*A_p{1});

temp=X{2}'*X{2}+rho2*eye(size(X{2}'*X{2}));

A{2}= pinv(temp)*(X{2}'*Y_p{2}+rho2*A_p{2});

AA=TV_A3(A_p{3}',Y_p{3}',X{3}',rhoA,opts5);

A{3}=AA';

%% update Y

Y{1} = (X{1}*A{1}+rho3*Y_p{1})/(1+rho1); Y1 = Fold(Y{1}', Nway, 1);

Y{2} = (X{2}*A{2}+rho3*Y_p{2})/(1+rho2); Y2 = Fold(Y{2}', Nway, 2);

Y{3} = (X{3}*A{3}+rho3*Y_p{3})/(1+rho3); Y3 = Fold(Y{3}', Nway, 3);

Y_tensor = alpha(1)*Y1+alpha(2)*Y2+alpha(3)*Y3;

Y_tensor(known) = Y_tensor0(known);

Res = norm(Y_tensor(:)-Y_tensorP(:))/norm(Y_tensorP(:));

Out.Res = [Out.Res,Res];

Y_tensorP = Y_tensor;

%% update Rank

resold1=resrank1;

resold2=resrank2;

resold3=resrank3;

YY1=Fold(A{1}'*X{1}',Nway,1); YY1(known)=Y_tensor0(known);

YY2=Fold(A{2}'*X{2}',Nway,2); YY2(known)=Y_tensor0(known);

YY3=Fold(A{3}'*X{3}',Nway,3); YY3(known)=Y_tensor0(known);

resrank1=norm(Y_tensor0(:)-YY1(:));

resrank2=norm(Y_tensor0(:)-YY2(:));

resrank3=norm(Y_tensor0(:)-YY3(:));

ifrank1 = abs(1-resrank1/resold1);

ifrank2 = abs(1-resrank2/resold2);

ifrank3 = abs(1-resrank3/resold3);

nowrank=[size(A{1},1),size(A{2},1),size(A{3},1)];

fprintf('Iteration: %i ',k);

fprintf('nowrank:');

fprintf(' %.0f',nowrank);

fprintf(' ');

fprintf('RelCha: %.6f ',Res);

fprintf('\n');

if ifrank1<0.01 && nowrank(1)<max_rank(1)

[A{1},X{1}]=rank_inc_adaptive(A{1},X{1},rank_inc);

end

if ifrank2<0.01 && nowrank(2)<max_rank(2)

[A{2},X{2}]=rank_inc_adaptive(A{2},X{2},rank_inc);

end

if ifrank3<0.01 && nowrank(3)<max_rank(3)

[A{3},X{3}]=rank_inc_adaptive(A{3},X{3},rank_inc);

end

%% check stopping criterion

if Res<tol

break

end

Y{1} = Unfold(Y_tensor,Nway,1); Y{1} = Y{1}';

Y{2} = Unfold(Y_tensor,Nway,2); Y{2} = Y{2}';

Y{3} = Unfold(Y_tensor,Nway,3); Y{3} = Y{3}';

Y_p=Y;X_p=X;A_p=A;

end

end

function [A,X]=rank_inc_adaptive(A,X,rank_inc)

% increase the estimated rank

for ii = 1:rank_inc

rdnx = rand(size(X,1),1);

rdna = rand(1,size(A,2));

X = [X,rdnx];

A = [A;rdna];

end

end
```