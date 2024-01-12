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
