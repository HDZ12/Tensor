# ImgSeq2Video
```Matlab
function ImgSeq2Video(preFileName, formatFile, nS, nE, isColor)
% transform the image sequence into the .mat format
if isColor == 0
    for i = nS:nE
        fileName = [preFileName, int2str(i), formatFile];
        data(:,:,i-nS+1) = double(imread(fileName));
    end
    save(['preFileName', 'gray.mat'], 'data');
end

if isColor == 1
    for i = nS:nE
        fileName = [preFileName, int2str(i), formatFile];
        I = double(imread(fileName));
        dataR(:,:,i-nS+1) = I(:,:,1);
        dataG(:,:,i-nS+1) = I(:,:,2);
        dataB(:,:,i-nS+1) = I(:,:,3);
    end
    save([preFileName, 'R.mat'], 'dataR');
    save([preFileName, 'G.mat'], 'dataG');
    save([preFileName, 'B.mat'], 'dataB');
end
```
这段代码的目的是将一系列图像转换为 .mat 格式。它接受五个参数：preFileName（文件名前缀），formatFile（文件格式），nS 和 nE（定义要处理的图像序列的开始和结束索引），以及 isColor（一个标志，指示图像是否为彩色）。

如果 isColor 为0，表示图像是灰度的，那么代码会读取每个图像，将其转换为双精度数组，并将其存储在一个三维数组 data 中。然后，它会将 data 保存到一个名为 ‘preFileNamegray.mat’ 的文件中。

如果 isColor 为1，表示图像是彩色的，**那么代码会分别读取每个图像的红色**、绿色和蓝色通道，并将它们存储在三个三维数组 dataR、dataG 和 dataB 中。然后，它会将这三个数组分别保存到 ‘preFileNameR.mat’，‘preFileNameG.mat’ 和 ‘preFileNameB.mat’ 这三个文件中。
# CreateMovie
```Matlab
function CreateMovie(T, fileName)
[height, width, layer] = size(T{1});
aviobj = avifile(fileName, 'fps', 30, 'compression', 'None');
for i = 1:layer
    img(:,:,1) = T{1}(:,:,i);
    img(:,:,2) = T{2}(:,:,i);
    img(:,:,3) = T{3}(:,:,i);
    aviobj = addframe(aviobj, img);
end
aviobj = close(aviobj);
```
这段代码的目的是将一系列图像转换为 .avi 格式的电影。它接受两个参数：T（一个包含图像序列的单元数组）和 fileName（输出文件的名称）。

这段代码的关键步骤如下：

获取第一帧图像的尺寸。
创建一个新的 .avi 文件，帧率为30，不使用压缩。
对于 T 中的每一帧图像，将其添加到 .avi 文件中。
关闭 .avi 文件，完成电影的创建。
这段代码应该能够正确地执行，只要提供的图像文件存在并且可以被 imread 函数读取
# CreateImageSeq
```Matlab
function CreateImageSeq(inData, filePath, format)
[height, width, layer] = size(inData{1});
for n = 1:layer
    Img(:,:,1) = inData{1}(:,:,n);
    Img(:,:,2) = inData{2}(:,:,n);
    Img(:,:,3) = inData{3}(:,:,n);
    fileName = [filePath, int2str(n), format];
    imwrite(Img, fileName);
end
```
这段代码的目的是将一系列图像保存为文件。它接受三个参数：inData（一个包含图像序列的单元数组），filePath（输出文件的路径）和 format（输出文件的格式）。

这段代码的关键步骤如下：

获取第一帧图像的尺寸。
对于 inData 中的每一帧图像，将其保存为一个文件。文件名由 filePath、帧的索引和 format 组成。
这段代码应该能够正确地执行，只要提供的图像文件存在并且可以被 imwrite 函数写入。
# MySVDtau,MySVD,MySVDs,Turncate
这部分代码均为截断操作。具体参考`Truncate`

# SiLRTC

$$min(X, M1, M2, M3,... Mn):(\gamma1||X_{(1)}-M1||^2 + \gamma2||X_{(2)}-T2||^2 + \gamma3||X_{(3)}-T3||^2 + ...)/2+\alpha1||M1||_{\*}+\alpha2||M2||\_{\*}+\alpha3||M3||\_{\*}+....$$

$$s.t. X_\Omega = T_\Omega$$
**输入参数：**

- `T`：带有缺失值的输入张量。
- `Omega`：表示`T`中观测条目的二进制矩阵（观测到的为1，缺失为0）。
- `alpha`：低秩项的正则化参数。
- `gamma`：张量每个模式的正则化参数。
- `maxIter`：最大迭代次数。
- `epsilon`：收敛阈值。
- `X`：已完成张量的初始猜测（可选）。
```Matlab
function [X, errList] = SiLRTC(T, Omega, alpha, gamma, maxIter, epsilon, X)

if nargin < 7

X = T;

X(logical(1-Omega)) = mean(T(Omega));

end

errList = zeros(maxIter, 1);

normT = norm(T(:));

%L = errList;

dim = size(T);

M = cell(ndims(T), 1);

gammasum = sum(gamma);

tau = alpha./ gamma;

%normT = norm(T(:));

for k = 1:maxIter

if mod(k, 20) == 0

fprintf('SiLRTC: iterations = %d difference=%f\n', k, errList(k-1));

end

Xsum = 0;

for i = 1:ndims(T)

M{i} = Fold(Pro2TraceNorm(Unfold(X, dim, i), tau(i)), dim, i);

Xsum = Xsum + gamma(i) * M{i};

end

Xlast = X;

X = Xsum / gammasum;

X(Omega) = T(Omega);

errList(k) = norm(X(:)-Xlast(:)) / normT;

if (errList(k) < epsilon)

errList = errList(1:k);

break;

end

%L(k) = norm(X(:)-T(:)) / normT;

end

fprintf('SiLRTC ends: total iterations = %d difference=%f\n\n', k, errList(k));
```
**具体信息请参考我的另一个仓库：SCI-programming**
# FaLRTC
$$min_{X} : \Psi(X)$$

$$s.t. : X_\Omega = M_\Omega$$

$$\Psi(X) = max_{Z_{i(i)} <= 1}: <X, \sum_i Y_i> - 0.5 \mu_i \|Y_i\|_F^2$$
函数的输入参数包括：

- `M`：需要进行补全的张量
- `Omega`：表示已知元素位置的矩阵
- `alpha`：用于调整低秩正则项的权重
- `mu`：用于调整数据一致性项的权重
- `L`：Lipschitz常数的初始值
- `C`：用于调整Lipschitz常数的参数
- `maxIter`：最大迭代次数
- `epsilon`：收敛阈值
- `X`：初始化的张量（如果没有提供，则使用`M`的已知元素和平均值进行初始化）

函数的输出参数包括：

- `Y`：补全后的张量
- `errList`：每次迭代的误差列表

函数的主要步骤包括：

- 初始化参数和变量
- 迭代更新`Y`，直到满足收敛条件或达到最大迭代次数
- 
```Matlab
function [Y, errList] = FaLRTC(M, Omega, alpha, mu, L, C, maxIter, epsilon, X)

% initialization
if nargin < 9
    X= M;
    X(logical(1-Omega)) = mean(M(Omega));
end

Y = X;
Z = X;
B = 0;

N = ndims(M);
dim = size(M);
Gx = zeros(dim);

%funValueList = zeros(K+1,1);
%funValueList(1) = inf;

%%%%%%
normM = norm(M(:));
%errList(1) = norm(Y(:)-M(:)) / normM;
errList = zeros(maxIter, 1);
%%%%%%

% a2m = alpha.^2 ./ mu;
% ma = mu ./ alpha;
% 
% fylast = inf;

Lmax = 10*sum(1./mu);

tmp = zeros(1, N);
for i = 1:N
    % tmp(i) = max(SingularValue(Unfold(X, dim, i))) * alpha(i) * 0.4;
    tmp(i) = max(SingularValue(Unfold(X, dim, i))) * alpha(i) * 0.3;
end
P = 1.15;
flatNum = 15;
slope = (tmp - mu) / (1-(maxIter-flatNum)^(-P));
offset = (mu*(maxIter-flatNum)^P - tmp) / ((maxIter-flatNum)^P-1); 
% offset = (mu*K^P - tmp) / ((K-flatNum)^P-1); 

mu0 = mu;
for k = 1:maxIter
    if mod(k, 20) == 0
         fprintf('FaLRTC: iterations = %d   difference=%f\n', k, errList(k-1));
    end
    
    % update mu
    mu = max(slope / (k^P) +offset, mu0);
    % mu = mu0;
    % mu = mu0 * (K-k+1)^2
    
    a2m = alpha.^2 ./ mu;
    ma = mu ./ alpha;
    
    Ylast = Y;
    %%  test L
    %L = L*C;
    while true
        b = (1+sqrt(1+4*L*B)) / (2*L);
        X = b/(B+b) * Z + B/(B+b) * Ylast;
        % compute f'(x) namely "Gx" and f(x) namely "fx"
        Gx = Gx * 0;
        fx = 0;
        for i = 1 : N
            [temp, sigma2] = Truncate(Unfold(X, dim, i), ma(i));
            temp = Fold(temp, dim, i);
            Gx = Gx + a2m(i) * temp;
            fx = fx + a2m(i)*(sum(sigma2) - sum(max(sqrt(sigma2)-ma(i), 0).^2));
        end
        Gx(Omega) = 0;

        % compute f(Ytest) namely fy
        Y = X - Gx / L;
        fy = 0;
        for i = 1 : N
            [sigma] = SingularValue(Unfold(Y, dim, i));
            fy = fy + a2m(i)*(sum(sigma.^2) - sum(max(sigma-ma(i), 0).^2));
        end
        % test if L(fx-fy) > \|Gx\|^2
        if (fx - fy)*L < sum(Gx(:).^2)
            if L > Lmax
                k
                % funValueList = funValueList(2:k);
                % disp('Exceed the Maximum Lipschitiz Constant');
                fprintf('FaLRTC: iterations = %d   difference=%f\n Exceed the Maximum Lipschitiz Constan\n\n', k, errList(k-1));
                errList = errList(1:k);
                return;
                % break;
            end
            L = L/C;
        else
             break;
        end
    end
    
    errList(k) =  norm(Y(:)-Ylast(:)) / normM;
    if errList(k) < epsilon
        break;
    end 
    %% update Z, Y, and B
    Z = Z - b*Gx;
    B = B+b;
end
errList = errList(1:k);
fprintf('FaLRTC ends: total iterations = %d   difference=%f\n\n', k, errList(k));
```
**具体参考另一仓库**
# HaLRTC
$$min_X: \sum_i \alpha_i \|X_{i(i)}\|_*$$

$$s.t.:  X_\Omega = T_\Omega$$

函数的输入参数包括：

- `T`：需要进行补全的张量
- `Omega`：表示已知元素位置的矩阵
- `alpha`：用于调整低秩正则项的权重
- `beta`：用于调整数据一致性项的权重
- `maxIter`：最大迭代次数
- `epsilon`：收敛阈值
- `X`：初始化的张量（如果没有提供，则使用`T`的已知元素和平均值进行初始化）

函数的输出参数包括：

- `X`：补全后的张量
- `errList`：每次迭代的误差列表

函数的主要步骤包括：

- 初始化参数和变量
- 迭代更新`X`，直到满足收敛条件或达到最大迭代次数
- 返回补全后的张量和误差列表
```Matlab
function [X, errList] = HaLRTC(T, Omega, alpha, beta, maxIter, epsilon, X)

if nargin < 7
    X = T;
    X(logical(1-Omega)) = mean(T(Omega));
    % X(logical(1-Omega)) = 10;
end

errList = zeros(maxIter, 1);
dim = size(T);
Y = cell(ndims(T), 1);
M = Y;

normT = norm(T(:));
for i = 1:ndims(T)
    Y{i} = X;
    M{i} = zeros(dim);
end

Msum = zeros(dim);
Ysum = zeros(dim);
for k = 1: maxIter
    if mod(k, 20) == 0
        fprintf('HaLRTC: iterations = %d   difference=%f\n', k, errList(k-1));
    end
    beta = beta * 1.05;
    
    % update Y
    Msum = 0*Msum;
    Ysum = 0*Ysum;
    for i = 1:ndims(T)
        Y{i} = Fold(Pro2TraceNorm(Unfold(X-M{i}/beta, dim, i), alpha(i)/beta), dim, i);
        Msum = Msum + M{i};
        Ysum = Ysum + Y{i};
    end
    
    % update X
    %X(logical(1-Omega)) = ( Msum(logical(1-Omega)) + beta*Ysum(logical(1-Omega)) ) / (ndims(T)*beta);
    lastX = X;
    X = (Msum + beta*Ysum) / (ndims(T)*beta);
    X(Omega) = T(Omega);
    
    % update M
    for i = 1:ndims(T)
        M{i} = M{i} + beta*(Y{i} - X);
    end
    
    % compute the error
    errList(k) = norm(X(:)-lastX(:)) / normT;
    if errList(k) < epsilon
        break;
    end
end

errList = errList(1:k);
fprintf('FaLRTC ends: total iterations = %d   difference=%f\n\n', k, errList(k));
```
**具体参考另一仓库**