# Tensor
对于张量补全相关算法的matlab整理和复现
## Unfold
$$unfold_k(X) := X^{(k)} \in R^{I_k \times (I_1...I_{k-1}I_{k+1}...I_n)}$$
==一个用于展开张量的 MATLAB 函数。函数的输入参数是：==
- X：需要展开的张量。
- dim：张量 X 的维度。
- i：指定要展开的模式。\
函数的输出是展开后的矩阵 X。\
函数的工作原理是首先使用 shiftdim 函数将张量 X 沿着第 i 个维度进行旋转，然后使用 reshape 函数将旋转后的张量重塑为一个二维矩阵。这个二维矩阵的行数等于第 i 个维度的大小，列数等于其他维度的大小的乘积。
```Matlab
function [X] = Unfold( X, dim, i )
X = reshape(shiftdim(X,i-1), dim(i), []);
```
**这个函数首先使用 `shiftdim` 函数将张量 `X` 沿着第 `i` 个维度进行旋转，然后使用 `reshape` 函数将旋转后的张量重塑为一个二维矩阵。这个二维矩阵的行数等于第 `i` 个维度的大小，列数等于其他维度的大小的乘积**
## Fold
$$fold_k(X^{(k)}) := X$$
这个 Fold 函数的工作原理是，首先使用 circshift 函数将 dim 向左旋转 i-1 个位置，然后使用 reshape 函数将输入张量 X 重塑为 dim 指定的形状。最后，使用 shiftdim 函数将重塑后的张量 X 沿着第 length(dim)+1-i 个维度进行旋转
```Matlab
function [X] = Fold(X, dim, i)
dim = circshift(dim, [1-i, 1-i]);
X = shiftdim(reshape(X, dim), length(dim)+1-i);
```
**与`unfold`操作相对一个是折叠一个是展开**
# SigularValue
这个 SingularValue 函数的目的是计算输入矩阵 A 的奇异值。奇异值是描述矩阵作用于某些向量的标量，都是描述向量模长变化幅度的数值
函数的工作原理如下：

首先，函数获取输入矩阵 A 的大小 m 和 n。

如果 2\*m < n，函数计算 A 的自乘 A*A'，然后对结果进行奇异值分解 svd(AAT)，并取出奇异值 V1。然后，函数计算 V 的平方根，并返回结果。

如果 m > 2*n，函数计算 A 的转置自乘 A'*A，然后对结果进行奇异值分解 svd(AAT)，并取出奇异值 V。然后，函数计算 V 的平方根，并返回结果。

如果上述两个条件都不满足，函数直接对 A 进行奇异值分解 svd(A)，并取出奇异值 V1。然后，函数返回 V 的对角线元素。

这个函数的目的是提供一种高效的方法来计算矩阵的奇异值，特别是对于非常大或者非常小的矩阵
```Matlab
function V = SingularValue(A)
[m, n] = size(A);
if 2*m < n
    AAT = A*A';
    [S, V, D] = svd(AAT);
    V = sqrt(diag(V));
    return;
end
if m > 2*n
    AAT = A'*A;
    [S, V, D] = svd(AAT);
    V = sqrt(diag(V));
    return;
end
[S,V,D] = svd(A);
V = diag(V);
```
# Truncate
$$T_{\tau}(X) = U\Sigma_{\overline{\tau}}V^T$$

 $$\Sigma_{\overline{\tau}} = \text{diag}(\min(\sigma_i, \tau))$$
这个运算符的作用是将矩阵X的奇异值分解中的大奇异值截断为τ，从而得到一个对角线元素不超过τ的矩阵。
```Matlab
function [X, Sigma2] = Truncate(Z, tau)
% [Y, n, Sigma2] = Pro2TraceNorm(Z, tau);
% X = Z-Y;

%% new
[m, n] = size(Z);
if 2*m < n
    AAT = Z*Z';
    [S, Sigma2, D] = svd(AAT);
    Sigma2 = diag(Sigma2);
    V = sqrt(Sigma2);
    tol = max(size(Z)) * eps(max(V));
    n = sum(V > max(tol, tau));
    mid = max(V(1:n)-tau, 0) ./ V(1:n) ;
    X = (eye(m)-S(:, 1:n) * diag(mid) * S(:, 1:n)') * Z;
    return;
end
if m > 2*n
    [X, Sigma2] = Truncate(Z', tau);
    X = X';
    return;
end
[S,V,D] = svd(Z, 0);
Sigma2 = diag(V).^2;
n = sum(diag(V) > tau);
X = Z - S(:, 1:n) * diag(diag(V(1:n,1:n))-tau) * D(:, 1:n)';
```
```Matlab
function [X, n, Sigma2] = Pro2TraceNorm(Z, tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min: 1/2*||Z-X||^2 + ||X||_tr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [S, V, D, Sigma2] = MySVDtau(Z, tau);
% V = max(diag(V) - tau, 0);
% n = sum(V > 0);
% X = S(:, 1:n) * diag(V(1:n)) * D(:, 1:n)';

%% new
[m, n] = size(Z);
if 2*m < n
    AAT = Z*Z';
    [S, Sigma2, D] = svd(AAT);
    Sigma2 = diag(Sigma2);
    V = sqrt(Sigma2);
    tol = max(size(Z)) * eps(max(V));
    n = sum(V > max(tol, tau));
    mid = max(V(1:n)-tau, 0) ./ V(1:n) ;
    X = S(:, 1:n) * diag(mid) * S(:, 1:n)' * Z;
    return;
end
if m > 2*n
    [X, n, Sigma2] = Pro2TraceNorm(Z', tau);
    X = X';
    return;
end
[S,V,D] = svd(Z);
Sigma2 = diag(V).^2;
n = sum(diag(V) > tau);
X = S(:, 1:n) * max(V(1:n,1:n)-tau, 0) * D(:, 1:n)';
```
这两段代码的主要区别在于如何处理奇异值。在第一段代码中，所有小于 tau 的奇异值都被设置为0。而在第二段代码中，所有的奇异值都减去了 tau，然后再将所有小于0的值设置为0。这意味着第二段代码实际上是在对所有的奇异值应用了一个阈值函数。

另一个区别是在第二段代码中，X 是通过 S、V 和 D 计算得到的，而在第一段代码中，X 是通过 S 和 Z 计算得到的。这可能会导致得到的 X 有所不同
# 


