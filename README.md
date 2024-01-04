# Tensor
对于张量补全相关算法的matlab复现
## Unfold
$$unfold_k(X) := X^{(k)} \in R^{I_k \times (I_1...I_{k-1}I_{k+1}...I_n)}$$
一个用于展开张量的 MATLAB 函数。函数的输入参数是：
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
