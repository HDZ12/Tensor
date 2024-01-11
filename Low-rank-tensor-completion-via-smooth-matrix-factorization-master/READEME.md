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
