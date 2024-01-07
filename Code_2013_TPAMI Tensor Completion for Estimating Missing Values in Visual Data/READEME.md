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
$$min(X, M1, M2, M3,... Mn): (\gamma1||X_{(1)}-M1||^2 + \gamma2||X_{(2)}-T2||^2 + \gamma3||X_{(3)}-T3||^2 + ...)/2 + \alpha1||M1||_* + \alpha2||M2||_* + \alpha3||M3||_* + ....$$
$$s.t. X_\Omega = T_\Omega$$