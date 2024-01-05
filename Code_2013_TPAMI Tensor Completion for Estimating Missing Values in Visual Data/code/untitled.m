% 读取图像
originalImage = imread('testImg.png');

% 获取图像的大小
[rows, cols, ~] = size(originalImage);

% 计算要删除的像素范围
removeRows = round(rows * 0.3); % 删除百分之三十的行
removeCols = round(cols * 0.3); % 删除百分之三十的列

% 计算要删除的像素范围的上下左右边界
top = round((rows - removeRows) / 2);
bottom = top + removeRows;
left = round((cols - removeCols) / 2);
right = left + removeCols;

% 将要删除的像素值设置为0
originalImage(top:bottom, left:right, :) = 0;

% 保存编辑后的图像
imwrite(originalImage, 'modifiedImg.png');

% 显示编辑后的图像
imshow(originalImage);
