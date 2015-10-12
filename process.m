close all
clear all

data = bfopen('B.nd2');
series = data{1, 1};%??

x_size = size(series{1,1},1);%series{1,1} first cell(picture)
y_size = size(series{1,1},2);
z_size = size(series,1);%series includes 128 cells(picture)

volume = zeros(x_size,y_size,z_size);%build a container

for z = 1:z_size
    plane = series{z,1};
    volume(:,:,z) = plane;
end

volume = volume(:,:,1:4:end); % 1:green; 2:green-yellow; 3:nir; 4:sum

% save_nii(make_nii(volume,[1 1 1], [0 0 0], 4),'B.nii.gz')

