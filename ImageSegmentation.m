function varargout = ImageSegmentation(varargin)%%do not edit
% IMAGESEGMENTATION MATLAB code for ImageSegmentation.fig
%      IMAGESEGMENTATION, by itself, creates a new IMAGESEGMENTATION or raises the existing
%      singleton*.
%
%      H = IMAGESEGMENTATION returns the handle to a new IMAGESEGMENTATION or the handle to
%      the existing singleton*.
%
%      IMAGESEGMENTATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGESEGMENTATION.M with the given input arguments.
%
%      IMAGESEGMENTATION('Property','Value',...) creates a new IMAGESEGMENTATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImageSegmentation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImageSegmentation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ImageSegmentation

% Last Modified by GUIDE v2.5 19-Jan-2018 01:19:43

% set(gcf,'units','normalized')%?????,??????????;

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ImageSegmentation_OpeningFcn, ...
    'gui_OutputFcn',  @ImageSegmentation_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ImageSegmentation is made visible.
function ImageSegmentation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImageSegmentation (see VARARGIN)
%set (handles.axes,'Xtick','off','Ytick','off');
%set (handles.axes5,'Xtick','off','Ytick','off');
%set(handles.axes,'visible','off');
% axes(handles.axes)%select axes
% set(handles.axes,'visible','off');
% axes(handles.axes)%select axes
% set(handles.axes,'visible','off');
% Choose default command line output for ImageSegmentation

% ha=axes('units','normalized','position',[0 0 1 1]);
% uistack(ha,'down')
% II=imread('med-x.jpg');
% image(II)
% colormap gray
% set(ha,'handlevisibility','off','visible','off');

handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
%str='Please load a file';
%set(handles.edit3,'string',str);
set(handles.slice,'string','slice');

% --- Outputs from this function are returned to the command line.
%output in matlab command window, do need to edit
function varargout = ImageSegmentation_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function slice_CreateFcn(hObject, eventdata, handles)

% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

%set (handles.axes,'Xtick','off','Ytick','off');
%set (handles.axes5,'Xtick','off','Ytick','off');
%set(handles.axes,'visible','off');

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
%read a nd2 file,get the operation file 'evalue'
function openfile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to openfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%str='Opening...Please wait.';
%set(handles.edit3,'string',str);
[~,PathName] = uigetfile('*.tif','Select the file');
channel = '1';
if strcmp(channel, 'all')
    Files = dir(fullfile(PathName,'*.tif'));
elseif strcmp(channel,  'part')
    Files1 = dir(fullfile(PathName,['*','.tif']));
    prompt = {'Start','End'};
    dlg_title = 'Input';
    num_lines = 1;
    defaultans = {'1','100'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    Files = Files1([str2num(answer{1,1}):str2num(answer{2,1})],1);
else
    Files = dir(fullfile(PathName,['*',num2str(channel),'.tif']));
end
place = [PathName,Files(1).name];
tif = imread(place); % z = 156 slices, x = 437 or tif(:,1,1), y = 1379 or tif(1,:,1)
sliceSize = length(Files(:,1));
Data = ones(length(tif(:,1,1)),length(tif(1,:,1)),sliceSize);
for k = 1:sliceSize
    place = [PathName,Files(k,1).name];
    tif = imread(place);
    Data(:,:,k) = tif(:,:,find(mean(mean(tif,1),2)));
end

Crop = Data;

x_size = size(Crop,1);%series{1,1} first cell(picture)
y_size = size(Crop,2);
z_size = size(Crop,3);%series includes 128 cells(picture)
global volumeAll;
volumeAll = zeros(x_size,y_size,z_size);%build a container

for z = 1:z_size
    plane = Crop(:,:,z);
    volumeAll(:,:,z) = plane;
end
global volume;
volume = volumeAll(:,:,:); % 1:green; 2:green-yellow; 3:nir; 4:sum% save_nii(make_nii(volume,[1 1 1], [0 0 0], 4),'B.nii.gz')
sz=size(volume);
set(handles.slice,'string', {1:1:sz(3)});%set the selection of popupmenu by number 1 to sz(3)
global evalue;
evalue=max(volume,[],3);% get the maxvalue of each pixel from all slice
%%cut the maximum part of the image for the effect of histogram equalization
maxvalue=max(max(max(evalue)));
A=find(evalue>maxvalue/4);
evalue(A)=maxvalue/4;
evalue=medfilt2(evalue);%median filter, remove the white noise
evalue_uint8=uint8(floor(evalue/16));
evalue_hist=histeq(evalue_uint8,256);%histogram equalization the image 'evalue', make the image easier to read
evalue=im2double(evalue_hist);  % change the vaule of evalue to (0,1)
axes(handles.axes);%the left axes
imshow(evalue,[]);

%%image enhancement by average the dataset
% sum_all=ones(sz(3),1)/sz(3);
%
% evalue=reshape(reshape(volume,[],sz(3))*sum_all,[sz(1:2) 1]);%???????????????????????????????????
% evalue1=floor(evalue/16);%?????????????????????????????????????????????????????
% evalue2=uint8(evalue1);
% evalue3=histeq(evalue2,256);
% evalue=im2double(evalue3);
% axes(handles.axes);
% imshow(evalue,[]);

%str='Finshed.';
%set(handles.edit3,'string',str);

%%initialize some global value
global segOutline; % A 2-D 0-1 matrix stand for the result of segmentation, the edge consist of pixels which value are 1
segOutline=zeros(sz(1),sz(2));
segOutline=im2double(segOutline);

global cflag;%true already cut
cflag=0;

global sflag;%true already segment
sflag=0;

global rect;%the array about the coordinate of vertex of cutting rectangle  [y1,y2,x1,x2]
rect=[1,sz(2),1,sz(1)];%be careful about the order

global savevalue;%container saves evalue, used in function of restart
savevalue=[];


% --- Executes on selection change in popupmenu2.
function slice_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global volume;
global segOutline;
global cflag;
global sflag;
global rect;
slct=get(handles.slice,'value');
slice=volume(:,:,slct);%get a slice for volume dataset
if cflag==1%means image have been cut
    slicecut=slice(rect(1):rect(2),rect(3):rect(4));
    slice=slicecut;
end
axes(handles.axes5);
k=max(max(slice));
if sflag==1
    slice=slice+segOutline*k;
end
imshow(slice,[0,k]);


% --- Executes on button press in selectpoint.
function selectpoint_Callback(hObject, eventdata, handles)
% hObject    handle to selectpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%str='Choose points in the round';
%set(handles.edit3,'string',str);

global label;
global seed;
global evalue;
global myseed;  %used for plot
myseed= [];
seed=[];
label=[];
[X,Y]=size(evalue);

% if ~isempty(label)
%     label = [0 label 1];
%     axes(handles.axes);
%     for i = 1:2
%         [x_p,y_p]=ginput(1);
%         x_p=floor(x_p);
%         y_p=floor(y_p);
%         hold on
%         if i == 1
%             plot(x_p,y_p,'r.','MarkerSize',5);
%             seed=[sub2ind([X Y],y_p,x_p),seed];
%             myseed = [myseed [0;x_p;y_p]];
%         else
%             plot(x_p,y_p,'g.','MarkerSize',5);
%             seed=[seed,sub2ind([X Y],y_p,x_p)];
%             myseed = [myseed [1;x_p;y_p]];
%         end
%     end
%     hold off
% else
%     i=0;
%     n=4;%you can change the num of seeds
%     myseed= [];
%     seed=[];
%     label=[0];
%     N=ones(1,n-1);
%     label=[label N];%first version, fixed 4 seeds and 2 types of labels.
%     %str='Choose points on the round';
%     %set(handles.edit3,'string',str);
%     axes(handles.axes);
%     while(i<n)
%         [x_p,y_p]=ginput(1);
%         x_p=floor(x_p);
%         y_p=floor(y_p);
%         hold on
%         if i ==0
%             plot(x_p,y_p,'r.','MarkerSize',5);
%             myseed = [myseed [0;x_p;y_p]];
%         else
%             plot(x_p,y_p,'g.','MarkerSize',5);
%             myseed = [myseed [1;x_p;y_p]];
%         end
%         seed=[seed,sub2ind([X Y],y_p,x_p)];
%         i=i+1;
%         %     pause(1);
%     end
%     %hold off
% end

%numberOfPoints = inputdlg({'number of inside points','number of outside points'});
numberOfPoints = {get(handles.numPoints, 'String'),get(handles.numPoints, 'String')};
%set the inside points 0, and put them in the front

label = zeros(1,str2num(cell2mat(numberOfPoints(1))));
label = [label,ones(1,str2num(cell2mat(numberOfPoints(2))))];

for i = 1:length(label)
    [x_p,y_p]=ginput(1);
    x_p=floor(x_p);
    y_p=floor(y_p);
    hold on
    if i <= str2num(cell2mat(numberOfPoints(1)))
        plot(x_p,y_p,'r.','MarkerSize',5);
        myseed = [myseed [0;x_p;y_p]];
    else
        
        plot(x_p,y_p,'g.','MarkerSize',5);
        myseed = [myseed [1;x_p;y_p]];
    end
    hold off
    seed=[seed,sub2ind([X Y],y_p,x_p)];
end

global volume;
slct=get(handles.slice,'value');
slice=volume(:,:,slct);
global cflag;
global rect;
if cflag==1
    slicecut=slice(rect(1):rect(2),rect(3):rect(4));
    slice=slicecut;
end

[mask,probabilities] = random_walker(evalue,seed,label);
%mask=reshape(mask,X,Y);
%%inflation
% D=[0 1 0,1 1 1,0 1 0];
% mask=imerode(mask,D);
% mask=imclose(mask,D);

global segOutline;

% [fx,fy]=gradient(mask);
% x_mask=find(fx);
% y_mask=find(fy);
% segOutline=zeros(X,Y);
% segOutline(x_mask)=1;
% segOutline(y_mask)=1;

I=bwmorph(mask,'remove');
I=~I;
%figure;imshow(I);
segOutline = I;
segOutline(1,:)=1;
segOutline(end,:)=1;
segOutline(:,1)=1;
segOutline(:,end)=1;
segOutline = ~segOutline;

global numofpixels;
global newsegOutline;
%numofpixels=str2num(cell2mat(numberOfPoints(3)));
%newsegOutline = Expand(segOutline,numofpixels);
axes(handles.axes);


hold on


outcome1=evalue+segOutline*max(max(evalue));
imshow(outcome1,[]);         %outline in original picture


for i = 1:length(label)
    if myseed(1,i)==0
        plot(myseed(2,i),myseed(3,i),'r.','MarkerSize',5);
    else
        plot(myseed(2,i),myseed(3,i),'g.','MarkerSize',5);
    end
end
hold off

axes(handles.axes5);
outcome2=slice+segOutline*max(max(slice));
imshow(outcome2,[]);        %outline in second picture
global sflag;
sflag=1;

% figure;
% subplot(3,3,1);
% imshow(segOutline);title('segOutline');
% subplot(3,3,2);
% D = bwdist(segOutline,'euclidean');
% imagesc(D);title('imshow D');
% subplot(3,3,3);
% subimage(mat2gray(D));
% imcontour(D);
% title('D imcontour');
% subplot(3,3,4);
% contour(D);title('D contour');
% subplot(3,3,5);
% imagesc(D);title('D imagesc');
% subplot(3,3,6);
% [px,py] = gradient(D);
% [c, h] = contour(D, [3, 3]);  %temporaryly set v = 3;
% title('D gradient');
% %level1 = 3, level2 = 3
% % C = [level1 x1 x2 x3 ... level2 x2 x2 x3 ...;
% %      pairs1 y1 y2 y3 ... pairs2 y2 y2 y3 ...]
% subplot(3,3,7);
% [C, H] = contour(D, [5, 1]);
% hold on, quiver(px,py), hold off
% a1 = round(c);
% %figure; plot(a1);
% a2 = a1(:,2:c(2,1)+1);
% a3 = zeros(1,c(2,1));
% a8 = zeros(max(a1(:)));
% for i = 1:length(a1)
%     a8(a1(1,i),a1(2,i))=1;
% end
% %figure;imshow(a8);
% % set the distance 5
% for i = 1: c(2,1)
% %   a3(i)=sqrt(px(a2(1,i),a2(2,i))^2+py(a2(1,i),a2(2,i))^2);
% end

% D = bwdist(segOutline,'euclidean');
% [px,py] = gradient(D);
% for j = 1: 30
%     figure(10);
%     axis equal;                     % equal axes
%     close(figure(10));
%     a1 = round(c);
%     a2 = a1(:,2:c(2,1)+1);
%     gray = zeros(length(segOutline(1,:)),length(segOutline(:,1)));
%     vector = 0;
%     for i = 1:c(2,1)
%         positionX = vector*px(a2(1,i),a2(2,i))+a2(1,i);
%         positionY = vector*py(a2(1,i),a2(2,i))+a2(2,i);
%         %         positionX = a2(1,i);
%         %         positionY = a2(2,i);
%         gray (round(positionX), round(positionY)) = 1;
%     end
%     gray = flipud(gray);
%     gray = rot90(gray,3);
%     temp = evalue.*gray;
%     nonzerovalue = temp(find(temp));
%     nonzerovalue = nonzerovalue';
%     disp(['mean = ',num2str(mean(nonzerovalue)),' v =',num2str(j)]);
%     %figure;
%     %plot(nonzerovalue); 
%     if j == 1
%         numpixel = length(nonzerovalue);
%         finalmatrix = nonzerovalue;
%     else
%         newmatrix = 1:numpixel;
%         newmatrix = round(newmatrix*length(nonzerovalue)/numpixel);
%         newnonzerovalue = nonzerovalue(newmatrix);
%         finalmatrix = [finalmatrix;newnonzerovalue];
%     end
%     
%     %disp(num(find(temp)));
% end
% figure;
% imshow(finalmatrix);
% for i = 1:c(2,1)
%     gray(i) = evalue(round(positionX), round(positionY));
% end
% plot(gray);title(vector);


% D = bwdist(segOutline,'euclidean');
% [px,py] = gradient(D);  % for example px(2,1) = (D(3,1)+D(1,1))/2
% figure;imagesc(D);
% %axis([1,108,90,1]);
% figure;contour(D);
% vector = 5;
% graymapIn = zeros(length(segOutline(1,:)),length(segOutline(:,1)));
% graymapOut = graymapIn;
% 
% [x,y]=find(segOutline);
% for i = 1:length(find(segOutline))    
%     positionX = vector*px(x(i),y(i))+x(i);
%     positionY = vector*py(x(i),y(i))+y(i);
%     graymapOut(round(positionX), round(positionY)) = 1;
%     positionX = (-1)*vector*px(x(i),y(i))+x(i);
%     positionY = (-1)*vector*py(x(i),y(i))+y(i);
%     graymapIn(round(positionX), round(positionY)) = 1;
% end
% figure;imshow(graymapOut);
% figure;imshow(graymapIn);
%  
% for i = 1:length(find(segOutline))    
%     posX(i) = px(x(i),y(i));
%     posY(i) = py(x(i),y(i));
% end
% figure;imshow(neighborIn);
% figure;imshow(segOutline);
% figure;imshow(neighborOut);


% --- Executes on button press in average. Line
function average_Callback(hObject, eventdata, handles)
% hObject    handle to average (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global segOutline;

global volumeAll;
global cflag;
global rect;
global newsegOutline;
Crop = volumeAll;
num_hist = 10; % the slice for histmatch
x = 9;
outHist = HistMatch(Crop,num_hist);%figure;imshow(outHist(:,:,x),[])
Col = ColumnMatch(outHist);%figure;imshow(Col(:,:,x),[])
filter_size =5; % the size of filer in the preprosessing to denoise
filter_num =1; % how many times the filter is applied
bi_threshold = 0.1; %the threshhold of binarization, larger means fewer area left
Bi = MeanBi(Col,bi_threshold,filter_num,filter_size);%figure;imshow(Bi(:,:,x),[])
BiVessel = bwareaopen(Bi,40,8);
se = strel('disk',7);
I = imdilate(BiVessel,se);
vessel = imerode(I,se);
save('vessel.mat','vessel')
axes(handles.axes);
imshow(vessel(:,:,x),[])



% --- Executes on button press in back.
function back_Callback(hObject, eventdata, handles)
% hObject    handle to back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global volume;
% global pflag;
global segOutline;
global rect
global cflag
global sflag;
slct=get(handles.slice,'value');
if slct>1
    slct=slct-1;
end
set(handles.slice,'value',slct);

slice=volume(:,:,slct);
if cflag==1
    slicecut=slice(rect(1):rect(2),rect(3):rect(4));
    slice=slicecut;
end
axes(handles.axes5);
k=max(max(slice));
if sflag==1
    slice=slice+segOutline*k;
end
imshow(slice,[0,k]);

% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global volume;
global segOutline;
global rect;
global cflag;
global sflag;
slct=get(handles.slice,'value');
sz=size(volume);
if slct<sz(3)
    slct=slct+1;
end
set(handles.slice,'value',slct);
slice=volume(:,:,slct);
if cflag==1
    slicecut=slice(rect(1):rect(2),rect(3):rect(4));
    slice=slicecut;
end
axes(handles.axes5);
k=max(max(slice));
if sflag==1
    slice=slice+segOutline*k;
end
imshow(slice,[0,k]);

% --- Executes on button press in cut.
function cut_Callback(hObject, eventdata, handles)
% hObject    handle to cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%---select and cut--
global volume
global evalue;
global rect;
%global slice;
%str='Choose a rectangle with two point';
%set(handles.edit3,'string',str);
axes(handles.axes);
[rx_1,ry_1]=ginput(1);
rx_1=floor(rx_1);ry_1=floor(ry_1);
hold on
plot(rx_1,ry_1,'r.','MarkerSize',5);
[rx_2,ry_2]=ginput(1);
rx_2=floor(rx_2);ry_2=floor(ry_2);
plot(rx_2,ry_2,'r.','MarkerSize',5);
hold off
rect=[ry_1,ry_2,rx_1,rx_2];
global savevalue;
savevalue={savevalue;evalue};
evaluecut=evalue(rect(1):rect(2),rect(3):rect(4));
evalue=evaluecut;
evalue = imadjust(evalue);
imshow(evalue,[]); %????????
%imagesc(evalue);
slct=get(handles.slice,'value');
slice=volume(:,:,slct);
slicecut=slice(rect(1):rect(2),rect(3):rect(4));

%slice=slicecut;
axes(handles.axes5);
%imshow(slice,[]);      %???????
imagesc(slicecut);
global cflag;%ture:cflag=1,operate on cut iamge
cflag=1;
global seed;
global label;
seed = [];
label = [];

% --- Executes on button press in restart.
function restart_Callback(hObject, eventdata, handles)
% hObject    handle to restart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global volume;
global evalue;
global savevalue;
if ~isempty(savevalue)
    evalue=savevalue{2};%forms some kind of stack
end
axes(handles.axes);
imshow(evalue,[]);
slct=get(handles.slice,'value');
slice=volume(:,:,slct);
axes(handles.axes5);
imshow(slice,[]);

sz=size(volume);
global segOutline;
segOutline=zeros(sz(1),sz(2));
segOutline=im2double(segOutline);

global cflag;
cflag=0;

global rect;
rect=[1,sz(2),1,sz(1)];

global sflag;
sflag=0;

savevalue=savevalue(1);
global label;
global seed;
label =[0];
seed=[];
set(handles.ContrastValue, 'String', 'contrast'); % 0 means black and 1 means white
set(handles.contrastIn, 'String', 'valueinside');
set(handles.contrastOut, 'String', 'valueoutside');
set(handles.text6,'string', 'value');


% --- Executes during object creation, after setting all properties.
function axes_CreateFcn(hObject, eventdata, handles)
set(hObject,'xTick',[]);
set(hObject,'ytick',[]);
set(hObject,'box','on');
% hObject    handle to axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes


% --- Executes during object creation, after setting all properties.
function axes5_CreateFcn(hObject, eventdata, handles)
set(hObject,'xTick',[]);
set(hObject,'ytick',[]);
set(hObject,'box','on');
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes5


% --- Executes on button press in expand.  Sphere
function expand_Callback(hObject, eventdata, handles)
% hObject    handle to expand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 global label;
% global seed;
 global myseed;  %used for plot
% label = label(2:end-1);
% seed = seed(2:end-1);
% myseed = myseed(:,1:end-2);
global segOutline;
global numofpixels;
global newsegOutline;
global evalue;
global volumeAll
Crop = volumeAll;
num_hist = 50;
outHist = HistMatch(Crop,num_hist);
%Col = ColumnMatch(outHist);
figure;imshow(outHist(:,:,1),[]);

filter_size =3;
filter_num =3;
bi_threshold = 0.45; % most important, larger means fewer left
Bi = MeanBi(outHist,bi_threshold,filter_num,filter_size);
figure;imshow(Bi(:,:,1));

area_size = 20; % delete the area whose size is smaller than area_size
se_size = 3;
cell = Cell(Bi,area_size,se_size);
save('microglia.mat','cell')
axes(handles.axes);
imshow(cell,[])



% --- Executes during object creation, after setting all properties.
function ExpandNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ExpandNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',3);


% --- Executes on button press in Contrast.
function Contrast_Callback(hObject, eventdata, handles)
% hObject    handle to Contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global volumeAll;
global cflag;
global rect;
global segOutline;
volumeContrast = volumeAll(:,:,1:4:end);
inside = str2num(get(handles.valueInside,'String'));
outside = str2num(get(handles.valueOutside,'String'));
neighborIn = Expand(segOutline,inside)-Expand(segOutline,inside+1);
neighborOut = Expand(segOutline,outside)-Expand(segOutline,outside-1);
if cflag==1
    slice=volumeContrast(:,:,double(get(handles.slice, 'Value')));
    slicecut=slice(rect(1):rect(2),rect(3):rect(4));
    slice=slicecut;
end
valueIn=sum(neighborIn.*slice)/sum(neighborIn);
valueOut=sum(neighborOut.*slice)/sum(neighborOut);
contrast = (valueOut-valueIn)/valueIn;
set(handles.ContrastValue, 'String', num2str(contrast)); % 0 means black and 1 means white
set(handles.contrastIn, 'String', num2str(valueIn));
set(handles.contrastOut, 'String', num2str(valueOut));


% expand the segOutline to several pixels
function [newsegOutline] = Expand(segOutline,numofpixels)
%msgbox('a');
%disp(segOutline);
[m,n]=find(segOutline);
center = [round(mean(m)),round(mean(n))];
newsegOutline=segOutline;

for i = 1:length(m(:))
    if m(i)>=center(1)
        if n(i)>=center(2)
            if numofpixels > 0
                newsegOutline(m(i):m(i)+numofpixels,n(i):n(i)+numofpixels)=1;
            else
                newsegOutline(m(i)+numofpixels:m(i),n(i)+numofpixels:n(i))=1;
            end
        else
            if numofpixels > 0
                newsegOutline(m(i):m(i)+numofpixels,n(i)-numofpixels:n(i))=1;
            else
                newsegOutline(m(i)+numofpixels:m(i),n(i):n(i)-numofpixels)=1;
            end
        end
    else
        if n(i)>=center(2)
            if numofpixels > 0
                newsegOutline(m(i)-numofpixels:m(i),n(i):n(i)+numofpixels)=1;
            else
                newsegOutline(m(i):m(i)-numofpixels,n(i)+numofpixels:n(i))=1;
            end
        else
            if numofpixels > 0
                newsegOutline(m(i)-numofpixels:m(i),n(i)-numofpixels:n(i))=1;
            else
                newsegOutline(m(i):m(i)-numofpixels,n(i):n(i)-numofpixels)=1;
            end
        end
    end
end
newsegOutline = newsegOutline-segOutline;



function ExpandNum_Callback(hObject, eventdata, handles)
% hObject    handle to ExpandNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ExpandNum as text
%        str2double(get(hObject,'String')) returns contents of ExpandNum as a double



function numPoints_Callback(hObject, eventdata, handles)
% hObject    handle to numPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numPoints as text
%        str2double(get(hObject,'String')) returns contents of numPoints as a double


% --- Executes during object creation, after setting all properties.
function numPoints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',5);



function valueInside_Callback(hObject, eventdata, handles)
% hObject    handle to valueInside (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of valueInside as text
%        str2double(get(hObject,'String')) returns contents of valueInside as a double


% --- Executes during object creation, after setting all properties.
function valueInside_CreateFcn(hObject, eventdata, handles)
% hObject    handle to valueInside (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function valueOutside_Callback(hObject, eventdata, handles)
% hObject    handle to valueOutside (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of valueOutside as text
%        str2double(get(hObject,'String')) returns contents of valueOutside as a double


% --- Executes during object creation, after setting all properties.
function valueOutside_CreateFcn(hObject, eventdata, handles)
% hObject    handle to valueOutside (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Cyliner.
function Cyliner_Callback(hObject, eventdata, handles)
% hObject    handle to Cyliner (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global volumeAll
Crop = volumeAll;
num_hist = 30; % choose the cleareat one
outHist = HistMatch(Crop,num_hist);
%Col = ColumnMatch(outHist);
figure;imshow(Crop(:,:,1),[])
filter_size =2;
filter_num =6;
bi_threshold = 0.3; % most important, larger means fewer left
Bi = MeanBi(outHist,bi_threshold,filter_num,filter_size);
figure;imshow(Bi(:,:,1),[])

se = strel('disk',6);
temp = imerode(Bi,se);
myelin = imdilate(temp,se);
save('myelin.mat','myelin')
axes(handles.axes);
imshow(myelin,[])
