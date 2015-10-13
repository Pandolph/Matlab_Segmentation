function varargout = ImageSegmentation(varargin)
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

% Last Modified by GUIDE v2.5 10-Oct-2015 09:51:38

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
function ImageSegmentation_OpeningFcn(hObject, eventdata, handles, varargin)%初始化，暂不需要
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImageSegmentation (see VARARGIN)

% Choose default command line output for ImageSegmentation
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
str='Please load a file';
set(handles.edit3,'string',str);
set(handles.popupmenu1,'string','slice');

% --- Outputs from this function are returned to the command line.
function varargout = ImageSegmentation_OutputFcn(hObject, eventdata, handles) %输出到matlab界面，不需要
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function savepicture_Callback(hObject, eventdata, handles)%保存键
% hObject    handle to savepicture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function openfile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to openfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str='Opening···Please wait.';
set(handles.edit3,'string',str);
[filename, pathname] = uigetfile('*.nd2','select the .nd2 file');
data = bfopen([ pathname,filename]);
%data = bfopen('B.nd2');
series = data{1, 1};%??

x_size = size(series{1,1},1);%series{1,1} first cell(picture)
y_size = size(series{1,1},2);
global z_size;
z_size = size(series,1);%series includes 128 cells(picture)
global volume;
volume = zeros(x_size,y_size,z_size);%build a container

for z = 1:z_size
    plane = series{z,1};
    volume(:,:,z) = plane;
end
volume = volume(:,:,1:4:end); % 1:green; 2:green-yellow; 3:nir; 4:sum% save_nii(make_nii(volume,[1 1 1], [0 0 0], 4),'B.nii.gz')
sz=size(volume);
set(handles.popupmenu1,'string', {1:1:sz(3)});

global evalue;
evalue=max(volume,[],3);
% maxvalue=max(max(max(evalue)));
% A=find(evalue>maxvalue/4);
% evalue(A)=maxvalue/4;
evalue=medfilt2(evalue);
% axes(handles.axes);
% imshow(evalue,[]);

evalue_uint8=uint8(floor(evalue/16));
evalue_hist=histeq(evalue_uint8,256);
evalue_hist1=evalue_hist;
axes(handles.axes);
imshow(evalue_hist1);
evalue=im2double(evalue_hist);


%---均值加强
%sum_all=ones(sz(3),1)/sz(3);
% global average;
% average=reshape(reshape(volume,[],sz(3))*sum_all,[sz(1:2) 1]);%取均值，图像加强。
% average1=floor(average/16);%改变取值范围以适应直方图均衡
% average2=uint8(average1);
% average3=histeq(average2,256);
% axes(handles.axes);
% imshow(average3,[]);


%---图像切割（限制边界范围）
% str='Choose a rectangle with two point';
% set(handles.edit3,'string',str);
% [rx_1,ry_1]=ginput(1);
% [rx_2,ry_2]=ginput(1);
% rx_1=floor(rx_1);rx_2=floor(rx_2);ry_1=floor(ry_1);ry_2=floor(ry_2);
% rectangle=zeros(x_size,y_size);
% rectangle(ry_1:ry_2,rx_1:rx_2)=1;
% global rectangle1;
% rectangle1=rectangle;
% rectangle=uint8(rectangle);
% average3=average3.*rectangle;

str='Choose a point in a round';
set(handles.edit3,'string',str);


global segOutline;
segOutline=zeros(sz(1),sz(2));
segOutline=im2double(segOutline);




% --- Executes on selection change in popupmenu2.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global volume;
% global pflag;
global segOutline;
global slct;
slct=get(handles.popupmenu1,'value');
slice=volume(:,:,slct);
axes(handles.axes5);
outcome=slice+segOutline*max(max(slice));
imshow(outcome,[]);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton.
function pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global evalue;
[X,Y]=size(evalue);
%[x_p1,y_p1]=ginput(1);
%x_p1=floor(x_p1);
%y_p1=floor(y_p1);
%[x_p2,y_p2]=ginput(1);x_p2=floor(x_p2);y_p2=floor(y_p2);

%i=1;
%[x_p,y_p]=ginput(1);
%x_p=floor(x_p);
%y_p=floor(y_p);
%seed=sub2ind([X Y],y_p,x_p);

i=0;
seed=[];
% global pflag;
% pflag=0;



str='Choose points on the round';set(handles.edit3,'string',str);
axes(handles.axes);
while(i<2)   
    [x_p,y_p]=ginput(1);
    x_p=floor(x_p);
    y_p=floor(y_p);
    hold on
    plot(x_p,y_p,'g.','MarkerSize',5);
    seed=[seed,sub2ind([X Y],y_p,x_p)];
    i=i+1;
%     if pflag==1
%         break;
%     end 
    pause(1);   
end
hold off

%global rectangle1;
%average4=average.*rectangle1;
%[mask,probabilities] = random_walker(average4,seed,[1:1:i]);


global volume;
slct=get(handles.popupmenu1,'value');
slice=volume(:,:,slct);
%slice1=max(max(slice))-slice;

[mask,probabilities] = random_walker(evalue,seed,[1:1:i]);

%[mask,probabilities] = random_walker(evalue,[sub2ind([X Y],y_p1,x_p1), ...
%    sub2ind([X Y],y_p2,x_p2)],[1,2]);
mask=reshape(mask,X,Y);
% D=[0 1 0,1 1 1,0 1 0];
% mask=imdilate(mask,D);
[fx,fy]=gradient(mask);
x_mask=find(fx);
y_mask=find(fy);
global segOutline;
segOutline=zeros(X,Y);
segOutline(x_mask)=1;
segOutline(y_mask)=1;%获取轮廓
axes(handles.axes);
outcome1=evalue+segOutline*max(max(evalue));
imshow(outcome1,[]);
axes(handles.axes5);
outcome2=slice+segOutline*max(max(slice));
imshow(outcome2,[]);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global segOutline;
global volume;
slct=get(handles.popupmenu1,'value');
slice=volume(:,:,slct);
num=sum(segOutline);
value=sum(segOutline.*slice)/num;
set(handles.text6,'string', num2str(value));


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global volume;
% global pflag;
global segOutline;
global slct;
slct=get(handles.popupmenu1,'value');
if slct>1
 slct=slct-1;
end
set(handles.popupmenu1,'value',slct);
slice=volume(:,:,slct);
axes(handles.axes5);
outcome=slice+segOutline*max(max(slice));
imshow(outcome,[]);    




% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global volume;
% global pflag;
global segOutline;
global slct;
slct=get(handles.popupmenu1,'value');
sz=size(volume);
if slct<sz(3)
slct=slct+1;
end
set(handles.popupmenu1,'value',slct);
slice=volume(:,:,slct);
axes(handles.axes5);
outcome=slice+segOutline*max(max(slice)); 
imshow(outcome,[]);
