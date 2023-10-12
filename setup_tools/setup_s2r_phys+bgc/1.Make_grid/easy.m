function varargout = easy(varargin)
% 
%  EASY is a GUI for the Easy Grid package
%  execute by typing: easy
%
%  (c) 2008, Jeroen Molemaker,  UCLA
%
%
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @easy_OpeningFcn, ...
                   'gui_OutputFcn',  @easy_OutputFcn, ...
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

% --- Executes just before easy is made visible.
function easy_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for easy
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using easy.
if strcmp(get(hObject,'Visible'),'off')

   nx     = 100;
   ny     = 100;
   size_x = 1000e3;
   size_y = 1000e3;
   rotate =    0;
   tra_lon=    5;
   tra_lat= 52.5;
 if exist('easy_settings.mat')
  load('easy_settings')
 end
 set(handles.edit1,'String',num2str(size_x/1e3));   %% size_x (km)
 set(handles.edit2,'String',num2str(size_y/1e3));   %% size_y (km)
 set(handles.edit3,'String',num2str(rotate));  %% Rotation
 set(handles.edit4,'String',num2str(tra_lon)); %% Lon Center
 set(handles.edit5,'String',num2str(tra_lat)); %% Lat Center
 set(handles.edit6,'String',num2str(nx));      %% nx
 set(handles.edit7,'String',num2str(ny));      %% ny

 size_x = str2double(get(handles.edit1,'String'))*1e3; %%
 size_y = str2double(get(handles.edit2,'String'))*1e3; %%
 rotate = str2double(get(handles.edit3,'String'));     %%
 tra_lon= str2double(get(handles.edit4,'String'));     %%
 tra_lat= str2double(get(handles.edit5,'String'));     %%
 nx     = str2num(get(handles.edit6,'String'));
 ny     = str2num(get(handles.edit7,'String'));

   [lon,lat,pm,pn,ang,lone,late] = easy_grid(nx,ny,size_x,size_y,tra_lon,tra_lat,rotate);
   out_lon = [lone(2,2:end-1) lone(2:end-1,end-1)' lone(end-1,end-1:-1:2) lone(end-1:-1:2,2)'];
   out_lat = [late(2,2:end-1) late(2:end-1,end-1)' late(end-1,end-1:-1:2) late(end-1:-1:2,2)'];
   r_earth = 6371315.;
   colormap(jet(256))
   dlon = max(size_x,size_y)/r_earth;
   if dlon*180/pi > 60
    lo0 = min(min(lon))*180/pi - 2;
    lo1 = max(max(lon))*180/pi + 2;
    la0 = min(min(lat))*180/pi - 2;
    la1 = max(max(lat))*180/pi + 2;
    dlon*180/pi
    lo0
    lo1
    disp('mercator')
    m_proj ('mercator','longitude',[lo0 lo1],'latitude',[la0 la1]);
   else
    m_proj('Gnomonic','lon',tra_lon,'lat',tra_lat,'rad',1.0*dlon*180/pi,'rec','on')
   end

   load coast_global
  % load /batavia/nmolem/OBSERV/COAST/coast_pac; clon = clon + 360;
   m_plot(out_lon*180/pi,out_lat*180/pi,'r')
   m_grid
   hold on
   m_coast('patch',[.7 .7 .7],'edgecolor','k');
   m_gshhs_i('patch',[.7 .7 .7],'edgecolor','k');
%  m_plot(clon-180, -clat,'r')
%  m_plot(clon, clat,'r')
%  m_gshhs_i('save','testing.mat');
%  m_usercoast('testing.mat','patch',[.5 .5 .5],'edgecolor','r')
%  m_plot(out_lon*180/pi,out_lat*180/pi,'r')
   hold off
end

% UIWAIT makes easy wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function varargout = easy_OutputFcn(hObject, eventdata, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;

function pushbutton1_Callback(hObject, eventdata, handles)
% --- Executes on button press in pushbutton1.
axes(handles.axes1);
cla;

 load coast_global
% load /batavia/nmolem/OBSERV/COAST/coast_pac; clon = clon + 360;
 popup_sel_index = get(handles.popupmenu1, 'Value');
 size_x = str2double(get(handles.edit1,'String'))*1e3; %% Width 
 size_y = str2double(get(handles.edit2,'String'))*1e3; 
 rotate = str2double(get(handles.edit3,'String'));
 if abs(rotate) > 45
   rotate = rotate/abs(rotate)*45;
   set(handles.edit3,'String',num2str( rotate));   %% Rotation
 end
 
 tra_lon= str2double(get(handles.edit4,'String'));
 tra_lat= str2double(get(handles.edit5,'String'));
 nx     = str2num(get(handles.edit6,'String'));
 ny     = str2num(get(handles.edit7,'String')); 
 save('easy_settings','nx','ny','size_x','size_y','rotate','tra_lon','tra_lat')
 [lon,lat,pm,pn,ang,lone,late] = easy_grid(nx,ny,size_x,size_y,tra_lon,tra_lat,rotate);


 colormap(jet(256))
 r_earth = 6371315.;
 dlon = max(size_x,size_y)/r_earth;
  if dlon*180/pi > 60
    %% TODO: move the lons smarter!
%   lon(lon>0)  = lon(lon>0)  - 2*pi;
%   lone(lone>0)= lone(lone>0)- 2*pi;
    lo0 = min(min(lon))*180/pi - 4;
    lo1 = max(max(lon))*180/pi + 4;
    la0 = min(min(lat))*180/pi - 4;
    la1 = max(max(lat))*180/pi + 4;
    disp('mercator')
    m_proj ('mercator','longitude',[lo0 lo1],'latitude',[la0 la1]);
  else
   m_proj('Gnomonic','lon',tra_lon,'lat',tra_lat,'rad',1.0*dlon*180/pi,'rec','on')
  end

  out_lon = [lone(2,2:end-1) lone(2:end-1,end-1)' lone(end-1,end-1:-1:2) lone(end-1:-1:2,2)'];
  out_lat = [late(2,2:end-1) late(2:end-1,end-1)' late(end-1,end-1:-1:2) late(end-1:-1:2,2)'];
%  out_lon = [out_lon 0.5*squeeze(lon(2:end-1,end)+lon(2:end-1,end-1))'];
%  out_lat = [out_lat 0.5*squeeze(lat(2:end-1,end)+lat(2:end-1,end-1))'];
%  out_lon(out_lon>0) = out_lon(out_lon>0) - 2*pi;
   

switch popup_sel_index
    case 1
      m_plot(out_lon*180/pi,out_lat*180/pi,'r')
      m_grid
      hold on
      m_gshhs_i('patch',[.7 .7 .7],'edgecolor','k')
      m_plot(out_lon*180/pi,out_lat*180/pi,'r')
%     m_coast('patch',[.7 .7 .7],'edgecolor','k');
%     m_plot(clon,clat,'k')
%     m_plot(clon-180, -clat,'r')
%     m_plot(out_lon*180/pi,out_lat*180/pi,'r')
      hold off
    case 2
%     m_coast('patch',[.7 .7 .7],'edgecolor','k');
      m_plot(lon*180/pi,lat*180/pi,'.g')
      m_grid
      hold on
%     m_coast('patch',[.7 .7 .7],'edgecolor','k');
      m_gshhs_i('patch',[.7 .7 .7],'edgecolor','k')
      m_plot(lone(2:end-1,2:end-1)*180/pi,late(2:end-1,2:end-1)*180/pi,'.b')
      m_plot(out_lon*180/pi,out_lat*180/pi,'r')
%     m_plot(clon,clat,'k')
      hold off
    case 3

%%    TODO:  be a bit smarter about extending the topo data set.
%%    TODO:  Also, make sure to take a sub-set of the data only.
      toponame = 'etopo5.nc'; 
      x = ncread(toponame,'topo_lon');
      y = ncread(toponame,'topo_lat');
      d = single(ncread(toponame,'topo'));
      x(find(x<0)) = x(find(x<0)) +360;
      xm = x-360;
      x = [xm' x']';
      d = [d' d'];
      di = interp2(x,y,d,lon*180/pi,lat*180/pi);
      di(di>0) = 0.;

      die  = zeros(ny+3,nx+3);
      die(1:end-1,1:end-1)  = di;

%     m_pcolor(lon*180/pi,lat*180/pi,di);shading flat;colorbar
      m_pcolor(lone*180/pi,late*180/pi,die);shading flat;colorbar
      m_grid
      hold on
      m_coast('patch',[.7 .7 .7],'edgecolor','k');
      m_plot(out_lon*180/pi,out_lat*180/pi,'r')
%     m_plot(lon*180/pi,lat*180/pi,'.g')
%     m_plot(lone*180/pi,late*180/pi,'.b')
%     m_plot(out_lon*180/pi,out_lat*180/pi,'r')
      hold off
    case 4
      m_pcolor(lon*180/pi,lat*180/pi,pm);shading flat;colorbar
      m_grid
      m_coast('patch',[.7 .7 .7],'edgecolor','k');
    case 5
      m_pcolor(lon*180/pi,lat*180/pi,pn);shading flat;colorbar
      m_grid
      m_coast('patch',[.7 .7 .7],'edgecolor','k');
    case 6
      m_pcolor(lon*180/pi,lat*180/pi,ang*180/pi);shading flat;colorbar
      m_grid
      m_coast('patch',[.7 .7 .7],'edgecolor','k');
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)

file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)

printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)

selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --------------------------------------------------------------------
function popupmenu1_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function popupmenu1_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

set(hObject, 'String', {'outline', 'grid', 'topo', 'pm', 'pn', 'angle'});


% --------------------------------------------------------------------
function edit1_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function edit1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function edit2_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function edit2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function edit3_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function edit3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function edit4_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function edit4_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function edit5_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function edit5_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function pushbutton4_Callback(hObject, eventdata, handles)
% --- Executes on button press in pushbutton4.

   size_x = str2double(get(handles.edit1,'String'))*1e3; %% Width 
   size_y = str2double(get(handles.edit2,'String'))*1e3; 
   rotate = str2double(get(handles.edit3,'String'));
   if abs(rotate) > 45
     rotate = rotate/abs(rotate)*45;
     set(handles.edit3,'String',num2str( rotate));   %% Rotation
   end
   tra_lon= str2double(get(handles.edit4,'String'));
   tra_lat= str2double(get(handles.edit5,'String'));
   nx     = str2num(get(handles.edit6,'String'));
   ny     = str2num(get(handles.edit7,'String')); 
   [lon,lat,pm,pn,ang,lone,late] = easy_grid(nx,ny,size_x,size_y,tra_lon,tra_lat,rotate);
   %% Getting topo
   toponame = 'etopo5.nc'; 
   x = ncread(toponame,'topo_lon');
   y = ncread(toponame,'topo_lat');
   d = double(ncread(toponame,'topo'));
   x(find(x<0)) = x(find(x<0)) +360;
   xm = x-360;
   x = [xm' x']';
   d = [d' d'];
   di = interp2(x,y,d,lon*180/pi,lat*180/pi);
   %% Make the grid file
   make_grid(nx,ny,lon,lat,pn,pm,di,ang,size_x,size_y,rotate,tra_lon,tra_lat,lone,late);


function edit6_Callback(hObject, eventdata, handles)

function edit6_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit7_Callback(hObject, eventdata, handles)

function edit7_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


