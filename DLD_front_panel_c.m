function varargout = DLD_front_panel_c(varargin)
% DLD_FRONT_PANEL_C M-file for DLD_front_panel_c.fig
%      DLD_FRONT_PANEL_C, by itself, creates a new DLD_FRONT_PANEL_C or raises the existing
%      singleton*.
%
%      H = DLD_FRONT_PANEL_C returns the handle to a new DLD_FRONT_PANEL_C or the handle to
%      the existing singleton*.
%
%      DLD_FRONT_PANEL_C('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DLD_FRONT_PANEL_C.M with the given input arguments.
%
%      DLD_FRONT_PANEL_C('Property','Value',...) creates a new
%      DLD_FRONT_PANEL_C or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DLD_front_panel_c_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DLD_front_panel_c_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DLD_front_panel_c

% Last Modified by GUIDE v2.5 26-Nov-2013 15:07:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DLD_front_panel_c_OpeningFcn, ...
                   'gui_OutputFcn',  @DLD_front_panel_c_OutputFcn, ...
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
end

% --- Executes just before DLD_front_panel_c is made visible.
function DLD_front_panel_c_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DLD_front_panel_c (see VARARGIN)

% Choose default command line output for DLD_front_panel_c
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DLD_front_panel_c wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%[a, MSGID] = lastwarn()
%warning('off', MATLAB:MKDIR:DirectoryExists)
warning('off', 'MKDIR:DirectoryExists')

handles.txy_data_all = [];
handles.txy_data_plot = [];
handles.can_replot = 1;

set(handles.fit_TF_average, 'Enable', 'off');
set(handles.fit_TF_temp_check, 'Enable', 'off');
set(handles.fit_spatial_average, 'Enable', 'off');
set(handles.gpu_checkbox,'Enable', 'off')

set(handles.plot_button,'Enable','off');

handles = repopulate_fp('config123.txt',hObject,handles,1);

get(handles.lock_windows_check,'Value');

if(get(handles.lock_windows_check,'Value'))
   handles.is_window_locked = 1;
else
   handles.is_window_locked = 0;
end

%handles.is_window_locked = 0;

%refresh(DLD_front_panel_c);
guidata(hObject, handles);
end












% --- Executes on button press in calculate_button.
function calculate_button_Callback(hObject, eventdata, handles)
% hObject    handle to calculate_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = read_fp_inputs(hObject,handles);
write_config(hObject,handles);

set(DLD_front_panel_c, 'HandleVisibility', 'off');

close all

fprintf('\n \n \n \n \n')
disp('=================================================')
disp('                beginning conversion')
disp('=================================================')
if handles.x64_status == '0'
    disp('32 bit')
    system('DLD_convert_and_correlations_for_MATLAB.exe');
elseif handles.x64_status == '1'
    disp('64 bit')
    system('DLD_convert_and_correlations_for_MATLAB_x64.exe');    
end

% CHECKING CORNERS OF DLD DATA
channels = dlmread([handles.filename_input,'_channels.txt']);
if (min(channels(1:4,1))/max(channels(1:4,1)) < 0.8 || channels(5,1)/min(channels(1:4,1)) > 0.1)
    set(handles.channels_text,'ForegroundColor','r');
else
    set(handles.channels_text,'ForegroundColor','k');
end
set(handles.channels_text,'String',['[Ch0, Ch1, Ch2, Ch3, Rubbish] = [',num2str(channels(1,1)),',',num2str(channels(2,1)),',',num2str(channels(3,1)),',',num2str(channels(4,1)),',',num2str(channels(5,1)),']']);


% plot original NVF unwindowed from C++
if (str2num(handles.nvf) == 1)
    
    ncf_cpp = dlmread([handles.filename_input,'_nvf.txt']);
    
    Mean_counts_per_file = mean(ncf_cpp(:,1));
    Std_counts_per_file = std(ncf_cpp(:,1));
    
    fig_1_note1 = {['Average hits per file = ',num2str(Mean_counts_per_file)]};
    fig_1_note2 = {['St. dev. hits per file = ',num2str(Std_counts_per_file)]};
    
    nvf_xaxis = [str2num(handles.first_file_input):str2num(handles.last_file_input)];
    
    if (max(size(ncf_cpp) == max(size(nvf_xaxis))))
        figure(1)
        NVF_plot = plot(nvf_xaxis,ncf_cpp);
        xlabel('File number')
        ylabel('Number of counts')
        title('Hits per file, no windows')
        uicontrol('Style', 'text','String', fig_1_note1,'Units','normalized','Position', [0.2 .11 0.28 0.035],'BackgroundColor','w');
        uicontrol('Style', 'text','String', fig_1_note2,'Units','normalized','Position', [0.6 .11 0.28 0.035],'BackgroundColor','w');
    else
        disp('NVF sizes not correct')
        size(nvf_xaxis)
        size(ncf_cpp)
    end
end


% plot original TOF unwindowed from C++
if (str2num(handles.tof) == 1)    
    tof_data = [];
    fileID = fopen([handles.filename_input,'_TOF.txt']);
    tof_data_cell = textscan(fileID, '%u32 ','delimiter', ';');
    tof_data = tof_data_cell{1};
    fclose(fileID);
    
    tof_data_double = double(tof_data);
        
    num_files_to_plot = (str2num(handles.last_file_input) - str2num(handles.first_file_input)+1);
       
    tof_xaxis = [str2num(handles.t_min) + 1e-5:1e-5:str2num(handles.t_max)];
    
    close(figure(2))
    
    if (max(size(tof_data)) == max(size(tof_xaxis)))
        figure(2)
        %plot(tof_xaxis,tof_data./num_files_to_plot./10)
         plot(tof_xaxis,tof_data_double./num_files_to_plot/10)
        xlabel('t(s)')
        ylabel('Average flux (MHz). Saturates at 1 MHz.')
        title('All files TOF at 10us resolution, no windowing')
    else
        disp('TOF window not same as original file from c++')
    end
end


% plot original spatial image unwindowed from C++
if (str2num(handles.spatial) == 1)
    spatial_image = [];
    fileID = fopen([handles.filename_input,'_spatial.txt']);
    spatial_image_cell = textscan(fileID, '%u32 ','delimiter', ',');
    spatial_image(:,1) = spatial_image_cell{1};
    fclose(fileID);
    
    spatial_image_square = vec2mat(spatial_image,607);
    spatial_image_square = spatial_image_square'/(str2num(handles.last_file_input) - str2num(handles.first_file_input)+1);
    figure(4)
    axis = [-40:80/606:40];
    spatial_plot_fig = pcolor(axis,axis,spatial_image_square);
    set(spatial_plot_fig,'EdgeColor','none')
    xlabel('x position (mm)');
    ylabel('y position (mm)');
    title('All files spatial image at DLD resolution (132 micron per bin), no windowing') 
    colormap gray
    % colorbar 
end

%handles = window_and_fit(hObject,handles,str2num(handles.nvf)); % the one means do each file separately
handles = window_and_fit(hObject,handles,0);

%set(handles.plot_button,'Enable','On');

%refresh(DLD_front_panel_c);
guidata(hObject, handles);
end





function y = fit_2d_gaussian(Param,binned_output)

y = 0;
size_vec = size(binned_output);
max_y = size_vec(2);
max_x = size_vec(1);

for xpos = 1:max_x
    for ypos = 1:max_y
        y = y + (Param(1)*exp(-0.5*((xpos-Param(5)).^2./(Param(2)).^2))*exp(-0.5*((ypos-Param(6)).^2./(Param(3)).^2)) + Param(4) - binned_output(xpos,ypos)).^2;
    end
end
end











% --- Executes on button press in plot_button.
function plot_button_Callback(hObject, eventdata, handles)
% hObject    handle to plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

first_file_input = get(handles.first_file_edit,'String');   %first_file_input
last_file_input = get(handles.last_file_edit,'String');    %last_file_input
filename_input = get(handles.filepath_edit,'String');

if (str2num(first_file_input) ~= handles.first_file_convert || str2num(last_file_input) ~= handles.last_file_convert || handles.can_replot == 0 || get(handles.nvf_check,'Value'))
    handles = window_and_fit(hObject,handles,1);  
else
    handles = window_and_fit(hObject,handles,0);
end

guidata(hObject, handles);
end



% --- Executes on button press in browse_button.
function browse_button_Callback(hObject, eventdata, handles)
% hObject    handle to browse_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,filepath,click] = uigetfile({'*.txt;*.bin;','Text or Binary files';},'Choose any raw file from data set to load',handles.filename_browse);
if (click == 1)
    filename = fliplr([filepath,filename]);
    [~,filename] = strtok(filename,'.');
    [~,filename] = strtok(filename(2:end),[char(32:47) char(58:64) char(65:90) char(91:96) char(97:122) char(123:126)]);
    set(handles.filepath_edit,'String',fliplr(filename));
    
    [~,filename] = strtok(filename,'\');
    handles.filename_browse = fliplr(filename(2:end));
end

guidata(hObject, handles);

end

% --- Executes on button press in load_old_button.
function load_old_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_old_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.filename_reload

[filename] = uigetdir(handles.filename_reload,'Choose folder from previously calculated settings to load');

if (filename ~= 0)
    handles.filename_reload = filename;
    
    if exist([filename,'\spatial_rebinned.fig'], 'file')
        open([filename,'\spatial_rebinned.fig'])
    end
    if exist([filename,'\NVF_plot.fig'], 'file')
        open([filename,'\NVF_plot.fig'])
    end
    if exist([filename,'\TOF_rebinned.fig'], 'file')
        open([filename,'\TOF_rebinned.fig'])
    end
    
    copyfile([filename,'\config123.txt'],'config123.txt');
    
    repopulate_fp('config123.txt',hObject,handles,0);
    
    handles.can_replot = 0;
end

guidata(hObject, handles);

end




function newhandles = repopulate_fp(filename,hObject,handles,redoreload)

% INITIALISE VALUES WITH PREVIOUS SAVED CONFIG FILE
fileID = fopen(filename,'r');

default_inputs = textscan(fileID,'%s', 'delimiter',';');
%celldisp(default_inputs)

filename_input_fixed = strrep(default_inputs{1}{1}, '\\', '\');
set(handles.filepath_edit,'String',filename_input_fixed);
set(handles.first_file_edit,'String',default_inputs{1}{2});
set(handles.last_file_edit,'String',default_inputs{1}{3});
set(handles.num_CPU_edit,'String',default_inputs{1}{4});
if (default_inputs{1}{5} == '1')
    set(handles.gpu_checkbox,'Value',1)
else
    set(handles.gpu_checkbox,'Value',0)
end
%skip - convert flag
set(handles.t_min_edit,'String',default_inputs{1}{7});
set(handles.t_max_edit,'String',default_inputs{1}{8});
set(handles.x_min_edit,'String',default_inputs{1}{9});
set(handles.x_max_edit,'String',default_inputs{1}{10});
set(handles.y_min_edit,'String',default_inputs{1}{11});
set(handles.y_max_edit,'String',default_inputs{1}{12});

if (default_inputs{1}{13} == '1')
    set(handles.dt_radio,'Value',1);
elseif (default_inputs{1}{13} == '2')
    set(handles.dx_radio,'Value',1);
elseif (default_inputs{1}{13} == '3')
    set(handles.dy_radio,'Value',1);
end

if (default_inputs{1}{14} == '0')
    set(handles.g0_radio,'Value', 1);
elseif (default_inputs{1}{14} == '2')
    set(handles.g2_radio,'Value',1);
elseif(default_inputs{1}{14} == '3')
    set(handles.g3_radio,'Value',1);
end

set(handles.t_bin_edit,'String',default_inputs{1}{15});
set(handles.x_bin_edit,'String',default_inputs{1}{16});
set(handles.y_bin_edit,'String',default_inputs{1}{17});
set(handles.number_bins_edit,'String',default_inputs{1}{18});   

if (default_inputs{1}{19} == '1')
    set(handles.tof_check,'Value',1)
else
    set(handles.tof_check,'Value',0)
end

if (default_inputs{1}{20} == '1')
    set(handles.spatial_check,'Value',1)
else
    set(handles.spatial_check,'Value',0)
end

if (default_inputs{1}{21} == '1')
    set(handles.nvf_check,'Value',1)
else
    set(handles.nvf_check,'Value',0)
end

if (default_inputs{1}{22} == '1')
    set(handles.fit_remove_TOF_check,'Value',1)
else
    set(handles.fit_remove_TOF_check,'Value',0)
end

if (default_inputs{1}{23} == '1')
    set(handles.fit_T_check,'Value',1)
else
    set(handles.fit_T_check,'Value',0)
end

set(handles.files_per_chunk_edit,'String',default_inputs{1}{24});
set(handles.temp_window_edit,'String',default_inputs{1}{25});
set(handles.num_window_edit,'String',default_inputs{1}{26});
set(handles.t_bin_plot_edit,'String',default_inputs{1}{27});
set(handles.bins_per_pixel_edit,'String',default_inputs{1}{28});
set(handles.dead_time_edit,'String',num2str(str2num(default_inputs{1}{29})/1e-9*25e-12));

if (default_inputs{1}{30} == '1')
    set(handles.exclusion_check,'Value',1)
else
    set(handles.exclusion_check,'Value',0)
end

set(handles.t_min_exc_edit,'String',default_inputs{1}{31});
set(handles.t_max_exc_edit,'String',default_inputs{1}{32});
set(handles.x_min_exc_edit,'String',default_inputs{1}{33});
set(handles.x_max_exc_edit,'String',default_inputs{1}{34});
set(handles.y_min_exc_edit,'String',default_inputs{1}{35});
set(handles.y_max_exc_edit,'String',default_inputs{1}{36});
handles.filename_browse = default_inputs{1}{37};
if (redoreload == 1)
    handles.filename_reload = default_inputs{1}{38};
else
    handles.filename_reload = '';
end

if (default_inputs{1}{39} == '1')
    set(handles.fit_TF_temp_check,'Value',1)
else
    set(handles.fit_TF_temp_check,'Value',0)
end

set(handles.f_x_edit,'String',default_inputs{1}{40});
set(handles.f_y_edit,'String',default_inputs{1}{41});
set(handles.f_z_edit,'String',default_inputs{1}{42});
set(handles.temp_guess_tf_edit,'String',default_inputs{1}{43});
set(handles.tf_radius_guess_edit,'String',default_inputs{1}{44});
set(handles.cond_percent_guess_edit,'String',default_inputs{1}{45});

if (default_inputs{1}{46} == '1')
    set(handles.force_recalc_txy_check,'Value',1)
else
    set(handles.force_recalc_txy_check,'Value',0)
end

set(handles.num_files_denom_edit,'String',default_inputs{1}{48});

if (default_inputs{1}{49} == '1')
    set(handles.sort_by_temp_check,'Value',1)
else
    set(handles.sort_by_temp_check,'Value',0)
end
% default_inputs{1}{47} is ready_for_gn

if (default_inputs{1}{50} == '1')
    set(handles.fit_remove_x0_check,'Value',1)
else
    set(handles.fit_remove_x0_check,'Value',0)
end

if (default_inputs{1}{51} == '1')
    set(handles.fit_remove_y0_check,'Value',1)
else
    set(handles.fit_remove_y0_check,'Value',0)
end

if (default_inputs{1}{52} == '1')
    set(handles.x64_check,'Value',1)
else
    set(handles.x64_check,'Value',0)
end

if (default_inputs{1}{53} == '1')
    set(handles.fit_thermal_average_check,'Value',1)
else
    set(handles.fit_thermal_average_check,'Value',0)
end

if (default_inputs{1}{54} == '1')
    set(handles.custom_file_list_check,'Value',1)
else
    set(handles.custom_file_list_check,'Value',0)
end

set(handles.custom_file_list_string_edit,'String',default_inputs{1}{55});
set(handles.rotation_edit,'String',default_inputs{1}{56});


newhandles = handles;

fclose(fileID);

guidata(hObject, handles);

end



function newhandles = read_fp_inputs(hObject,handles)

% GET VALUES FROM USER INPUT

handles.filename_input = get(handles.filepath_edit,'String');
handles.filename_input_new = strrep(handles.filename_input, '\', '\\\\');   %for some reason need to double the double so c++ is happy

handles.first_file_input = get(handles.first_file_edit,'String');   %first_file_input
handles.first_file_convert = str2num(handles.first_file_input);
handles.last_file_input = get(handles.last_file_edit,'String');   %last_file_input
handles.last_file_convert = str2num(handles.last_file_input);

handles.number_CPU = get(handles.num_CPU_edit,'String');
if (get(handles.gpu_checkbox,'Value')) %plot TOF profile
    handles.GPU_flag = num2str(1);
else
    handles.GPU_flag = num2str(0);
end

handles.t_min = get(handles.t_min_edit,'String');
handles.t_max = get(handles.t_max_edit,'String');
handles.x_min = get(handles.x_min_edit,'String');
handles.x_max = get(handles.x_max_edit,'String');
handles.y_min = get(handles.y_min_edit,'String');
handles.y_max = get(handles.y_max_edit,'String');

if (get(handles.g0_radio,'Value') == 1)
	order = 0;
elseif (get(handles.g2_radio,'Value') == 1)
	order = 2;
elseif (get(handles.g3_radio,'Value') == 1)
	order = 3;
end

if (get(handles.dt_radio,'Value') == 1)
	axis = 1;
elseif (get(handles.dx_radio,'Value') == 1)
	axis = 2;
elseif (get(handles.dy_radio,'Value') == 1)
	axis = 3;
end

handles.order_s = num2str(order);
handles.axis_s = num2str(axis);

handles.t_bin = get(handles.t_bin_edit,'String');
handles.x_bin = get(handles.x_bin_edit,'String');
handles.y_bin = get(handles.y_bin_edit,'String');
handles.num_bins = get(handles.number_bins_edit,'String');

if (get(handles.tof_check,'Value')) %plot TOF profile
    handles.tof = num2str(1);
else
    handles.tof = num2str(0);
end

if (get(handles.spatial_check,'Value')) %plot spatial profile
    handles.spatial = num2str(1);
else
    handles.spatial = num2str(0);
end

if (get(handles.nvf_check,'Value')) %plot number versus file
    handles.nvf = num2str(1);
else
    handles.nvf = num2str(0);
end

if (get(handles.fit_remove_TOF_check,'Value')) %fit and remove t0
    handles.fit_remove_TOF = num2str(1);
else
    handles.fit_remove_TOF = num2str(0);
end

if (get(handles.fit_T_check,'Value')) %fit and remove t0
    handles.fit_T = num2str(1);
else
    handles.fit_T = num2str(0);
end

if (get(handles.fit_thermal_average_check,'Value')) %fit and remove t0
    handles.fit_thermal_average = num2str(1);
else
    handles.fit_thermal_average = num2str(0);
end

handles.files_per_chunk = get(handles.files_per_chunk_edit,'String');
handles.temp_window = get(handles.temp_window_edit,'String');
handles.num_window = get(handles.num_window_edit,'String');

handles.t_bin_plot = get(handles.t_bin_plot_edit,'String');
handles.bins_per_pixel = get(handles.bins_per_pixel_edit,'String');
handles.bin_size_disp = str2num(handles.bins_per_pixel)*0.132;
%set(handles.pixel_size_txt,'String',num2str(bin_size_disp));

handles.dead_time = get(handles.dead_time_edit,'String');

if (get(handles.exclusion_check,'Value')) % use exclusion zone
    handles.exclusion_active =  num2str(1);
else
    handles.exclusion_active = num2str(0);
end

handles.t_min_exc = get(handles.t_min_exc_edit,'String');
handles.t_max_exc = get(handles.t_max_exc_edit,'String');
handles.x_min_exc = get(handles.x_min_exc_edit,'String');
handles.x_max_exc = get(handles.x_max_exc_edit,'String');
handles.y_min_exc = get(handles.y_min_exc_edit,'String');
handles.y_max_exc = get(handles.y_max_exc_edit,'String');

if (get(handles.force_recalc_txy_check,'Value'))
    handles.force_recalc_txy = num2str(1);
else
    handles.force_recalc_txy = num2str(0);
end

handles.convert_flag = 1;  %convert_flag

%fileInfo = dir([filename_input,'_txy_forc_AGM_',first_file_input,'.txt']);

if (exist([handles.filename_input,'_txy_forc_AGM_',handles.first_file_input,'.txt'],'file'))
    if (exist([handles.filename_input,'_txy_forc_AGM_',handles.last_file_input,'.txt'],'file'))
        if (handles.force_recalc_txy == '0')
            disp('no need to convert')
            handles.convert_flag = 0;
        end
    end
end

handles.convert_flag_s = num2str(handles.convert_flag);

handles.dead_time_str = num2str(str2num(handles.dead_time)*1e-9/25e-12);

if (get(handles.fit_TF_temp_check,'Value')) % use exclusion zone
    handles.TF_fit =  num2str(1);
else
    handles.TF_fit = num2str(0);
end

handles.f_x = get(handles.f_x_edit,'String');
handles.f_y = get(handles.f_y_edit,'String');
handles.f_z = get(handles.f_z_edit,'String');
handles.temp_guess_tf = get(handles.temp_guess_tf_edit,'String');
handles.tf_radius_guess = get(handles.tf_radius_guess_edit,'String');
handles.cond_percent_guess = get(handles.cond_percent_guess_edit,'String');

if (get(handles.force_recalc_txy_check,'Value')) % override flag to skip convert
    handles.force_recalc_txy =  num2str(1);
else
    handles.force_recalc_txy = num2str(0);
end

handles.ready_for_gn = '0';

handles.num_files_denom = get(handles.num_files_denom_edit,'String');

if (get(handles.sort_by_temp_check,'Value')) % override flag to skip convert
    handles.sort_by_temp =  num2str(1);
else
    handles.sort_by_temp = num2str(0);
end

if (get(handles.fit_remove_x0_check,'Value')) % override flag to skip convert
    handles.fit_remove_x0 =  num2str(1);
else
    handles.fit_remove_x0 = num2str(0);
end

if (get(handles.fit_remove_y0_check,'Value')) % override flag to skip convert
    handles.fit_remove_y0 =  num2str(1);
else
    handles.fit_remove_y0 = num2str(0);
end

if (get(handles.x64_check,'Value'))
    handles.x64_status =  num2str(1);
else
    handles.x64_status = num2str(0);
end

if (get(handles.custom_file_list_check,'Value'))
    handles.custom_file_list =  num2str(1);
else
    handles.custom_file_list = num2str(0);
end

handles.custom_file_list_string = get(handles.custom_file_list_string_edit,'String');
handles.rotation_string = get(handles.rotation_edit,'String');

newhandles = handles;
guidata(hObject, handles);

end


function write_config(hObject,handles)

% NOW WRITE CONFIG FILE READ BY C++ PROGRAM

fileID = fopen('config123.txt','wt');
fprintf(fileID,[handles.filename_input_new,'; \n']);
fprintf(fileID,[handles.first_file_input,'; \n']);
fprintf(fileID,[handles.last_file_input,'; \n']);
fprintf(fileID,[handles.number_CPU,'; \n']);
fprintf(fileID,[handles.GPU_flag,'; \n']);
fprintf(fileID,[handles.convert_flag_s,'; \n']);
fprintf(fileID,[handles.t_min,'; \n']);
fprintf(fileID,[handles.t_max,'; \n']);
fprintf(fileID,[handles.x_min,'; \n']);
fprintf(fileID,[handles.x_max,'; \n']);
fprintf(fileID,[handles.y_min,'; \n']);
fprintf(fileID,[handles.y_max,'; \n']);
fprintf(fileID,[handles.axis_s,'; \n']);
fprintf(fileID,[handles.order_s,'; \n']);
fprintf(fileID,[handles.t_bin,'; \n']);
fprintf(fileID,[handles.x_bin,'; \n']);
fprintf(fileID,[handles.y_bin,'; \n']);
fprintf(fileID,[handles.num_bins,'; \n']);
fprintf(fileID,[handles.tof,'; \n']);
fprintf(fileID,[handles.spatial,'; \n']);
fprintf(fileID,[handles.nvf,'; \n']);
fprintf(fileID,[handles.fit_remove_TOF,'; \n']);
fprintf(fileID,[handles.fit_T,'; \n']);
fprintf(fileID,[handles.files_per_chunk,'; \n']);
fprintf(fileID,[handles.temp_window,'; \n']);
fprintf(fileID,[handles.num_window,'; \n']);
fprintf(fileID,[handles.t_bin_plot,'; \n']);
fprintf(fileID,[handles.bins_per_pixel,'; \n']);
fprintf(fileID,[handles.dead_time_str,'; \n']);
fprintf(fileID,[handles.exclusion_active,'; \n']);
fprintf(fileID,[handles.t_min_exc,'; \n']);
fprintf(fileID,[handles.t_max_exc,'; \n']);
fprintf(fileID,[handles.x_min_exc,'; \n']);
fprintf(fileID,[handles.x_max_exc,'; \n']);
fprintf(fileID,[handles.y_min_exc,'; \n']);
fprintf(fileID,[handles.y_max_exc,'; \n']);
handles.filename_browse_write = strrep(handles.filename_browse, '\', '\\');
fprintf(fileID,[handles.filename_browse_write,'; \n']);
handles.filename_reload_write = strrep(handles.filename_reload, '\', '\\');
fprintf(fileID,[handles.filename_reload_write,'; \n']);
fprintf(fileID,[handles.TF_fit,'; \n']);
fprintf(fileID,[handles.f_x,'; \n']);
fprintf(fileID,[handles.f_y,'; \n']);
fprintf(fileID,[handles.f_z,'; \n']);
fprintf(fileID,[handles.temp_guess_tf,'; \n']);
fprintf(fileID,[handles.tf_radius_guess,'; \n']);
fprintf(fileID,[handles.cond_percent_guess,'; \n']);
fprintf(fileID,[handles.force_recalc_txy,'; \n']);
fprintf(fileID,[handles.ready_for_gn,'; \n']);
fprintf(fileID,[handles.num_files_denom,'; \n']);
fprintf(fileID,[handles.sort_by_temp,'; \n']);
fprintf(fileID,[handles.fit_remove_x0,'; \n']);
fprintf(fileID,[handles.fit_remove_y0,'; \n']);
fprintf(fileID,[handles.x64_status,'; \n']);
fprintf(fileID,[handles.fit_thermal_average,'; \n']);
fprintf(fileID,[handles.custom_file_list,'; \n']);
fprintf(fileID,[handles.custom_file_list_string,'; \n']);
fprintf(fileID,[handles.rotation_string,'; \n']);
fclose(fileID);

end







function [t0,Temperature,error,cond_percent,fit_out] = fit_Temp_t0(handles,time_vec,TOF,TF_fit,show_fits,T_guess,TF_radius_guess,cond_percent_guess)

t0 = -999;
Temperature = -999;
error = -999;
cond_percent = 0;
fit_out = [];

max_countrate = max(TOF);

TOF = TOF./max_countrate;
[~,peak] = max(TOF);
length_TOF = max(size(TOF));

time_of_peak = time_vec(peak);

time_vec = time_vec - time_of_peak;

if (TF_fit == '0')
    
    for i =1:length_TOF
        if TOF(i)> exp(-1)
            fw1=time_vec(i);
            j=i;
            break
        end
    end
    
    for i =j+3:length_TOF
        if TOF(i)< exp(-1)
            fw2=time_vec(i);
            break
        end
    end
    
    FW1_e=(abs(fw1)+abs(fw2)); % approximate FWHM

    guess(1)=1;
    guess(2)=time_vec(peak);
    guess(3)=FW1_e;
    guess(4)=0;
    
    [fit_thermal_out,error_thermal] = fminsearch(@(Param)thermal_TOF_dist(Param,time_vec,TOF),guess,optimset('MaxFunEvals',30000,'MaxIter',20000));     
    
    error = error_thermal;
    
    v0=9.8*fit_thermal_out(3)/2;
    Temperature=6.64e-27*v0^2/(2*1.38e-23);

    t0 = fit_thermal_out(2) + time_of_peak;
    
    if (show_fits == 1)
        figure(12345)
        y=fit_thermal_out(1)*exp(-((2*(time_vec-fit_thermal_out(2)))/fit_thermal_out(3)).^2)+fit_thermal_out(4);
        plot(time_vec,TOF,'b.',time_vec,y,'r')
        xlabel('Time')
    end
    
    fit_thermal_out(1) = fit_thermal_out(1)*max_countrate;
    
    fit_out = fit_thermal_out;
    
elseif (TF_fit == '1')

    disp('TF fit not yet working')
    
    return
    
    disp('BEC')

    guess(1)=0;   % background
    guess(2)=T_guess;    % T guess
    guess(3)=1 - cond_percent_guess/100; % thermal amplitude
    guess(4)=cond_percent_guess/100; % condensate amplitude
    guess(5)=TF_radius_guess;  % TF user input radius guess
    guess(6)=0;   % t0
    vel=4.0768; % 9.8m/s/s * 0.416s
    wz=2*pi*str2num(handles.f_z);
    wr=2*pi*sqrt(str2num(handles.f_x)*str2num(handles.f_y));
    nmax = 200; % number of thermal states to sum over (AGT sets to 200)
    tau = 0.416;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% thermal
    
    z0=sqrt(((2*1.38e-23*guess(2))/(6.64e-27*wr^2))*(1+(wr^2)*0.416^2));
    thermal_state_vec=1:1:nmax;
    time_vec=time_vec+guess(6);
    z_vec=vel*time_vec; %convert flight into distance
    
    thermID = zeros(length(time_vec),1);
    for i=1:length(time_vec)
        thermID(i)=guess(3)*sum(((exp(-(z_vec(i)^2/z0^2))).^thermal_state_vec)./thermal_state_vec.^(5/2));
    end
    
    [TF_thermal_fit,Temp_fit_TF,error_thermal,fit_thermal_out] = TF_thermal(time_vec,thermID);
    
    ID = thermID;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% condensate

for i=1:length(time_vec)
    if time_vec(i)>-fit_thermal_out(5) && time_vec(i)<fit_thermal_out(5)
        ID(i)=ID(i)+(fit_thermal_out(4)*((1-z_vec(i)^2/(fit_thermal_out(5)*vel)^2)^2));
    end
end
condID=ID-thermID;
y=ID+fit_thermal_out(1);

tot=sum(condID)+sum(thermID);
T=fit_thermal_out(2);
thermfrac=sum(thermID)/tot;
condfrac=sum(condID)/tot;
mu=(0.5*6.64e-27*wr^2*(Temp_fit_TF*vel)^2/(1+(wr^2*0.416^2)))/(6.626e-34);

%%%%%%%%%%%%%%%%%%%% working out N0 from Chem. pot

wbar=(wz*wr^2)^(1/3);
aho=sqrt(hbar/(Hemass*wbar));

cp=0.5*6.64e-27*(wr^2/(1+wr^2*.416^2))*(Temp_fit_TF*vel)^2;
N01=((2*cp/(1.055e-34*wbar))^(5/2))*aho/(15*7.512e-9);

%%%%%%%%%%%%%%%%%%%% working out thermal fraction from critical temp

Ntherm=1.202*(1.38e-23*T/(1.055e-34*wbar))^3;
total=Ntherm/thermfrac;
N02=total*condfrac;
a=((2*cp/(hbar*wbar))^(5/2))*aho/(15*N02);

%     fitt=fminsearch('TC_and_below_agm',fp,OPTIONS,time_vec,TOF,vel,tau,wz,nmax);
%     [fit_BEC_out,error_BEC] = fminsearch(@(Param)TC_and_below_agm(Param,time_vec,TOF),guess,optimset('MaxFunEvals',30000,'MaxIter',20000));     
% 
%     %[fit_thermal_out,error_thermal] = fminsearch(@(Param)thermal_TOF_dist(Param,time_vec,TOF),guess,optimset('MaxFunEvals',30000,'MaxIter',20000));     
% 
%     
%     t0 = fp(6);
%     cond_percent = fp(4);
%     Temperature = fp(2);

Temp_fit_TF
mu
N01
N02

end
end















function [r0,countrate_fit,r_width,error] = fit_1D_spatial(handles,r_vec,r_profile)

max_countrate = max(r_profile);

r_profile = r_profile./max_countrate;
[~,peak] = max(r_profile);
length_r_profile = max(size(r_profile));
r_of_peak = r_vec(peak);
r_vec = r_vec - r_of_peak;

j=0;

for i =1:length_r_profile
    if r_profile(i)> exp(-1)
        fw1=r_vec(i);
        j=i;
        break
    end
end

for i =j+3:length_r_profile
    if r_profile(i)< exp(-1)
        fw2=r_vec(i);
        break
    end
end

FW1_e=(abs(fw1)+abs(fw2)); % approximate FWHM

guess(1)=1;
guess(2)=r_vec(peak);
guess(3)=FW1_e;
guess(4)=0;

[fit_thermal_out,error] = fminsearch(@(Param)thermal_1d_spatial_dist(Param,r_vec,r_profile),guess,optimset('MaxFunEvals',30000,'MaxIter',20000));

r_width = fit_thermal_out(3);
r0 = fit_thermal_out(2) + r_of_peak;
countrate_fit = fit_thermal_out(1)*max_countrate;
end





















function [TF_thermal_fit,Temp_fit_TF,error_thermal,fit_thermal_out]=TF_thermal(time_vec,TOF)

max_countrate = max(TOF);

TOF = TOF./max_countrate;
[~,peak] = max(TOF);
length_TOF = max(size(TOF));

time_vec = time_vec - time_vec(peak);
    
    for i =1:length_TOF
        if TOF(i)> exp(-1)
            fw1=time_vec(i);
            j=i;
            break
        end
    end
    
    for i =j+3:length_TOF
        if TOF(i)< exp(-1)
            fw2=time_vec(i);
            break
        end
    end
    
    FW1_e=(abs(fw1)+abs(fw2)); % approximate FWHM

    guess(1)=1;
    guess(2)=time_vec(peak);
    guess(3)=FW1_e;
    guess(4)=0;
    
    [fit_thermal_out,error_thermal] = fminsearch(@(Param)thermal_TOF_dist(Param,time_vec,TOF),guess,optimset('MaxFunEvals',30000,'MaxIter',20000));     
    
    %error = error_thermal;
    
    v0=9.8*fit_thermal_out(3)/2;
    Temp_fit_TF=6.64e-27*v0^2/(2*1.38e-23);
    TF_thermal_fit=fit_thermal_out(1)*exp(-((2*(time_vec-fit_thermal_out(2)))/fit_thermal_out(3)).^2)+fit_thermal_out(4);

end

function error = thermal_TOF_dist(Param,time_vec,TOF)
val = Param(1)*exp(-((2*(time_vec-Param(2)))/Param(3)).^2)+Param(4)-TOF';
error = sum(val.^2);
end

function error = thermal_1d_spatial_dist(Param,time_vec,TOF)
val = Param(1)*exp(-((time_vec-Param(2)).^2)/(2*((Param(3)).^2)))+Param(4)-TOF';
error = sum(val.^2);
end







function [newhandles] = window_and_fit(hObject,handles,keepeachfile)

handles.txy_data_plot = [];

first_file_input = get(handles.first_file_edit,'String');   %first_file_input
last_file_input = get(handles.last_file_edit,'String');    %last_file_input
filename_input = get(handles.filepath_edit,'String');

handles.t_min_val = str2num(get(handles.t_min_edit,'String'));
handles.t_max_val = str2num(get(handles.t_max_edit,'String'));
handles.x_min_val = str2num(get(handles.x_min_edit,'String'))*1e-3;
handles.x_max_val = str2num(get(handles.x_max_edit,'String'))*1e-3;
handles.y_min_val = str2num(get(handles.y_min_edit,'String'))*1e-3;
handles.y_max_val = str2num(get(handles.y_max_edit,'String'))*1e-3;

handles.t_min_exc_val = str2num(get(handles.t_min_exc_edit,'String'));
handles.t_max_exc_val = str2num(get(handles.t_max_exc_edit,'String'));
handles.x_min_exc_val = str2num(get(handles.x_min_exc_edit,'String'))*1e-3;
handles.x_max_exc_val = str2num(get(handles.x_max_exc_edit,'String'))*1e-3;
handles.y_min_exc_val = str2num(get(handles.y_min_exc_edit,'String'))*1e-3;
handles.y_max_exc_val = str2num(get(handles.y_max_exc_edit,'String'))*1e-3;

handles.t_bin_plot_val = str2num(get(handles.t_bin_plot_edit,'String'))*1e-3;

TOF_axis = handles.t_min_val:handles.t_bin_plot_val:handles.t_max_val;
TOF = zeros(max(size(TOF_axis)),1);
handles.TOF_all = zeros(max(size(TOF_axis)),1);
handles.TOF_centred_all = zeros(max(size(TOF_axis)),1);

bins_per_pixel_number = str2num(get(handles.bins_per_pixel_edit,'String'));
spatial_bin_size = 0.00013158*bins_per_pixel_number;
x_axis = handles.x_min_val:spatial_bin_size:handles.x_max_val;
y_axis = handles.y_min_val:spatial_bin_size:handles.y_max_val;
spatial_image = zeros(max(size(x_axis)),max(size(y_axis)));
handles.spatial_all = zeros(max(size(x_axis)),max(size(y_axis)));
x_profile_per_file = zeros(max(size(x_axis)),1);
y_profile_per_file = zeros(max(size(y_axis)),1);

size_nvf = str2num(handles.last_file_input) - str2num(handles.first_file_input) + 1;
handles.nvf_vec = zeros(size_nvf,1);
temperature_vec = zeros(size_nvf,1);
temperature_error_vec = zeros(size_nvf,1);
t0_vec = ones(size_nvf,1);
x0_vec = ones(size_nvf,1);
y0_vec = ones(size_nvf,1);
condensate_percent_vec = ones(size_nvf,1);

TOF_centred_all = zeros(max(size(TOF_axis)),1);
x_profile_centred_all = zeros(max(size(x_axis)),1);
y_profile_centred_all = zeros(max(size(y_axis)),1);

handles.txy_data_all = [];

file_list_to_use = str2num(handles.first_file_input):str2num(handles.last_file_input);
number_files_to_plot = str2num(handles.last_file_input) - str2num(handles.first_file_input) + 1;

if  (get(handles.custom_file_list_check,'Value') || get(handles.temp_window_edit,'String') ~= '0' || get(handles.num_window_edit,'String') ~= '0')
    file_list_to_use = dlmread([handles.filename_input,'_file_list.txt']);
    
    %%% get rid of duplicates
    file_list_to_use = unique(file_list_to_use);
    
    number_files_to_plot = max(size(file_list_to_use));
end



    
%for file = str2num(handles.first_file_input):str2num(handles.last_file_input)
% for file_dummy = 1:number_files_to_plot
for file_dummy = 1:number_files_to_plot
    
    file = file_list_to_use(file_dummy);
    
    filename_with_ext = [handles.filename_input,'_txy_forc_AGM_',num2str(file),'.txt'];   
    fileInfo = dir(filename_with_ext);
       
    if ~exist(filename_with_ext,'file')
        disp(['file missing ', filename_with_ext]);
       continue 
    end
    
    if (fileInfo.bytes > 1)
        
        skip_flag = 0;
        
        if (handles.is_window_locked == 0 || exist([handles.filename_input,'_txy_forc_AGM_correlations_fwrite_t',num2str(file),'.bin']) == 0)
            data_this_file = [];
            fileID = fopen(filename_with_ext);
            data_this_file_cell = textscan(fileID, '%f %f %f ','delimiter', ',');
            size_cell_1 = size(data_this_file_cell{1});
            size_cell_2 = size(data_this_file_cell{2});
            size_cell_3 = size(data_this_file_cell{3});
            
            if (size_cell_1(1) ~= size_cell_2(1) || size_cell_1(1) ~= size_cell_3(1))
                disp(['TXY data file size mismatch on file ',num2str(file)])
                skip_flag = 1;
                continue
            end

            %[file,fileInfo.bytes,size(data_this_file_cell{1}),size(data_this_file_cell{2}),size(data_this_file_cell{3})]
            
            data_this_file(:,1) = data_this_file_cell{1};
            data_this_file(:,2) = data_this_file_cell{2};
            data_this_file(:,3) = data_this_file_cell{3};
            fclose(fileID);
            
            % keep all unwindowed data in this matrix
            if (keepeachfile == 1)
                handles.txy_data_all = cat(1, handles.txy_data_all, data_this_file);
            end
            
            rowstokeep1 = data_this_file(:,1) < handles.t_max_val;
            rowstokeep2 = data_this_file(:,1) > handles.t_min_val;
            rowstokeep3 = data_this_file(:,2) < handles.x_max_val;
            rowstokeep4 = data_this_file(:,2) > handles.x_min_val;
            rowstokeep5 = data_this_file(:,3) < handles.y_max_val;
            rowstokeep6 = data_this_file(:,3) > handles.y_min_val;
            rowstokeep_all = (rowstokeep1&rowstokeep2&rowstokeep3&rowstokeep4&rowstokeep5&rowstokeep6);
            data_this_file = data_this_file(rowstokeep_all,:);
            
            if (get(handles.exclusion_check,'Value')) % use exclusion zone
                
                if ((size(handles.t_min_exc_val) == 1))
                    if ((size(handles.t_max_exc_val) == 1))
                        rowstokeep = (data_this_file(:,1) < handles.t_min_exc_val);
                        rowstokeep2 = (data_this_file(:,1) > handles.t_max_exc_val);
                        rowstokeep3 = logical(rowstokeep+rowstokeep2);
                        data_this_file = data_this_file(rowstokeep3,:);
                    end
                end
                
                %rowstokeep_xy = 0;
                rowstokeep3 = 0;
                rowstokeep4 = 0;
                
                if ((size(handles.x_min_exc_val) == 1))
                    if ((size(handles.x_max_exc_val) == 1))
                        rowstokeep = (data_this_file(:,2) < handles.x_min_exc_val);
                        rowstokeep2 = (data_this_file(:,2) > handles.x_max_exc_val);
                        rowstokeep3 = logical(rowstokeep+rowstokeep2);
                        %handles.txy_data_plot = handles.txy_data_plot(rowstokeep3,:);
                    end
                end
                
                if ((size(handles.y_min_exc_val) == 1))
                    if ((size(handles.y_max_exc_val) == 1))
                        rowstokeep = (data_this_file(:,3) < handles.y_min_exc_val);
                        rowstokeep2 = (data_this_file(:,3) > handles.y_max_exc_val);
                        rowstokeep4 = logical(rowstokeep+rowstokeep2);
                        %handles.txy_data_plot = handles.txy_data_plot(rowstokeep3,:);
                    end
                end
                
                if and(max(size(rowstokeep4))>1,max(size(rowstokeep3))>1)
                    rowstokeep_xy = logical(rowstokeep3+rowstokeep4);
                    data_this_file = data_this_file(rowstokeep_xy,:);
                elseif max(size(rowstokeep4))>1
                    data_this_file = data_this_file(rowstokeep4,:);
                elseif max(size(rowstokeep3))>1
                    data_this_file = data_this_file(rowstokeep3,:);
                end
            end
        end
        
        if skip_flag == 1
            continue
        end
        
        if (handles.is_window_locked == 1 && exist([handles.filename_input,'_txy_forc_AGM_correlations_fwrite_t',num2str(file),'.bin']) == 2)
            fid = fopen([handles.filename_input,'_txy_forc_AGM_correlations_fwrite_t',num2str(file),'.bin'], 'r');
            data_this_file(:,1) = fread(fid, 'double');
            fclose(fid);
            
            fid = fopen([handles.filename_input,'_txy_forc_AGM_correlations_fwrite_x',num2str(file),'.bin'], 'r');
            data_this_file(:,2) = fread(fid, 'double');
            fclose(fid);
            
            fid = fopen([handles.filename_input,'_txy_forc_AGM_correlations_fwrite_y',num2str(file),'.bin'], 'r');
            data_this_file(:,3) = fread(fid, 'double');
            fclose(fid);
        end
      
        %FITTING

        size_data_file = size(data_this_file);
        if size_data_file(1) > 0
            for hit = 1:(size(data_this_file,1))
                spatial_image(ceil((data_this_file(hit,2) - handles.x_min_val)/spatial_bin_size),ceil((data_this_file(hit,3) - handles.y_min_val)/spatial_bin_size)) = spatial_image(ceil((data_this_file(hit,2) - handles.x_min_val)/spatial_bin_size),ceil((data_this_file(hit,3) - handles.y_min_val)/spatial_bin_size)) + 1;
                TOF(ceil((data_this_file(hit,1) - handles.t_min_val)/handles.t_bin_plot_val)) = TOF(ceil((data_this_file(hit,1) - handles.t_min_val)/handles.t_bin_plot_val)) + 1;
                x_profile_per_file(ceil((data_this_file(hit,2) - handles.x_min_val)/spatial_bin_size)) = x_profile_per_file(ceil((data_this_file(hit,2) - handles.x_min_val)/spatial_bin_size)) + 1;
                y_profile_per_file(ceil((data_this_file(hit,3) - handles.y_min_val)/spatial_bin_size)) = y_profile_per_file(ceil((data_this_file(hit,3) - handles.y_min_val)/spatial_bin_size)) + 1;
            end
        end

        % Fit TOF - Temp and t0
        if (get(handles.fit_remove_TOF_check,'Value') || get(handles.fit_T_check,'Value'))
            
            handles.temp_guess_tf = get(handles.temp_guess_tf_edit,'String');
            handles.tf_radius_guess = get(handles.tf_radius_guess_edit,'String');
            handles.cond_percent_guess = get(handles.cond_percent_guess_edit,'String');
            
            showfits = 0;
            [t0,Temperature,error,cond_percent,fit_out] = fit_Temp_t0(handles,TOF_axis,TOF,handles.TF_fit,showfits,str2num(handles.temp_guess_tf)*1e-9,str2num(handles.tf_radius_guess)*1e-3,str2num(handles.cond_percent_guess));
            
            if (get(handles.fit_remove_TOF_check,'Value'))
                t0_vec(file - str2num(handles.first_file_input) + 1) = t0;
                data_this_file(:,1) = data_this_file(:,1) - t0 + handles.t_min_val;
                %filename_with_ext_write = [handles.filename_input,'_txy_forc_AGM_correlations',num2str(file),'.txt'];
                
                t0_in_window = t0 - handles.t_min_val;
                
                t0_bin = ceil((t0_in_window)/handles.t_bin_plot_val);
                size_TOF = max(size(TOF_axis));
                TOF_centred = zeros(size_TOF,1);
                
                min_t_centred = max(t0_bin - ceil(size_TOF/2),1);
                max_t_centred = min(t0_bin + ceil(size_TOF/2),size_TOF);
                
                offset_centroid = round(max(abs(min_t_centred-1)/2,abs(size_TOF - max_t_centred))/2);
                
                min_t_centred_for_target = max(min_t_centred - offset_centroid ,1);
                max_t_centred_for_target = min(max_t_centred - offset_centroid ,size_TOF);
                
                if max_t_centred - min_t_centred ~= max_t_centred_for_target -  min_t_centred_for_target
                    if min_t_centred < min_t_centred_for_target
                        max_t_centred = max_t_centred - (min_t_centred_for_target - min_t_centred);
                    end
                    if max_t_centred < max_t_centred_for_target
                        min_t_centred = min_t_centred + (max_t_centred_for_target - max_t_centred);
                    end
                    if min_t_centred > min_t_centred_for_target
                        min_t_centred_for_target = min_t_centred_for_target + (min_t_centred_for_target - min_t_centred);
                    end
                    if max_t_centred > max_t_centred_for_target
                        max_t_centred_for_target = max_t_centred_for_target - (max_t_centred_for_target - max_t_centred);
                    end
                    
                end
                TOF_centred(min_t_centred_for_target:max_t_centred_for_target) = TOF(min_t_centred:max_t_centred);
                TOF_centred_all = TOF_centred_all + TOF_centred;
            end
            
            if (get(handles.fit_T_check,'Value'))
                temperature_vec(file - str2num(handles.first_file_input) + 1) = Temperature;
                temperature_error_vec(file - str2num(handles.first_file_input) + 1) = error;
                condensate_percent_vec(file - str2num(handles.first_file_input) + 1) = cond_percent;
            end
        end
        
        if (get(handles.fit_remove_x0_check,'Value'))
            [x0,countrate_fit,x_width,error] = fit_1D_spatial(handles,x_axis,x_profile_per_file);

            
            %             r0
            %             countrate_fit
            %             r_width
%             figure(456981)
%             hold on
%             plot(x_axis,x_profile_per_file,'-b')
%             fit_x_profile = countrate_fit*exp(-((x_axis-x0).^2)/(2*((x_width).^2)));
%             plot(x_axis,fit_x_profile,'-k')
%             hold off
%             x0
%             pause
%             close
            
            %x0_vec(file - str2num(handles.first_file_input) + 1) = x0 + handles.x_min_val;
            
             size_x_profile = max(size(x_axis));
             bin_size_x = (handles.x_max_val- handles.x_min_val)/size_x_profile;
            
            x0_vec(file - str2num(handles.first_file_input) + 1) = x0;
            data_this_file(:,2) = data_this_file(:,2) - x0 + handles.x_min_val;
            
            x0_in_window = x0 - handles.x_min_val;
            
            x0_bin = ceil((x0_in_window)/bin_size_x);
           
            x_profile_centred = zeros(size_x_profile,1);
            
            min_x_centred = max(x0_bin - ceil(size_x_profile/2),1);
            max_x_centred = min(x0_bin + ceil(size_x_profile/2),size_x_profile);
            
            offset_centroid = round(max(abs(min_x_centred-1)/2,abs(size_x_profile - max_x_centred))/2);
            
            min_x_centred_for_target = max(min_x_centred - offset_centroid ,1);
            max_x_centred_for_target = min(max_x_centred - offset_centroid ,size_x_profile);
            
            if max_x_centred - min_x_centred ~= max_x_centred_for_target -  min_x_centred_for_target
                if min_x_centred < min_x_centred_for_target
                    max_x_centred = max_x_centred - (min_x_centred_for_target - min_x_centred);
                end
                if max_x_centred < max_x_centred_for_target
                    min_x_centred = min_x_centred + (max_x_centred_for_target - max_x_centred);
                end
                if min_x_centred > min_x_centred_for_target
                    min_x_centred_for_target = min_x_centred_for_target + (min_x_centred_for_target - min_x_centred);
                end
                if max_x_centred > max_x_centred_for_target
                    max_x_centred_for_target = max_x_centred_for_target - (max_x_centred_for_target - max_x_centred);
                end
                
            end
            x_profile_centred(min_x_centred_for_target:max_x_centred_for_target) = x_profile_per_file(min_x_centred:max_x_centred);
            x_profile_centred_all = x_profile_centred_all + x_profile_centred;
        end
        
        if (get(handles.fit_remove_y0_check,'Value'))
            [y0,countrate_fit,y_width,error] = fit_1D_spatial(handles,y_axis,y_profile_per_file);

            size_y_profile = max(size(y_axis));
            bin_size_y = (handles.y_max_val- handles.y_min_val)/size_y_profile;
            
            y0_vec(file - str2num(handles.first_file_input) + 1) = y0;
            data_this_file(:,2) = data_this_file(:,2) - y0 + handles.y_min_val;
            
            y0_in_window = y0 - handles.y_min_val;
            
            y0_bin = ceil((y0_in_window)/bin_size_y);
           
            y_profile_centred = zeros(size_y_profile,1);
            
            min_y_centred = max(y0_bin - ceil(size_y_profile/2),1);
            max_y_centred = min(y0_bin + ceil(size_y_profile/2),size_y_profile);
            
            offset_centroid = round(max(abs(min_y_centred-1)/2,abs(size_y_profile - max_y_centred))/2);
            
            min_y_centred_for_target = max(min_y_centred - offset_centroid ,1);
            max_y_centred_for_target = min(max_y_centred - offset_centroid ,size_y_profile);
            
            if max_y_centred - min_y_centred ~= max_y_centred_for_target -  min_y_centred_for_target
                if min_y_centred < min_y_centred_for_target
                    max_y_centred = max_y_centred - (min_y_centred_for_target - min_y_centred);
                end
                if max_y_centred < max_y_centred_for_target
                    min_y_centred = min_y_centred + (max_y_centred_for_target - max_y_centred);
                end
                if min_y_centred > min_y_centred_for_target
                    min_y_centred_for_target = min_y_centred_for_target + (min_y_centred_for_target - min_y_centred);
                end
                if max_y_centred > max_y_centred_for_target
                    max_y_centred_for_target = max_y_centred_for_target - (max_y_centred_for_target - max_y_centred);
                end
                
            end
            y_profile_centred(min_y_centred_for_target:max_y_centred_for_target) = y_profile_per_file(min_y_centred:max_y_centred);
            y_profile_centred_all = y_profile_centred_all + y_profile_centred;
        end
        
        
        
        handles.TOF_all = handles.TOF_all + TOF;
        handles.spatial_all = handles.spatial_all + spatial_image;
        
        % keep all windowed data in this matrix
        if (keepeachfile == 1)
            handles.txy_data_plot = cat(1, handles.txy_data_plot, data_this_file);
        end
        
        if (get(handles.g2_radio,'Value') || get(handles.g3_radio,'Value') || get(handles.fit_remove_TOF_check,'Value') || get(handles.fit_remove_x0_check,'Value')|| get(handles.fit_remove_y0_check,'Value'))         
            if (handles.is_window_locked == 0)
                filename_with_ext_write = [filename_input,'_txy_forc_AGM_correlations',num2str(file),'.txt'];

                fid = fopen([handles.filename_input,'_txy_forc_AGM_correlations_fwrite_t',num2str(file),'.bin'], 'w');
                fwrite(fid, data_this_file(:,1), 'double');
                fclose(fid);

                fid = fopen([handles.filename_input,'_txy_forc_AGM_correlations_fwrite_x',num2str(file),'.bin'], 'w');
                fwrite(fid, data_this_file(:,2), 'double');
                fclose(fid);

                fid = fopen([handles.filename_input,'_txy_forc_AGM_correlations_fwrite_y',num2str(file),'.bin'], 'w');
                fwrite(fid, data_this_file(:,3), 'double');
                fclose(fid);
            end

        end
        
        size_of_data_file = size(data_this_file);
        handles.nvf_vec(file - str2num(handles.first_file_input) + 1) = size_of_data_file(1);
    else
        handles.nvf_vec(file - str2num(handles.first_file_input) + 1) = 0;
    end
    
    TOF = zeros(max(size(TOF_axis)),1);
    x_profile_per_file = zeros(max(size(x_axis)),1);
    y_profile_per_file = zeros(max(size(y_axis)),1);
end

if (keepeachfile == 1)   
    rowstokeep1 = handles.txy_data_all(:,1) < handles.t_max_val;
    rowstokeep2 = handles.txy_data_all(:,1) > handles.t_min_val;
    rowstokeep3 = handles.txy_data_all(:,2) < handles.x_max_val;
    rowstokeep4 = handles.txy_data_all(:,2) > handles.x_min_val;
    rowstokeep5 = handles.txy_data_all(:,3) < handles.y_max_val;
    rowstokeep6 = handles.txy_data_all(:,3) > handles.y_min_val;
    rowstokeep = (rowstokeep1&rowstokeep2&rowstokeep3&rowstokeep4&rowstokeep5&rowstokeep6);
    handles.txy_data_plot = handles.txy_data_all(rowstokeep,:);
    
    if (get(handles.exclusion_check,'Value')) % use exclusion zone

        if ((size(handles.t_min_exc_val) == 1))
            if ((size(handles.t_max_exc_val) == 1))
                rowstokeep = (handles.txy_data_plot(:,1) < handles.t_min_exc_val);
                rowstokeep2 = (handles.txy_data_plot(:,1) > handles.t_max_exc_val);
                rowstokeep3 = logical(rowstokeep+rowstokeep2);
                handles.txy_data_plot = handles.txy_data_plot(rowstokeep3,:);
            end
        end
        
        rowstokeep3 = 0;
        rowstokeep4 = 0;
        
        if ((size(handles.x_min_exc_val) == 1))
            if ((size(handles.x_max_exc_val) == 1))
                rowstokeep = (handles.txy_data_plot(:,2) < handles.x_min_exc_val);
                rowstokeep2 = (handles.txy_data_plot(:,2) > handles.x_max_exc_val);
                rowstokeep3 = logical(rowstokeep+rowstokeep2);
                %handles.txy_data_plot = handles.txy_data_plot(rowstokeep3,:);
            end
        end
        
        if ((size(handles.y_min_exc_val) == 1))
            if ((size(handles.y_max_exc_val) == 1))
                rowstokeep = (handles.txy_data_plot(:,3) < handles.y_min_exc_val);
                rowstokeep2 = (handles.txy_data_plot(:,3) > handles.y_max_exc_val);
                rowstokeep4 = logical(rowstokeep+rowstokeep2);
                %handles.txy_data_plot = handles.txy_data_plot(rowstokeep3,:);
            end
        end
               
        if and(max(size(rowstokeep4))>1,max(size(rowstokeep3))>1)
            rowstokeep_xy = logical(rowstokeep3+rowstokeep4);
            handles.txy_data_plot = handles.txy_data_plot(rowstokeep_xy,:);
        elseif max(size(rowstokeep4))>1
            handles.txy_data_plot = handles.txy_data_plot(rowstokeep4,:);
        elseif max(size(rowstokeep3))>1
            handles.txy_data_plot = handles.txy_data_plot(rowstokeep3,:);
        end
    end 
    
    %FITTING
    
    size_data_file = size(handles.txy_data_plot);
    if size_data_file(1) > 0
        for hit = 1:(size(handles.txy_data_plot,1))
            spatial_image(ceil((handles.txy_data_plot(hit,2) - handles.x_min_val)/spatial_bin_size),ceil((handles.txy_data_plot(hit,3) - handles.y_min_val)/spatial_bin_size)) = spatial_image(ceil((handles.txy_data_plot(hit,2) - handles.x_min_val)/spatial_bin_size),ceil((handles.txy_data_plot(hit,3) - handles.y_min_val)/spatial_bin_size)) + 1;
            TOF(ceil((handles.txy_data_plot(hit,1) - handles.t_min_val)/handles.t_bin_plot_val)) = TOF(ceil((handles.txy_data_plot(hit,1) - handles.t_min_val)/handles.t_bin_plot_val)) + 1;
        end
    end
    
    % Fit TOF - Temp and t0
    if (get(handles.fit_remove_TOF_check,'Value') || get(handles.fit_T_check,'Value'))
        
        handles.temp_guess_tf = get(handles.temp_guess_tf_edit,'String');
        handles.tf_radius_guess = get(handles.tf_radius_guess_edit,'String');
        handles.cond_percent_guess = get(handles.cond_percent_guess_edit,'String');
        
        showfits = 0;
        [t0,Temperature,error,cond_percent,fit_out] = fit_Temp_t0(handles,TOF_axis,TOF,handles.TF_fit,showfits,str2num(handles.temp_guess_tf)*1e-9,str2num(handles.tf_radius_guess)*1e-3,str2num(handles.cond_percent_guess));
        
        if (get(handles.fit_remove_TOF_check,'Value'))
            t0_vec(file - str2num(handles.first_file_input) + 1) = t0;
            handles.txy_data_plot(:,1) = handles.txy_data_plot(:,1) - t0 + handles.t_min_val;
            filename_with_ext_write = [handles.filename_input,'_txy_forc_AGM_correlations',num2str(file),'.txt'];
            dlmwrite(filename_with_ext_write,handles.txy_data_plot,'delimiter', ';', 'newline', 'pc');
            t0_bin = ceil((t0)/handles.t_bin_plot_val);
            size_TOF = max(size(TOF_axis));
            TOF_centred = zeros(size_TOF,1);
            min_t_centred = max(t0_bin - ceil(size_TOF/2),1);
            max_t_centred = min(t0_bin + ceil(size_TOF/2),size_TOF);
            
            TOF_centred(min_t_centred + ceil(size_TOF/2 - t0_bin) + 1:max_t_centred + ceil(size_TOF/2 - t0_bin) + 1) = TOF(min_t_centred:max_t_centred);
            handles.TOF_centred_all = handles.TOF_centred_all + TOF_centred;
        end
        
        if (get(handles.fit_T_check,'Value'))
            temperature_vec(file - str2num(handles.first_file_input) + 1) = Temperature;
            temperature_error_vec(file - str2num(handles.first_file_input) + 1) = error;
            condensate_percent_vec(file - str2num(handles.first_file_input) + 1) = cond_percent;
        end
    end
end




% Window by counts
Mean_counts_per_file = mean(handles.nvf_vec(:,1));
Std_counts_per_file = std(handles.nvf_vec(:,1));
nvf_xaxis = [str2num(handles.first_file_input):str2num(handles.last_file_input)];

rowstokeep = handles.nvf_vec(:,1) < Mean_counts_per_file*(1+str2num(handles.num_window)/100);
nvf_selected = handles.nvf_vec(rowstokeep,1);
%temperature_vec =  temperature_vec(rowstokeep,1);

rowstokeep = nvf_selected(:,1) > Mean_counts_per_file*(1-str2num(handles.num_window)/100);
nvf_selected = nvf_selected(rowstokeep,1);
%temperature_vec =  temperature_vec(rowstokeep,1);

if (get(handles.fit_T_check,'Value'))
    figure(6666)
    hold on
    T_plot = plot(temperature_vec);
    %rectangle('Position',[min(nvf_xaxis),Mean_counts_per_file*(1+handles.num_window/100),max(size(nvf_xaxis)),Mean_counts_per_file*(2*handles.num_window/100)],'FaceColor',[.1,.1,.1]);
    xlabel('File number')
    ylabel('T')
    title('Temperature per file')  
%     uicontrol('Style', 'text','String', fig_1_note1,'Units','normalized','Position', [0.4 0.2 0.4 0.05],'BackgroundColor','w');
%     uicontrol('Style', 'text','String', fig_1_note2,'Units','normalized','Position', [0.4 0.15 0.4 0.05],'BackgroundColor','w');
    hold off 
end


if (get(handles.sort_by_temp_check,'Value'))
    nvf_and_temp = [nvf_selected temperature_vec];
    size(nvf_selected)
    size(temperature_vec)
    size(nvf_and_temp)
    disp('ok what?');
end

if (str2num(handles.nvf) == 1)

    fig_1_note1 = {['Average hits per file = ',num2str(Mean_counts_per_file)]};
    fig_1_note2 = {['St. dev. hits per file = ',num2str(Std_counts_per_file)]};

    figure(6)
    hold on
    NVF_plot = plot(nvf_xaxis,handles.nvf_vec);
    %rectangle('Position',[min(nvf_xaxis),Mean_counts_per_file*(1+handles.num_window/100),max(size(nvf_xaxis)),Mean_counts_per_file*(2*handles.num_window/100)],'FaceColor',[.1,.1,.1]);
    xlabel('File number')
    ylabel('Number of counts')
    title('Hits per file with windowing')  
    uicontrol('Style', 'text','String', fig_1_note1,'Units','normalized','Position', [0.2 .11 0.28 0.035],'BackgroundColor','w');
    uicontrol('Style', 'text','String', fig_1_note2,'Units','normalized','Position', [0.6 .11 0.28 0.035],'BackgroundColor','w');
    hold off
end

if (str2num(handles.tof) == 1)
    bin_factor = floor(str2num(handles.t_bin_plot)*1e-3/1e-5);
    figure(7)
    hold on
    %tof_xaxis_binned = [str2num(handles.t_min)+ 1e-5*bin_factor:1e-5*bin_factor:str2num(handles.t_max)]; 
    tof_xaxis_binned = [handles.t_min_val:handles.t_bin_plot_val:handles.t_max_val]; 
    TOF_fig = plot(tof_xaxis_binned,handles.TOF_all'./(str2num(handles.last_file_input) - str2num(handles.first_file_input)+1)*1e-6/(1e-5*bin_factor));
    %size(tof_xaxis_binned)
    %size(handles.TOF_all'./(str2num(handles.last_file_input) - str2num(handles.first_file_input)+1)*1e-6/(1e-5*bin_factor))
    xlabel('t(s)')
    ylabel(['Average flux (MHz), bin size = ',num2str(1e-5*bin_factor),' s. Saturates at 1 MHz.'])
    title('All files TOF at chosen resolution with windowing')
    if(get(handles.fit_thermal_average_check,'Value'))
        showfits = 0;
        [t0_ave_nofit,Temperature_ave_nofit,error_ave_nofit,cond_percent_ave_nofit,fit_out_ave_nofit] = fit_Temp_t0(handles,TOF_axis,handles.TOF_all,handles.TF_fit,showfits,str2num(handles.temp_guess_tf)*1e-9,str2num(handles.tf_radius_guess)*1e-3,str2num(handles.cond_percent_guess));
        temp_fit_to_centred_average_nofit = ['T = ',num2str(Temperature_ave_nofit*1e9),' nK'];
        thermal_fit_nofit=(fit_out_ave_nofit(1)*exp(-((2*(tof_xaxis_binned-t0_ave_nofit))/fit_out_ave_nofit(3)).^2)+fit_out_ave_nofit(4));
        thermal_fit_nofit = thermal_fit_nofit./max(thermal_fit_nofit)*max(handles.TOF_all'./(str2num(handles.last_file_input) - str2num(handles.first_file_input)+1));
        plot(tof_xaxis_binned,thermal_fit_nofit*1e-6/(1e-5*bin_factor),'-r')
        uicontrol('Style', 'text','String', temp_fit_to_centred_average_nofit,'Units','normalized','Position', [0.5 0.8 0.4 0.05],'BackgroundColor','w');
        set(handles.thermal_temp_disp,'String',num2str(Temperature_ave_nofit*1e9));
    end
    hold off
    
    if (get(handles.fit_remove_TOF_check,'Value'))

        figure(700)
        %tof_xaxis_binned = [str2num(handles.t_min)+ 1e-5*bin_factor:1e-5*bin_factor:str2num(handles.t_max)]; 
        tof_xaxis_binned = [str2num(handles.t_min):1e-5*bin_factor:str2num(handles.t_max)];
        hold on
        TOF_centred_fig = plot(tof_xaxis_binned,TOF_centred_all'./(str2num(handles.last_file_input) - str2num(handles.first_file_input)+1),'-k');
       % TOF_fig = plot(tof_xaxis_binned,handles.TOF_all'./(str2num(handles.last_file_input) - str2num(handles.first_file_input)+1));
        
        
        if(get(handles.fit_thermal_average_check,'Value'))
            [t0_ave,Temperature_ave,error_ave,cond_percent_ave,fit_out_ave] = fit_Temp_t0(handles,TOF_axis,TOF_centred_all,handles.TF_fit,showfits,str2num(handles.temp_guess_tf)*1e-9,str2num(handles.tf_radius_guess)*1e-3,str2num(handles.cond_percent_guess));
            set(handles.thermal_temp_fit_disp,'String',num2str(Temperature_ave*1e9));
            set(handles.cond_perc_disp,'String',num2str(cond_percent_ave));
            thermal_fit=(fit_out_ave(1)*exp(-((2*(tof_xaxis_binned-t0_ave))/fit_out_ave(3)).^2)+fit_out_ave(4));
            thermal_fit = thermal_fit./max(thermal_fit)*max(TOF_centred_all'./(str2num(handles.last_file_input) - str2num(handles.first_file_input)+1));
            plot(tof_xaxis_binned,thermal_fit,'-r')
            temp_fit_to_centred_average = ['T = ',num2str(Temperature_ave*1e9),' nK'];
            uicontrol('Style', 'text','String', temp_fit_to_centred_average,'Units','normalized','Position', [0.5 0.8 0.4 0.05],'BackgroundColor','w');
        end
        
        
        xlabel('t(s)')
        ylabel(['Average flux, bin size = ',num2str(1e-5*bin_factor),' s'])
        title('All files TOF at chosen resolution with windowing and centred') 
        hold off
    end
end

if (str2num(handles.spatial) == 1)
    figure(8)
    %axis_rebin_x = [handles.x_min_val:80/606*bins_per_pixel_number:handles.x_max_val];
    %axis_rebin_y = [handles.y_min_val:80/606*bins_per_pixel_number:handles.y_max_val];
    axis_rebin_x = [handles.x_min_val:(handles.x_max_val - handles.x_min_val)/(size(handles.spatial_all,1)-1):handles.x_max_val]';
    axis_rebin_y = [handles.y_min_val:(handles.y_max_val - handles.y_min_val)/(size(handles.spatial_all,2)-1):handles.y_max_val]';
    
%     size(axis_rebin_x)
%     size(axis_rebin_y)
%     size(spatial_all')
    
    spatial_plot_fig_rebinned = pcolor(axis_rebin_x,axis_rebin_y,handles.spatial_all');
    %spatial_plot_fig_rebinned = pcolor(spatial_all');
    set(spatial_plot_fig_rebinned,'EdgeColor','none')
    xlabel({'x position (mm)'});
    ylabel({'y position (mm)'});
    title(['All files spatial image at ',num2str(handles.bin_size_disp),' mm resolution with windowing'])
    colormap gray
    %colorbar
    
    figure(880)
    hold on
    x_profile = sum(handles.spatial_all',1);
    x_profile_plot = plot(axis_rebin_x,x_profile);
    xlabel({'x position (mm)'});
    ylabel({'Density'});
    title(['X profile, all files at ',num2str(handles.bin_size_disp),' mm resolution with windowing'])
    if (get(handles.fit_remove_x0_check,'Value'))
        [x0,countrate_fit,x_width,error] = fit_1D_spatial(handles,axis_rebin_x,x_profile);
        fit_x_profile = countrate_fit*exp(-((x_axis-x0).^2)/(2*((x_width).^2)));
        plot(x_axis,fit_x_profile,'-k')
        x_fit_to_centred_average = ['x width = ',num2str(x_width*1e3),' mm'];
        uicontrol('Style', 'text','String', x_fit_to_centred_average,'Units','normalized','Position', [0.7 0.8 0.2 0.05],'BackgroundColor','w');
    end
    hold off
    
%     x0
%     sum(x_axis/1000.*fit_x_profile/countrate_fit)
    
    if (get(handles.fit_remove_x0_check,'Value'))
        figure(881)
        hold on
        x_profile_cented_plot = plot(x_axis,x_profile_centred_all);
        xlabel({'x position (mm)'});
        ylabel({'Density'});
        title(['X profile, all files at ',num2str(handles.bin_size_disp),' mm resolution with windowing and centred'])
        [x0,countrate_fit,x_width,error] = fit_1D_spatial(handles,x_axis,x_profile_centred_all);
        fit_x_profile = countrate_fit*exp(-((x_axis-x0).^2)/(2*((x_width).^2)));
        plot(x_axis,fit_x_profile,'-k')
        x_fit_to_centred_average = ['x width = ',num2str(x_width*1e3),' mm'];
        uicontrol('Style', 'text','String', x_fit_to_centred_average,'Units','normalized','Position', [0.7 0.8 0.2 0.05],'BackgroundColor','w');
        set(handles.spatial_fit_x_disp,'String',num2str(x_width*1e3));
        hold off
        
%         x0
%         sum(x_axis.*fit_x_profile)/countrate_fit/1000
    end

%     figure(890)
%     y_profile = sum(handles.spatial_all',2);
%     y_profile_plot = plot(axis_rebin_y,y_profile);
%     xlabel({'y position (mm)'});
%     ylabel({'Density'});
%     title(['Y profile, all files at ',num2str(handles.bin_size_disp),' mm resolution with windowing'])
    
   figure(890)
    hold on
    y_profile = sum(handles.spatial_all',2);
    y_profile_plot = plot(axis_rebin_y,y_profile);
    xlabel({'y position (mm)'});
    ylabel({'Density'});
    title(['Y profile, all files at ',num2str(handles.bin_size_disp),' mm resolution with windowing'])
    if (get(handles.fit_remove_y0_check,'Value'))
        [y0,countrate_fit,y_width,error] = fit_1D_spatial(handles,axis_rebin_y,y_profile');
        fit_y_profile = countrate_fit*exp(-((y_axis-y0).^2)/(2*((y_width).^2)));
        plot(y_axis,fit_y_profile,'-k')
        y_fit_to_centred_average = ['y width = ',num2str(y_width*1e3),' mm'];
        uicontrol('Style', 'text','String', y_fit_to_centred_average,'Units','normalized','Position', [0.7 0.8 0.2 0.05],'BackgroundColor','w');
    end
    hold off
    
%     x0
%     sum(x_axis/1000.*fit_x_profile/countrate_fit)
    
    if (get(handles.fit_remove_y0_check,'Value'))
        figure(891)
        hold on
        y_profile_cented_plot = plot(y_axis,y_profile_centred_all);
        xlabel({'y position (mm)'});
        ylabel({'Density'});
        title(['Y profile, all files at ',num2str(handles.bin_size_disp),' mm resolution with windowing and centred'])
        [y0,countrate_fit,y_width,error] = fit_1D_spatial(handles,y_axis,y_profile_centred_all);
        fit_y_profile = countrate_fit*exp(-((y_axis-y0).^2)/(2*((y_width).^2)));
        plot(y_axis,fit_y_profile,'-k')
        y_fit_to_centred_average = ['y width = ',num2str(y_width*1e3),' mm'];
        uicontrol('Style', 'text','String', y_fit_to_centred_average,'Units','normalized','Position', [0.7 0.8 0.2 0.05],'BackgroundColor','w');
        set(handles.spatial_fit_y_disp,'String',num2str(y_width*1e3));
        hold off
        
%         x0
%         sum(x_axis.*fit_x_profile)/countrate_fit/1000
    end
    
    
    
    
    
    %%%%FIX!!!
    if (get(handles.fit_spatial_average,'Value'))
           
        binned_output = handles.spatial_all';
                
        %guess = [1, radius_fit_guess_num, radius_fit_guess_num, 0, y_fit_guess_num, x_fit_guess_num];
        guess = [1, 1e-3, 1e-3, 0, 0, 0];
        
        guess(2) = 5*20/bins_per_pixel_number;
        guess(3) = 5*20/bins_per_pixel_number;
        guess(5) = ceil(max(size(binned_output))/2);
        guess(6) = ceil(max(size(binned_output))/2);
        brightest_pixel = max(max(binned_output));
        binned_output = binned_output./brightest_pixel;
        
        [Param_out,~] = fminsearch(@(Param)fit_2d_gaussian(Param,binned_output),guess);
        fitted_profile = zeros(max(size(binned_output)),max(size(binned_output)));
        for xpos = 1:max(size(binned_output))
            for ypos = 1:max(size(binned_output))
                fitted_profile(xpos,ypos) = Param_out(1)*exp(-0.5*((xpos-Param_out(5))^2/(Param_out(2))^2))*exp(-0.5*((ypos-Param_out(6))^2/(Param_out(3))^2)) + Param_out(4);
            end
        end
        
        figure(1e6)
        spatial_plot_fig_rebinned = pcolor(axis_rebin,axis_rebin,fitted_profile');
        set(spatial_plot_fig_rebinned,'EdgeColor','none')
        xlabel({'x position (mm)'});
        ylabel({'y position (mm)'});
        title(['Fitted Thermal Profile to all files spatial image at ',num2str(handles.bin_size_disp),' mm resolution with windowing'])
        colormap gray
        colorbar
        
        pixel_size = 13.158e-6 * bins_per_pixel_number;
        set(handles.spatial_fit_x_disp,'String',num2str(Param_out(2)*1e3*pixel_size));
        set(handles.spatial_fit_y_disp,'String',num2str(Param_out(3)*1e3*pixel_size));
    end
    
end


if (get(handles.save_check,'Value')) % use exclusion zone
    fit_TOF_string = '';
    if (handles.fit_remove_TOF == '1')
        fit_TOF_string = '_fit_TOF';
    end
    
    fit_x0_string = '';
    if (handles.fit_remove_x0 == '1')
        fit_x0_string = '_fit_x0';
    end
    
    fit_y0_string = '';
    if (handles.fit_remove_y0 == '1')
        fit_y0_string = '_fit_y0';
    end
    
    temp_window_string = '';
    if (get(handles.temp_window_edit,'String') ~= '0')
        temp_window_string = ['temp_window',get(handles.temp_window_edit,'String')];
    end
    
    num_window_string = '';
    if (get(handles.num_window_edit,'String') ~= '0')
        num_window_string = ['num_window',get(handles.num_window_edit,'String')];
    end
    
    exclusion_string = '';
    if (handles.exclusion_active == '1')
        
        t_min_exc_val = str2num(get(handles.t_min_exc_edit,'String'));
        t_max_exc_val = str2num(get(handles.t_max_exc_edit,'String'));
        x_min_exc_val = str2num(get(handles.x_min_exc_edit,'String'))*1e-3;
        x_max_exc_val = str2num(get(handles.x_max_exc_edit,'String'))*1e-3;
        y_min_exc_val = str2num(get(handles.y_min_exc_edit,'String'))*1e-3;
        y_max_exc_val = str2num(get(handles.y_max_exc_edit,'String'))*1e-3;
        
        if ((size(t_min_exc_val) == 1))
            if ((size(t_max_exc_val) == 1))
                exclusion_string = [exclusion_string,'t_ex',num2str(t_min_exc_val),'to',num2str(t_max_exc_val)];
            end
        end
        
        if ((size(x_min_exc_val) == 1))
            if ((size(x_max_exc_val) == 1))
                exclusion_string = [exclusion_string,'x_ex',num2str(x_min_exc_val),'to',num2str(x_max_exc_val)];
            end
        end
        
        if ((size(y_min_exc_val) == 1))
            if ((size(y_max_exc_val) == 1))
                exclusion_string = [exclusion_string,'y_ex',num2str(y_min_exc_val),'to',num2str(y_max_exc_val)];
            end
        end
        
    end
    
    custom_file_list_message = '';    
    if (get(handles.custom_file_list_check,'Value'))
        handles.custom_file_list_string = get(handles.custom_file_list_string_edit,'String');
        custom_file_list_message = ['_',handles.custom_file_list_string];
    end

    new_dir_name = ['outputs\file',handles.first_file_input,'to',handles.last_file_input,'t',handles.t_min,'to',handles.t_max,'x',handles.x_min,'to',handles.x_max,'y',handles.y_min,'to',handles.y_max,'axis_',handles.axis_s,'_tbin',handles.t_bin,'xbin',handles.x_bin,'ybin',handles.y_bin,fit_TOF_string,fit_x0_string,fit_y0_string,temp_window_string,num_window_string,exclusion_string,custom_file_list_message];
    if (~exist([handles.filename_input,'\',new_dir_name],'dir'))
        mkdir(handles.filename_input,new_dir_name);
    end
    copyfile([handles.filename_input,'_channels.txt'],[handles.filename_input,'\',new_dir_name,'\channels.txt']);
    copyfile([handles.filename_input,'_nvf.txt'],[handles.filename_input,'\',new_dir_name,'\nvf.txt']);
    copyfile([handles.filename_input,'_TOF.txt'],[handles.filename_input,'\',new_dir_name,'\TOF.txt']);
    copyfile([handles.filename_input,'_spatial.txt'],[handles.filename_input,'\',new_dir_name,'\spatial.txt']);
    copyfile('config123.txt',[handles.filename_input,'\',new_dir_name,'\config123.txt']);
    if(exist('NVF_plot'))
        saveas(NVF_plot,[handles.filename_input,'\',new_dir_name,'\NVF_plot.fig']);
    end
    if(exist('TOF_fig'))
        saveas(TOF_fig,[handles.filename_input,'\',new_dir_name,'\TOF_rebinned.fig']);
    end
    if(exist('TOF_centred_fig'))
        saveas(TOF_centred_fig,[handles.filename_input,'\',new_dir_name,'\TOF_centred_rebinned.fig']);
    end
    if(exist('spatial_plot_fig_rebinned'))
        saveas(spatial_plot_fig_rebinned,[handles.filename_input,'\',new_dir_name,'\spatial_rebinned.fig']);
    end
    if(exist('x_profile_cented_plot'))
        saveas(x_profile_cented_plot,[handles.filename_input,'\',new_dir_name,'\x_profile_cented_plot.fig']);
    end
    if(exist('y_profile_cented_plot'))
        saveas(y_profile_cented_plot,[handles.filename_input,'\',new_dir_name,'\y_profile_cented_plot.fig']);
    end
    
    
    dlmwrite([handles.filename_input,'\',new_dir_name,'\nvf_windowed.txt'],(handles.nvf_vec)','precision', '%9.0f','delimiter', ';', 'newline', 'pc');
    dlmwrite([handles.filename_input,'_nvf_windowed.txt'],(handles.nvf_vec),'precision', '%9.0f','delimiter', ';', 'newline', 'pc');
    if (get(handles.fit_remove_TOF_check,'Value'))
        dlmwrite([handles.filename_input,'\',new_dir_name,'\t0.txt'],(t0_vec)','delimiter', ';', 'newline', 'pc');
        dlmwrite([handles.filename_input,'_t0.txt'],t0_vec','delimiter', ';', 'newline', 'pc');
    end
    if (get(handles.fit_remove_x0_check,'Value'))
        dlmwrite([handles.filename_input,'\',new_dir_name,'\x0.txt'],(x0_vec)','delimiter', ';', 'newline', 'pc');
        dlmwrite([handles.filename_input,'_x0.txt'],x0_vec','delimiter', ';', 'newline', 'pc');
    end
    if (get(handles.fit_remove_y0_check,'Value'))
        dlmwrite([handles.filename_input,'\',new_dir_name,'\y0.txt'],(y0_vec)','delimiter', ';', 'newline', 'pc');
        dlmwrite([handles.filename_input,'_y0.txt'],y0_vec','delimiter', ';', 'newline', 'pc');
    end   
    if (get(handles.fit_T_check,'Value'))
        dlmwrite([handles.filename_input,'\',new_dir_name,'\temperature.txt'],temperature_vec','delimiter', ';', 'newline', 'pc');
        dlmwrite([handles.filename_input,'\',new_dir_name,'\temperature_error.txt'],temperature_error_vec','delimiter', ';', 'newline', 'pc');
        dlmwrite([handles.filename_input,'\',new_dir_name,'\condensate_percent.txt'],condensate_percent_vec','delimiter', ';', 'newline', 'pc');
        dlmwrite([handles.filename_input,'_temperature.txt'],temperature_vec','delimiter', ';', 'newline', 'pc');
        dlmwrite([handles.filename_input,'_temperature_error.txt'],temperature_error_vec','delimiter', ';', 'newline', 'pc');
        dlmwrite([handles.filename_input,'_condensate_percent.txt'],condensate_percent_vec','delimiter', ';', 'newline', 'pc');
    end
       
    if (get(handles.g2_radio,'Value') || get(handles.g3_radio,'Value'))
        
        file_list = str2num(handles.first_file_input):str2num(handles.last_file_input);

        if (str2num(get(handles.temp_window_edit,'String')) ~= 0)
            %disp('got here 123')
            mean_temp = mean(temperature_vec);
            rows_to_throw = abs(temperature_vec - mean_temp) > str2num(handles.temp_window)*mean_temp/100;
            file_list(rows_to_throw) = [];  
        end
        if (str2num(get(handles.num_window_edit,'String')) ~= 0)
            %disp('got here 1234')
            mean_num = mean(handles.nvf_vec);
            rows_to_throw = abs(handles.nvf_vec - mean_num) > str2num(handles.num_window)*mean_num/100;
            file_list(rows_to_throw) = [];  
        end    
        
        if (get(handles.custom_file_list_check,'Value'))
            disp('Will be using custom file list, not any computed file lists. Uncheck custom file list to window by number or temperature from front panel calculations.')
             copyfile([handles.filename_input,'_file_list.txt'],[handles.filename_input,'\',new_dir_name,'\file_list.txt']);
        else    
            dlmwrite([handles.filename_input,'\',new_dir_name,'\file_list.txt'],file_list','precision', '%6.0f','delimiter', ';', 'newline', 'pc');
            dlmwrite([handles.filename_input,'_file_list.txt'],file_list','precision', '%6.0f','delimiter', ';', 'newline', 'pc');        
        end

        handles.ready_for_gn = '1';
        write_config(hObject,handles);

        fprintf('\n')
        disp('             ******************************')
        disp('             *  beginning g2 calculation  *')
        disp('             ******************************')
        if handles.x64_status == '0'
            system('DLD_convert_and_correlations_for_MATLAB.exe');
        elseif handles.x64_status == '1'
            system('DLD_convert_and_correlations_for_MATLAB_x64.exe');
        end
       
        handles.ready_for_gn = '0';
        write_config(hObject,handles);
        
        if(exist([handles.filename_input,'_g2_0.txt']))
            copyfile([handles.filename_input,'_g2_0.txt'],[handles.filename_input,'\',new_dir_name,'\g2_0.txt']);
            copyfile([handles.filename_input,'_g2_1.txt'],[handles.filename_input,'\',new_dir_name,'\g2_1.txt']);
            copyfile([handles.filename_input,'_g2_2.txt'],[handles.filename_input,'\',new_dir_name,'\g2_2.txt']);
            copyfile([handles.filename_input,'_g2_3.txt'],[handles.filename_input,'\',new_dir_name,'\g2_3.txt']);
            copyfile([handles.filename_input,'_g2_m1.txt'],[handles.filename_input,'\',new_dir_name,'\g2_m1.txt']);
            copyfile([handles.filename_input,'_g2_m2.txt'],[handles.filename_input,'\',new_dir_name,'\g2_m2.txt']);
            copyfile([handles.filename_input,'_g2_m3.txt'],[handles.filename_input,'\',new_dir_name,'\g2_m3.txt']);
            copyfile([handles.filename_input,'_g2_den_0.txt'],[handles.filename_input,'\',new_dir_name,'\g2_den_0.txt']);
            copyfile([handles.filename_input,'_g2_den_1.txt'],[handles.filename_input,'\',new_dir_name,'\g2_den_1.txt']);
            copyfile([handles.filename_input,'_g2_den_2.txt'],[handles.filename_input,'\',new_dir_name,'\g2_den_2.txt']);
            copyfile([handles.filename_input,'_g2_den_3.txt'],[handles.filename_input,'\',new_dir_name,'\g2_den_3.txt']);
            copyfile([handles.filename_input,'_g2_den_m1.txt'],[handles.filename_input,'\',new_dir_name,'\g2_den_m1.txt']);
            copyfile([handles.filename_input,'_g2_den_m2.txt'],[handles.filename_input,'\',new_dir_name,'\g2_den_m2.txt']);
            copyfile([handles.filename_input,'_g2_den_m3.txt'],[handles.filename_input,'\',new_dir_name,'\g2_den_m3.txt']);
        end
        if((exist([handles.filename_input,'_g3_0.txt'])) && get(handles.g3_radio,'Value'))
            copyfile([handles.filename_input,'_g3_0.txt'],[handles.filename_input,'\',new_dir_name,'\g3_0.txt']);
            copyfile([handles.filename_input,'_g3_1.txt'],[handles.filename_input,'\',new_dir_name,'\g3_1.txt']);
            copyfile([handles.filename_input,'_g3_2.txt'],[handles.filename_input,'\',new_dir_name,'\g3_2.txt']);
            copyfile([handles.filename_input,'_g3_3.txt'],[handles.filename_input,'\',new_dir_name,'\g3_3.txt']);
            copyfile([handles.filename_input,'_g3_m1.txt'],[handles.filename_input,'\',new_dir_name,'\g3_m1.txt']);
            copyfile([handles.filename_input,'_g3_m2.txt'],[handles.filename_input,'\',new_dir_name,'\g3_m2.txt']);
            copyfile([handles.filename_input,'_g3_m3.txt'],[handles.filename_input,'\',new_dir_name,'\g3_m3.txt']);
            copyfile([handles.filename_input,'_g3_den_0.txt'],[handles.filename_input,'\',new_dir_name,'\g3_den_0.txt']);
            copyfile([handles.filename_input,'_g3_den_1.txt'],[handles.filename_input,'\',new_dir_name,'\g3_den_1.txt']);
            copyfile([handles.filename_input,'_g3_den_2.txt'],[handles.filename_input,'\',new_dir_name,'\g3_den_2.txt']);
            copyfile([handles.filename_input,'_g3_den_3.txt'],[handles.filename_input,'\',new_dir_name,'\g3_den_3.txt']);
            copyfile([handles.filename_input,'_g3_den_m1.txt'],[handles.filename_input,'\',new_dir_name,'\g3_den_m1.txt']);
            copyfile([handles.filename_input,'_g3_den_m2.txt'],[handles.filename_input,'\',new_dir_name,'\g3_den_m2.txt']);
            copyfile([handles.filename_input,'_g3_den_m3.txt'],[handles.filename_input,'\',new_dir_name,'\g3_den_m3.txt']);
        end
        
        g2_0 = dlmread([handles.filename_input,'_g2_0.txt']);
        g2_1 = dlmread([handles.filename_input,'_g2_1.txt']);
        g2_2 = dlmread([handles.filename_input,'_g2_2.txt']);
        g2_3 = dlmread([handles.filename_input,'_g2_3.txt']);
        g2_m1 = dlmread([handles.filename_input,'_g2_m1.txt']);
        g2_m2 = dlmread([handles.filename_input,'_g2_m2.txt']);
        g2_m3 = dlmread([handles.filename_input,'_g2_m3.txt']);
        g2_den_0 = dlmread([handles.filename_input,'_g2_den_0.txt']);
        g2_den_1 = dlmread([handles.filename_input,'_g2_den_1.txt']);
        g2_den_2 = dlmread([handles.filename_input,'_g2_den_2.txt']);
        g2_den_3 = dlmread([handles.filename_input,'_g2_den_3.txt']);
        g2_den_m1 = dlmread([handles.filename_input,'_g2_den_m1.txt']);
        g2_den_m2 = dlmread([handles.filename_input,'_g2_den_m2.txt']);
        g2_den_m3 = dlmread([handles.filename_input,'_g2_den_m3.txt']);
        
        if (get(handles.dt_radio,'Value') == 1)
            axis = 1;
        elseif (get(handles.dx_radio,'Value') == 1)
            axis = 2;
        elseif (get(handles.dy_radio,'Value') == 1)
            axis = 3;
        end

        handles.t_bin = str2num(get(handles.t_bin_edit,'String'));
        handles.x_bin = str2num(get(handles.x_bin_edit,'String'));
        handles.y_bin = str2num(get(handles.y_bin_edit,'String'));
        handles.num_bins = str2num(get(handles.number_bins_edit,'String'));
        
        switch(axis)
            case 1
                hist_vec = handles.t_bin*[1:handles.num_bins];
                x_title = 'dt (ms)';
                y_title = 'g(2)(t)';
                y_title_g3 = 'g(3)(t)';
                title_1 = ' mm in x, and ';
                title_2 = ' mm in y';
            case 2
                hist_vec = handles.x_bin*[1:handles.num_bins];
                x_title = 'dx (mm)';
                y_title = 'g(2)(x)';
                y_title_g3 = 'g(3)(x)';
                title_1 = ' ms in t, and ';
                title_2 = ' mm in y';
            case 3
                hist_vec = handles.y_bin*[1:handles.num_bins];
                x_title = 'dy (mm)';
                y_title = 'g(2)(y)';
                y_title_g3 = 'g(3)(y)';
                title_1 = ' ms in t, and ';
                title_2 = ' mm in x';
        end
            
        handles.first_file = str2num(get(handles.first_file_edit,'String'));
        handles.last_file = str2num(get(handles.last_file_edit,'String'));
        num_files = handles.last_file - handles.first_file + 1;
        
        num_files_in_denom = str2num(get(handles.num_files_denom_edit,'String'));
        
        normalisation_factor = 1 + var(handles.nvf_vec)/(mean(handles.nvf_vec)^2)
        
        g2_0_norm = g2_0./(g2_den_0 - g2_0)*(num_files-1)/(num_files/num_files_in_denom)/normalisation_factor;
        g2_1_norm = g2_1./(g2_den_1 - g2_1)*(num_files-1)/(num_files/num_files_in_denom)/normalisation_factor;
        g2_2_norm = g2_2./(g2_den_2 - g2_2)*(num_files-1)/(num_files/num_files_in_denom)/normalisation_factor;
        g2_3_norm = g2_3./(g2_den_3 - g2_3)*(num_files-1)/(num_files/num_files_in_denom)/normalisation_factor;
        g2_m1_norm = g2_m1./(g2_den_m1 - g2_m1)*(num_files-1)/(num_files/num_files_in_denom)/normalisation_factor;
        g2_m2_norm = g2_m2./(g2_den_m2 - g2_m2)*(num_files-1)/(num_files/num_files_in_denom)/normalisation_factor;
        g2_m3_norm = g2_m3./(g2_den_m3 - g2_m3)*(num_files-1)/(num_files/num_files_in_denom)/normalisation_factor;
        
       size_g2 = max(size(g2_0_norm));
       
        g2_0_norm = g2_0_norm ./ mean(g2_0_norm(ceil(size_g2/2):size_g2));
        g2_1_norm = g2_1_norm ./ mean(g2_1_norm(ceil(size_g2/2):size_g2));
        g2_2_norm = g2_2_norm ./ mean(g2_2_norm(ceil(size_g2/2):size_g2));
        g2_3_norm = g2_3_norm ./ mean(g2_3_norm(ceil(size_g2/2):size_g2));
        g2_m1_norm = g2_m1_norm ./ mean(g2_m1_norm(ceil(size_g2/2):size_g2));
        g2_m2_norm = g2_m2_norm ./ mean(g2_m2_norm(ceil(size_g2/2):size_g2));
        g2_m3_norm = g2_m3_norm ./ mean(g2_m3_norm(ceil(size_g2/2):size_g2));
        
        
        % FIT RESULTANT G2
        
        % time is in units of ms
        
        A_g = max(g2_0_norm)-1;
        his_vec_centre = ceil(max(size(hist_vec))/2);
        sigma2_g = hist_vec(his_vec_centre);   %units are ms
        guess = [A_g,sigma2_g];
        [Param_out,~]=fminsearch(@(Param)fit_g2_properly(Param,g2_0_norm',hist_vec),guess);
        g2_fit_0 = 1+Param_out(1)*exp(-((hist_vec).^2)/(2*(Param_out(2))^2));
        
        A_g = max(g2_1_norm)-1;
        guess = [A_g,sigma2_g];
        [Param_out1,~]=fminsearch(@(Param)fit_g2_properly(Param,g2_1_norm',hist_vec),guess);
        g2_fit_1 = 1+Param_out1(1)*exp(-((hist_vec).^2)/(2*(Param_out1(2))^2));
        
        A_g = max(g2_2_norm)-1;
        guess = [A_g,sigma2_g];
        [Param_out2,~]=fminsearch(@(Param)fit_g2_properly(Param,g2_2_norm',hist_vec),guess);
        g2_fit_2 = 1+Param_out2(1)*exp(-((hist_vec).^2)/(2*(Param_out2(2))^2));
        
        A_g = max(g2_3_norm)-1;
        guess = [A_g,sigma2_g];
        [Param_out3,~]=fminsearch(@(Param)fit_g2_properly(Param,g2_3_norm',hist_vec),guess);
        g2_fit_3 = 1+Param_out3(1)*exp(-((hist_vec).^2)/(2*(Param_out3(2))^2));
        
        A_g = max(g2_m1_norm)-1;
        guess = [A_g,sigma2_g];
        [Param_outm1,~]=fminsearch(@(Param)fit_g2_properly(Param,g2_m1_norm',hist_vec),guess);
        g2_fit_m1 = 1+Param_outm1(1)*exp(-((hist_vec).^2)/(2*(Param_outm1(2))^2));
        
        A_g = max(g2_m2_norm)-1;
        guess = [A_g,sigma2_g];
        [Param_outm2,~]=fminsearch(@(Param)fit_g2_properly(Param,g2_m2_norm',hist_vec),guess);
        g2_fit_m2 = 1+Param_outm2(1)*exp(-((hist_vec).^2)/(2*(Param_outm2(2))^2));
        
        A_g = max(g2_m3_norm)-1;
        guess = [A_g,sigma2_g];
        [Param_outm3,~]=fminsearch(@(Param)fit_g2_properly(Param,g2_m3_norm',hist_vec),guess);
        g2_fit_m3 = 1+Param_outm3(1)*exp(-((hist_vec).^2)/(2*(Param_outm3(2))^2));
        
        
        
        
        % and plot g2
        
        figure(99901)
        hold on
        g2_0_norm_plot = plot(hist_vec,g2_0_norm,'b-');
        plot(hist_vec,g2_fit_0,'-r');
        xlabel({x_title});
        ylabel({y_title});
        title(['Bins are ',num2str(handles.x_bin),title_1,num2str(handles.y_bin),title_2])
        g2_0_message = ['g(2)(0) = ',num2str(Param_out(1)+1),10,'corr length = ',num2str(Param_out(2))];
        uicontrol('Style', 'text','String', g2_0_message,'Units','normalized','Position', [0.7 .8 .2 0.1],'BackgroundColor','w');
        %[left bottom width height]
        hold off
        saveas(g2_0_norm_plot,[handles.filename_input,'\',new_dir_name,'\g2_0_norm_plot.fig']);
        
        figure(99902)
        hold on
        g2_1_norm_plot = plot(hist_vec,g2_1_norm);
        plot(hist_vec,g2_fit_1,'-r');
        xlabel({x_title});
        ylabel({y_title});
        title(['Bins are ',num2str(handles.x_bin*1.2),title_1,num2str(handles.y_bin*1.2),title_2])
        g2_1_message = ['g(2)(0) = ',num2str(Param_out1(1)+1),10,'corr length = ',num2str(Param_out1(2))];
        uicontrol('Style', 'text','String', g2_1_message,'Units','normalized','Position', [0.7 .8 .2 0.1],'BackgroundColor','w');
        hold off
        saveas(g2_1_norm_plot,[handles.filename_input,'\',new_dir_name,'\g2_1_norm_plot.fig']);
        
        figure(99903)
        hold on
        g2_2_norm_plot = plot(hist_vec,g2_2_norm);
        plot(hist_vec,g2_fit_2,'-r');
        xlabel({x_title});
        ylabel({y_title});
        title(['Bins are ',num2str(handles.x_bin*1.4),title_1,num2str(handles.y_bin*1.4),title_2])
        g2_2_message = ['g(2)(0) = ',num2str(Param_out2(1)+1),10,'corr length = ',num2str(Param_out2(2))];
        uicontrol('Style', 'text','String', g2_2_message,'Units','normalized','Position', [0.7 .8 .2 0.1],'BackgroundColor','w');
        hold off
        saveas(g2_2_norm_plot,[handles.filename_input,'\',new_dir_name,'\g2_2_norm_plot.fig']);
        
        figure(99904)
        hold on
        g2_3_norm_plot = plot(hist_vec,g2_3_norm);
        plot(hist_vec,g2_fit_3,'-r');
        xlabel({x_title});
        ylabel({y_title});
        title(['Bins are ',num2str(handles.x_bin*1.6),title_1,num2str(handles.y_bin*1.6),title_2])
        g2_3_message = ['g(2)(0) = ',num2str(Param_out3(1)+1),10,'corr length = ',num2str(Param_out3(2))];
        uicontrol('Style', 'text','String', g2_3_message,'Units','normalized','Position', [0.7 .8 .2 0.1],'BackgroundColor','w');
        hold off
        saveas(g2_3_norm_plot,[handles.filename_input,'\',new_dir_name,'\g2_3_norm_plot.fig']);
        
        figure(99905)
        hold on
        g2_m1_norm_plot = plot(hist_vec,g2_m1_norm);
        plot(hist_vec,g2_fit_m1,'-r');
        xlabel({x_title});
        ylabel({y_title});
        title(['Bins are ',num2str(handles.x_bin*0.8),title_1,num2str(handles.y_bin*0.8),title_2])
        g2_m1_message = ['g(2)(0) = ',num2str(Param_outm1(1)+1),10,'corr length = ',num2str(Param_outm1(2))];
        uicontrol('Style', 'text','String', g2_m1_message,'Units','normalized','Position', [0.7 .8 .2 0.1],'BackgroundColor','w');
        hold off
        saveas(g2_m1_norm_plot,[handles.filename_input,'\',new_dir_name,'\g2_m1_norm_plot.fig']);
        
        figure(99906)
        hold on
        g2_m2_norm_plot = plot(hist_vec,g2_m2_norm);
        plot(hist_vec,g2_fit_m2,'-r');
        xlabel({x_title});
        ylabel({y_title});
        title(['Bins are ',num2str(handles.x_bin*0.6),title_1,num2str(handles.y_bin*0.6),title_2])
        g2_m2_message = ['g(2)(0) = ',num2str(Param_outm2(1)+1),10,'corr length = ',num2str(Param_outm2(2))];
        uicontrol('Style', 'text','String', g2_m2_message,'Units','normalized','Position', [0.7 .8 .2 0.1],'BackgroundColor','w');
        hold off
        saveas(g2_m2_norm_plot,[handles.filename_input,'\',new_dir_name,'\g2_m2_norm_plot.fig']);
        
        figure(99907)
        hold on
        g2_m3_norm_plot = plot(hist_vec,g2_m3_norm);
        plot(hist_vec,g2_fit_m3,'-r');
        xlabel({x_title});
        ylabel({y_title});
        title(['Bins are ',num2str(handles.x_bin*0.4),title_1,num2str(handles.y_bin*0.4),title_2])
        g2_m3_message = ['g(2)(0) = ',num2str(Param_outm3(1)+1),10,'corr length = ',num2str(Param_outm3(2))];
        uicontrol('Style', 'text','String', g2_m3_message,'Units','normalized','Position', [0.7 .8 .2 0.1],'BackgroundColor','w');
        hold off
        saveas(g2_m3_norm_plot,[handles.filename_input,'\',new_dir_name,'\g2_m3_norm_plot.fig']);
        
     
        if(exist([handles.filename_input,'_g3_0.txt'])  && get(handles.g3_radio,'Value') )
            g3_0 = dlmread([handles.filename_input,'_g3_0.txt']);
            g3_1 = dlmread([handles.filename_input,'_g3_1.txt']);
            g3_2 = dlmread([handles.filename_input,'_g3_2.txt']);
            g3_3 = dlmread([handles.filename_input,'_g3_3.txt']);
            g3_m1 = dlmread([handles.filename_input,'_g3_m1.txt']);
            g3_m2 = dlmread([handles.filename_input,'_g3_m2.txt']);
            g3_m3 = dlmread([handles.filename_input,'_g3_m3.txt']);
            g3_den_0 = dlmread([handles.filename_input,'_g3_den_0.txt']);
            g3_den_1 = dlmread([handles.filename_input,'_g3_den_1.txt']);
            g3_den_2 = dlmread([handles.filename_input,'_g3_den_2.txt']);
            g3_den_3 = dlmread([handles.filename_input,'_g3_den_3.txt']);
            g3_den_m1 = dlmread([handles.filename_input,'_g3_den_m1.txt']);
            g3_den_m2 = dlmread([handles.filename_input,'_g3_den_m2.txt']);
            g3_den_m3 = dlmread([handles.filename_input,'_g3_den_m3.txt']);
            
            g3_0_norm = g3_0./(g3_den_0 - g3_0)*((num_files)^2 - 1)/(((num_files/num_files_in_denom))^2);
            g3_1_norm = g3_1./(g3_den_1 - g3_1)*((num_files)^2 - 1)/(((num_files/num_files_in_denom))^2);
            g3_2_norm = g3_2./(g3_den_2 - g3_2)*((num_files)^2 - 1)/(((num_files/num_files_in_denom))^2);
            g3_3_norm = g3_3./(g3_den_3 - g3_3)*((num_files)^2 - 1)/(((num_files/num_files_in_denom))^2);
            g3_m1_norm = g3_m1./(g3_den_m1 - g3_m1)*((num_files)^2 - 1)/(((num_files/num_files_in_denom))^2);
            g3_m2_norm = g3_m2./(g3_den_m2 - g3_m2)*((num_files)^2 - 1)/(((num_files/num_files_in_denom))^2);
            g3_m3_norm = g3_m3./(g3_den_m3 - g3_m3)*((num_files)^2 - 1)/(((num_files/num_files_in_denom))^2);
            
            figure(999901)
            g3_0_norm_plot = surf(hist_vec,hist_vec,g3_0_norm);
            xlabel({x_title});
            ylabel({x_title});
            zlabel({y_title_g3});
            title(['Bins are ',num2str(handles.x_bin),title_1,num2str(handles.y_bin),title_2])
            saveas(g3_0_norm_plot,[handles.filename_input,'\',new_dir_name,'\g3_0_norm_plot.fig']);
            
            figure(999902)
            g3_1_norm_plot = surf(hist_vec,hist_vec,g3_1_norm);
            xlabel({x_title});
            ylabel({x_title});
            zlabel({y_title_g3});
            title(['Bins are ',num2str(handles.x_bin*1.2),title_1,num2str(handles.y_bin*1.2),title_2])
            saveas(g3_1_norm_plot,[handles.filename_input,'\',new_dir_name,'\g3_1_norm_plot.fig']);
            
            figure(999903)
            g3_2_norm_plot = surf(hist_vec,hist_vec,g3_2_norm);
            xlabel({x_title});
            ylabel({x_title});
            zlabel({y_title_g3});
            title(['Bins are ',num2str(handles.x_bin*1.4),title_1,num2str(handles.y_bin*1.4),title_2])
            saveas(g3_2_norm_plot,[handles.filename_input,'\',new_dir_name,'\g3_2_norm_plot.fig']);
            
            figure(999904)
            g3_3_norm_plot = surf(hist_vec,hist_vec,g3_3_norm);
            xlabel({x_title});
            ylabel({x_title});
            zlabel({y_title_g3});
            title(['Bins are ',num2str(handles.x_bin*1.6),title_1,num2str(handles.y_bin*1.6),title_2])
            saveas(g3_3_norm_plot,[handles.filename_input,'\',new_dir_name,'\g3_3_norm_plot.fig']);
            
            figure(999905)
            g3_m1_norm_plot = surf(hist_vec,hist_vec,g3_m1_norm);
            xlabel({x_title});
            ylabel({x_title});
            zlabel({y_title_g3});
            title(['Bins are ',num2str(handles.x_bin*0.8),title_1,num2str(handles.y_bin*0.8),title_2])
            saveas(g3_m1_norm_plot,[handles.filename_input,'\',new_dir_name,'\g3_m1_norm_plot.fig']);
            
            figure(999906)
            g3_m2_norm_plot = surf(hist_vec,hist_vec,g3_m2_norm);
            xlabel({x_title});
            ylabel({x_title});
            zlabel({y_title_g3});
            title(['Bins are ',num2str(handles.x_bin*0.6),title_1,num2str(handles.y_bin*0.6),title_2])
            saveas(g3_m2_norm_plot,[handles.filename_input,'\',new_dir_name,'\g3_m2_norm_plot.fig']);
            
            figure(999907)
            g3_m3_norm_plot = surf(hist_vec,hist_vec,g3_m3_norm);
            xlabel({x_title});
            ylabel({x_title});
            zlabel({y_title_g3});
            title(['Bins are ',num2str(handles.x_bin*0.4),title_1,num2str(handles.y_bin*0.4),title_2])
            saveas(g3_m3_norm_plot,[handles.filename_input,'\',new_dir_name,'\g3_m3_norm_plot.fig']);
        end
    end
end

newhandles = handles;

end

function w=fit_g2_properly(Param,g2_vec_norm,independent_vec)
%fit_g2_properly(Param,plotrange,g2_vec_norm,time_bin,independent_vec)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

one_vec = ones(1,max(size(independent_vec)));
w = sum((one_vec+Param(1).*exp(-((independent_vec).^2)/(2*(Param(2))^2)) - g2_vec_norm).^2);

% y=0;
% for count = 1:plotrange
%     y = y + (1+Param(1)*exp(-((count*time_bin)^2)/(2*(time_bin*Param(2))^2)) - g2_vec_norm(count))^2;
% end
% 
% w-y

end
























































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FRONT PANEL STUFF %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






function first_file_edit_Callback(hObject, eventdata, handles)
% hObject    handle to first_file_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of first_file_edit as text
%        str2double(get(hObject,'String')) returns contents of first_file_edit as a double
end

% --- Executes during object creation, after setting all properties.
function first_file_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to first_file_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function last_file_edit_Callback(hObject, eventdata, handles)
% hObject    handle to last_file_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of last_file_edit as text
%        str2double(get(hObject,'String')) returns contents of last_file_edit as a double
end

% --- Executes during object creation, after setting all properties.
function last_file_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to last_file_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function t_min_edit_Callback(hObject, eventdata, handles)
% hObject    handle to t_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_min_edit as text
%        str2double(get(hObject,'String')) returns contents of t_min_edit as a double
end

% --- Executes during object creation, after setting all properties.
function t_min_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function t_max_edit_Callback(hObject, eventdata, handles)
% hObject    handle to t_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_max_edit as text
%        str2double(get(hObject,'String')) returns contents of t_max_edit as a double
end

% --- Executes during object creation, after setting all properties.
function t_max_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function x_min_edit_Callback(hObject, eventdata, handles)
% hObject    handle to x_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_min_edit as text
%        str2double(get(hObject,'String')) returns contents of x_min_edit as a double
end

% --- Executes during object creation, after setting all properties.
function x_min_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function x_max_edit_Callback(hObject, eventdata, handles)
% hObject    handle to x_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_max_edit as text
%        str2double(get(hObject,'String')) returns contents of x_max_edit as a double
end

% --- Executes during object creation, after setting all properties.
function x_max_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function y_min_edit_Callback(hObject, eventdata, handles)
% hObject    handle to y_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_min_edit as text
%        str2double(get(hObject,'String')) returns contents of y_min_edit as a double
end

% --- Executes during object creation, after setting all properties.
function y_min_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function y_max_edit_Callback(hObject, eventdata, handles)
% hObject    handle to y_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_max_edit as text
%        str2double(get(hObject,'String')) returns contents of y_max_edit as a double
end

% --- Executes during object creation, after setting all properties.
function y_max_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function x_bin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to x_bin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_bin_edit as text
%        str2double(get(hObject,'String')) returns contents of x_bin_edit as a double
end

% --- Executes during object creation, after setting all properties.
function x_bin_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_bin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function t_bin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to t_bin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_bin_edit as text
%        str2double(get(hObject,'String')) returns contents of t_bin_edit as a double
end

% --- Executes during object creation, after setting all properties.
function t_bin_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_bin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function y_bin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to y_bin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_bin_edit as text
%        str2double(get(hObject,'String')) returns contents of y_bin_edit as a double
end

% --- Executes during object creation, after setting all properties.
function y_bin_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_bin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function number_bins_edit_Callback(hObject, eventdata, handles)
% hObject    handle to number_bins_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_bins_edit as text
%        str2double(get(hObject,'String')) returns contents of number_bins_edit as a double
end

% --- Executes during object creation, after setting all properties.
function number_bins_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_bins_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function files_per_chunk_edit_Callback(hObject, eventdata, handles)
% hObject    handle to files_per_chunk_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of files_per_chunk_edit as text
%        str2double(get(hObject,'String')) returns contents of files_per_chunk_edit as a double
end

% --- Executes during object creation, after setting all properties.
function files_per_chunk_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to files_per_chunk_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function temp_window_edit_Callback(hObject, eventdata, handles)
% hObject    handle to temp_window_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of temp_window_edit as text
%        str2double(get(hObject,'String')) returns contents of temp_window_edit as a double
end

% --- Executes during object creation, after setting all properties.
function temp_window_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to temp_window_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function num_window_edit_Callback(hObject, eventdata, handles)
% hObject    handle to num_window_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_window_edit as text
%        str2double(get(hObject,'String')) returns contents of num_window_edit as a double
end

% --- Executes during object creation, after setting all properties.
function num_window_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_window_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in tof_check.
function tof_check_Callback(hObject, eventdata, handles)
% hObject    handle to tof_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tof_check
end

% --- Executes on button press in spatial_check.
function spatial_check_Callback(hObject, eventdata, handles)
% hObject    handle to spatial_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spatial_check
end

% --- Executes on button press in corr_check.
function corr_check_Callback(hObject, eventdata, handles)
% hObject    handle to corr_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of corr_check
end

% --- Executes on button press in nvf_check.
function nvf_check_Callback(hObject, eventdata, handles)
% hObject    handle to nvf_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nvf_check
end

% --- Executes on button press in fit_remove_TOF_check.
function fit_remove_TOF_check_Callback(hObject, eventdata, handles)
% hObject    handle to fit_remove_TOF_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fit_remove_TOF_check
end

% --- Executes on button press in fit_T_check.
function fit_T_check_Callback(hObject, eventdata, handles)
% hObject    handle to fit_T_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fit_T_check
end


function t_bin_plot_edit_Callback(hObject, eventdata, handles)
% hObject    handle to t_bin_plot_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_bin_plot_edit as text
%        str2double(get(hObject,'String')) returns contents of t_bin_plot_edit as a double
end

% --- Executes during object creation, after setting all properties.
function t_bin_plot_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_bin_plot_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function bins_per_pixel_edit_Callback(hObject, eventdata, handles)
% hObject    handle to bins_per_pixel_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bins_per_pixel_edit as text
%        str2double(get(hObject,'String')) returns contents of bins_per_pixel_edit as a double
end

% --- Executes during object creation, after setting all properties.
function bins_per_pixel_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bins_per_pixel_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function dead_time_edit_Callback(hObject, eventdata, handles)
% hObject    handle to dead_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dead_time_edit as text
%        str2double(get(hObject,'String')) returns contents of dead_time_edit as a double
end

% --- Executes during object creation, after setting all properties.
function dead_time_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dead_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function t_min_exc_edit_Callback(hObject, eventdata, handles)
% hObject    handle to t_min_exc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_min_exc_edit as text
%        str2double(get(hObject,'String')) returns contents of t_min_exc_edit as a double
end

% --- Executes during object creation, after setting all properties.
function t_min_exc_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_min_exc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in exclusion_check.
function exclusion_check_Callback(hObject, eventdata, handles)
% hObject    handle to exclusion_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of exclusion_check
end


function t_max_exc_edit_Callback(hObject, eventdata, handles)
% hObject    handle to t_max_exc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_max_exc_edit as text
%        str2double(get(hObject,'String')) returns contents of t_max_exc_edit as a double
end

% --- Executes during object creation, after setting all properties.
function t_max_exc_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_max_exc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function x_min_exc_edit_Callback(hObject, eventdata, handles)
% hObject    handle to x_min_exc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_min_exc_edit as text
%        str2double(get(hObject,'String')) returns contents of x_min_exc_edit as a double
end

% --- Executes during object creation, after setting all properties.
function x_min_exc_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_min_exc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function x_max_exc_edit_Callback(hObject, eventdata, handles)
% hObject    handle to x_max_exc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_max_exc_edit as text
%        str2double(get(hObject,'String')) returns contents of x_max_exc_edit as a double
end

% --- Executes during object creation, after setting all properties.
function x_max_exc_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_max_exc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function y_min_exc_edit_Callback(hObject, eventdata, handles)
% hObject    handle to y_min_exc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_min_exc_edit as text
%        str2double(get(hObject,'String')) returns contents of y_min_exc_edit as a double
end

% --- Executes during object creation, after setting all properties.
function y_min_exc_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_min_exc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function y_max_exc_edit_Callback(hObject, eventdata, handles)
% hObject    handle to y_max_exc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_max_exc_edit as text
%        str2double(get(hObject,'String')) returns contents of y_max_exc_edit as a double
end

% --- Executes during object creation, after setting all properties.
function y_max_exc_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_max_exc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Outputs from this function are returned to the command line.
function varargout = DLD_front_panel_c_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


function filepath_edit_Callback(hObject, eventdata, handles)
% hObject    handle to filepath_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filepath_edit as text
%        str2double(get(hObject,'String')) returns contents of filepath_edit as a double
end

% --- Executes during object creation, after setting all properties.
function filepath_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filepath_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in fit_thermal_average_check.
function fit_thermal_average_check_Callback(hObject, eventdata, handles)
% hObject    handle to fit_thermal_average_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fit_thermal_average_check
end

% --- Executes on button press in fit_spatial_average.
function fit_spatial_average_Callback(hObject, eventdata, handles)
% hObject    handle to fit_spatial_average (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fit_spatial_average
end

% --- Executes on button press in fit_TF_average.
function fit_TF_average_Callback(hObject, eventdata, handles)
% hObject    handle to fit_TF_average (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fit_TF_average
end


function num_CPU_edit_Callback(hObject, eventdata, handles)
% hObject    handle to num_CPU_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_CPU_edit as text
%        str2double(get(hObject,'String')) returns contents of num_CPU_edit as a double
end

% --- Executes during object creation, after setting all properties.
function num_CPU_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_CPU_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in gpu_checkbox.
function gpu_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to gpu_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gpu_checkbox
end

% --- Executes on button press in save_check.
function save_check_Callback(hObject, eventdata, handles)
% hObject    handle to save_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_check
end

% --- Executes on button press in fit_TF_temp_check.
function fit_TF_temp_check_Callback(hObject, eventdata, handles)
% hObject    handle to fit_TF_temp_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fit_TF_temp_check

end

function f_x_edit_Callback(hObject, eventdata, handles)
% hObject    handle to f_x_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f_x_edit as text
%        str2double(get(hObject,'String')) returns contents of f_x_edit as a double

end

% --- Executes during object creation, after setting all properties.
function f_x_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f_x_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function f_y_edit_Callback(hObject, eventdata, handles)
% hObject    handle to f_y_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f_y_edit as text
%        str2double(get(hObject,'String')) returns contents of f_y_edit as a double

end

% --- Executes during object creation, after setting all properties.
function f_y_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f_y_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function f_z_edit_Callback(hObject, eventdata, handles)
% hObject    handle to f_z_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f_z_edit as text
%        str2double(get(hObject,'String')) returns contents of f_z_edit as a double

end

% --- Executes during object creation, after setting all properties.
function f_z_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f_z_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end



function temp_guess_tf_edit_Callback(hObject, eventdata, handles)
% hObject    handle to temp_guess_tf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of temp_guess_tf_edit as text
%        str2double(get(hObject,'String')) returns contents of temp_guess_tf_edit as a double
end

% --- Executes during object creation, after setting all properties.
function temp_guess_tf_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to temp_guess_tf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function tf_radius_guess_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tf_radius_guess_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tf_radius_guess_edit as text
%        str2double(get(hObject,'String')) returns contents of tf_radius_guess_edit as a double
end

% --- Executes during object creation, after setting all properties.
function tf_radius_guess_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tf_radius_guess_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function cond_percent_guess_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cond_percent_guess_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cond_percent_guess_edit as text
%        str2double(get(hObject,'String')) returns contents of cond_percent_guess_edit as a double
end

% --- Executes during object creation, after setting all properties.
function cond_percent_guess_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cond_percent_guess_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in force_recalc_txy_check.
function force_recalc_txy_check_Callback(hObject, eventdata, handles)
% hObject    handle to force_recalc_txy_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of force_recalc_txy_check
end



function num_files_denom_edit_Callback(hObject, eventdata, handles)
% hObject    handle to num_files_denom_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_files_denom_edit as text
%        str2double(get(hObject,'String')) returns contents of num_files_denom_edit as a double
end

% --- Executes during object creation, after setting all properties.
function num_files_denom_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_files_denom_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in sort_by_temp_check.
function sort_by_temp_check_Callback(hObject, eventdata, handles)
% hObject    handle to sort_by_temp_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sort_by_temp_check
end


% --- Executes on button press in lock_windows_check.
function lock_windows_check_Callback(hObject, eventdata, handles)
% hObject    handle to lock_windows_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lock_windows_check

if (get(hObject,'Value'))
    set(handles.first_file_edit, 'Enable', 'off');
    set(handles.last_file_edit, 'Enable', 'off');
    set(handles.t_min_edit, 'Enable', 'off');
    set(handles.t_max_edit, 'Enable', 'off');
    set(handles.x_min_edit, 'Enable', 'off');
    set(handles.x_max_edit, 'Enable', 'off');
    set(handles.y_min_edit, 'Enable', 'off');
    set(handles.y_max_edit, 'Enable', 'off');
    set(handles.t_min_exc_edit, 'Enable', 'off');
    set(handles.t_max_exc_edit, 'Enable', 'off');
    set(handles.x_min_exc_edit, 'Enable', 'off');
    set(handles.x_max_exc_edit, 'Enable', 'off');
    set(handles.y_min_exc_edit, 'Enable', 'off');
    set(handles.y_max_exc_edit, 'Enable', 'off');
    set(handles.save_check, 'Enable', 'off');
    set(handles.exclusion_check, 'Enable', 'off');
    set(handles.lock_window_string, 'String', 'will NOT recompute .bin files for c++');
    handles.is_window_locked = 1;
else
    set(handles.first_file_edit, 'Enable', 'on');
    set(handles.last_file_edit, 'Enable', 'on');
    set(handles.t_min_edit, 'Enable', 'on');
    set(handles.t_max_edit, 'Enable', 'on');
    set(handles.x_min_edit, 'Enable', 'on');
    set(handles.x_max_edit, 'Enable', 'on');
    set(handles.y_min_edit, 'Enable', 'on');
    set(handles.y_max_edit, 'Enable', 'on');
    set(handles.t_min_exc_edit, 'Enable', 'on');
    set(handles.t_max_exc_edit, 'Enable', 'on');
    set(handles.x_min_exc_edit, 'Enable', 'on');
    set(handles.x_max_exc_edit, 'Enable', 'on');
    set(handles.y_min_exc_edit, 'Enable', 'on');
    set(handles.y_max_exc_edit, 'Enable', 'on');
    set(handles.save_check, 'Enable', 'on');
    set(handles.exclusion_check, 'Enable', 'on');    
    set(handles.lock_window_string, 'String', 'WILL recompute .bin files for c++');
    handles.is_window_locked = 0;
end

guidata(hObject, handles);

end


% --- Executes on button press in clear_plot_windows_button.
function clear_plot_windows_button_Callback(hObject, eventdata, handles)
% hObject    handle to clear_plot_windows_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(DLD_front_panel_c, 'HandleVisibility', 'off');

close all

set(DLD_front_panel_c, 'HandleVisibility', 'on');

end


% --- Executes on button press in fit_remove_x0_check.
function fit_remove_x0_check_Callback(hObject, eventdata, handles)
% hObject    handle to fit_remove_x0_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fit_remove_x0_check

end

% --- Executes on button press in fit_remove_y0_check.
function fit_remove_y0_check_Callback(hObject, eventdata, handles)
% hObject    handle to fit_remove_y0_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fit_remove_y0_check

end


% --- Executes on button press in properties_button.
function properties_button_Callback(hObject, eventdata, handles)
% hObject    handle to properties_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = helpdlg(sprintf('Conversion will output files in TXY format (*_txy_forc_AGM_#.txt), same as older Matlab converter (DLD_raw_txy.m). This output does not depend on windowing etc.  If computed once, will not recompute unless "force recalculate TXY" chosen.  \n\nFiles output for correlations (*_correlations_fwrite_[t,x,y]#.bin) depend on window settings and will be rewritten every time, unless "lock windowing" selected. \n\nDead time for DLD to TXY conversion: 5 ns (constant).\nUser controlled dead time is for correlations only. \n\nRotation of 0.61 radians automatically applied to align DLD with trap axes. \n\nUsing 4 corners for conversion only, no reconstruction from 3 corners or MCP pulse.\n\nIf calculation fails, ensure that correct source files are selected (i.e. raw files direct from DLD, not "txy_forc" or .bin files).\n\nTo use this on a different computer, copy contents of "DLD_front_panel_c" Matlab folder. If conversion does not occur (.exe does not run), may need to install MS Visual Studio Express C++ 2012 (v11.0) or later for multicore libraries. Perhaps Visual C++ Redistributable for Visual Studio 2012 works, not tested. Compiled with both Win32 and x64 as separate .exe apps, so use appropriate version. Other packages such as Armadillo do not appear to be necessary.\n\nIf the channel information lights up in red, either one channel is at least 20 percent different to another channel, or rubbish counts are at least 10 percent of total channel counts.\n\n"Load old settings" is a bit broken, avoid using it.\n\nFitting is dodgy at the moment, currently working on fixing that.\n\nEnsure that correct processor architecture is chosen. Do not use x64 for an x86 (32 bit) computer/OS. \n\nTo use a custom file list, place a file named [folder name]_file_list.txt (e.g. 2012_08_24_MEW_file_list.txt) in the same folder as the raw files. Must be a column vector of numbers, delimited by newlines.'),'Properties of DLD_front_panel_c');

end


% --- Executes on button press in x64_check.
function x64_check_Callback(hObject, eventdata, handles)
% hObject    handle to x64_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of x64_check

end


% --- Executes on button press in custom_file_list_check.
function custom_file_list_check_Callback(hObject, eventdata, handles)
% hObject    handle to custom_file_list_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of custom_file_list_check
end


function custom_file_list_string_edit_Callback(hObject, eventdata, handles)
% hObject    handle to custom_file_list_string_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of custom_file_list_string_edit as text
%        str2double(get(hObject,'String')) returns contents of custom_file_list_string_edit as a double
end

% --- Executes during object creation, after setting all properties.
function custom_file_list_string_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to custom_file_list_string_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function rotation_edit_Callback(hObject, eventdata, handles)
% hObject    handle to rotation_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rotation_edit as text
%        str2double(get(hObject,'String')) returns contents of rotation_edit as a double
end

% --- Executes during object creation, after setting all properties.
function rotation_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rotation_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end