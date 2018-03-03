function varargout = GPLreview(varargin)
% GPLREVIEW MATLAB code for GPLreview.fig
%      GPLREVIEW, by itself, creates a new GPLREVIEW or raises the existing
%      singleton*.
%
%      H = GPLREVIEW returns the handle to a new GPLREVIEW or the handle to
%      the existing singleton*.
%
%      GPLREVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GPLREVIEW.M with the given input arguments.
%
%      GPLREVIEW('Property','Value',...) creates a new GPLREVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GPLreview_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GPLreview_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GPLreview

% Last Modified by GUIDE v2.5 27-Jan-2014 12:06:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name', mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GPLreview_OpeningFcn, ...
    'gui_OutputFcn',  @GPLreview_OutputFcn, ...
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


% --- Executes just before GPLreview is made visible.
function GPLreview_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GPLreview (see
 
handles.j = 1;
handles.marker_count = 0;
handles.dim_coords = 0;
handles.filter = 0;
handles.brightness = 0.3;
handles.NextFile = 0;

% Choose default command line output for GPLreview
handles.output = hObject;
set(handles.figure1,'KeyPressFcn',@initialize_buttons);

% Compute starting values
handles.FFTLVal = str2double(get(handles.FFTL,'String'));
handles.SampleFreqVal = str2double(get(handles.sample_freq,'String'));
handles.StartFreqVal = str2double(get(handles.start_freq,'String'));
handles.EndFreqVal = str2double(get(handles.end_freq,'String'));
handles.OverlapVal = str2double(get(handles.overlap,'String'));
handles.PlotLengthVal = str2double(get(handles.plot_length,'String'));
handles.PlotLengthVal = str2double(get(handles.plot_length,'String'));
handles.MarkerNumberVal = str2double(get(handles.marker_number,'String'));
handles.BufferVal = str2double(get(handles.buffer,'String'));
handles.Whiten = 0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GPLreview wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GPLreview_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function FFTL_Callback(hObject, ~, handles)
% hObject    handle to FFTL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.FFTLVal = str2double(get(handles.FFTL,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function FFTL_CreateFcn(hObject, ~, ~)
% hObject    handle to FFTL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function sample_freq_Callback(hObject, ~, handles)
% hObject    handle to sample_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SampleFreqVal = str2double(get(handles.sample_freq,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function sample_freq_CreateFcn(hObject, ~, handles)
% hObject    handle to sample_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function start_freq_Callback(hObject, ~, handles)
% hObject    handle to start_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.StartFreqVal = str2double(get(handles.start_freq,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function start_freq_CreateFcn(hObject, ~, handles)
% hObject    handle to start_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function end_freq_Callback(hObject, ~, handles)
% hObject    handle to end_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.EndFreqVal = str2double(get(handles.end_freq,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function end_freq_CreateFcn(hObject, ~, handles)
% hObject    handle to end_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function plot_length_Callback(hObject, ~, handles)
% hObject    handle to plot_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.PlotLengthVal = str2double(get(handles.plot_length,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function plot_length_CreateFcn(hObject, ~, handles)
% hObject    handle to plot_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function overlap_Callback(hObject, ~, handles)
% hObject    handle to overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.OverlapVal = str2double(get(handles.overlap,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function overlap_CreateFcn(hObject, ~, handles)
% hObject    handle to overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_spec.
function plot_spec_Callback(hObject, ~, handles)
% hObject    handle to plot_spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'WaveFile')
    error('Please identify audio data folder.')
end
if ~isfield(handles,'AudioData') || isempty(handles.AudioData)
    fprintf('Loading audio data\n for this window.')
    pushbutton2_Callback(hObject, 1, handles)
elseif isempty(handles.bt)
    disp('No detections in this file. Skipping.')
    handles.NextFile = 1;
    pushbutton2_Callback(hObject,[],handles)
else
    disp('Redrawing spectrogram.')
    draw_spectogram(handles);
end
guidata(hObject, handles);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % Move forward button '>'
if(handles.NextFile == 1) || handles.j==size(handles.bt,2)
    disp('Advancing to next detection file')
    % If we need to advance to next file, 
    bt = [];
    while isempty(bt) && handles.CurrentDetectionFileIdx<length(handles.DetectionFileList)
        % Load the next detection file and identify the associated audio file
        handles.CurrentDetectionFileIdx = handles.CurrentDetectionFileIdx+1;
        handles.DetectionFile = handles.DetectionFileList(handles.CurrentDetectionFileIdx).name;
        detectionFileStem = strrep(handles.DetectionFile,'.mat','');
        handles.CurrentWavIdx = find(~cellfun(@isempty,strfind({handles.WaveFileList.name},...
            detectionFileStem)));
        if isempty(handles.CurrentWavIdx)
            error('No audio file found to match detection file %s\n',handles.DetectionFile)
        else
            % storing file name of new wav file
            handles.WaveFile = fullfile(handles.WaveFileList(handles.CurrentWavIdx).pathName,...
                handles.WaveFileList(handles.CurrentWavIdx).name);
            % load new detection file
            fprintf('Opening %s\n',handles.DetectionFile)
            fprintf('Associated audio file is %s\n',handles.WaveFile)
            load(fullfile(handles.DetectionFilePath,handles.DetectionFile));
            if isempty(bt)
                fprintf('No detections in this file. Skipping.\n')
            else
                fprintf('This file contains %.f detections\n', size(bt,1))
            end
        end
    end
    handles.bt = bt;
    handles.j = 1;
    handles.ViewStart = 1;
    handles.marker_count = 0;
else
    disp('Moving forward in file.')
end
audioData = []; handles.AudioData=[]; handles.markers=0; length_index=0;

bt = handles.bt;
audioInfo = audioinfo(handles.WaveFile);
audioSize = audioInfo.TotalSamples;
buffer = handles.BufferVal*handles.SampleFreqVal;
handles.ViewStart = handles.j;

fprintf('Window starting at detection %0.0f\n',handles.ViewStart)
while(length(handles.AudioData) < handles.SampleFreqVal*handles.PlotLengthVal)...
        && (handles.j <= size(bt,1))
    % data=wavread(handles.WaveFile,[bt(handles.j,1)-buffer,bt(handles.j,2)+buffer]);

    % Read in the audio data that goes with this detection (currently
    % doesn't check if buffer reads into header).
    audioData = audioread(handles.WaveFile,[bt(handles.j,1)-buffer,...
        min(bt(handles.j,2) + buffer,audioSize)]);
    
    handles.AudioData=[handles.AudioData;audioData];
    length_index = length_index+length(audioData);
    handles.markers = [handles.markers,length_index];
    handles.j = handles.j+1;
end
handles.ViewEnd = handles.j-1;
fprintf('Window ending at detection %.0f\n',handles.ViewEnd)
set(handles.percent_completed,'String',(handles.ViewEnd)/size(handles.bt,1)*100)

handles.marker_count = handles.marker_count + length(handles.markers);

handles.plot_length_prev = get(handles.plot_length,'string');
if(handles.ViewEnd >= length(bt))
    set(handles.plot_length,'String',length(handles.AudioData)/...
        handles.SampleFreqVal)
end
guidata(hObject, handles);
plot_spec_Callback(hObject, eventdata, handles)

if(handles.ViewEnd>=size(bt,1))
    handles.NextFile = 1;
    fprintf('Reached end of this detection file.\n')
else
    handles.NextFile =0;
end
guidata(hObject, handles);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % Scroll backward button '<'
disp('Moving backward.')
if handles.ViewStart<=1
    % need to open the previous file
    if handles.CurrentDetectionFileIdx>1
        bt = [];
        while isempty(bt)
            handles.CurrentDetectionFileIdx = handles.CurrentDetectionFileIdx-1;
            handles.DetectionFile = handles.DetectionFileList(handles.CurrentDetectionFileIdx).name;
            fprintf('Opening previous file: %s\n', handles.DetectionFile)
            detectionFileStem = strrep(handles.DetectionFile,'.mat','');
            handles.CurrentWavIdx = find(~cellfun(@isempty,strfind({handles.WaveFileList.name},...
                detectionFileStem)));
            if isempty(handles.CurrentWavIdx)
                warning('No audio file found to match detection file %s\n',handles.DetectionFile)
            else
                % storing file name of new wav file
                handles.WaveFile = fullfile(handles.WaveFileList(handles.CurrentWavIdx).pathName,...
                    handles.WaveFileList(handles.CurrentWavIdx).name);
                % load new detection file
                fprintf('Opening %s\n',handles.DetectionFile)
                load(fullfile(handles.DetectionFilePath,handles.DetectionFile));
                if isempty(bt)
                    fprintf('No detections in this file, backing up further.\n')
                else
                    fprintf('This file contains %.f detections\n', size(bt,1))
                end
            end
        end
        % whos bt
        handles.bt = bt;
        
        handles.j = size(bt,1);
        handles.marker_count=0;
        %handles.reverse_vector=0;
        %handles.reverse_counter=0;
        handles.ViewStart = handles.j;
        handles.ViewEnd = handles.j;

    else
        fprintf('There are no earlier files in this folder\n')
    end
else
    handles.j = handles.ViewStart;
    handles.ViewEnd = handles.j-1;
    viewStart = handles.ViewStart-1;
end
%     guidata(hObject, handles);
audioInfo = audioinfo(handles.WaveFile);
audioSize = audioInfo.TotalSamples;

handles.AudioData = [];
buffer = handles.BufferVal*handles.SampleFreqVal;
handles.markers = 0;

fprintf('Window ending at detection %0.0f\n', handles.ViewEnd)
while(length(handles.AudioData) < handles.SampleFreqVal*handles.PlotLengthVal)...
        && (viewStart > 0)||(viewStart == 1)
    % Read in the audio data that goes with this detection (currently
    % doesn't check if buffer reads into header).
    audioData = audioread(handles.WaveFile,[handles.bt(viewStart,1)-buffer,...
        min(handles.bt(viewStart,2) + buffer,audioSize)]);
    
    % concatenate data in reverse
    handles.AudioData = [audioData;handles.AudioData];
    
    reverseLengthIndex = length(audioData);
    handles.markers = [0,handles.markers+reverseLengthIndex];
    viewStart = viewStart-1;
end
handles.ViewStart = viewStart+1;
fprintf('Window starting at detection %0.0f\n',handles.ViewStart)
set(handles.percent_completed,'String',handles.ViewEnd/size(handles.bt,1)*100)
handles.marker_count = handles.marker_count+length(handles.markers);

guidata(hObject, handles);
plot_spec_Callback(hObject, eventdata, handles)
guidata(hObject, handles);



% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.brightness = get(hObject,'Value');
guidata(hObject,handles);
plot_spec_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in deliminate_calls.
function deliminate_calls_Callback(hObject, eventdata, handles)
% hObject    handle to deliminate_calls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of deliminate_calls
handles.deliminate_calls = get(hObject,'Value');
guidata(hObject,handles);


% --- Executes on button press in mark.
function mark_Callback(hObject, eventdata, handles)
% hObject    handle to mark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles  structure with handles and user data (see GUIDATA)
% % Mark selected picks as true
coordinates = ginput(2);
handles.value_mark = handles.MarkerNumberVal;

[~, k2] = find(handles.markers/handles.OverlapVal > coordinates(1,1));
[~, k4] = find(handles.markers/handles.OverlapVal <= coordinates(2,1));

k5 = intersect(k2-1,k4);

handles.bt(handles.ViewStart+k5-1,3) = handles.value_mark;

bt = handles.bt;
save(strcat([handles.DetectionFilePath,handles.DetectionFile]), 'bt','-append')
guidata(hObject,handles);

plot_spec_Callback(hObject, eventdata, handles)


% --- Executes on button press in whiten.
function whiten_Callback(hObject, eventdata, handles)
% hObject    handle to whiten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Whiten = get(hObject,'Value');
guidata(hObject,handles);


% --- Executes on button press in play_audio.
function play_audio_Callback(hObject, eventdata, handles)
% hObject    handle to play_audio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

coordinates = ginput(2);
round(coordinates);

[~, k2] = find(handles.markers/handles.OverlapVal > coordinates(1,1));
[~, k4] = find(handles.markers/handles.OverlapVal <= coordinates(2,1));

%added so that Amanda can see exact time of audio that's playing:
k5 = intersect(k2-1,k4);
datestr(handles.bt(handles.ViewStart+k5(1)-1,4))
datestr(handles.bt(handles.ViewStart+k5(end),4))

if(handles.filter == 1)
    % Apply a band pass filter
    lower_freq = coordinates(1,2)/handles.FFTLVal*handles.SampleFreqVal;
    upper_freq = coordinates(2,2)/handles.FFTLVal*handles.SampleFreqVal;
    
    fs = handles.SampleFreqVal; % sampling rate
    
    F = [lower_freq-50, lower_freq, upper_freq, upper_freq+50];  % band limits
    A = [0 1 0];                % band type: 0='stop', 1='pass'
    dev = [0.0001, 10^(0.1/20)-1, 0.0001]; % ripple/attenuation spec
    [M,Wn,beta,typ] = kaiserord(F,A,dev,fs);  % window parameters
    b = fir1(M,Wn,typ,kaiser(M+1,beta),'noscale'); % filter design
    DATA = filter(b,1,handles.AudioData(coordinates(1,1)*...
        handles.OverlapVal:coordinates(2,1)...
        *handles.OverlapVal));
    
    handles.dim_coords = [floor(coordinates(1,2)),floor(coordinates(2,2))];
    
    plot_spec_Callback(hObject, eventdata, handles)
    
    if(handles.speedup==1)
        soundsc(DATA,5*handles.SampleFreqVal);
    else
        soundsc(DATA,handles.SampleFreqVal);
    end
    handles.dim_coords=0;
    plot_spec_Callback(hObject, eventdata, handles)
    
else
    
    if(handles.speedup==1)
        soundsc(handles.AudioData(coordinates(1,1)*handles.OverlapVal:coordinates(2,1)*handles.OverlapVal),5*handles.SampleFreqVal);
    else
        soundsc(handles.AudioData(coordinates(1,1)*handles.OverlapVal:coordinates(2,1)*handles.OverlapVal),handles.SampleFreqVal);
    end
end

guidata(hObject,handles);


% --- Executes on button press in filter.
function filter_Callback(hObject, eventdata, handles)
% hObject    handle to filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filter
handles.filter = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in mark_yes.
function mark_yes_Callback(hObject, eventdata, handles)
% hObject    handle to mark_yes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % Mark everything as true
fprintf('Labeling %0.0f detections as TRUE.\n',...
    handles.ViewEnd-handles.ViewStart)
handles.bt(handles.ViewStart:handles.ViewEnd,3) = 1;

bt = handles.bt;
save(strcat([handles.DetectionFilePath,handles.DetectionFile]), 'bt','-append')
guidata(hObject,handles);
% automatically advances? Why not just draw it?
pushbutton2_Callback(hObject, eventdata, handles)


% --- Executes on button press in mark_no.
function mark_no_Callback(hObject, eventdata, handles)
% hObject    handle to mark_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fprintf('Labeling %0.0f detections as FALSE.\n',...
    handles.ViewEnd-handles.ViewStart)
handles.bt(handles.ViewStart:handles.ViewEnd,3)=0;

bt = handles.bt;
save(strcat([handles.DetectionFilePath,handles.DetectionFile]), 'bt','-append')
guidata(hObject,handles);

% plot_spec_Callback(hObject, eventdata, handles)
pushbutton2_Callback(hObject, eventdata, handles)

% --- Executes on key press with focus on filter and none of its controls.
function filter_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to filter (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


function initialize_buttons(src, eventdata, handles, hObject)
%this function takes in two inputs by default

%src is the gui figure
%evnt is the keypress information

%this line brings the handles structures into the local workspace
%now we can use handles.cats in this subfunction!

handles = guidata(src);
hObject = handles.output;
%switch evnt.Key
switch eventdata.Key
    case 'leftarrow'
        pushbutton3_Callback(hObject, eventdata, handles)
    case 'rightarrow'
        pushbutton2_Callback(hObject, eventdata, handles)
    case 'y'
        mark_yes_Callback(hObject, eventdata, handles)
    case 'n'
        mark_no_Callback(hObject, eventdata, handles)
    case 'a'
        pushbutton10_Callback(hObject, eventdata, handles)
    case 'm'
        mark_Callback(hObject, eventdata, handles)
        
    case 'escape'
end;


% --------------------------------------------------------------------
function uipushtool3_ClickedCallback(hObject, eventdata, handles)
% Choose folder of detection files to be reviewed, and specify start file
% hObject    handle to uipushtool3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName,FilterIndex] = uigetfile('.mat',...
    'Select starting detection file');

if FilterIndex >0
    fprintf('Starting with detection file %s\n',FileName)
    addpath(PathName);
    handles.DetectionFileList = dir(strcat(PathName,'\*.mat')); % get list of files to process, 
    % assumes they are returned sorted by name, and therefore time, but
    % is that always true?
    handles.DetectionFile = FileName; % This is the starting file
    % get index of starting file
    handles.CurrentDetectionFileIdx = find(strcmp({(handles.DetectionFileList.name)}, FileName));
    handles.DetectionFilePath = PathName;
    
    load([PathName,handles.DetectionFile]);
    handles.bt = bt;
    handles.FileName = FileName;
    guidata(hObject,handles);
else
    error('No detection file selected. \n')
end


% --------------------------------------------------------------------
function uipushtool5_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PathName = uigetdir('','Select base folder of audio files');

% make sure detection file has been selected
if ~isfield(handles,'DetectionFile')
    error('Please select a detection file first')
end

% Search folder and subfolders for audio files
if PathName~=0
    dirList = dir(PathName); %are there subfolders?
    dirList = dirList(3:end);
    wavEnding = strfind({dirList.name},'.wav');
    wavFileFlag = ~cellfun(@isempty,wavEnding);
    fileList = [];
    if sum(wavFileFlag)>0
        fileList = dirList(~cellfun(@isempty,wavEnding));
        for iFile0 = 1:length(fileList)
            fileList(iFile0).pathName = fullfile(PathName);
        end
    end
    for iDir = 1:length(dirList)
        if dirList(iDir).isdir
            subDirList = dir(fullfile(PathName,dirList(iDir).name));
            wavEnding = strfind({subDirList.name},'.wav');
            newFileSet = subDirList(~cellfun(@isempty,wavEnding));
            for iFile = 1:length(newFileSet)
                newFileSet(iFile).pathName = fullfile(PathName,dirList(iDir).name);
            end
            fileList = [fileList;newFileSet];
        end
    end
    %addpath(PathName);
    handles.WaveFileList = fileList;
    targetFileName = strrep(handles.DetectionFile,'.mat','');
    findUscores = strfind(targetFileName,'_');
    matchNameTF = strfind({fileList(:).name},targetFileName);
    audioIdx = find(~cellfun(@isempty,matchNameTF));
    if isempty(audioIdx)
        error('No audio file found to match detection file %s\n',handles.DetectionFile)
    else
        handles.WaveFile = fullfile(fileList(audioIdx).pathName,fileList(audioIdx).name);
        handles.CurrentWavIdx = audioIdx;        
        fprintf('Found matching wav file %s\n',handles.WaveFile)

    end
    guidata(hObject,handles);
else
    fprintf('No audio file folder selected.\n')
end
%%

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fprintf('Labeling %0.0f detections as airguns.\n',...
    handles.ViewEnd-handles.ViewStart)
handles.bt(handles.ViewStart:handles.ViewEnd,3)=2;

bt = handles.bt;
save(strcat([handles.DetectionFilePath,handles.DetectionFile]), 'bt','-append')
guidata(hObject,handles);

pushbutton2_Callback(hObject, eventdata, handles)


% --- Executes on button press in speedup.
function speedup_Callback(hObject, eventdata, handles)
% hObject    handle to speedup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.speedup = get(hObject,'Value');
guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of speedup


function marker_number_Callback(hObject, eventdata, handles)
% hObject    handle to marker_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.MarkerNumberVal = str2double(get(handles.marker_number,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function marker_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to marker_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.marker_number=get(hObject,'Value');
handles.marker_number

guidata(hObject,handles);


function buffer_Callback(hObject, eventdata, handles)
% hObject    handle to buffer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.BufferVal = str2double(get(handles.buffer,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function buffer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to buffer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mark_subset_no.
function mark_subset_no_Callback(hObject, eventdata, handles)
% hObject    handle to mark_subset_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

coordinates = ginput(2);
handles.value_mark = str2double(get(handles.marker_number,'String'));

% identify which detections were picked
[~, k2] = find(handles.markers/handles.OverlapVal > coordinates(1,1));
[~, k4] = find(handles.markers/handles.OverlapVal <= coordinates(2,1));

k5 = intersect(k2-1,k4);
% Flag them as zeros
fprintf('Flagging %.0f detection(s) as false.\n',length(k5))
handles.bt(handles.ViewStart+k5-1,3) = 0;

bt = handles.bt;
save(strcat([handles.DetectionFilePath,handles.DetectionFile]), 'bt','-append')
guidata(hObject,handles);

plot_spec_Callback(hObject, eventdata, handles)