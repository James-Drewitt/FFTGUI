function varargout = FFTGUI(varargin)
% FFTGUI M-file for FFTGUI.fig
%      FFTGUI, Fast Fourier Transform GUI
%      James W E Drewitt
%      james.drewitt@bristol.ac.uk
%
%      H = FFTGUI returns the handle to a new FFTGUI or the handle to
%      the existing singleton*.
%
%      FFTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FFTGUI.M with the given input arguments.
%
%      FFTGUI('Property','Value',...) creates a new FFTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FFTGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FFTGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FFTGUI

% Last Modified by GUIDE v2.5 30-Aug-2011 18:15:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FFTGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @FFTGUI_OutputFcn, ...
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

% --- Executes just before FFTGUI is made visible.
function FFTGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FFTGUI (see VARARGIN)

% Choose default command line output for FFTGUI
handles.output = hObject;

set(handles.DoFFT,'enable','off')
set(handles.set_lowr,'enable','off')
set(handles.DoBFFT,'enable','off')
set(handles.Spline_fit,'enable','off')
set(handles.Lorch_window,'enable','off')
set(handles.CosineWindow,'enable','off')
set(handles.Reset,'enable','off')

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes FFTGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = FFTGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Open_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear SQ gr;
axes(handles.axes1)
sqdat=uigetfile('*.dat', 'Select S(Q) file');
global SQ
SQ=importdata(sqdat);
plot(SQ(:,1),SQ(:,2));
xlabel('Scattering vector, Q (1/A)');
ylabel('Structure factor, S(Q)');
global SQ2
SQ2=SQ;
set(handles.DoFFT,'enable','on')
set(handles.Reset,'enable','on')
set(handles.set_lowr,'enable','off')
set(handles.DoBFFT,'enable','off')
set(handles.Spline_fit,'enable','on')
set(handles.Lorch_window,'enable','on')
set(handles.CosineWindow,'enable','on')

% --------------------------------------------------------------------
function Fourier_Callback(hObject, eventdata, handles)
% hObject    handle to Fourier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function DoFFT_Callback(hObject, eventdata, handles)
% hObject    handle to DoFFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Fast Fourier Transform of S(Q)
global SQ2
global gr
global gr2
global rho
axes(handles.axes1)
[gr]=FastFT(SQ2,rho,13,27);
        plot(gr(:,1),gr(:,2));
        xlim([0 8]);
        xlabel('Distance r (A)');
        ylabel('Pair distribution function G(r)');
set(handles.DoFFT,'enable','off')
set(handles.Reset,'enable','on')
set(handles.DoBFFT,'enable','on')
set(handles.set_lowr,'enable','on')
set(handles.Spline_fit,'enable','off')
set(handles.Lorch_window,'enable','off')
set(handles.CosineWindow,'enable','off')
dlmwrite('Gr.dat',gr,'delimiter','\t','precision','%.4f')
gr2=gr;

% --------------------------------------------------------------------
function DoBFFT_Callback(hObject, eventdata, handles)
% hObject    handle to DoBFFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Back transform: from G(r) to BT(Q)
global SQ
global gr2
global rho
gr2(:,2)=gr2(:,2)-1;
axes(handles.axes1)
[BT]=BFFT(gr2,rho,13);
        plot(SQ(:,1),SQ(:,2),BT(:,1),BT(:,2));
        xlim([0,25])
        xlabel('Scattering vector, Q (1/A)');
        ylabel('Structure factor, S(Q)');
set(handles.DoFFT,'enable','on')
set(handles.Reset,'enable','on')
set(handles.DoBFFT,'enable','off')
set(handles.set_lowr,'enable','off')
set(handles.Spline_fit,'enable','on')
set(handles.Lorch_window,'enable','on')
set(handles.CosineWindow,'enable','on')

% --------------------------------------------------------------------
function menu_Fit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Menu_Mod_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_Mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Lorch_window_Callback(hObject, eventdata, handles)
% hObject    handle to Lorch_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global SQ
global SQ2
SQ2=SQ;
axes(handles.axes1)
Qmax=(length(SQ2(:,1))-1)*(SQ2(2,1)-SQ2(1,1))
SQ2(:,2)=SQ2(:,2).*((sin((pi.*SQ2(:,1))./Qmax))./((pi.*SQ2(:,1))./Qmax));
plot(SQ(:,1),SQ(:,2),SQ2(:,1),SQ2(:,2));
xlabel('Scattering vector, Q (1/A)');
ylabel('Structure factor, S(Q)');

% --------------------------------------------------------------------
function CosineWindow_Callback(hObject, eventdata, handles)
% hObject    handle to CosineWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
global SQ
global SQ2
SQ2=SQ;
axes(handles.axes1)
qstep=SQ2(2,1)-SQ2(1,1);
cosmax=(length(SQ2(:,1))-1)*qstep
prompt={'First Q value for cosine window function: '};
title='Cosine Window Function';
default_ans={'10'};
cospar=inputdlg(prompt,title,1,default_ans);
cosmin=str2num(cospar{1});
n=(cosmax-cosmin)/qstep;
x=1:n;
Win(x)=0+((pi*(x/n))/2);
for j=1+(cosmin/qstep):(cosmax/qstep);
    SQ2(j,2)=SQ2(j,2)*cos(Win(j-(cosmin/qstep)));
end
plot(SQ(:,1),SQ(:,2),SQ2(:,1),SQ2(:,2));
xlabel('Scattering vector, Q (1/A)');
ylabel('Structure factor, S(Q)');

% --------------------------------------------------------------------
function Spline_fit_Callback(hObject, eventdata, handles)
% hObject    handle to Spline_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Make a spline fit to the S(Q)
global SQ
global SQ2
axes(handles.axes1)
x=SQ(:,1);
y=SQ(:,2);
prompt={'Qmax: ', 'Qstep: ', 'Number of knots: '};
title='Set fit parameters';
default_ans={'23','0.05','150'};
fitpar=inputdlg(prompt,title,1,default_ans);
Qmax=str2num(fitpar{1});
Qstep=str2num(fitpar{2});
knots=str2num(fitpar{3});
x2=(0:Qstep:Qmax);
PP=splinefit(x,y,knots);
y2=ppval(PP,x2);
plot(x,y,x2,y2);
xlabel('Scattering vector, Q (1/A)');
ylabel('Structure factor, S(Q)');
OUT(:,1)=x2(:);
OUT(:,2)=y2(:);
SQ2=OUT;

% --------------------------------------------------------------------
function Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global SQ
global SQ2
axes(handles.axes1)
plot(SQ(:,1),SQ(:,2));
xlabel('Scattering vector, Q (1/A)');
ylabel('Structure factor, S(Q)');
set(handles.DoFFT,'enable','on')
set(handles.Reset,'enable','on')
set(handles.DoBFFT,'enable','off')
set(handles.set_lowr,'enable','off')
set(handles.Spline_fit,'enable','on')
set(handles.Lorch_window,'enable','on')
set(handles.CosineWindow,'enable','on')
SQ2=SQ;

% --------------------------------------------------------------------
function set_lowr_Callback(hObject, eventdata, handles)
% hObject    handle to set_lowr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gr
global gr2
gr2=gr;
prompt={'No. low-r data points: ', 'G(r=0) limit: '};
title='Set low-r';
default_ans={'24','0'};
lowr=inputdlg(prompt,title,1,default_ans);
Num=str2num(lowr{1});
G0=str2num(lowr{2});
cut=1:Num;
gr2(cut,2)=G0;
plot(gr2(:,1),gr2(:,2),gr(cut,1),gr(cut,2));
xlim([0 8]);
        xlabel('Distance r (A)');
        ylabel('Pair distribution function G(r)');

function set_rho_Callback(hObject, eventdata, handles)
% hObject    handle to set_rho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_rho as text
%        str2double(get(hObject,'String')) returns contents of set_rho as a double
global rho
rho=str2num(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function set_rho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_rho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global rho
rho=str2num(get(hObject,'String'));


% --------------------------------------------------------------------
function Analyse_Callback(hObject, eventdata, handles)
% hObject    handle to Analyse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function CN_Callback(hObject, eventdata, handles)
% hObject    handle to CN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gr
prompt={'low r: ', 'high r: ', 'xlim: ', 'ylim:'};
title='Set r-limits';
default_ans={'51','71','1','2.5'};
fitpar=inputdlg(prompt,title,1,default_ans);
lowr=str2num(fitpar{1});
hir=str2num(fitpar{2});
xlim=str2num(fitpar{3});
ylim=str2num(fitpar{4});
N=size(gr,1);
x=1:lowr;
y=hir:N;
figure(1)
grCN(:,1)=gr(:,1);
grCN(:,2)=gr(:,2);
grCN(x,2)=0;
grCN(y,2)=0;
subplot(2,1,1);plot(gr(:,1),gr(:,2),grCN(:,1),grCN(:,2));
axis([xlim ylim -1 3])
subplot(2,1,2);plot(gr(:,1),gr(:,2));
axis([xlim ylim -1 3])
Ap=Ipeak(grCN);
Ap(71,2)

% --------------------------------------------------------------------
function EXIT_Callback(hObject, eventdata, handles)
% hObject    handle to EXIT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear all
clear global
close all force