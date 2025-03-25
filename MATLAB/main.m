clc
clearvars 
close all


%% Connect to Arduino
port = instrhwinfo('serial');                                               %Grab information about the available serial ports.
if isempty(port)                                                            %If no serial ports were found on this computer...
    error(['ERROR IN CONNECT_StressTest: There are no available serial '...
        'ports on this computer.']);                                        %Show an error.
end
port = port.SerialPorts                                                     %Pair down the list of serial ports to only those available.
busyports = instrfind;                                                      %Grab all the ports currently in use.
if ~isempty(busyports)                                                      %If there are any ports currently in use...
    busyports = {busyports.Port};                                           %Make a list of their port addresses.
    if iscell(busyports{1})                                                 %If there's more than one busy port.
        busyports = busyports{1}';                                          %Kick out the extraneous port name parts.
    end
else                                                                        %Otherwise...
    busyports = {};                                                         %Make an empty list for comparisons.
end

uiheight = 2;                                                               %Set the height for all buttons.
temp = [10,length(port)*(uiheight+0.1)+0.1];                                %Set the width and height of the port selection figure.
set(0,'units','centimeters');                                               %Set the screensize units to centimeters.
pos = get(0,'ScreenSize');                                                  %Grab the screensize.
pos = [pos(3)/2-temp(1)/2,pos(4)/2-temp(2)/2,temp(1),temp(2)];              %Scale a figure position relative to the screensize.
fig1 = figure('units','centimeters',...
    'Position',pos,...
    'resize','off',...
    'MenuBar','none',...
    'name','Select A Serial Port',...
    'numbertitle','off');                                                   %Set the properties of the figure.
        
for i = 1:length(port);
    if any(strcmpi(port{i}, busyports))
        txt = [' (' port{i} '): busy (reset?)'];
    else
        txt = [' (' port{i} '): available'];
    end
            
           
uicontrol(fig1,'style','pushbutton',...
    'string',txt,...
    'units','centimeters',...
    'position',[0.1 temp(2)-i*(uiheight+0.1) 9.8 uiheight],...
    'fontweight','bold',...
    'fontsize',14,...
    'callback',...
    ['guidata(gcbf,' num2str(i) '); uiresume(gcbf);']);                     %Make a button for the port showing that it is busy.
end
        
uiwait(fig1);
if ishandle(fig1)
    i = guidata(fig1);
    port = port{i};
    close(fig1);
    checkreset = 0;
else
    port = [];
end        

listbox = [];                                                               %Create a variable to hold a listbox handle.

if ~isempty(port) && checkreset && any(strcmpi(port,busyports))             %If that serial port is busy...
    i = questdlg(['Serial port ''' port ''' is busy. Reset and use '...
        'this port?'],['Reset ''' port '''?'],'Reset','Cancel','Reset');    %Ask the user if they want to reset the busy port.
    if strcmpi(i,'Cancel')                                                  %If the user selected "Cancel"...
        port = [];                                                          %...set the selected port to empty.
    end
end


if any(strcmpi(port,busyports))                                             %If the specified port is already busy...
    i = find(strcmpi(port,busyports));                                      %Find the index of the specified ports in the list of all busy ports.
    temp = instrfind;                                                       %Grab all the open serial connections.
    fclose(temp(i));                                                        %Close the busy serial connection.
    delete(temp(i));                                                        %Delete the existing serial connection.
end

serialcon = serial(port,'baudrate',9600);                                   %Set up the serial connection on the specified port.

try                                                                         %Try to open the serial port for communication.
    fopen(serialcon);                                                       %Open the serial port.
catch err                                                                   %If no connection could be made to the serial port...
    delete(serialcon);                                                      %...delete the serial object...
    error(['ERROR IN CONNECT: Could not open a serial '...
        'connection on port ''' port ''' (' err.message ').']);             %Show an error.
end

%% Sampling frequency
fs = 360;
ts = 1/fs;

%% Train and test records
recs_train = ["100", "105", "106"]; % 100, 105, 106, 209, 220
recs_test = ["209", "220"]; % 100, 105, 106, 209, 220

%% Train and save SVM classifier
SVMModel = train_svm_classifier(recs_train);

%% Create continous connection to app
quitProg = "n";

while quitProg == "n"
    
    %% Select test signal
    t_select = 20; % time in seconds
    ecg = generate_test_signal(t_select, recs_test);

    %% Prep signal for plotting
    positiveSignal = ecg+abs(min(ecg));
    rangeOfSignal = max(positiveSignal);
    positiveSignal = positiveSignal/(rangeOfSignal/255)

    %% Plot 10 seconds of signal on App
    for i = 1:12:3600
        fwrite(serialcon,positiveSignal(i:i+12),"uint8");
        pause(0.05)
    end
    
    %% Wait for button press
    decision = 0;
    t0 = clock;
    while decision == 0
        if serialcon.BytesAvailable >= 1
            decision = fread(serialcon);
        end
        if etime(clock, t0) > 30
            break
        end
    end

    %% Return answer to app
    fwrite(serialcon,decision,"uint8");
    
    %% Classify test signal if button is pressed
    n_abnorm = 0;
    n_total = 0;
    
    if decision == 1
        [n_abnorm, n_total] = classify_beats(ecg, SVMModel);
        fwrite(serialcon,"Number of abnormal beats detected: "+n_abnorm+"                Total beats analyzed: "+n_total,"uint8");
    end
    
    %% Wait a bit before starting another test session
    t0 = clock;
    while etime(clock, t0) < 30

    end
    
    
end






