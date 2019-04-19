
format long
close all
clc

% STUDY OF THE INTERLINK BETWEEN SMALL SATELLITES IN A CONSTELLATION
% Author: David Puig Puig
% Tutor: Miquel Sureda Anfres
% ESEIAAT - UPC

% Visual contact for two satellites analysis

% Reminder: All times are in UTC

%% Menu module

% Introduction and information

input_tle_list = {'Author: David Puig', 'Tutor: Miquel Sureda', 'ESEIAAT - UPC'};
[indx,tf] = listdlg('ListString',input_tle_list,'Name','InterLink','PromptString','This tool is used to analyse visibility windows in satellite constellations',...
                    'SelectionMode','single','ListSize',[500,300],'OKString','Next','CancelString','Quit');

if tf == 0
    disp('User selected Quit');
    return
end

%% Input Parameters module

% Input Celestial Object System

input_tle_list = {'Earth', 'Other'};
[indx,tf] = listdlg('ListString',input_tle_list,'Name','Celestial Object System','PromptString','Select your analysis system:',...
                    'SelectionMode','single','ListSize',[500,300],'OKString','Next','CancelString','Quit');

if tf == 0
    disp('User selected Quit');
    return
end

% Common System Parameters

global mu;
k = 2*pi;      % Factor from [rev/s] to [rad/s]

if indx == 1
    % Earth System parameters
    body_radius = 6.378e6;                          % Radius of the primary body [m]
    extra_radius = 0;                               % Extra radius for the primary body [m]
    S = body_radius + extra_radius;                 % Magnitude of the rise-set vector [m]
    mu = 3.986004418e14;                            % Standard gravitational parameter [m^3/s^2]
else
    % Other system parameters
    prompt = {'Primary body radius [m]:', 'Extra radius (atmosphere and other effects) [m]:', 'Mu parameter [m^3/s^2]:'};
    dlgtitle = 'Celestial Object System';
    dims = [1 70; 1 70; 1 70];
    system_answer = inputdlg(prompt,dlgtitle,dims);
    body_radius = str2double(system_answer{1});     % Radius of the primary body [m]
    extra_radius = str2double(system_answer{2});    % Extra radius for the primary body [m]
    S = body_radius + extra_radius;                 % Magnitude of the rise-set vector [m]
    mu = str2double(system_answer{1});              % Standard gravitational parameter [m^3/s^2]
end

% TLE input menu
input_tle_list = {'Examples', 'From .txt file (without blank lines between set)', 'Paste'};
[indx,tf] = listdlg('ListString',input_tle_list,'Name','Two Line Element (TLE)','PromptString','Select a TLE input mode:',...
                    'SelectionMode','single','ListSize',[500,300],'OKString','Next','CancelString','Quit');

if tf == 0
    disp('User selected Quit');
    return
end

if indx == 1
    input_examples_list = {'EGYPTSAT 1', 'TRMM', 'GOES 3', 'NOAA 3', 'NAVSTAR 46'};
    [indx,tf] = listdlg('ListString',input_examples_list,'Name','Two Line Element (TLE)','PromptString','Select two or more TLE to analyse:',...
                        'SelectionMode','multiple','ListSize',[500,300],'OKString','Run','CancelString','Quit');
    
    % Hard-Coded TLE as input examples
    possible_example_answers = {{'EGYPTSAT 1                                                           ';
                                '1 31117U 07012A 08142.74302347 .00000033 00000-0 13654-4 0 2585      ';
                                '2 31117 098.0526 218.7638 0007144 061.2019 298.9894 14.69887657 58828'};
                                {'TRMM                                                                 ';
                                '1 25063U 97074A 08141.84184490 .00002948 00000-0 41919-4 0 7792      ';
                                '2 25063 034.9668 053.5865 0001034 271.1427 088.9226 15.55875272598945'};
                                {'GOES 3                                                               ';
                                '1 10953U 78062A 08140.64132336 -.00000110 00000-0 10000-3 0 1137     ';
                                '2 10953 014.2164 003.1968 0001795 336.4858 023.4617 01.00280027 62724'};
                                {'NOAA 3                                                               ';
                                '1 06920U 73086A 08141.92603915 -.00000030 00000-0 +10000-3 0 00067   ';
                                '2 06920 101.7584 171.9430 0006223 187.3360 172.7614 12.40289355563642'};
                                {'NAVSTAR 46                                                           ';
                                '1 25933U 99055A 08142.14123352 .00000019 00000-0 10000-3 0 00126     ';
                                '2 25933 051.0650 222.9439 0079044 032.8625 327.6958 02.00568102 63184'};
                                {'EGYPTSAT 1                                                           ';
                                '1 31117U 07012A 08142.74302347 .00000033 00000-0 13654-4 0 2585      ';
                                '2 31117 098.0526 218.7638 0007144 061.2019 298.9894 14.69887657 58828'}};
              
    if tf == 0
        disp('User selected Quit');
        return
    end
    
    if size(indx) == 1
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        msgbox('A minimum of two TLE set are needed to compute visibility','Error',CreateStruct);
        return
    else
        % TLE variables extraction
        selected_example_answers = cell(1,1);
        count=1;
        for i=1:size(indx,2)  
            for j=1:3
            selected_example_answers{1}{count,1} = possible_example_answers{indx(i),1}{j};
            count=count+1;
            end
        end
        
        % Find the total number of satellites in the file
        num_satellites = size(indx,2);

        % Initialize array
        sat_id_line = zeros(1,num_satellites);
        line_count = 1;
        for i=1:num_satellites
            % Take every 3rd line
            sat_id_line(i) = line_count;
            txt_data = textscan(selected_example_answers{1}{line_count,1},'%s %s %s %s %s %s %s %s %s');

            OrbitData.ID(i) = txt_data{1};
            if isempty(txt_data{2})
                OrbitData.designation(i) = {''};
            else
                OrbitData.designation(i) = txt_data{2};
            end
            
            if isempty(txt_data{3})
                OrbitData.PRN(i) = {''};
            else
                OrbitData.PRN(i) = txt_data{3};
            end

            % Jump to the next Satellite Name line
            line_count  = line_count + 3;
        end

        % Find the two lines corresponding to the spacecraft in question
        for j=1:length(sat_id_line)

            % Find the first line of the first satellite
            index = sat_id_line(j);
            txt_data_second = textscan(selected_example_answers{1}{index+2,1},'%s %s %s %s %s %s %s %s %s');

            % Translate two line element data into obital elements
            OrbitData.i(j)     = str2double(txt_data_second{1,3})*(pi/180); % Inclination [deg] to [rad/s]
            OrbitData.RAAN(j)  = str2double(txt_data_second{1,4})*(pi/180); % Right ascention of the ascending node[deg] to [rad/s]
            OrbitData.omega(j) = str2double(txt_data_second{1,6})*(pi/180); % Argument of the periapsis [deg] to [rad/s]
            OrbitData.M(j)     = str2double(txt_data_second{1,7})*(pi/180); % Mean anomaly [deg] to [rad/s]
            n                  = str2double(txt_data_second{1,8});          % Unperturbed mean motion [rev/day]
            OrbitData.n(j)     = n*2*pi/24/60/60;                           % Unperturbed mean motion [rad/s]
            OrbitData.a(j)     = ( mu / OrbitData.n(j)^2 )^(1/3);           % Semi-major axis [m]
            OrbitData.e(j)     = str2double(txt_data_second{1,5})*1e-7;     % Eccentricity [unitless]

            % Compute the UTC date / time
            txt_data_first    = textscan(selected_example_answers{1}{index+1,1},'%s %s %s %s %s %s %s %s %s');
            temp2             = txt_data_first{1,4};
            yy                = str2double(temp2{1}(1:2));
            yyyy              = 2000 + yy;
            start             = datenum([yyyy 0 0 0 0 0]);
            secs              = str2double(temp2{1}(3:length(temp2{1})))*24*3600;
            date1             = datevec(addtodate(start,floor(secs),'second'));
            remainder         = [0 0 0 0 0 mod(secs,1)];
            OrbitData.date{j} = datestr(date1+remainder,'dd-mmm-yyyy HH:MM:SS.FFF');

            % Compute ballistic coefficient in SI units
            temp3 = txt_data_first{1,7};
            if length(temp3{1}) == 7
                base  = str2double(temp3{1}(1:5));
                expo  = str2double(temp3{1}(6:7));
            elseif length(temp3{1}) == 8
                base  = str2double(temp3{1}(2:6));
                expo  = str2double(temp3{1}(7:8));
            else
                fprintf('Error in ballistic coefficient calculation\n')
                CreateStruct.Interpreter = 'tex';
                CreateStruct.WindowStyle = 'modal';
                msgbox('Error in ballistic coefficient calculation\n','Error',CreateStruct);            
                error('End program')
            end

            Bstar = base*10^expo;
            OrbitData.BC(j) = 1/12.741621/Bstar;                            % Ballistic coefficient [kg/m^2]

        end
    end
    
elseif indx == 2

    [file,path] = uigetfile('*.txt');

    if isequal(file,0)
        disp('User selected Cancel');
        return
    
    else
        disp(['User selected ', fullfile(path,file)]);
        % TLE file name and variables extraction
        fid_input = fopen(fullfile(path,file));
        txt_data = textscan(fid_input,'%s %s %s %s %s %s %s %s %s');

        % Find the total number of satellites in the file
        num_satellites = length(txt_data{1})/3;

        % Initialize array
        sat_id_line = zeros(1,num_satellites);
        line_count = 1;
        for i=1:num_satellites
            % Take every 3rd line
            sat_id_line(i) = line_count;

            OrbitData.ID(i) = txt_data{1}(line_count);
            if isempty(txt_data{2}(line_count))
                OrbitData.designation(i) = {''};
            else
                OrbitData.designation(i) = txt_data{2}(line_count);
            end
            
            if isempty(txt_data{3}(line_count))
                OrbitData.PRN(i) = {''};
            else
                OrbitData.PRN(i) = txt_data{3}(line_count);
            end
            
            % Jump to the next Satellite Name line
            line_count  = line_count + 3;
        end

        % Find the two lines corresponding to the spacecraft in question
        for j=1:length(sat_id_line)
            
            % Find the first line of the first satellite
            index = sat_id_line(j);

            % Translate two line element data into obital elements
            OrbitData.i(j)     = str2double(txt_data{1,3}{index+2})*(pi/180);           % Inclination [deg] to [rad/s]
            OrbitData.RAAN(j)  = str2double(txt_data{1,4}{index+2})*(pi/180);           % Right ascention of the ascending node[deg] to [rad/s]
            OrbitData.omega(j) = str2double(txt_data{1,6}{index+2})*(pi/180);           % Argument of the periapsis [deg] to [rad/s]
            OrbitData.M(j)     = str2double(txt_data{1,7}{index+2})*(pi/180);           % Mean anomaly [deg] to [rad/s]
            n                  = str2double(txt_data{1,8}{index+2});                    % Unperturbed mean motion [rev/day]
            OrbitData.n(j)     = n*2*pi/24/60/60;                                       % Unperturbed mean motion [rad/s]
            OrbitData.a(j)     = ( mu / OrbitData.n(j)^2 )^(1/3);                                    % Semi-major axis [m]
            OrbitData.e(j)     = str2double(txt_data{1,5}{index+2})*1e-7;               % Eccentricity [unitless]

            % Compute the UTC date / time
            temp2             = txt_data{1,4}{index+1};
            yy                = str2double(temp2(1:2));
            yyyy              = 2000 + yy;
            start             = datenum([yyyy 0 0 0 0 0]);
            secs              = str2double(temp2(3:length(temp2)))*24*3600;
            date1             = datevec(addtodate(start,floor(secs),'second'));
            remainder         = [0 0 0 0 0 mod(secs,1)];
            OrbitData.date{j} = datestr(date1+remainder,'dd-mmm-yyyy HH:MM:SS.FFF');

            % Compute ballistic coefficient in SI units
            temp3 = txt_data{1,7}{index+1};
            if length(temp3) == 7
                base  = str2double(temp3(1:5));
                expo  = str2double(temp3(6:7));
            elseif length(temp3) == 8
                base  = str2double(temp3(2:6));
                expo  = str2double(temp3(7:8));
            else
                fprintf('Error in ballistic coefficient calculation\n')
                CreateStruct.Interpreter = 'tex';
                CreateStruct.WindowStyle = 'modal';
                msgbox('Error in ballistic coefficient calculation\n','Error',CreateStruct);
                error('end program')
                return
            end

            Bstar = base*10^expo;
            OrbitData.BC(j) = 1/12.741621/Bstar; % [kg/m^2]
            
        end
    end

elseif indx == 3
    
    prompt = 'How many TLE do you want to analyse?';
    dlgtitle = 'Paste TLE';
    dims = [1 69];
    tle_num_answer = inputdlg(prompt,dlgtitle,dims);
    
    if isempty(tle_num_answer) 
        disp('User selected Cancel');
        return
    end
    
    number_of_tle = str2double(tle_num_answer);
    
    if number_of_tle < 2
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        msgbox('A minimum of two TLE set are needed to compute visibility','Error',CreateStruct);
        return
    end
    
    prompt = sprintf('Enter %d sets of TLE without blank lines between set:', number_of_tle);
    dlgtitle = 'Paste TLE';
    dims = [3*number_of_tle 69];
    tle_pasted_answer = inputdlg(prompt,dlgtitle,dims);
    
    if isempty(tle_pasted_answer)
        disp('User selected Cancel');
        return
    end
    
    % TLE variables extraction

    % Find the total number of satellites in the file
    num_satellites = str2double(tle_num_answer);

    % Initialize array
    sat_id_line = zeros(1,num_satellites);
    line_count = 1;
    for i=1:num_satellites
        % Take every 3rd line
        sat_id_line(i) = line_count;
        txt_data = textscan(tle_pasted_answer{1}(line_count,1:69),'%s %s %s %s %s %s %s %s %s');
         
        OrbitData.ID(i) = txt_data{1};
            if isempty(txt_data{2})
                OrbitData.designation(i) = {''};
            else
                OrbitData.designation(i) = txt_data{2};
            end
            
            if isempty(txt_data{3})
                OrbitData.PRN(i) = {''};
            else
                OrbitData.PRN(i) = txt_data{3};
            end
        
        % Jump to the next Satellite Name line
        line_count  = line_count + 3; 
    end

    % Find the two lines corresponding to the spacecraft in question
    for j=1:length(sat_id_line)
        
        % Find the first line of the first satellite
        index = sat_id_line(j);
        txt_data_second = textscan(tle_pasted_answer{1}(index+2,1:69),'%s %s %s %s %s %s %s %s %s');
        
        % Translate two line element data into obital elements
        OrbitData.i(j)     = str2double(txt_data_second{1,3})*(pi/180);     % Inclination [deg] to [rad/s]
        OrbitData.RAAN(j)  = str2double(txt_data_second{1,4})*(pi/180);     % Right ascention of the ascending node[deg] to [rad/s]
        OrbitData.omega(j) = str2double(txt_data_second{1,6})*(pi/180);     % Argument of the periapsis [deg] to [rad/s]
        OrbitData.M(j)     = str2double(txt_data_second{1,7})*(pi/180);     % Mean anomaly [deg] to [rad/s]
        n                  = str2double(txt_data_second{1,8});              % Unperturbed mean motion [rev/day]
        OrbitData.n(j)     = n*2*pi/24/60/60;                               % Unperturbed mean motion [rad/s]
        OrbitData.a(j)     = ( mu / OrbitData.n(j)^2 )^(1/3);               % Semi-major axis [m]
        OrbitData.e(j)     = str2double(txt_data_second{1,5})*1e-7;         % Eccentricity [unitless]

        % Compute the UTC date / time
        txt_data_first = textscan(tle_pasted_answer{1}(index+1,1:69),'%s %s %s %s %s %s %s %s %s');
        temp2             = txt_data_first{1,4};
        yy                = str2double(temp2{1}(1:2));
        yyyy              = 2000 + yy;
        start             = datenum([yyyy 0 0 0 0 0]);
        secs              = str2double(temp2{1}(3:length(temp2{1})))*24*3600;
        date1             = datevec(addtodate(start,floor(secs),'second'));
        remainder         = [0 0 0 0 0 mod(secs,1)];
        OrbitData.date{j} = datestr(date1+remainder,'dd-mmm-yyyy HH:MM:SS.FFF');

        % Compute ballistic coefficient in SI units
        temp3 = txt_data_first{1,7};
        if length(temp3{1}) == 7
            base  = str2double(temp3{1}(1:5));
            expo  = str2double(temp3{1}(6:7));
        elseif length(temp3{1}) == 8
            base  = str2double(temp3{1}(2:6));
            expo  = str2double(temp3{1}(7:8));
        else
            fprintf('Error in ballistic coefficient calculation\n')
            CreateStruct.Interpreter = 'tex';
            CreateStruct.WindowStyle = 'modal';
            msgbox('Error in ballistic coefficient calculation\n','Error',CreateStruct);            
            error('End program')
        end
        
        Bstar = base*10^expo;
        OrbitData.BC(j) = 1/12.741621/Bstar; % [kg/m^2]
        
    end
    
end

% Simulation Parameters menu

input_simulation_list = {'From Now to Tomorrow (24h simulation)', 'Other'};
[indx,tf] = listdlg('ListString',input_simulation_list,'Name','Simulation Time','PromptString','Select a time for your analysis:',...
                    'SelectionMode','single','ListSize',[500,300],'OKString','Next','CancelString','Quit');

if tf == 0
    disp('User selected Quit');
    return
end

if indx == 1
    % Simulation Parameters
    start_time = '22-May-2008 12:00:00';
    %start_time = datetime('now', 'TimeZone', 'UTC');
    
    start_time_unix = posixtime(datetime(start_time));
    fprintf('Conversion of the simulation start time: %s is %d in Unix time\n', start_time, start_time_unix);                       % Command window print
    start_time_to_log = sprintf('Conversion of the simulation start time: %s is %d in Unix time', start_time, start_time_unix);
    t = start_time_unix;                                                                                                            % Start simulation time in Unix time [s]
    
    end_time = '23-May-2008 00:00:00';
    %end_time = datetime('now', 'TimeZone', 'UTC')+days(1);
    
    end_time_unix = posixtime(datetime(end_time));
    fprintf('Conversion of the simulation end time: %s is %d in Unix time\n', end_time, end_time_unix);                             % Command window print
    end_time_to_log = sprintf('Conversion of the simulation end time: %s is %d in Unix time', end_time, end_time_unix);
    t_end = end_time_unix;                                                                                                          % End of simulation time in Unix time [s]
    
else
    prompt = {'Simulation start:', 'Simulation end:'};
    dlgtitle = 'Simulation Time. Example: 22-Jan-2019 13:22:22';
    dims = [1 70; 1 70];
    simulation_answer = inputdlg(prompt,dlgtitle,dims);
    start_time = simulation_answer{1};
    end_time = simulation_answer{2};
  
    start_time_unix = posixtime(datetime(start_time));
    fprintf('Conversion of the simulation start time: %s is %d in Unix time\n', start_time, start_time_unix);                       % Command window print
    start_time_to_log = sprintf('Conversion of the simulation start time: %s is %d in Unix time', start_time, start_time_unix);
    fprintf(fid_log, '%s: %s\n\n', datestr(datetime('now', 'TimeZone', 'UTC')), start_time_to_log);                                                             % Appending simulation end time to log file
    t = start_time_unix;                                                                                                            % Start simulation time in Unix time [s]
    end_time_unix = posixtime(datetime(end_time));
    fprintf('Conversion of the simulation end time: %s is %d in Unix time\n', end_time, end_time_unix);                             % Command window print
    end_time_to_log = sprintf('Conversion of the simulation end time: %s is %d in Unix time', end_time, end_time_unix);
    t_end = end_time_unix;                                                                                                          % End of simulation time in Unix time [s]
end

time_divisons = 500; %4320
increment = (end_time_unix-start_time_unix)/time_divisons;                                                                          % Time increment [s]
num_steps = round(((end_time_unix-start_time_unix)/increment)+1);                                                                   % Number of time steps

% Satellite orbit parameters

for i=1:num_satellites
    OrbitData.epoch(i) = posixtime(datetime(char(OrbitData.date(i))));
    OrbitData.T(i) = OrbitData.epoch(i)-OrbitData.M(i)/OrbitData.n(i);
    OrbitData.q(i) = OrbitData.a(i)*(1-OrbitData.e(i));
end

num_pairs=0;
for sat1=1:num_satellites-1
    for sat2=sat1+1:num_satellites    
        num_pairs = num_pairs +1; 
    end
end

% Preallocated variables
for i=1:num_satellites
    n = zeros(1, num_satellites);                                           % Unperturbed mean motion [rev/day]
    M = zeros(1, num_satellites);                                           % Mean anomaly [degrees]
    Fn = zeros(1, num_satellites);                                          % Eccentric anomaly from Kepler's Equation for hyperbolic orbit (n) [degrees or rad]
    Fn1 = zeros(1, num_satellites);                                         % Eccentric anomaly from Kepler's Equation for hyperbolic orbit (n+1) [degrees or rad] 
    f = zeros(1, num_satellites);                                           % True Anomaly [degrees]
    A = zeros(1, num_satellites);                                           % Barker's Equation parameter [degrees]
    B = zeros(1, num_satellites);                                           % Barker's Equation parameter
    C = zeros(1, num_satellites);                                           % Barker's Equation parameter
    En = zeros(1, num_satellites);                                          % Eccentric anomaly from Kepler's Equation (n) [degrees or rad]
    En1 = zeros(1, num_satellites);                                         % Eccentric anomaly from Kepler's Equation (n+1) [degrees or rad]
    Px = zeros(1, num_satellites);                                          % First component of the unit orientation vector (dynamical center-periapsis) [m] 
    Py = zeros(1, num_satellites);                                          % Second component of the unit orientation vector (dynamical center-periapsis) [m] 
    Pz = zeros(1, num_satellites);                                          % Third component of the unit orientation vector (dynamical center-periapsis) [m] 
    Qx = zeros(1, num_satellites);                                          % First component of the unit orientation vector (advanced to P by a right angle in the motion direction) [m] 
    Qy = zeros(1, num_satellites);                                          % Second component of the unit orientation vector (advanced to P by a right angle in the motion direction) [m] 
    Qz = zeros(1, num_satellites);                                          % Third component of the unit orientation vector (advanced to P by a right angle in the motion direction) [m] 
    r = zeros(1, num_satellites);                                           % Magnitude of the vector from the center of the primary body to the satellite [m]
    xi = zeros(1, num_satellites);                                          % Component of r_vector in the periapsis line [m]
    eta = zeros(1, num_satellites);                                         % Component of r_vector in the descending node line [m]
    r_fullvector = [0 0 0];                                                 % Intermediate vector to store the pair of r_vectors [m]
    r_vector = zeros(num_satellites, 3);                                    % Vector from the center of the primary body to the satellite [m]
    parameter = zeros(1, num_satellites);                                   % Semi-parameter of the orbit [m]
    Rsimple1 = 0;                                                           % Visibility parameter [m]
    Rsimple2 = 0;                                                           % Visibility parameter [m]
    Rcomplex = zeros(num_steps, num_pairs);                                 % Visibility parameter [m]
    Rangle = 0;                                                             % Visibility parameter [m]
    Rv = 0;                                                                 % Distance from earth to satellite-satellite line
end

%% Log file module

prompt = 'Name this analysis: Log file format will be yyyymmddHHMMSS-Name.txt';
dlgtitle = 'Log file name';
dims = [1 69];
name_log = inputdlg(prompt,dlgtitle,dims);

if isempty(name_log)
    disp('User selected Cancel');
    return
end

date_log = datestr(now,'yyyymmddHHMMSS');
full_name_log = sprintf('%s-%s.txt',date_log,name_log{1});

fopen(fullfile([pwd, '/logs'], full_name_log), 'w'); % Create log file overwriting old one
fid_log = fopen(fullfile([pwd, '/logs'], full_name_log), 'a'); % Setting log file to append mode
if fid_log == -1
  error('Cannot open log file.');
end

fprintf(fid_log, '%s: %s\n\n', datestr(datetime('now', 'TimeZone', 'UTC')), start_time_to_log);         % Appending simulation start time to log file
fprintf(fid_log, '%s: %s\n\n', datestr(datetime('now', 'TimeZone', 'UTC')), end_time_to_log);           % Appending simulation end time to log file
disp('TLE data collected:');
disp(OrbitData);                                                            % Print TLE parameters in command window

% Log file is closed with "fclose" function once the algorithm is ended

%% 3D Visuals module

% Add path to the Earth plotting function and TLE data preparation
addpath([pwd, '\Plot Earth']);
addpath([pwd, '\TLE Plotter']);

% Simulation Start Date
simStart = start_time;

% Compute sidereal time
GMST = utc2gmst(datevec(simStart)); % [rad]

% Create a time vector
tSim = linspace(start_time_unix, end_time_unix, num_steps);

% Allocate space
RSave = NaN(length(tSim), 3, num_satellites);

%% Algorithm

tic; % Runtime start

num_pairs = 0;

for sat1=1:num_satellites-1
    
    for sat2=sat1+1:num_satellites
        
        num_pairs= num_pairs + 1;

        step_count=1;

        for t=t:increment:t_end % Simulation time and time discretization

            for i=sat1:sat2

                % Step 1 - Finding unperturbed mean motion
                if OrbitData.e(i) > 1
                    n(i) = k*sqrt(mu/-OrbitData.a(i)^3);
                elseif OrbitData.e(i) == 1
                    n(i) = k*sqrt(mu/(2*OrbitData.q(i)^3));
                elseif OrbitData.e(i) < 1 && OrbitData.e(i) >= 0
                    % Mean motion method 1
                    % n(i) = k*sqrt(mu/OrbitData.a(i)^3);
                    % Mean motion method 2
                    n(i) = OrbitData.n(i);
                else
                    error('Eccentricity cannot be a negative value')
                end

                % Step 2 - Solving Mean Anomaly
                % Mean anomaly method 1
                % M(i) = n(i)*(t-OrbitData.T(i));
                % Mean anomaly method 2
                M(i) = OrbitData.M(i) + n(i)*(t-start_time_unix);

                % Step 3 - Finding true anomaly 
                if OrbitData.e(i) > 1
                    % Iteration method 1
                    Fn(i) = 6*M(i);
                    error = 1;
                    while error > 1e-8
                        Fn1(i) = Fn(i) + (M(i)-OrbitData.e(i)*sinh(Fn(i))+Fn(i))/(OrbitData.e(i)*cosh(Fn(i))-1);
                        error = abs(Fn1(i)-Fn(i));
                        Fn(i) = Fn1(i);
                    end

                    % Iteration method 2 
                    % TODO

                    f(i) = atan((-sinh(Fn(i))*sqrt(OrbitData.e(i)^2-1))/(cosh(Fn(i))-OrbitData.e(i)));

                elseif OrbitData.e(i) == 1
                    A(i) = (3/2)*M(i);
                    B(i) = (sqrt(A(i)^2+1)+A(i))^(1/3);
                    C(i) = B(i)-1/B(i);
                    f(i) = 2*atan(C(i));

                elseif OrbitData.e(i) < 1 && OrbitData.e(i) >= 0
                    % Iteration method 1
        %             En(i) = M(i);
        %             error = 1;
        %             while error > 1e-8
        %                 En1(i) = En(i) + (M(i)-OrbitData.e(i)*sin(En(i))-En(i))/(1-OrbitData.e(i)*cos(En(i)));
        %                 error = abs(En1(i)-En(i));
        %                 En(i) = En1(i);
        %             end

                    % Iteration method 2
                    if M(i) < pi % careful with negatives
                        Einicial = M(i) + OrbitData.e(i)/2;
                    else 
                        Einicial = M(i) - OrbitData.e(i)/2;
                    end
                    E = Einicial;
                    TOL = 1;
                    while TOL > 1e-8
                        fdee = E - OrbitData.e(i)*sin(E)-M(i);
                        fprimadee = 1-OrbitData.e(i)*cos(E);
                        TOL = abs(fdee/fprimadee);
                        En(i)=E;
                        E=E-fdee/fprimadee;
                    end

                    f(i) = atan((sin(En(i))*sqrt(1-OrbitData.e(i)^2))/(cos(En(i))-OrbitData.e(i))); % TODO

                else
                    error('Eccentricity cannot be a negative value')
                end

                % Step 4 - Finding primary body center to satellite distance
                r(i) = (1+OrbitData.e(i))*OrbitData.q(i)/(1+OrbitData.e(i)*cos(f(i)));

                % Step 5 - Finding standard orientation vectors
                Px(i) = cos(OrbitData.omega(i))*cos(OrbitData.RAAN(i))-sin(OrbitData.omega(i))*sin(OrbitData.RAAN(i))*cos(OrbitData.i(i));
                Py(i) = cos(OrbitData.omega(i))*sin(OrbitData.RAAN(i))+sin(OrbitData.omega(i))*cos(OrbitData.RAAN(i))*cos(OrbitData.i(i));
                Pz(i) = sin(OrbitData.omega(i))*sin(OrbitData.i(i));
                Qx(i) = -sin(OrbitData.omega(i))*cos(OrbitData.RAAN(i))+cos(OrbitData.omega(i))*sin(OrbitData.RAAN(i))*cos(OrbitData.i(i));
                Qy(i) = -sin(OrbitData.omega(i))*sin(OrbitData.RAAN(i))+cos(OrbitData.omega(i))*cos(OrbitData.RAAN(i))*cos(OrbitData.i(i));
                Qz(i) = cos(OrbitData.omega(i))*sin(OrbitData.i(i));

                % Step 6 - Finding components of the primary body center to satellite vector in the orbital plane
                xi(i) = r(i)*cos(f(i));
                eta(i) = r(i)*sin(f(i));

                % Step 7 - Finding primary body center to satellite vector
                r_fullvector = xi(i)*[Px(i) Py(i) Pz(i)] + eta(i)*[Qx(i) Qy(i) Qz(i)];
                for j=1:3
                    r_vector(i,j) = r_fullvector(j);
                end

                % Step 8 - Finding Parameter or Semi-parameter
                parameter(i) = OrbitData.a(i)*(1-OrbitData.e(i)^2);

                % Step 9 - Transformation for 3D visuals
                % Adjust RAAN such that we are consisten with Earth's current
                % orientation. This is a conversion to Longitude of the
                % Ascending Node (LAN). 
                RAAN2 = OrbitData.RAAN(i) - GMST;

                % Convert to ECI and save the data.
                [X,~] = COE2RV(OrbitData.a(i), OrbitData.e(i), OrbitData.i(i), RAAN2, OrbitData.omega(i), M(i));
                RSave(step_count,:,i) = X';  
            end

            % Step 10 - Solving visibility equation

            P1 = [Px(sat1) Py(sat1) Pz(sat1)];
            P2 = [Px(sat2) Py(sat2) Pz(sat2)];
            Q1 = [Qx(sat1) Qy(sat1) Qz(sat1)];
            Q2 = [Qx(sat2) Qy(sat2) Qz(sat2)];
            A1 = dot(P1,P2);
            A2 = dot(Q1,P2);
            A3 = dot(P1,Q2);
            A4 = dot(Q1,Q2);

            gamma1 = asin(A2/sqrt(A1^2+A2^2));
            psi1 = asin(A4/sqrt(A3^2+A4^2));

            gamma2 = acos(A1/sqrt(A1^2+A2^2));
            psi2 = acos(A3/sqrt(A3^2+A4^2));

            sin_gamma = A2/sqrt(A1^2+A2^2);
            cos_gamma = A1/sqrt(A1^2+A2^2);
            sin_psi = A4/sqrt(A3^2+A4^2);
            cos_psi = A3/sqrt(A3^2+A4^2);

            D1 = sqrt(A1^2+A2^2);
            D2 = sqrt(A3^2+A4^2);

            % r1dotr2calc = r_vector(1,1)*r_vector(2,1) + r_vector(1,2)*r_vector(2,2) + r_vector(1,3)*r_vector(2,3) 
            % r1dotr2simple = (parameter(1)*parameter(2)/((1+OrbitData.e(1)*cos(f(1)))*(1+OrbitData.e(2)*cos(f(2)))))*(A1*cos(f(1))*cos(f(2))+A3*cos(f(1))*sin(f(2))+A2*sin(f(1))*cos(f(2))+A4*sin(f(1))*sin(f(2)));
            r1dotr2complex = (parameter(sat1)*parameter(sat2)/((1+OrbitData.e(sat1)*cos(f(sat1)))*(1+OrbitData.e(sat2)*cos(f(sat2)))))* ...
                             (D1*cos(f(sat2))*(cos_gamma*cos(f(sat1))+sin_gamma*sin(f(sat1)))+D2*sin(f(sat2))*(cos_psi*cos(f(sat1))+sin_psi*sin(f(sat1))));
            % Rsimple1 = r1dotr2simple^2 - r(2)^2*r(1)^2 + (r(2)^2 + r(1)^2)*S^2 - 2*S^2*(r1dotr2simple);
            % Rsimple2 = r1dotr2complex^2 - r(2)^2*r(1)^2 + (r(2)^2 + r(1)^2)*S^2 - 2*S^2*(r1dotr2complex);
            Rcomplex(step_count, num_pairs) = parameter(sat1)^2 * parameter(sat2)^2 * ( D1*cos(f(sat2))*(cos_gamma*cos(f(sat1))+sin_gamma*sin(f(sat1))) + ...
                D2*sin(f(sat2))*(cos_psi*cos(f(sat1))+sin_psi*sin(f(sat1))) )^2 - parameter(sat1)^2*parameter(sat2)^2 + S^2*( parameter(sat1)^2* ...
                (1+OrbitData.e(sat2)*cos(f(sat2)))^2 + parameter(sat2)^2*(1+OrbitData.e(sat1)*cos(f(sat1)))^2 ) - 2*S^2*parameter(sat1)*parameter(sat2)* ...
                ( D1*cos(f(sat2))* ( cos_gamma*cos(f(sat1))+sin_gamma*sin(f(sat1)) ) + D2*sin(f(sat2))* ( cos_psi*cos(f(sat1))+sin_psi*sin(f(sat1)) ) ) * ...
                (1+OrbitData.e(sat1)*cos(f(sat1))) * (1+OrbitData.e(sat2)*cos(f(sat2)));
            % Rangle = parameter(1)^2*parameter(2)^2*(D1*cos(f(2))*cos(gamma-f(1))+D2*sin(f(2))*cos(psi-f(1)))^2-parameter(1)^2*parameter(2)^2+S^2*(parameter(1)^2*(1+OrbitData.e(2)*cos(f(2)))^2+parameter(2)^2*(1+OrbitData.e(1)*cos(f(1)))^2)-2*S^2*parameter(1)*parameter(2)*(D1*cos(f(2))*cos(gamma-f(1))+D2*sin(f(2))*cos(psi-f(1)))*(1+OrbitData.e(1)*cos(f(1)))*(1+OrbitData.e(2)*cos(f(2)));
            Rv = sqrt((r(sat1)^2 * r(sat2)^2 - r1dotr2complex^2)/(r(sat1)^2 + r(sat2)^2 - 2*r1dotr2complex)) - body_radius;
            % Rv_Tot = (r(1)^2*r(2)^2-r1dotr2complex^2)/(r(1)^2 + r(2)^2-2*r1dotr2complex)
            % Rv_Numerador = r(1)^2 * r(2)^2 - r1dotr2complex^2
            % Rv_Denominador = r(1)^2 + r(2)^2 - 2*r1dotr2complex

            % Step 9: Print Results for the given epoch time 
            pair_result = 'The result for %s%s and %s%s at %s is %d ';
            visibility = '--- Direct line of sight';
            non_visibility= '--- Non-visibility';

            result_to_log = sprintf(pair_result, OrbitData.ID{sat1}, OrbitData.designation{sat1}, OrbitData.ID{sat2}, OrbitData.designation{sat2}, datetime(t, 'ConvertFrom', 'posixtime'), Rcomplex(step_count, num_pairs));
            fprintf(result_to_log); % Command window print

            if Rcomplex(step_count, num_pairs) < 0
                disp(visibility); % Command window print
                fprintf(fid_log, '%s: %s%s\n\n', datestr(datetime('now', 'TimeZone', 'UTC')), result_to_log, visibility); % Appending visibility analysis result to log file
            else
                disp(non_visibility); % Command window print
                fprintf(fid_log, '%s: %s%s\n\n', datestr(datetime('now', 'TimeZone', 'UTC')), result_to_log, non_visibility); % Appending visibility analysis result to log file
            end

            step_count = step_count + 1;

        end

        t = start_time_unix;
    
    end

end

toc; % Runtime end

%% Plot the orbit

% Static plot
% Plot the Earth
% If you want a color Earth, use 'neomap', 'BlueMarble'
% If you want a black and white Earth, use 'neomap', 'BlueMarble_bw'
% A smaller sample step gives a finer resolution Earth
disp('Opening plot module...')

% Simualtion Unix time vector converted to DateTimes strings inside a cell
tSim_strings = {step_count-1};
for t=1:step_count-1
    tSim_strings{t} = datestr(datetime(tSim(t),'ConvertFrom','posixtime'));
end

plot_list = {'Static Plot', 'Live Plots (See color legend in 1 vs 1 plot to identify satellites. It may take a lot of time)'};
[indx,tf] = listdlg('ListString',plot_list,'Name','3D Plot','PromptString','Select a plot mode:','SelectionMode','single','ListSize',[600,300],'OKString','Plot','CancelString','Quit');

if tf == 0
    disp('User selected Quit');
    return
end

colors = lines(num_satellites);

if indx == 2
    h1 = plotearth('neomap', 'BlueMarble_bw', 'SampleStep', 1);
    for i=1:num_satellites
        plot3(RSave(:,1,i) / body_radius, RSave(:,2,i) / body_radius, RSave(:,3,i) / body_radius,...
              'color', colors(i,:), 'LineWidth', 1, 'DisplayName', strcat(OrbitData.ID{i}, OrbitData.designation{i}, ' - Orbit'))
        plot3(RSave(1,1,i) / body_radius, RSave(1,2,i) / body_radius, RSave(1,3,i) / body_radius,...
              '.', 'color', colors(i,:), 'MarkerSize', 10, 'DisplayName', strcat(OrbitData.ID{i}, OrbitData.designation{i}, ' - Starting Point'))
    end
    lgd2 = legend('AutoUpdate', 'off');

    % Live 3D plot
    num_pairs = 0;
    h2 = plotearth('neomap', 'BlueMarble_bw', 'SampleStep', 1);
    for sat1=1:num_satellites-1
        
        for sat2=sat1+1:num_satellites
            num_pairs = num_pairs + 1;
            hold on
            
            for t=1:step_count-1
 
                    for i=sat1:sat2
                        if Rcomplex(t, num_pairs) < 0 && i == sat1
                            curve = animatedline('LineWidth',2,'color', [100, 255, 110] / 255, 'DisplayName', 'Visibility', 'HandleVisibility', 'on'); % Green color
                        elseif Rcomplex(t, num_pairs) >= 0 && i == sat1
                            curve = animatedline('LineWidth',2,'color', [225, 90, 90] / 255, 'DisplayName', 'Non-visibility', 'HandleVisibility', 'on'); % Red color
                        elseif Rcomplex(t, num_pairs) < 0 && i == sat2
                            curve = animatedline('LineWidth',2,'color', [100, 255, 110] / 255, 'DisplayName', 'Visibility', 'HandleVisibility', 'off'); % Green color
                        elseif Rcomplex(t, num_pairs) >= 0 && i == sat2
                            curve = animatedline('LineWidth',2,'color', [225, 90, 90] / 255, 'DisplayName', 'Non-visibility', 'HandleVisibility', 'off'); % Red color
                        end
                        if i == sat1
                            head1 = scatter3(RSave(t,1,sat1) / body_radius, RSave(t,2,sat1) / body_radius, RSave(t,3,sat1) / body_radius, 'filled', 'MarkerFaceColor', colors(sat1,:),... 
                                        'DisplayName', strcat(OrbitData.ID{sat1}, OrbitData.designation{sat1}), 'HandleVisibility', 'on');
                        elseif i == sat2
                            head2 = scatter3(RSave(t,1,sat2) / body_radius, RSave(t,2,sat2) / body_radius, RSave(t,3,sat2) / body_radius, 'filled', 'MarkerFaceColor', colors(sat2,:),... 
                                    'DisplayName', strcat(OrbitData.ID{sat2}, OrbitData.designation{sat2}), 'HandleVisibility', 'on');
                        end
                        addpoints(curve, RSave(1:t,1,i) / body_radius, RSave(1:t,2,i) / body_radius, RSave(1:t,3,i) / body_radius);
                        drawnow;
                        pause(0.01)
                    end
                    
                    lgd = legend(tSim_strings{t});
                    lgd.FontSize = 15;
                    delete(head1);
                    delete(head2);
            end

        end

    end
    
else
    % Static plot
    h1 = plotearth('neomap', 'BlueMarble_bw', 'SampleStep', 1);
    for i=1:num_satellites
        plot3(RSave(:,1,i) / body_radius, RSave(:,2,i) / body_radius, RSave(:,3,i) / body_radius,...
              'color', colors(i,:), 'LineWidth', 1, 'DisplayName', strcat(OrbitData.ID{i}, OrbitData.designation{i}, ' - Orbit'))
        plot3(RSave(1,1,i) / body_radius, RSave(1,2,i) / body_radius, RSave(1,3,i) / body_radius,...
              '.', 'color', colors(i,:), 'MarkerSize', 10, 'DisplayName', strcat(OrbitData.ID{i}, OrbitData.designation{i}, ' - Starting Point'))
    end
    lgd2 = legend();
end

disp('Program ended successfully')

%% CSV output file module

fid_csv = fopen(fullfile([pwd, '\Data Output File'],'InterlinkData.csv'), 'a');
toadd = (1:4);
dlmwrite(fullfile([pwd, '\Data Output File'],'InterlinkData.csv'),toadd,'-append','delimiter',';');
fclose(fid_csv);
