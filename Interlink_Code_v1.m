tic; % Runtime start

format long
close all
clc

fopen(fullfile('C:\Users\david\Desktop\Uni\TFG\Matlab David\logs', 'log_file.txt'), 'w'); % Create log file overwriting old one
fid = fopen(fullfile('C:\Users\david\Desktop\Uni\TFG\Matlab David\logs', 'log_file.txt'), 'a'); % Setting log file to append mode
if fid == -1
  error('Cannot open log file.');
end

% STUDY OF THE INTERLINK BETWEEN SMALL SATELLITES IN A CONSTELLATION
% Author: David Puig Puig

% Visual contact for two satellites analysis

% TLE of a pair of Satellites

% EGYPTSAT 1
% 1 31117U 07012A 08142.74302347 .00000033 00000-0 13654-4 0 2585
% 2 31117 098.0526 218.7638 0007144 061.2019 298.9894 14.69887657 58828

% TRMM
% 1 25063U 97074A 08141.84184490 .00002948 00000-0 41919-4 0 7792
% 2 25063 034.9668 053.5865 0001034 271.1427 088.9226 15.55875272598945

%% Input parameters

% Simulation Parameters
start_time = '22-May-2008 12:00:00';
start_time_unix = posixtime(datetime(start_time));
fprintf('Conversion of the simulation start time: %s is %d in Unix time\n', start_time, start_time_unix); % Command window print
start_time_to_log = sprintf('Conversion of the simulation start time: %s is %d in Unix time', start_time, start_time_unix);
fprintf(fid, '%s: %s\n\n', datestr(now, 0), start_time_to_log); % Appending simulation end time to log file
t = start_time_unix; % Start simulation time in Unix time [s]
end_time = '23-May-2008 00:00:00';
end_time_unix = posixtime(datetime(end_time));
fprintf('Conversion of the simulation end time: %s is %d in Unix time\n', end_time, end_time_unix); % Command window print
end_time_to_log = sprintf('Conversion of the simulation end time: %s is %d in Unix time', end_time, end_time_unix);
fprintf(fid, '%s: %s\n\n', datestr(now, 0), end_time_to_log); % Appending simulation end time to log file
t_end = end_time_unix; % End of simulation time in Unix time [s]
increment = 10; % Time increment [s]

% System parameters (Earth)
body_radius = 6.378e6; % Radius of the primary body [m]
extra_radius = 10000; % Extra radius for the primary body [m]
S = body_radius + extra_radius; % Magnitude of the rise-set vector [m]
k = 2*pi; % Factor from [rev/s] to [rad/s]
mu = 3.986004418e14; % Standard gravitational parameter [m^3/s^2]

% Satellite orbit parameters
inclination = [98.0526*pi/180 34.9668*pi/180]; % Inclination [degrees] converted to [rad]
argument_periapsis = [61.2019*pi/180 271.1427*pi/180]; % Argument of the periapsis [degrees] converted to [rad]
longitude_ascending_node = [218.7638*pi/180 53.5865*pi/180]; % Longitude of the ascending node [degrees] converted to [rad]
mean_anomaly_tle = [298.9894*pi/180 88.9226*pi/180]; % Mean anomaly extracted from TLE [degrees] converted to [rad]
mean_motion_tle = [14.69887657*2*pi/86400 15.558752725*2*pi/86400]; % Unperturbed mean motion extracted from TLE [rev/day] converted to [rad/s]
epoch_year_tle = [2008 2008]; % Epoch year extracted from TLE [s]
epoch_day_tle = [142.74302347 141.84184490]; % Epoch day extracted from TLE [s] 
hours1 = 0.74302347*24;
minutes1 = abs(hours1-fix(hours1))*60;
seconds1 = abs(minutes1-fix(minutes1))*60;
hours2 = 0.84184490*24;
minutes2 = abs(hours2-fix(hours2))*60;
seconds2 = abs(minutes2-fix(minutes2))*60;
epoch = [posixtime(datetime('22-May-2008 17:49:57.2278')) posixtime(datetime('21-May-2008 20:12:15.3994'))]; % Epoch from Unix time [s] 
T = [epoch(1)-mean_anomaly_tle(1)/mean_motion_tle(1) epoch(2)-mean_anomaly_tle(2)/mean_motion_tle(2)]; % Time of perifocal passage [s]
semimajor_axis = [mu^(1/3)/(mean_motion_tle(1))^(2/3) mu^(1/3)/(mean_motion_tle(2))^(2/3)]; % Semi-major axis [m]
eccentricity = [0.0007144 0.0001034]; % Eccentricity [dimensionless]
periapsis_distance = [semimajor_axis(1)*(1-eccentricity(1)) semimajor_axis(2)*(1-eccentricity(2))]; % Periapsis Distance [m]

% Preallocated variables
n = [0 0]; % Unperturbed mean motion [rev/day]
M = [0 0]; % Mean anomaly [degrees]
Fn = [0 0]; % Eccentric anomaly from Kepler's Equation for hyperbolic orbit (n) [degrees or rad]
Fn1 = [0 0]; % Eccentric anomaly from Kepler's Equation for hyperbolic orbit (n+1) [degrees or rad] 
f = [0 0]; % True Anomaly [degrees]
A = [0 0]; % Barker's Equation parameter [degrees]
B = [0 0]; % Barker's Equation parameter
C = [0 0]; % Barker's Equation parameter
En = [0 0]; % Eccentric anomaly from Kepler's Equation (n) [degrees or rad]
En1 = [0 0]; % Eccentric anomaly from Kepler's Equation (n+1) [degrees or rad]
Px = [0 0]; % First component of the unit orientation vector (dynamical center-periapsis) [m] 
Py = [0 0]; % Second component of the unit orientation vector (dynamical center-periapsis) [m] 
Pz = [0 0]; % Third component of the unit orientation vector (dynamical center-periapsis) [m] 
Qx = [0 0]; % First component of the unit orientation vector (advanced to P by a right angle in the motion direction) [m] 
Qy = [0 0]; % Second component of the unit orientation vector (advanced to P by a right angle in the motion direction) [m] 
Qz = [0 0]; % Third component of the unit orientation vector (advanced to P by a right angle in the motion direction) [m] 
r = [0 0]; % Magnitude of the vector from the center of the primary body to the satellite [m]
xi = [0 0]; % Component of r_vector in the periapsis line [m]
eta = [0 0]; % Component of r_vector in the descending node line [m]
r_fullvector = [0 0 0]; % Intermediate vector to store the pair of r_vectors [m]
r_vector = [0 0 0; 0 0 0]; % Vector from the center of the primary body to the satellite [m]
parameter = [0 0]; % Semi-parameter of the orbit [m]
R1 = 0; % Visibility parameter [m]
R2 = 0; % Visibility parameter [m]

%% Algorithm

for t=t:increment:t_end % Simulation time and time discretization

    for i=1:2

        % Step 1 - Finding unperturbed mean motion
        if eccentricity(i) > 1
            n(i) = k*sqrt(mu/-semimajor_axis(i)^3);
        elseif eccentricity(i) == 1
            n(i) = k*sqrt(mu/(2*periapsis_distance(i)^3));
        elseif eccentricity(i) < 1 && eccentricity(i) >= 0
            n(i) = k*sqrt(mu/semimajor_axis(i)^3);
            % n(i) = mean_motion_tle(i);
        else
            error('Eccentricity can''t be a negative value')
        end

        % Step 2 - Solving Mean Anomaly
        M(i) = n(i)*(t-T(i));

        % Step 3 - Finding true anomaly 
        if eccentricity(i) > 1
            Fn(i) = 6*M(i);
            error = 1;
            while error > 1e-8
                Fn1(i) = Fn(i) + (M(i)-eccentricity(i)*sinh(Fn(i))+Fn(i))/(eccentricity(i)*cosh(Fn(i))-1);
                error = abs(Fn1(i)-Fn(i));
                Fn(i) = Fn1(i);
            end
            f(i) = atan((-sinh(Fn(i))*sqrt(eccentricity(i)^2-1))/(cosh(Fn(i))-eccentricity(i)));
        elseif eccentricity(i) == 1
            A(i) = (3/2)*M(i);
            B(i) = (sqrt(A(i)^2+1)+A(i))^(1/3);
            C(i) = B(i)-1/B(i);
            f(i) = 2*atan(C(i));
        elseif eccentricity(i) < 1 && eccentricity(i) >= 0
            % En(i) = M(i);
            % error = 1;
            % while error > 1e-8
                % En1(i) = En(i) + (M(i)-eccentricity(i)*sin(En(i))-En(i))/(1-eccentricity(i)*cos(En(i)));
                % error = abs(En1(i)-En(i));
                % En(i) = En1(i);
                % end
            if M(i) < pi % careful with negatives
                Einicial = M(i) + eccentricity(i)/2;
            else 
                Einicial = M(i) - eccentricity(i)/2;
            end
            E = Einicial;
            TOL = 1;
            while TOL > 1e-8
                fdee = E - eccentricity(i)*sin(E)-M(i);
                fprimadee = 1-eccentricity(i)*cos(E);
                TOL = abs(fdee/fprimadee);
                En(i)=E;
                E=E-fdee/fprimadee;
            end
            f(i) = atan((sin(En(i))*sqrt(1-eccentricity(i)^2))/(cos(En(i))-eccentricity(i))); % TODO
        else
            error('Eccentricity can''t be a negative value')
        end

        % Step 4 - Finding primary body center to satellite distance
        r(i) = (1+eccentricity(i))*periapsis_distance(i)/(1+eccentricity(i)*cos(f(i)));

        % Step 5 - Finding standard orientation vectors
        Px(i) = cos(argument_periapsis(i))*cos(longitude_ascending_node(i))-sin(argument_periapsis(i))*sin(longitude_ascending_node(i))*cos(inclination(i));
        Py(i) = cos(argument_periapsis(i))*sin(longitude_ascending_node(i))+sin(argument_periapsis(i))*cos(longitude_ascending_node(i))*cos(inclination(i));
        Pz(i) = sin(argument_periapsis(i))*sin(inclination(i));
        Qx(i) = -sin(argument_periapsis(i))*cos(longitude_ascending_node(i))+cos(argument_periapsis(i))*sin(longitude_ascending_node(i))*cos(inclination(i));
        Qy(i) = -sin(argument_periapsis(i))*sin(longitude_ascending_node(i))+cos(argument_periapsis(i))*cos(longitude_ascending_node(i))*cos(inclination(i));
        Qz(i) = cos(argument_periapsis(i))*sin(inclination(i));

        % Step 6 - Finding components of the primary body center to satellite vector in the orbital plane
        xi(i) = r(i)*cos(f(i));
        eta(i) = r(i)*sin(f(i));

        % Step 7 - Finding primary body center to satellite vector
        r_fullvector = xi(i)*[Px(i) Py(i) Pz(i)] + eta(i)*[Qx(i) Qy(i) Qz(i)];
        for j=1:3
            r_vector(i,j) = r_fullvector(j);
        end
        
        % Step 8 - Finding Parameter or Semi-parameter
        parameter(i) = semimajor_axis(i)*(1-eccentricity(i)^2);
        
    end
    
    % Step 9 - Solving visibility equation

    P1 = [Px(1) Py(1) Pz(1)];
    P2 = [Px(2) Py(2) Pz(2)];
    Q1 = [Qx(1) Qy(1) Qz(1)];
    Q2 = [Qx(2) Qy(2) Qz(2)];
    A1 = dot(P1,P2);
    A2 = dot(Q1,P2);
    A3 = dot(P1,Q2);
    A4 = dot(Q1,Q2);

    gamma = asin(A2/sqrt(A1^2+A2^2));
    psi = asin(A4/sqrt(A3^2+A4^2));

    sin_gamma = A2/sqrt(A1^2+A2^2);
    cos_gamma = A1/sqrt(A1^2+A2^2);
    sin_psi = A4/sqrt(A3^2+A4^2);
    cos_psi = A3/sqrt(A3^2+A4^2);
    
    D1 = sqrt(A1^2+A2^2);
    D2 = sqrt(A3^2+A4^2);
    
    % cos(gamma-f(1)) == cos(gamma)*cos(f(1))+sin(gamma)*sin(f(1)) == ((A1/sqrt(A1^2+A2^2))*cos(f(1))+(A2/sqrt(A1^2+A2^2))*sin(f(1)))      
    % cos(psi-f(1))) == cos(psi)*cos(f(1))+sin(psi)*sin(f(1)) == ((A4/sqrt(A3^2+A4^2))*cos(f(1))+(A4/sqrt(A3^2+A4^2))*sin(f(1)))
    
    % R1 = parameter(1)^2*parameter(2)^2*(D1*cos(f(2))*cos(gamma-f(1))+D2*sin(f(2))*cos(psi-f(1)))^2-parameter(1)^2*parameter(2)^2+S^2*(parameter(1)^2*(1+eccentricity(2)*cos(f(2)))^2+parameter(2)^2*(1+eccentricity(1)*cos(f(1)))^2)-2*S^2*parameter(1)*parameter(2)*(D1*cos(f(2))*cos(gamma-f(1))+D2*sin(f(2))*cos(psi-f(1)))*(1+eccentricity(1)*cos(f(1)))*(1+eccentricity(2)*cos(f(2)));
    R2 = parameter(1)^2*parameter(2)^2*(D1*cos(f(2))*((A1/sqrt(A1^2+A2^2))*cos(f(1))+(A2/sqrt(A1^2+A2^2))*sin(f(1)))+D2*sin(f(2))*((A3/sqrt(A3^2+A4^2))*cos(f(1))+(A4/sqrt(A3^2+A4^2))*sin(f(1))))^2-parameter(1)^2*parameter(2)^2+S^2*(parameter(1)^2*(1+eccentricity(2)*cos(f(2)))^2+parameter(2)^2*(1+eccentricity(1)*cos(f(1)))^2)-2*S^2*parameter(1)*parameter(2)*(D1*cos(f(2))*((A1/sqrt(A1^2+A2^2))*cos(f(1))+(A2/sqrt(A1^2+A2^2))*sin(f(1)))+D2*sin(f(2))*((A3/sqrt(A3^2+A4^2))*cos(f(1))+(A4/sqrt(A3^2+A4^2))*sin(f(1))))*(1+eccentricity(1)*cos(f(1)))*(1+eccentricity(2)*cos(f(2)));
    
    % Step 9: Print Results for the given epoch time 
    pair_result = 'The result for the pair of satellites at %s is %d ';
    visibility = '--- Direct line of sight';
    non_visibility= '--- Non-visibility';

    result_to_log = sprintf(pair_result,datetime(t, 'ConvertFrom', 'posixtime'),R2);
    fprintf(result_to_log); % Command window print

    if R2 < 0
        disp(visibility); % Command window print
        fprintf(fid, '%s: %s%s\n\n', datestr(now, 0), result_to_log, visibility); % Appending visibility analysis result to log file
    else
        disp(non_visibility); % Command window print
        fprintf(fid, '%s: %s%s\n\n', datestr(now, 0), result_to_log, non_visibility); % Appending visibility analysis result to log file
    end

end

fclose(fid); % Closing log file

toc; % Runtime end