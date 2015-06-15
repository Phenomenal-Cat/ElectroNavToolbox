function EPS = TCP2EPS

%============================= EN_TCP2EPS.m ===============================
% Communicate with AlphaOmega Electrode Positioning System (EPS) software
% via a TCP ethernet connection.
%
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy © Copyleft 2015, GNU General Public License
%==========================================================================


EPS.IP = '156.40.249.83';                           % Set IP address of EPS PC
EPS.obj = tcpip(EPS.IP,'InputBufferSize',1024);     % Create TCP object
fopen(EPS.obj);                                     % Open connection
x = fread(EPS.obj, 10)                              % Read data

fclose(EPS.obj);                                    % Close connection
echotcpip('off');                                   


%============= 
% readTimeOut = 0.001;
readTimeOut = 2;                                 	% how long to wait for the socket info
sockcon=pnet('tcpsocket',4610);                     
pnet(sockcon,'setreadtimeout',readTimeOut);
con=pnet(sockcon,'tcplisten');
if con~=-1
    pnet(con,'setreadtimeout',readTimeOut);
end
pnet(con,'printf', reply);