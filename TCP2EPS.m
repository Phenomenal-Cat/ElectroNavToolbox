function EPS = ENT_TCP2EPS

%============================= ENT_TCP2EPS.m ==============================
% Communicate with AlphaOmega Electrode Positioning System (EPS) software
% via a TCP ethernet connection.
%
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy © Copyleft 2015, GNU General Public License
%==========================================================================

EPS.Port    = 8120;                                  	% Set port number of EPS
EPS.IP      = '156.40.249.83';                       	% Set IP address of EPS PC
readTimeOut = 2;                                       	% reading times out (seconds)

con         = pnet('tcpconnect',EPS.IP,EPS.Port);     	% Open TCP/IP connection
sockcon     = pnet('tcpsocket',EPS.Port);            	% Open socket
pnet(sockcon,'setreadtimeout',readTimeOut);             % Set a time out limit
a           = pnet(sockcon,'tcplisten');                % Check that EPS has connected to socket
% if a =-1
    str1     = pnet(con,'read',100,[],[],[],'noblock') 	% Read data sent from EPS
% end

pnet(con,'printf', 'a=PEER_INFO&TCP_APP_NAME=EPSClt&TCP_SERVICE_TYPE=CONTROLEPS&TCP_CONN_BLOCKED=-1');

str2     = pnet(con,'read')
if ~strcmpi(str2, 'STATUS=0')
    
end

% elements    = pnet(con,'write', data [,swapping])

stat        = pnet(con,'status');
[ip,port]   = pnet(con,'gethost');
pnet(con,'close')
pnet('closeall')

% %============= 
% % readTimeOut = 0.001;
% readTimeOut = 2;                                 	% how long to wait for the socket info
% sockcon=pnet('tcpsocket',4610);                     
% pnet(sockcon,'setreadtimeout',readTimeOut);
% con=pnet(sockcon,'tcplisten');
% if con~=-1
%     pnet(con,'setreadtimeout',readTimeOut);
% end
% pnet(con,'printf', reply);


% EPS.obj     = tcpip(EPS.IP,'InputBufferSize',1024);     % Create TCP object
% fopen(EPS.obj);                                         % Open connection
% x = fread(EPS.obj, 10)                              % Read data
% 
% fclose(EPS.obj);                                    % Close connection
% echotcpip('off');                                   