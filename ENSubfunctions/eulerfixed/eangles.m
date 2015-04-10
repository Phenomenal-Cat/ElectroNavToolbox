function eangles
	% EANGLES Graphical Display of Euler Angle Characteristics
	%
	% eangles() opens a graphical user interface which gives the user two
	% options. The first (forward) takes values for angles alpha, beta and gamma
	% (in radians) and gives the output of the 3x3 rotation matrix R.
	%
	% The second option (backward) takes the inputted 3x3 rotation matrix R and
	% gives the values (in radians) for the three euler angles alpha, beta and
	% gamma.
	%
	% In both cases, one of the twelve conventions for Euler Angles needs to be
	% given (for example zyz), and a small animation showing these rotations
	% (according to the specified convention) will be ran.
	%
	% Created on MATLAB R2007a (7.4)
	% Husam Aldahiyat, Jordan, 2009
	%
	
	% create main figure
	figure('units','normalized','position',[.2 .2 .7 .7],'color','w','menubar','non','numbertitle','off',...
		'name','Euler Angles')
	
	% create main axes
	am = axes('position',[.5 -.025 .55 .55]);
	view(3)
	ame = axes('position',[.5 .3 .55 .55]);
	view(3)
	amf = axes('position',[.5 .625 .55 .55]);
	view(3)
	
	% create an edit and text boxe for the rotation matrix
	eR = uicontrol('style','edit','max',2,'un','n','pos',[.05 .6 .35 .25],'fonts',8,'fontn','courier',...
		'backgroundcol',[1,1,1]);
	uicontrol('sty','te','un','n','pos',[.05 .85 .35 .07],'string','Rotation Matrix','fonts',20,'fontweig','b',...
		'backgroundcol','w')
	
	% three edits for the angles
	ea = uicontrol('sty','ed','un','n','pos',[.1 .45 .1 .1],'backgroundc','w','fonts',18,'string','20');
	eb = uicontrol('sty','ed','un','n','pos',[.1 .3 .1 .1],'backgroundc','w','fonts',18,'string','30');
	eg = uicontrol('sty','ed','un','n','pos',[.1 .15 .1 .1],'backgroundc','w','fonts',18,'string','60');
	
	% create three small axes to the left of each angle edit box,
	% then put a nice looking text in each axes
	axes('position',[0 .45 .1 .1]);
	axis off
	text('position',[0.3 .6],'string','\alpha','fontsize',50)
	
	axes('position',[0 .3 .1 .1]);
	axis off
	text('position',[0.3 .6],'string','\beta','fontsize',50)
	
	axes('position',[0 .15 .1 .1]);
	axis off
	text('position',[0.4 .7],'string','\gamma','fontsize',50)
	
	% create the button group at the bottom of the figure (backward/forward)
	hg = uibuttongroup('unit','no','pos',[.05 .025 .4 .1],'backgroundc','w');
	
	h1 = uicontrol('style','radiobutton','parent',hg,'un','n','pos',[0,0,.5,1],'background','w','string','Forward','fonts',20);
	h2 = uicontrol('style','radiobutton','parent',hg,'un','n','pos',[0.5,0,.5,1],'background','w','string','Backward','fonts',20);
	
	% pushbutton uicontrol
	uicontrol('un','n','pos',[.25 .15 .1 .1],'sty','pu','string','Go','fonts',20,'backgroundc','w','callback',@go1)
	
	% create convntion edit and text boxes
	uicontrol('un','n','pos',[.25 .5 .2 .1],'sty','te','string','Convention','fonts',20,'backgroundc','w')
	ce = uicontrol('un','n','pos',[.25 .45 .2 .1],'sty','e','string','xyz','fonts',20,'backgroundc','w');
	
	% revert back to main axes for future plots
	axes(am)
	rotate3d
	hold on
	
	% draw the three coloured frame lines
	Lx = line([0 1],[0,0],[0,0],'color',[1,0,0]);
	Ly = line([0 0],[0,1],[0,0],'color',[0,1,0]);
	Lz = line([0 0],[0,0],[0,1],'color',[0,0,1]);
	
	set([Lx,Ly,Lz],'linewidth',3)
	
	% draw the three base frame axis lines
	line([0 1.5],[0,0],[0,0],'color',[0,0,0],'linewidth',2);
	line([0 0],[0,1.5],[0,0],'color',[0,0,0],'linewidth',2);
	line([0 0],[0,0],[0,1.5],'color',[0,0,0],'linewidth',2);
	
	text('position',[1.3 0 .1],'string','x_0','fontw','b');
	text('position',[ 0 1.3 .1],'string','y_0','fontw','b');
	text('position',[.05 0.05 1.3],'string','z_0','fontw','b');
	
	
	% text for x,y and z for each frame line
	tX = text('position',[.7 0 .1],'string','x','fontw','b');
	tY = text('position',[ 0 .7 .1],'string','y','fontw','b');
	tZ = text('position',[.05 0.05 .7],'string','z','fontw','b');
	
	% create our sphere and plot it
	[xx,yy,zz]=sphere;
	
	pp = surf(xx.*.34,yy.*.34,zz.*.34);
	colormap bone
	
	% plot origin point
	plot3(0,0,0,'.','markersize',20,'color','k')
	
	% more axes options (changing axis properties)
	axis square
	axis equal
	axis off
	
	% wR stands for "what to rotate"
	wR = [tX,tZ,tY,Lx,Ly,Lz,pp];
	
	% pushbutton callback
	function go1(varargin)
		
		axes(am)
		% delete everything (start over)
		delete(wR)
		
		% create new lines and sphere
		[xx,yy,zz]=sphere;
		
		pp = surf(xx.*.34,yy.*.34,zz.*.34);
		colormap bone
		
		Lx = line([0 1],[0,0],[0,0],'color',[1,0,0]);
		Ly = line([0 0],[0,1],[0,0],'color',[0,1,0]);
		Lz = line([0 0],[0,0],[0,1],'color',[0,0,1]);
		
		set([Lx,Ly,Lz],'linewidth',3)
		
		tX = text('position',[.7 0 .1],'string','x','fontw','b');
		tY = text('position',[ 0 .7 .1],'string','y','fontw','b');
		tZ = text('position',[.05 0.05 .7],'string','z','fontw','b');
		
		wR = [tX,tZ,tY,Lx,Ly,Lz,pp];
		
		% get the value of the button group (forward or backward)
		v = get(hg,'selectedobject');
		
		% get value for convention
		conv = get(ce,'string');
		
		switch v
			
			% forward
			case h1
				
				% get angles values
				a = str2double(get(ea,'string'))/180*pi;
				b = str2double(get(eb,'string'))/180*pi;
				g = str2double(get(eg,'string'))/180*pi;
				
				% obtain rotation matrix
				s = eulfor(a,b,g,conv);
				
				% display rotation matrix
				set(eR,'string',num2str(s,5))
				
				% backward
			case h2
				
				% get rotation matrix and use eul() to obtain three angles
				s = eul(str2num(get(eR,'string')),conv); %#ok
				
				% display three angles
				set(ea,'string',num2str(s.a*180/pi));
				set(eb,'string',num2str(s.b*180/pi));
				set(eg,'string',num2str(s.g*180/pi));
				
		end
		
		% get Euler angles
		a = str2double(get(ea,'string'))/180*pi;
		b = str2double(get(eb,'string'))/180*pi;
		g = str2double(get(eg,'string'))/180*pi;
		
		% rotate first time
		for k = 20:20
			if strcmp(conv(1),'x')
				rotate(wR,[1,0,0],rad2deg(a),[0,0,0])
			elseif strcmp(conv(1),'y')
				rotate(wR,[0,1,0],rad2deg(a),[0,0,0])
			elseif strcmp(conv(1),'z')
				rotate(wR,[0,0,1],rad2deg(a),[0,0,0])
			end
			% 			pause(.05)
		end
		
		% get new vectors for x,y and z (for next rotation)
		pxx = eval(['get(L',conv(2),',','''xdata''',')']);
		pyy = eval(['get(L',conv(2),',','''ydata''',')']);
		pzz = eval(['get(L',conv(2),',','''zdata''',')']);
		
		
		% rotate second time
		for k = 20:20
			
			rotate(wR,[pxx(2)-pxx(1),pyy(2)-pyy(1),pzz(2)-pzz(1)],rad2deg(b),[0,0,0])
			
			% 			pause(.05)
		end
		
		pxx2 = eval(['get(L',conv(3),',','''xdata''',')']);
		pyy2 = eval(['get(L',conv(3),',','''ydata''',')']);
		pzz2 = eval(['get(L',conv(3),',','''zdata''',')']);
		
		% rotate third time
		for k = 20:20
			
			rotate(wR,[pxx2(2)-pxx2(1),pyy2(2)-pyy2(1),pzz2(2)-pzz2(1)],rad2deg(g),[0,0,0])
			
			% 			pause(.05)
		end
		
		axes(ame)
		% delete everything (start over)
		delete(wRe)
		
		% create new lines and sphere
		[xx,yy,zz]=sphere;
		
		pp = surf(xx.*.34,yy.*.34,zz.*.34);
		colormap bone
		
		Lx = line([0 1],[0,0],[0,0],'color',[1,0,0]);
		Ly = line([0 0],[0,1],[0,0],'color',[0,1,0]);
		Lz = line([0 0],[0,0],[0,1],'color',[0,0,1]);
		
		set([Lx,Ly,Lz],'linewidth',3)
		
		tX = text('position',[.7 0 .1],'string','x','fontw','b');
		tY = text('position',[ 0 .7 .1],'string','y','fontw','b');
		tZ = text('position',[.05 0.05 .7],'string','z','fontw','b');
		
		wRe = [tX,tZ,tY,Lx,Ly,Lz,pp];
		% 		wR = wRe;
		
		% get the value of the button group (forward or backward)
		v = get(hg,'selectedobject');
		
		% get value for convention
		conv = get(ce,'string');
		
		switch v
			
			% forward
			case h1
				
				% get angles values
				a = str2double(get(ea,'string'))/180*pi;
				b = str2double(get(eb,'string'))/180*pi;
				g = str2double(get(eg,'string'))/180*pi;
				
				% obtain rotation matrix
				s = eulfor(a,b,g,conv);
				
				% display rotation matrix
				set(eR,'string',num2str(s,5))
				
				% backward
			case h2
				
				% get rotation matrix and use eul() to obtain three angles
				s = eul(str2num(get(eR,'string')),conv); %#ok
				
				% display three angles
				set(ea,'string',num2str(s.a*180/pi));
				set(eb,'string',num2str(s.b*180/pi));
				set(eg,'string',num2str(s.g*180/pi));
				
		end
		
		% get Euler angles
		a = str2double(get(ea,'string'))/180*pi;
		b = str2double(get(eb,'string'))/180*pi;
		g = str2double(get(eg,'string'))/180*pi;
		
		% rotate first time
		for k = 1:20
			if strcmp(conv(1),'x')
				rotate(wRe,[1,0,0],rad2deg(a/20),[0,0,0])
			elseif strcmp(conv(1),'y')
				rotate(wRe,[0,1,0],rad2deg(a/20),[0,0,0])
			elseif strcmp(conv(1),'z')
				rotate(wRe,[0,0,1],rad2deg(a/20),[0,0,0])
			end
			pause(.05)
		end
		
		% get new vectors for x,y and z (for next rotation)
		pxx = eval(['get(L',conv(2),',','''xdata''',')']);
		pyy = eval(['get(L',conv(2),',','''ydata''',')']);
		pzz = eval(['get(L',conv(2),',','''zdata''',')']);
		
		
		% rotate second time
		for k = 1:20
			
			rotate(wRe,[pxx(2)-pxx(1),pyy(2)-pyy(1),pzz(2)-pzz(1)],rad2deg(b/20),[0,0,0])
			
			pause(.05)
		end
		
		pxx2 = eval(['get(L',conv(3),',','''xdata''',')']);
		pyy2 = eval(['get(L',conv(3),',','''ydata''',')']);
		pzz2 = eval(['get(L',conv(3),',','''zdata''',')']);
		
		% rotate third time
		for k = 1:20
			
			rotate(wRe,[pxx2(2)-pxx2(1),pyy2(2)-pyy2(1),pzz2(2)-pzz2(1)],rad2deg(g/20),[0,0,0])
			
			pause(.05)
		end
		
		axes(amf)
		% delete everything (start over)
				delete(wRf)
		
		% create new lines and sphere
		[xx,yy,zz]=sphere;
		
		pp = surf(xx.*.34,yy.*.34,zz.*.34);
		colormap bone
		
		Lx = line([0 1],[0,0],[0,0],'color',[1,0,0]);
		Ly = line([0 0],[0,1],[0,0],'color',[0,1,0]);
		Lz = line([0 0],[0,0],[0,1],'color',[0,0,1]);
		
		set([Lx,Ly,Lz],'linewidth',3)
		
		tX = text('position',[.7 0 .1],'string','x','fontw','b');
		tY = text('position',[ 0 .7 .1],'string','y','fontw','b');
		tZ = text('position',[.05 0.05 .7],'string','z','fontw','b');
		
		wRf = [tX,tZ,tY,Lx,Ly,Lz,pp];
		% 		wR = wRf;
		
		% get the value of the button group (forward or backward)
		v = get(hg,'selectedobject');
		
		% get value for convention
		conv = get(ce,'string');
		
		switch v
			
			% forward
			case h1
				
				% get angles values
				a = str2double(get(ea,'string'))/180*pi;
				b = str2double(get(eb,'string'))/180*pi;
				g = str2double(get(eg,'string'))/180*pi;
				
				% obtain rotation matrix
				s = eulfor(a,b,g,conv);
				
				% display rotation matrix
				set(eR,'string',num2str(s,5))
				
				% backward
			case h2
				
				% get rotation matrix and use eul() to obtain three angles
				s = eul(str2num(get(eR,'string')),conv); %#ok
				
				% display three angles
				set(ea,'string',num2str(s.a*180/pi));
				set(eb,'string',num2str(s.b*180/pi));
				set(eg,'string',num2str(s.g*180/pi));
				
		end
		
		% get Euler angles
		a = str2double(get(ea,'string'))/180*pi;
		b = str2double(get(eb,'string'))/180*pi;
		g = str2double(get(eg,'string'))/180*pi;
		
		% rotate first time
		for k = 1:20
			if strcmp(conv(3),'x')
				rotate(wRf,[1,0,0],rad2deg(g/20),[0,0,0])
			elseif strcmp(conv(3),'y')
				rotate(wRf,[0,1,0],rad2deg(g/20),[0,0,0])
			elseif strcmp(conv(3),'z')
				rotate(wRf,[0,0,1],rad2deg(g/20),[0,0,0])
			end
			pause(.05)
		end
		%
		% 		% get new vectors for x,y and z (for next rotation)
		% 		pxx = eval(['get(L',conv(2),',','''xdata''',')']);
		% 		pyy = eval(['get(L',conv(2),',','''ydata''',')']);
		% 		pzz = eval(['get(L',conv(2),',','''zdata''',')']);
		%
		%
		% rotate second time
		for k = 1:20
			
			if strcmp(conv(2),'x')
				rotate(wRf,[1,0,0],rad2deg(b/20),[0,0,0])
			elseif strcmp(conv(2),'y')
				rotate(wRf,[0,1,0],rad2deg(b/20),[0,0,0])
			elseif strcmp(conv(2),'z')
				rotate(wRf,[0,0,1],rad2deg(b/20),[0,0,0])
			end
			
			pause(.05)
		end
		%
		% 		pxx2 = eval(['get(L',conv(3),',','''xdata''',')']);
		% 		pyy2 = eval(['get(L',conv(3),',','''ydata''',')']);
		% 		pzz2 = eval(['get(L',conv(3),',','''zdata''',')']);
		
		% rotate third time
		for k = 1:20
			
			if strcmp(conv(1),'x')
				rotate(wRf,[1,0,0],rad2deg(a/20),[0,0,0])
			elseif strcmp(conv(1),'y')
				rotate(wRf,[0,1,0],rad2deg(a/20),[0,0,0])
			elseif strcmp(conv(1),'z')
				rotate(wRf,[0,0,1],rad2deg(a/20),[0,0,0])
			end
			
			pause(.05)
		end
		
	end
	axes(amf)
	rotate3d
	hold on
	
	% draw the three coloured frame lines
	Lx = line([0 1],[0,0],[0,0],'color',[1,0,0]);
	Ly = line([0 0],[0,1],[0,0],'color',[0,1,0]);
	Lz = line([0 0],[0,0],[0,1],'color',[0,0,1]);
	
	set([Lx,Ly,Lz],'linewidth',3)
	
	% draw the three base frame axis lines
	line([0 1.5],[0,0],[0,0],'color',[0,0,0],'linewidth',2);
	line([0 0],[0,1.5],[0,0],'color',[0,0,0],'linewidth',2);
	line([0 0],[0,0],[0,1.5],'color',[0,0,0],'linewidth',2);
	
	text('position',[1.3 0 .1],'string','x_0','fontw','b');
	text('position',[ 0 1.3 .1],'string','y_0','fontw','b');
	text('position',[.05 0.05 1.3],'string','z_0','fontw','b');
	
	
	% text for x,y and z for each frame line
	tX = text('position',[.7 0 .1],'string','x','fontw','b');
	tY = text('position',[ 0 .7 .1],'string','y','fontw','b');
	tZ = text('position',[.05 0.05 .7],'string','z','fontw','b');
	
	% create our sphere and plot it
	[xx,yy,zz]=sphere;
	
	pp = surf(xx.*.34,yy.*.34,zz.*.34);
	colormap bone
	
	% plot origin point
	plot3(0,0,0,'.','markersize',20,'color','k')
	
	% more axes options (changing axis properties)
	axis square
	axis equal
	axis off
	
	% wR stands for "what to rotate"
	wRf = [tX,tZ,tY,Lx,Ly,Lz,pp];
	
	axes(ame)
	rotate3d
	hold on
	
	% draw the three coloured frame lines
	Lx = line([0 1],[0,0],[0,0],'color',[1,0,0]);
	Ly = line([0 0],[0,1],[0,0],'color',[0,1,0]);
	Lz = line([0 0],[0,0],[0,1],'color',[0,0,1]);
	
	set([Lx,Ly,Lz],'linewidth',3)
	
	% draw the three base frame axis lines
	line([0 1.5],[0,0],[0,0],'color',[0,0,0],'linewidth',2);
	line([0 0],[0,1.5],[0,0],'color',[0,0,0],'linewidth',2);
	line([0 0],[0,0],[0,1.5],'color',[0,0,0],'linewidth',2);
	
	text('position',[1.3 0 .1],'string','x_0','fontw','b');
	text('position',[ 0 1.3 .1],'string','y_0','fontw','b');
	text('position',[.05 0.05 1.3],'string','z_0','fontw','b');
	
	
	% text for x,y and z for each frame line
	tX = text('position',[.7 0 .1],'string','x','fontw','b');
	tY = text('position',[ 0 .7 .1],'string','y','fontw','b');
	tZ = text('position',[.05 0.05 .7],'string','z','fontw','b');
	
	% create our sphere and plot it
	[xx,yy,zz]=sphere;
	
	pp = surf(xx.*.34,yy.*.34,zz.*.34);
	colormap bone
	
	% plot origin point
	plot3(0,0,0,'.','markersize',20,'color','k')
	
	% more axes options (changing axis properties)
	axis square
	axis equal
	axis off
	
	% wR stands for "what to rotate"
	wRe = [tX,tZ,tY,Lx,Ly,Lz,pp];
	
	uicontrol('style','text','units','normalized','pos',[.55 .5 .1 .05],'string','Euler','fontweight','bold',...
		'fontsize',15,'backgroundcolor',[1,1,1])
	uicontrol('style','text','units','normalized','pos',[.55 .8 .1 .05],'string','Fixed','fontweight','bold',...
		'fontsize',15,'backgroundcolor',[1,1,1])
	
end