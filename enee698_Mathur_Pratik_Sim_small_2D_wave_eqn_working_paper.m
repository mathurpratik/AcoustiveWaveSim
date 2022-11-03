%Name: Pratik Mathur
%Course: ENEE 698 - Project in Electrical Enginering
%Description: Simulates a viscous wave equation propagating towards a
%hydrophone sitting about 100 meters away from plasma source. A technique
%known as Finite Difference Method was used to implement the viscous wave
%equation 16 found in this article: 
% https://vtechworks.lib.vt.edu/bitstream/handle/10919/35983/Chapter2.pdf?sequence=4&isAllowed=y

L_x=199e-6; % distance covered by the simulation (m) 1000e-6
L_y=199e-6; % distance covered by the simulation (m) 1000e-6
end_time = 500e-9; %stop simulation after 120ns

dt=1e-10; % time-step for use by Finite Difference Methods
NT=end_time/dt;
v = .001/998; % kinematic viscosity 
c0 = 1480; % speed of sound in water

% create the distance mesh
%x_mesh = linspace(0,L,x_mesh_param);
dx = 2.7000e-06;
dy = 2.2000e-06;
x = 0:dx:499*dx;
y = 0:dy:499*dy;
[x_mesh, y_mesh] = meshgrid(x,y);

% invoke the simulation
states = simulate(x_mesh, y_mesh, dt,v,c0,NT);


% simulate () - takes as input the distance mesh, time step, viscocity, and
% sound speed in water coefficients and computes all the states for each
% time step until the wave propagates 100 meters (distance to a hydrophone)

function states = simulate(x_mesh, y_mesh,dt,v,c0,NT)     
    I=length(x_mesh);
    J=length(y_mesh);
    N = I*J;
    dx=x_mesh(1,2)-x_mesh(1,1);
    dy=y_mesh(2,1)-y_mesh(1,1);
    
    L_x = x_mesh(I,I);
    L_y = y_mesh(J,J);
    states = [];
    
    % define what a acoustic pulse looks like underwater
    % initially just making it a simple Guassian that starts around x=450
    % meters, and it will propagate to about x = 550 meters (total 100 m)
    
    
%     
    
    a = .545;
    E = 50e-6;
    c_1 = 5190; 
    c_2 = 25306; 
    c=1483;

    load('paper_wave.mat');
%     syms t;
%     f = ( 25/( 3*pi)* (E/998) )^0.2 * 0.4 * (1e9)^a*t^(-a)+c;
%     for t=1:200 % ns
%         tt = t/2;
%         r(t) = int(f,0,tt) ;
%         
%         u_s = shock_velocity(tt,E,1);
%          p_s(t) = c_1*998*u_s*(10^( (u_s-c0 )/c_2)-1 );
%     end
%     
%     r = double(r) * 1e-9;
%     save('paper_wave.mat','r','p_s')
    %u_s(1) = u_s(2);
%    
%     
%     
%     wave = zeros(size(x_mesh)); zeros(size(xmesh, size(ymesh);
%     
%     for b=1:length(x_mesh)
%         x = abs(x_mesh(b) -L/2 ); % put laser in the middle
%         
%         if (x >= max(r))
%             % assuming that if we have a value in this mesh that is bigger than the
%             % maximum value of r, then the pressure is 0 there
%             wave(b) = 0; 
%             continue;
%         end
%         
%         if (x <= min(r))
%             wave(b) = p_s(1);
%             continue;
%         end
%         
%         % do some interpolation
%         % give this value of x, find the values of r it is between
%         nearest = dsearchn(r',[x]);
%         
%         % nearest index could be above or below x
%         % check which it is
%         if (r(nearest) < x)
%             % r(nearest) < x < r(nearest+1)
%             delta = x - r(nearest);
%             wave(b) = (1-delta)*p_s(nearest) + delta*p_s(nearest+1);
%         else
%             % r(nearest-1) < x < r(nearest)
%             delta = x - r(nearest-1);
%             wave(b) = (1-delta)*p_s(nearest-1) + delta*p_s(nearest);
%         end
%         
%     end
    %wave = p_s;
    
    
    for a=1:length(x_mesh)
        for b=1:length(y_mesh)
            
            
            x = abs(x_mesh(a,b) -L_x/2 ); % put laser in the middle
            y = abs(y_mesh(a,b) -L_y/2 ); % put laser in the middle
            
            if (a == 250 && b == 250)
                disp('a=b=250');
            end
            rad = sqrt(x^2 + y^2);
            if (rad >= max(r))
                % assuming that if we have a value in this mesh that is bigger than the
                % maximum value of r, then the pressure is 0 there
                wave(a,b) = 0; 
                continue;
            end

            if (rad <= min(r))
                wave(a,b) = p_s(1);
                continue;
            end
    %         
            % do some interpolation
            % give this value of rad, find the values of r it is between
            nearest = dsearchn(r',[rad]);

            % nearest index could be above or below rad
            % check which it is
            if (r(nearest) < rad)
                % r(nearest) < rad < r(nearest+1)
                delta = rad - r(nearest);
                wave(a,b) = (1-delta)*p_s(nearest) + delta*p_s(nearest+1);
            else
                % r(nearest-1) < rad < r(nearest)
                delta = rad - r(nearest-1);
                wave(a,b) = (1-delta)*p_s(nearest-1) + delta*p_s(nearest);
            end
        end
    end
    %wave = exp(-(x_mesh-L_x/2).^2./(1e-11) -(y_mesh-L_y/2).^2./(1e-11)  ); 
    %wave = P(1:1000);
    % here I am duplicating the Guassian pulse because the way I am
    % performing these computations is keeping a vector of size 2N, where
    % the first N represents the current state the the last N represents
    % the old state.  So at each time step, we always have a current state
    % and old state.  This will help with performing matrix transformations
    % to compute a new state. 
    
    %wave = wave';
    
    wave_1D = convertFrom2Dto1D(wave);
    
    wave_1D = [wave_1D; wave_1D];
    % initially assuming the water is still, so "state_initial" starts out
    % just being all 0s. 
    state_initial = zeros(1,2 * I * J);
    %state_initial(1) = 0;
    %state_initial(length(state_initial)) = 0;
    
    % initializes the first state to all 0s (i.e. water is still and no
    % laser has been beamed yet)
    % state_first = firstIteration(state_initial, x_mesh,dt,v,c0 );
    
    % keep track of current state and previous state
    % state_initial = [state_first state_initial];
    
    % appending this to our states matrix ("states" will keep track of the
    % entire simulation for all time steps).
    %states = [states; state_initial'];
    states = state_initial';
    % we are now ready to perform new_state = D * current_state
    %D = transform_matrix(N,c0,dx,dt,v);
    D = transform_matrix_sparse_efficient_2D(I,J,c0,dx,dy,dt,v); 
    
    %E = eig(full(D));
    %max_eigenvalue = max(abs(E));
    
%     if (max_eigenvalue > 1)
%         
%         return;
%     end
    
    %if max_eigenvalue >= 1
    %    disp('max eigenvalue of transform matrix bigger or equal to 1. This will cause numerical instability. Quitting.')
    %    return;
    %end
    
    % plotting states in "real" time as new states are computed
    figure(1); clf
    
    %axis tight manual

    ax = gca;
    ax.NextPlot = 'replaceChildren';
    
    
    new_state = states(:,1);
    mkdir('test')
    wobj = VideoWriter('test1.avi')
    wobj.FrameRate = 10;                  % frames per second (video speed)
    open(wobj);     
    counter =0;
    for i=1:NT % NT frames of time steps
        % comput the new state from old state
        new_state = D * new_state;
        counter = counter + 1
         % I send two pulses (one at time step i=79 and another at i=2*79)
        % You might be wondering why I chose multiples of 79. This is
        % because I don't want the pulses to interfere with each other as
        % they get send out. Waiting 79 * dt time steps allows for the
        % first pulse to propagate a little before sending the second one.
        % This also makes it less confusing for the hydrophone to sense the
        % pulse. Notice there is no noise added (keeping it simple for
        % now).
        if i == 1%fix(NT/20) %|| i == fix(2*NT/20)
            new_state = new_state + wave_1D;
        end
        ax = gca();
        surf(ax,x_mesh,y_mesh,reshape(new_state(1:N), [I J])');
        if (i==1)
             
             %plot3(ax,[x_mesh x_mesh+100e-6], [y_mesh y_mesh+100e-6],reshape(states(1:N,i), [I J])');
             
             % h = plot(ax, x1, y1, x2, y2, 'Linewidth',3);
             set(ax, 'XLimMode', 'manual', 'YLimMode', 'manual', 'ZLimMode', 'manual');
        else
             %plot3(ax,[x_mesh x_mesh+100e-6], [y_mesh y_mesh+100e-6],reshape(states(1:N,i), [I J])');
            
        end
        
        axis tight
        axis([ -L_x L_x -L_y L_y 0 1e9])
        
        max_pressure = abs(max(new_state(1:N)));
        max_pres_idx = find (abs(new_state(1:N)) == max_pressure);
        [x_idx y_idx] = row_and_col(max_pres_idx(1),I,J);
        max_dist = sqrt( (x_mesh(x_idx,y_idx)-L_x/2)^2+ (y_mesh(x_idx,y_idx) -L_y/2 )^2 );
        
        text(0,0,['distance = ' num2str(max_dist) ' time = ' num2str(i*dt)])
        disp(['distance = ' num2str(max_dist) ' time = ' num2str(i*dt)])
        fname = ['test' filesep 'test' num2str(i)]; % full name of image
        print('-djpeg','-r200',fname)     % save image with '-r200' resolution
          currImage = imread([fname '.jpg']);       % read saved image
        frame = im2frame(currImage);              % convert image to frame
        writeVideo(wobj,frame);           % save frame into video
        
       
       
        
        %states = [states new_state]; 
    end
    
    close(wobj);
end

function p2 = firstIteration(p1, x_mesh, dt, v, c0)
    dx=x_mesh(2)-x_mesh(1);
    
    p2 = zeros(size(p1));
    
    for i=2:length(x_mesh)-1
        p_txx = 0;
        p_xx  = (p1(i+1)-2*p1(i)+p1(i-1)) / (dx^2);
        p2(i) = (dt)^2 * (4.0/3 * v * p_txx + c0^2*p_xx)+p1(i);
    end
    
    % forward difference (NO BOUNDARY conditions)
    p_txx = 0;
    p_xx = (p1(3)-2*p1(2)+p1(1)) / (dx^2);
    %p2[0] = (dt)**2 * (4.0/3 * v * p_txx + c0**2*p_xx)+p1[0]
    
    % taylor series
    %p2[0]=1/(1-1/dt+1/(2*(dt)**2))* ((1-1/dt+1/((dt)**2))*p2[1]-p2[2]/(2*(dt**2)))
    
    %p2[0] = p1[1] # NEUMANN 
    
    % if you comment out lines 34,37,39 you get Dirichlet Conditions
    
    % backward difference
    p_txx = 0;
    p_xx = (p1(length(p1))-2*p1(length(p1)-1)+p1(length(p1)-2)) / (dx^2);
    %p2[-1] = (dt)**2 * (4.0/3 * v * p_txx + c0**2*p_xx)+p1[-1]
    
    
end

function D = transform_matrix(N,c0,dx,dt,v)
    D1 = zeros(2*N,2*N);
    
    for i=2:(N-1)
        D1(i,i) = 2;
        D1(i,i+N) = -1;
    end
    
    D1(N+1:(2*N), 1:N) = eye(N);
    
    block = zeros(N,N);
     for i=2:(N-1)
        block(i,i)=-2;
        block(i,i-1)=1;
        block(i,i+1)=1;
     end
     
    D2 = zeros(2*N,2*N);
    D2(1:N,1:N)=block*c0^2*((dt/dx)^2);
    
    D3 = zeros(2*N,2*N);
    
    D3(1:N,1:N)=block;
    D3(1:N,N+1:(2*N)) = -block;
    
    
    D3 = D3 * (4.0/3)*v*dt/(dx^2);
    
    
    D = D1+D2+D3;
end
function D = transform_matrix_2D(I,J,c0,dx,dy,dt,v)
    
    N = I*J;
    D1 = zeros(2*N,2*N);
    
    for i=1:(I-2)
        for j=1:(J-2)
            ij1 = i*J+j+1; % row
            ij2 = i*J+j+1 + I*J; % col 
            D1(ij1,ij1) = 2;
            D1(ij1,ij2) = -1;
        end
    end
    
    D1(N+1:(2*N), 1:N) = eye(N);
    
%     block = zeros(N,N);
%      for i=2:(N-1)
%         block(i,i)=-2;
%         block(i,i-1)=1;
%         block(i,i+1)=1;
%      end
     
    D2A = zeros(2*N,2*N);
    D2B = zeros(2*N,2*N);
        for i=1:(I-2)
            for j=1:(J-2)
                ij1 = i*J+j+1; % row
                ij2 = (i-1)*J+j+1 ; % col1
                ij3 = (i)*J+j+1 ; % col2
                ij4 = (i+1)*J+j+1 ; % col3
                
                D2A(ij1,ij2) = 1;
                D2A(ij1,ij3) = -2;
                D2A(ij1,ij4) = 1;
                
                
                %%%%%%%%%%%%%%%%%%%
                ij5 = i*J+(j-1)+1;
                D2B(ij1,ij5) = 1;
                D2B(ij1,ij5+1) = -2;
                D2B(ij1,ij5+2) = 1;
                
            end
        end
    D2A = (c0*dt/dx)^2 * D2A;
    D2B = (c0*dt/dy)^2 * D2B;
    %D2(1:N,1:N)=block*c0^2*((dt/dx)^2);
    
    %D3 = zeros(2*N,2*N);
    
     D3A = zeros(2*N,2*N);
     D3B = zeros(2*N,2*N);
        for i=1:(I-2)
            for j=1:(J-2)
                ij1 = i*J+j+1; % row
                ij2 = (i-1)*J+j+1 ; % col1
                ij3 = (i)*J+j+1 ; % col2
                ij4 = (i+1)*J+j+1 ; % col3
                
                %ij2 > I*J
                
                
                if (ij4+N  <= 2*N)
                    D3A(ij1,ij2) = 1;
                    D3A(ij1,ij3) = -2;
                    D3A(ij1,ij4) = 1;
                
                
                    D3A(ij1,ij2+N) = -1;
                    D3A(ij1,ij3+N) = 2;
                    D3A(ij1,ij4+N) = -1;
                end
                %%%%%%%%%%%%%%%%%%%%%
                ij5 = i*J+(j-1)+1 ; % col1
                
                if (ij5+2+N  <= 2*N)
                    D3B(ij1,ij5) = 1;
                    D3B(ij1,ij5+1) = -2;
                    D3B(ij1,ij5+2) = 1;
                
                
                    D3B(ij1,ij5+N) = -1;
                    D3B(ij1,ij5+1+N) = 2;
                    D3B(ij1,ij5+2+N) = -1;
                end
            end
        end
        
        %D3B=D3A;
    %D3(1:N,1:N)=block;
    %D3(1:N,N+1:(2*N)) = -block;
    D3A = 4*v*dt / (dx^2) * D3A;
    D3B = 4*v*dt / (dy^2) * D3B;
    
    %D3 = D3 * (4.0/3)*v*dt/(dx^2);
    
    
    D = D1+D2A+D2B+D3A+D3B;
    disp('D');
end

function OneDimVec = convertFrom2Dto1D(D)
    
    
    [M N] = size(D);
    OneDimVec=zeros(M*N,1);
    counter = 0; 
    for m=1:M
        for n=1:N
            counter = counter + 1;
            OneDimVec(counter) = D(m,n);
        end
    end
end

function [i j] = row_and_col(idx,I,J)
    i=fix(idx/J)+1;
    j=idx - (i-1)*J;
end