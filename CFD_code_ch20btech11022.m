clc;clear;close all;

imax =32;
jmax = 32;
u_init =0.002;
re = 1000;
% re = u_init * H / nu
% u_init and H is adjusted so that Reynold's number is within laminar flow

L=1;
H=0.5;
nu = 1e-6;
e = 1e-6;% e to be used in Gauss seidal method

u0 = zeros(jmax,imax);
v0=zeros(jmax,imax);
u = zeros(jmax,imax);
v=zeros(jmax,imax);
w=zeros(jmax,imax);
w0 = zeros(jmax,imax);
psi0=zeros(jmax,imax);
psi=zeros(jmax,imax);

% u,v,w,psi correspond to (n+1)th time (New time points)
% u0,v0,w0,psi0 correspond to (n)th time (Old time points)

delx = L/(imax-1);
dely = H/(jmax-1);

y_axis = 0:dely:H;
tmax =10;
delt = 0.005;
nt = tmax/delt;

%top wall
u0(1,:)=0;
v0(1,:)=0;
psi0(1,:) = u_init*H;

%bottom wall
u0(jmax,:) = 0;
v0(jmax,:)=0;
psi0(jmax,:) = 0;

%left wall
u0(:,1) =u_init;
v0(:,1)=0;
psi0(:,1) = u0(:,1).* (H:-dely:0)';

%updating vorticity for top wall

for i=1:imax
    w0(1,i) = -1* ((u0(1,i)-u0(2,i))/dely);
end

w=w0;
u=u0;
v=v0;
psi=psi0;

t=0;
plotTimings = [0 nt/4 nt/2 3*(nt/4) nt];
if (any(plotTimings(:) == t))
    figure(1);
    plot(u(:, 1), y_axis);
    hold on;
    figure(2);
    plot(abs(v(:, 1)), y_axis);
    hold on;
end
while(t < nt)
    %% calculating w at new time step
    for i= 2:imax-1
        w(1, i) = 2 * (psi0(1, i) - psi0(2, i)) / dely^2;
        w(jmax, i) = 2 * (psi0(jmax, i) - psi0(jmax - 1, i)) / dely^2;
    end
    for i=2:imax-1
        for j=2:jmax-1
            temp = (1/re)*( (w0(j,i-1) -2*w0(j,i) + w0(j,i+1))/delx^2 + (w0(j-1,i) -2*w0(j,i) + w0(j+1,i))/dely^2 );
            w(j,i) = delt*( temp - u0(j,i)*( (w0(j,i+1) - w0(j,i-1)) / (2*delx) ) - v0(j,i)*( (w0(j-1,i) - w0(j+1,i)) / (2*dely) )) + w0(j,i);
        end
    end
    w(:,imax) = w(:,imax-1);
    %% psi at new time step and also implementing gauss seidal method
    satisfied = false;
    while (satisfied == false)
        for i=2:imax-1
            for j=2:jmax-1
                psi(j,i) = 0.5 * ((delx^2 * dely^2)/(delx^2 + dely^2)) * (w(j,i) + (psi0(j,i+1) + psi0(j,i-1))/delx^2 + (psi0(j+1,i) + psi0(j-1,i))/dely^2);
            end
        end
        psi(:,imax) = psi(:,imax-1);
        curr = true;
        for i=2:imax-1
            for j=2:jmax-1
                if(abs(psi(j,i)-psi0(j,i)) > e)
                    curr = false;
                end
            end
        end
        if(curr==true)
            satisfied = true;
        end
        psi0 = psi;
    end
    %% calculating u, v at new time step
    for i=2:imax-1
        for j=2:jmax-1
            u(j,i) = (psi(j-1,i) - psi(j+1,i))/(2*dely);
            v(j,i) = -1 * (psi(j,i+1) - psi(j,i-1))/(2*delx);
        end
    end
    u(:,imax) = u(:,imax-1);
    v(:,imax) = v(:,imax-1);
    %% Updating old time step info with the new time step's
    w0=w;
    u0=u;
    v0=v;
    t=t+1;
    %% Plotting u, v on y-axis
    if (any(plotTimings(:) == t))
        figure(1);
        plot(u(:, imax), y_axis);
        hold on;
        figure(2);
        plot(abs(v(:, jmax)), y_axis);
        hold on;
    end
end
figure(1);
title('Velocity along x-axis vs. height');
xlabel('u');
ylabel('h');
legend('0',num2str(nt * delt / 4),num2str(nt * delt / 2), num2str(3 * nt * delt / 4), num2str(nt * delt));
hold off;
figure(2);
title('Velocity along y-axis vs. height');
xlabel('v');
ylabel('h');
legend('0',num2str(nt * delt / 4),num2str(nt * delt / 2), num2str(3 * nt * delt / 4), num2str(nt * delt));
hold off;