figure(1)

clf;
clear;
colormap(jet);

model = 'Bonny_invivo';

%% Numerical parameters
L=3;  % domain length (µm)
W=0.5;  % domain width (µm)
depth=1; % domain depth (µm)
dx=0.1; % spacial step  (µm)
dt=0.002; % time step (s)
Xgrid=(0:dx:L)';
N=length(Xgrid);
nsteps=20*500;

%% Plot params
plotEvery=4*500;

time=[1:nsteps]*dt;
tmax=max(time);
Dcyttime = zeros(size(time));
Ecyttime = zeros(size(time));
Dmemtime = zeros(size(time));
DEmemtime = zeros(size(time));
Ememtime = zeros(size(time));

Dmemall = [];
Ememall = [];

%% Diffusion constants (Bonny 2013)
switch model
    case 'Bonny_invitro'
        % in vitro
        Dmem=0.3; % MinD membrane diffusion in µm^2/s
        Emem=1.8; % MinE membrane diffusion in µm^2/s
        Dcyt=50; % cyt/bulk diffusion in µm^2/s
    case 'Bonny_invivo'
        % in vivo
        Dmem=0.06; % MinD membrane diffusion in µm^2/s
        Emem=0.3; % MinE membrane diffusion in µm^2/s
        Dcyt=14; % cyt/bulk diffusion in µm^2/s
        
endswitch

%% A - diffusion matrix, no flux
Ident=speye(N);
A=-2*Ident;

for i = 1:N-1
  A(i,i+1)=1;
  A(i+1,i)=1;
endfor

A(1,1)=-1;
A(N,N)=-1;
A=A*dt/(2*dx^2);
A=sparse(A);

%% Initial conditions 
switch model
    case 'Bonny_invitro'
        % in vitro
        D=4.8e2*ones(N,1)*depth*W; % in units of 1/µm^3
        E=7.0e2*ones(N,1)*depth*W; % in units of 1/µm^3
    case 'Bonny_invivo'
        % in vivo
        D=2.2e3*ones(N,1)*depth*W; % in units of 1/µm^3
        E=1.5e3*ones(N,1)*depth*W; % in units of 1/µm^3
        
endswitch

u1=zeros(N,1); % concentration of membrane species 1 in units of molecules per µm^2
u2=zeros(N,1); 
u3=zeros(N,1); 

u1(1:round(N/3),:)=3*2.2e3*ones(round(N/3),1)*depth*W; % in units of 1/µm^3
D=zeros(N,1);

%load BonnyIC1m1c.m
u4=zeros(N,1); 
u5=zeros(N,1); 
u6=zeros(N,1); 
u7=zeros(N,1); 

%% LOOP
for k=1:nsteps

    RHS=reactions(model,D(:)/(depth),E(:)/(depth),u1(:),u2(:),u3(:),u4(:),u5(:),u6(:),u7(:));
    
    D(:) = (Ident-Dcyt*A) \ ( (Ident+Dcyt*A)*D(:) + dt*RHS(:,1) );
    E(:) = (Ident-Dcyt*A) \ ( (Ident+Dcyt*A)*E(:) + dt*RHS(:,2) );

    u1(:) = (Ident-Dmem*A) \ ( (Ident+Dmem*A)*u1(:) + dt*RHS(:,3) );
    u2(:) = (Ident-Dmem*A) \ ( (Ident+Dmem*A)*u2(:) + dt*RHS(:,4) );
    u3(:) = (Ident-Emem*A) \ ( (Ident+Emem*A)*u3(:) + dt*RHS(:,5) );
    %u4(:) = (Ident-Dmem*A) \ ( (Ident+Dmem*A)*u4(:) + dt*RHS(:,6) );
    %u5(:) = (Ident-Dmem*A) \ ( (Ident+Dmem*A)*u5(:) + dt*RHS(:,7) );
    %u6(:) = (Ident-Dmem*A) \ ( (Ident+Dmem*A)*u6(:) + dt*RHS(:,8) );
    %u7(:) = (Ident-Dmem*A) \ ( (Ident+Dmem*A)*u7(:) + dt*RHS(:,9) );

    Dcyttime(k)=D(1,1);
    Ecyttime(k)=E(1,1);
    Dmemtime(k)=u1(1,1);
    DEmemtime(k)=u2(1,1);
    Ememtime(k)=u3(1,1);

    Dmemall(k,:)=u1'+u2';
    Ememall(k,:)=u2'+u3';

    % Plot
    if mod(k,plotEvery)==0

        subplot(2,3,1);
        plot(Xgrid',D');
        title(['MinD (cyt)   t = ' int2str(k*dt) ' s']);

        subplot(2,3,2);
        plot(Xgrid',E');
        title('MinE (cyt)');

        subplot(2,3,4);
        plot(Xgrid',u1');
        title('MinD (mem)');

        subplot(2,3,5);
        plot(Xgrid',u2');
        title('MinDE (mem)');

        subplot(2,3,6);
        plot(Xgrid',u3');
        title('MinE (mem)');
        
        subplot(2,3,3);
        plot(time,Dcyttime,time,Ecyttime,time,Dmemtime,time,DEmemtime,time,Ememtime);
        legend('MinDc','MinEc','MinDm','MinDEm','MinEm');
        legend('boxoff');

        pause(0.001);

    endif

endfor

%% input for imagesc can have a max length of 2^14-1
for k=1:N
  Dmemplot(:,k)=interp1(time,Dmemall(:,k),linspace(dt,max(time),2^14-1));
  Ememplot(:,k)=interp1(time,Ememall(:,k),linspace(dt,max(time),2^14-1));
endfor

%% plot kymograph
maxD=max(max(Dmemall));
maxE=max(max(Ememall));
subplot(2,1,1);
imagesc([0 tmax],[0 L],Dmemplot',[0,maxD]);
xlabel('Time (s)');
ylabel('Position (\mu{m})');
title('Kymograph for MinD on Membrane');
axis('xy');

subplot(2,1,2);
imagesc([0 tmax],[0 L],Ememplot',[0,maxE]);
xlabel('Time (s)');
ylabel('Position (\mu{m})');
title('Kymograph for MinE on Membrane');
axis('xy');

%% save kymograph
orient("portrait");
fname=["PDE_" model "_length" int2str(L)];
%print([fname ".pdf"]);

%save BonnyIC1m1c.m D E u1 u2 u3;