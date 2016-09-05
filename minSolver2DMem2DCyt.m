figure(1)

clf;
clear;
colormap(jet);

model = 'Bonny_invivo';

%% Numerical parameters

%{
% in vitro dimensions
L=100;  % domain length
W=100;  % domain width
depth=90; % domain depth
dx=1;
dt=0.01;
%}

L=3;  % domain length
W=1;  % domain width
depth=1; % domain depth
dx=0.1;
Xgrid=(0:dx:L)';
Ygrid=(0:dx:W)';
N=length(Xgrid);
M=length(Ygrid);
NM=N*M;
nsteps=20*500;
dt=0.002;

%% Plot params
plotEvery=4*500;

time=[1:nsteps]*dt;
Dcyttime = zeros(size(time));
Ecyttime = zeros(size(time));
Dmemtime = zeros(size(time));
DEmemtime = zeros(size(time));
Ememtime = zeros(size(time));
Dmemtime2 = zeros(size(time));

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
        
end

%% A - diffusion matrix, no flux
Ablock11 = -3*diag(ones(1,N)) + diag(ones(1,N-1),1) + diag(ones(1,N-1),-1);
Ablock11(1,1) = -2;
Ablock11(N,N) = -2;
Ablock1=-4*diag(ones(1,N)) + diag(ones(1,N-1),1) + diag(ones(1,N-1),-1);
Ablock1(1,1) = -3;
Ablock1(N,N) = -3;
Ablock2=1*diag(ones(1,N));
AblockMN = -3*diag(ones(1,N)) + diag(ones(1,N-1),1) + diag(ones(1,N-1),-1);
AblockMN(N,N) = -2;
AblockMN(1,1) = -2;

A=zeros(NM,NM);
A(1:N,1:2*N)=[Ablock11 Ablock2];

for k=2:M-1
    A((k-1)*N+1:k*N,(k-2)*N+1:(k+1)*N) = [Ablock2 Ablock1 Ablock2];
end
A((M-1)*N+1:M*N,(M-2)*N+1:M*N)=[Ablock2 AblockMN];

%%
A=A*dt/(2*dx^2);
A=sparse(A);
Ident=speye(NM);

%% Initial conditions 
switch model
    case 'Bonny_invitro'
        % in vitro
        D=4.8e2*ones(N,M)*depth; % in units of 1/µm^3
        E=7.0e2*ones(N,M)*depth; % in units of 1/µm^3
    case 'Bonny_invivo'
        % in vivo
        D=2.2e3*ones(N,M)*depth; % in units of 1/µm^3
        E=1.5e3*ones(N,M)*depth; % in units of 1/µm^3
        
end

u1=zeros(N,M); % concentration of membrane species 1 in units of molecules per µm^2
%u1(1:round(N/2),:)=u1(1:round(N/2),:)+10000;
%D(1:round(N/2),:)=D(1:round(N/2),:)-10000;
%u1(round(N/2),round(M/2))=10000;

u1(1:round(N/3),:)=3*2.2e3*ones(round(N/3),M)*depth; % in units of 1/µm^3
D=zeros(N,M);

u2=zeros(N,M); 
%u2(:,1:round(M/2))=u2(:,1:round(M/2))+3000;
%D(:,1:round(M/2))=D(:,1:round(M/2))-3000;
%E(:,1:round(M/2))=E(:,1:round(M/2))-3000;

u3=zeros(N,M); 

%load BonnyIC2m2c.m
u4=zeros(N,M); 
u5=zeros(N,M); 
u6=zeros(N,M); 
u7=zeros(N,M); 

%% LOOP

for k=1:nsteps

    RHS=reactions(model,D(:)/depth,E(:)/depth,u1(:),u2(:),u3(:),u4(:),u5(:),u6(:),u7(:));
    
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
    Dmemtime2(k)=u1(N,M);
    DEmemtime(k)=u2(1,1);
    Ememtime(k)=u3(1,1);
    
    Dmemsum=mean(u1')+mean(u2');
    Ememsum=mean(u2')+mean(u3');
    Dmemsumtime(k,:)=Dmemsum;
    Ememsumtime(k,:)=Ememsum;

    % Plot
    if mod(k,plotEvery)==0
        maxD = max(max(u1+u2))+eps;
        maxE = max(max(u3+u2))+eps;
        maxDE = max(max(u2))+eps;
        maxDcyt = max(max(D))/depth+eps;
        maxEcyt = max(max(E))/depth+eps;

        subplot(2,3,1);
        imagesc(Xgrid',Ygrid',D'/depth,[0,maxDcyt]);
        axis equal;
        axis('xy');
        xlabel('x (µm)');
        ylabel('y (µm)');
        title(['MinD (cyt)   t = ' int2str(k*dt) ' s']);
        colorbar('southoutside');
        daspect([1 1 1])

        subplot(2,3,2);
        imagesc(Xgrid',Ygrid',E'/depth,[0,maxEcyt]);
        axis equal;
        axis('xy');
        xlabel('x (µm)');
        ylabel('y (µm)');
        title('MinE (cyt)');
        colorbar('southoutside');
        daspect([1 1 1])

        subplot(2,3,4);
        imagesc(Xgrid',Ygrid',u1'+u2',[0,maxD]);
        axis equal;
        axis('xy');
        xlabel('x (µm)');
        ylabel('y (µm)');
        title('MinD (mem)');
        colorbar('southoutside');
        daspect([1 1 1])

        subplot(2,3,5);
        imagesc(Xgrid',Ygrid',u2',[0,maxDE]);
        axis equal;
        axis('xy');
        xlabel('x (µm)');
        ylabel('y (µm)');
        title('MinDE (mem)');
        colorbar('peer', gca,'southoutside');
        daspect([1 1 1])

        subplot(2,3,6);
        imagesc(Xgrid',Ygrid',u3'+u2',[0,maxE]);
        axis equal;
        axis('xy');
        xlabel('x (µm)');
        ylabel('y (µm)');
        title('MinE (mem)');
        colorbar('southoutside');
        daspect([1 1 1])
        
        subplot(2,3,3);
        plot(time,Dcyttime,time,Ecyttime,time,Dmemtime,time,DEmemtime,time,Ememtime);
        legend('MinDc','MinEc','MinDm','MinDEm','MinEm');
        legend('boxoff');
        xlabel('time (s)');
        ylabel('concentration (µm^{-2})');

        pause(0.001);

    endif

endfor

%% plot
maxD=max(max(Dmemsumtime));
maxE=max(max(Ememsumtime));
subplot(2,1,1);
imagesc(time,Xgrid',Dmemsumtime',[0,maxD]);
xlabel('Time (s)');
ylabel('Position (µm)');
title('Kymograph for MinD on Membrane');
axis('xy');
%colorbar('eastoutside');

subplot(2,1,2);
imagesc(time,Xgrid',Ememsumtime',[0,maxE]);
xlabel('Time (s)');
ylabel('Position (µm)');
title('Kymograph for MinE on Membrane');
axis('xy');
%colorbar('eastoutside');

%save BonnyIC2m2c.m D E u1 u2 u3;