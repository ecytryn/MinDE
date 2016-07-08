figure(1)

clf;
clear;
colormap(gray);

model = 'Bonny in vitro';

%% Numerical parameters
L=100;  % domain length
W=100;  % domain width
depth=90; % domain depth
dx=5;
Xgrid=(0:dx:L)';
Ygrid=(0:dx:W)';
N=length(Xgrid);
M=length(Ygrid);
NM=N*M;
nsteps=2000;
dt=0.01;

%% Plot params
plotEvery=100;

time=[1:nsteps]*dt;
Dcyttime = zeros(size(time));
Ecyttime = zeros(size(time));
Dmemtime = zeros(size(time));
DEmemtime = zeros(size(time));
Ememtime = zeros(size(time));

%% Diffusion constants (Bonny 2013)
switch model
    case 'Bonny in vitro'
        % in vitro
        Dmem=0.3; % MinD membrane diffusion in µm^2/s
        Emem=1.8; % MinE membrane diffusion in µm^2/s
        Dcyt=50; % cyt/bulk diffusion in µm^2/s
    case 'Bonny in vivo'
        % in vivo
        Dmem=0.06; % MinD membrane diffusion in µm^2/s
        Emem=0.3; % MinE membrane diffusion in µm^2/s
        Dcyt=14; % cyt/bulk diffusion in µm^2/s
        
end

%% A - diffusion matric, no flux
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
    case 'Bonny in vitro'
        % in vitro
        D=4.8e2*ones(N,M)*depth; % in units of 1/µm^3
        E=7.0e2*ones(N,M)*depth; % in units of 1/µm^3
    case 'Bonny in vivo'
        % in vivo
        D=1.4e3*ones(N,M)*depth; % in units of 1/µm^3
        E=9.7e2*ones(N,M)*depth; % in units of 1/µm^3
        
end

u1=zeros(N,M); % concentration of membrane species 1 in units of molecules per µm^2
u1(1:round(N/2),1:end)=u1(1:round(N/2),1:end)+10000;
D(1:round(N/2),1:end)=D(1:round(N/2),1:end)-10000;
u2=zeros(N,M); 
u2(1:end,1:round(M/2))=u2(1:end,1:round(M/2))+3000;
D(1:end,1:round(M/2))=D(1:end,1:round(M/2))-3000;
E(1:end,1:round(M/2))=E(1:end,1:round(M/2))-3000;
u3=zeros(N,M); 

load BonnyIC.m
%u3(1:end,1:round(M/2))=u3(1:end,1:round(M/2))+200;
%E(1:end,1:round(M/2))=E(1:end,1:round(M/2))-200;
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
    u4(:) = (Ident-Dmem*A) \ ( (Ident+Dmem*A)*u4(:) + dt*RHS(:,6) );
    u5(:) = (Ident-Dmem*A) \ ( (Ident+Dmem*A)*u5(:) + dt*RHS(:,7) );
    u6(:) = (Ident-Dmem*A) \ ( (Ident+Dmem*A)*u6(:) + dt*RHS(:,8) );
    u7(:) = (Ident-Dmem*A) \ ( (Ident+Dmem*A)*u7(:) + dt*RHS(:,9) );

    Dcyttime(k)=D(1,1);
    Ecyttime(k)=E(1,1);
    Dmemtime(k)=u1(1,1);
    DEmemtime(k)=u2(1,1);
    Ememtime(k)=u3(1,1);

    % Plot
    if mod(k,plotEvery)==0

        subplot(2,3,1);
        imagesc(Xgrid',Ygrid',D',[0,45000]);
        axis equal;
        title(['MinD (cyt)   t = ' int2str(k*dt) ' s']);

        subplot(2,3,2);
        imagesc(Xgrid',Ygrid',E',[0,65000]);
        axis equal;
        title('MinE (cyt)');

        subplot(2,3,4);
        imagesc(Xgrid',Ygrid',u1',[0,12000]);
        axis equal;
        title('MinD (mem)');

        subplot(2,3,5);
        imagesc(Xgrid',Ygrid',u2',[0,5000]);
        axis equal;
        title('MinDE (mem)');

        subplot(2,3,6);
        imagesc(Xgrid',Ygrid',u3',[0,3000]);
        axis equal;
        title('MinE (mem)');
        
        subplot(2,3,3);
        plot(time,Dcyttime,time,Ecyttime,time,Dmemtime,time,DEmemtime,time,Ememtime);
        legend('MinDc','MinEc','MinDm','MinDEm','MinEm');
        %plot(time,Dmemtime,time,DEmemtime,time,Ememtime);
        %legend('MinDm','MinDEm','MinEm');
        legend('boxoff');

        pause(0.001);

    end
    

end

%save BonnyIC2.m D E u1 u2 u3;