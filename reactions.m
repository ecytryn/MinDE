function RHS=reactions(model,D,E,u1,u2,u3,u4,u5,u6,u7)

RHS=zeros(length(u1),9);

switch model
    case 'Bonny in vitro' % Bonny et al. 2013
        % Kinetic parameters
        Dmax=2.75e4; % µm^(-2)
        wD=5e-4;     % µm/s
        wdD=3.18e-3; % µm^3/s
        wE=1.36e-4;  % µm^3/s
        wed=4.9e-3;  % µm^2/s
        wdec=0.16;   % 1/s
        wdem=2.52;   % 1/s
        we=0.5;      %1/s
        
        RHS(:,1)= - D.*(wD + wdD*u1).*(Dmax-u1-u2)/Dmax + wdec*u2 +wdem*u2; % D in cyt/bulk
        RHS(:,2)= - wE*E.*u1 + wdec*u2 +we*u3; % E in cyt/bulk
        RHS(:,3)= D.*(wD + wdD*u1).*(Dmax-u1-u2)/Dmax - wE*E.*u1 - wed*u1.*u3 ; % D on membrane
        RHS(:,4)= wE*E.*u1 - wdec*u2 - wdem*u2 + wed*u1.*u3; % DE on membrane
        RHS(:,5)= wdem*u2 - wed*u1.*u3 - we*u3; % E on membrane
        
    case 'Bonny in vivo' % Bonny et al. 2013
        % Kinetic parameters
        Dmax=5.4e3; % µm^(-2)
        wD=0.1;     % µm/s
        wdD=8.8e-3; % µm^3/s
        wE=6.96e-5;  % µm^3/s
        wed=0.139;  % µm^2/s
        wdec=0.08;   % 1/s
        wdem=1.5;   % 1/s
        we=0.5;      %1/s
        
        RHS(:,1)= - D.*(wD + wdD*u1).*(Dmax-u1-u2)/Dmax + wdec*u2 +wdem*u2; % D in cyt/bulk
        RHS(:,2)= - wE*E.*u1 + wdec*u2 +we*u3; % E in cyt/bulk
        RHS(:,3)= D.*(wD + wdD*u1).*(Dmax-u1-u2)/Dmax - wE*E.*u1 - wed*u1.*u3 ; % D on membrane
        RHS(:,4)= wE*E.*u1 - wdec*u2 - wdem*u2 + wed*u1.*u3; % DE on membrane
        RHS(:,5)= wdem*u2 - wed*u1.*u3 - we*u3; % E on membrane
        
    case 'PetrasekMe2'
        
end