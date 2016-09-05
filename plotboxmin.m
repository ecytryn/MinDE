figure(1);

clf;
clear;
colormap(jet);

L=2;
model="bonny_invivo";
%model="bonny_invitro";
%model="petrasek_Me2";

z1=load ([model "_D_" int2str(L) "out.txt"]); % MinD
z2=load ([model "_E_" int2str(L) "out.txt"]); % MinE
z3=load ([model "_DE_" int2str(L) "out.txt"]); % MinDE

t=z1(:,1);
y1=z1(:,2:end)+z3(:,2:end); % MinD total
y2=z2(:,2:end)+z3(:,2:end); % MinE total
h=linspace(0,L,length(y1(1,:)))';
maxy1=max(max(y1(round(length(y1)/2):end,:)));
maxy2=max(max(y2(round(length(y2)/2):end,:)));
#maxy1=max(max(y1));
#maxy2=max(max(y2));
maxt=max(t);

%% plot kymograph
set(0,'DefaultAxesFontSize',12 );
set(0,'DefaultAxesFontName','Arial');
subplot(2,1,1);
imagesc(t',[0 L],y1',[0,maxy1]);
xlabel('Time (s)');
ylabel('Position (\mu{m})');
title('Kymograph for MinD on Membrane');
axis('xy');
%colorbar('eastoutside');

subplot(2,1,2);
imagesc(t',[0 L],y2',[0,maxy2]);
xlabel('Time (s)');
ylabel('Position (\mu{m})');
title('Kymograph for MinE on Membrane');
axis('xy');
%colorbar('eastoutside');

%% save plot
orient("portrait");
fname=[model "_length" int2str(L)];
%print([fname ".pdf"]);