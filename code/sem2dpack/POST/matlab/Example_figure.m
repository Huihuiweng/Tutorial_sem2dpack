% Modify the path of the output files
datadir='/home/weng/Works/Softwares/sem2dpack/EXAMPLES/Thermpres_SWF';

grid = sem2d_read_specgrid(datadir);
data = sem2d_read_fault('Flt01');

figure(1)
subplot(2,2,1)
% Slip rate image
imagesc((0:data.nt-1)*data.dt,data.x/1e3,data.v);
xlabel('Time (s)');
ylabel('Along strike distance (km)');
clim([0,20]);
gca = colorbar;
ylabel(gca,'Slip rate (m/s)');

subplot(2,2,2)
% Slip rate image
imagesc((0:data.nt-1)*data.dt,data.x/1e3,(data.st0+data.st)/1e6);
xlabel('Time (s)');
ylabel('Along strike distance (km)');
%clim([0,20]);
gca = colorbar;
ylabel(gca,'Shear stress (MPa)');

subplot(2,2,3)
p1=int32(data.nx/4.0);
p2=int32(data.nx/2.0);
p3=int32(3*data.nx/4.0);
plot((0:data.nt-1)*data.dt,data.v(p1,:));
hold on
plot((0:data.nt-1)*data.dt,data.v(p2,:));
plot((0:data.nt-1)*data.dt,data.v(p3,:));
hold off
xlim([0 inf]);
ylim([-inf inf]);
xlabel('Time (s)');
ylabel('Slip rate (m/s)');

subplot(2,2,4)
p1=int32(data.nx/4.0);
p2=int32(data.nx/2.0);
p3=int32(3*data.nx/4.0);
plot((0:data.nt-1)*data.dt,(data.st(p1,:)+data.st0(p1,:))/1e6);
hold on
plot((0:data.nt-1)*data.dt,(data.st(p2,:)+data.st0(p2,:))/1e6);
plot((0:data.nt-1)*data.dt,(data.st(p3,:)+data.st0(p3,:))/1e6);

hold off
xlim([0 inf]);
ylim([0 30]);
xlabel('Time (s)');
ylabel('Shear stress (MPa)');


%%  if tp is used
if ~isempty(data.P)
    figure(2)
    subplot(2,2,1)
    % slip rate image
    imagesc((0:data.nt-1)*data.dt,data.x/1e3,data.P/1e6);
    xlabel('Time (s)');
    ylabel('Along strike distance (km)');
    gca = colorbar;
    ylabel(gca,'Pore pressure (MPa)');

    subplot(2,2,2)
    % slip rate image
    imagesc((0:data.nt-1)*data.dt,data.x/1e3,data.T);
    xlabel('Time (s)');
    ylabel('Along strike distance (km)');
    gca = colorbar;
    ylabel(gca,'Temperature (K)');

    subplot(2,2,3)
    p1=int32(data.nx/4.0);
    p2=int32(data.nx/2.0);
    p3=int32(3*data.nx/4.0);
    plot((0:data.nt-1)*data.dt,data.P(p1,:));
    hold on
    plot((0:data.nt-1)*data.dt,data.P(p2,:));
    plot((0:data.nt-1)*data.dt,data.P(p3,:));
    hold off
    xlim([0 inf]);
    ylim([-inf inf]);
    xlabel('Time (s)');
    ylabel('Pore pressure (MPa)');

    subplot(2,2,4)
    p1=int32(data.nx/4.0);
    p2=int32(data.nx/2.0);
    p3=int32(3*data.nx/4.0);
    plot((0:data.nt-1)*data.dt,data.T(p1,:));
    hold on
    plot((0:data.nt-1)*data.dt,data.T(p2,:));
    plot((0:data.nt-1)*data.dt,data.T(p3,:));

    hold off
    xlim([0 inf]);
    xlabel('Time (s)');
    ylabel('Temperature (K)');
end