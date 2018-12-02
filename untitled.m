clear
clc

dyn=load('./realaxi-dyna.txt');
t = dyn(:,1);
nt = length(t);
dt = t(2)-t(1);
zzt = dyn(:,4);
dxt = dyn(:,5);
dzt = dyn(:,6);

dxw = fftshift(fft(dxt - mean(dxt)));
dzw = fftshift(fft(dzt - mean(dzt)));
zzw = fftshift(fft(zzt - mean(zzt)));
if ~mod(nt,2)
    w = 2*pi/(nt*dt)*(-nt/2:nt/2-1);
else
    w = 2*pi/(nt*dt)*(-(nt-1)/2:(nt-1)/2);
end

figure
subplot(211)
plot(t,dxt/sqrt(2),t,dzt)
xlabel('\omega_0 t')
ylabel('\Deltax/l_0')
subplot(212)
plot(w/sqrt(5),abs(dxw),w/sqrt(5),abs(dzw))
xlabel('$\omega/(\sqrt{5}\omega_0)$','interpreter','latex')
ylabel('\Deltax(\omega)')
xlim([-4 4])
grid on
