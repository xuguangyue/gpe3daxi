clear
clc

dyn=load('./amp=0.10_f=4.286/realaxi-dyna.txt');
t = dyn(:,1);
nt = length(t);
dt = t(2)-t(1);
zzt = dyn(:,4);
dxt = dyn(:,5);
dzt = dyn(:,6);

dxw = fftshift(fft(dxt - mean(dxt)));
zzw = fftshift(fft(zzt - mean(zzt)));
if ~mod(nt,2)
    w = 2*pi/(nt*dt)*(-nt/2:nt/2-1);
else
    w = 2*pi/(nt*dt)*(-(nt-1)/2:(nt-1)/2);
end

subplot(211)
plot(t,dxt/sqrt(2),t,dzt)
xlabel('\omega_0 t')
ylabel('\Deltax/l_0')
subplot(212)
plot(w/sqrt(5),abs(dxw))
xlabel('$\omega/(\sqrt{5}\omega_0)$','interpreter','latex')
ylabel('\Deltax(\omega)')
xlim([-3 3])
grid on
