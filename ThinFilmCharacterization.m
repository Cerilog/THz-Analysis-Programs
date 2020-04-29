%% 
% Iterative Calculation of Complex Refractive Index Extraction
% For thin films deposit on a thick substrate
% No echoes are allowed
% Remember to check n3 at line 65, the R.I. of Substrate
% Revised on 2020/04/10
%% Basic I/O
clear;
Dsub = fscanf(fopen('Reference-SIGaAs-1-1024-10um-cut.txt','r'),'%g %g',[2,inf]);
Dsam = fscanf(fopen('SIGaAs-RTA-1e14-1-1024-10um-cut.txt','r'),'%g %g',[2,inf]);
fileID = fopen('ThinFilm-RIX-RTA-1e14-1.txt','w');
d = 600e-9; % in m, for thin film
delta_d = -1e-6; % Thickness of (Sample - Substrate) in m
%% Frequency-Domain Analysis
t = Dsub(1,:)*1e-12;
Esub = Dsub(2,:);
Esam = Dsam(2,:);
N = length(t);
dt = t(2)-t(1);
fs = 1/dt;
f = linspace(-fs/2,fs/2,N);
fTHz = f/1e+12;
w = 2*pi*f;
Fsub = fftshift(fft(Esub));
Fsam = fftshift(fft(Esam));

figure(1);
plot(fTHz,db((abs(Fsub).^2)),fTHz,db((abs(Fsam).^2)),'Linewidth',0.9);
grid on
xlabel('Frequency(THz)');
ylabel('Spectral Power');
title('Spectrums');
legend('Substrate','Sample');
axis([0,5,-inf,inf]);

H = Fsam./Fsub;
modulus = log(abs(H));
argument = angle(Fsam./Fsub);
phaseshift = unwrap(argument);
%% Extrapolation to correct phase
p1 = find(fTHz>=0.5,1);
p2 = find(fTHz>=1.0,1);
m = (phaseshift(p2)-phaseshift(p1))/(fTHz(p2)-fTHz(p1));
b = phaseshift(p1) - m*fTHz(p1);
for i=1:N
    phaseshift(i) = phaseshift(i) - b;
end    
%% Frequency Domain Plots
figure(2);
sgtitle('Complex Transmittance');
subplot(2,1,1);
plot(fTHz,modulus,'Linewidth',0.9);
title('Modulus');
xlabel('Frequency(THz)');
ylabel('Transmittance');
axis([0,2.5,-inf,inf]);
subplot(2,1,2);
plot(fTHz,phaseshift,'Linewidth',0.9);
title('Phaseshift');
xlabel('Frequency(THz)');
ylabel('Phaseshift');
axis([0,2.5,-inf,inf]);
%% Model Description and Parameter Extraction
n1 = 1+0*i; % Refractive index of free space
n3 = 3.584+5e-5*i; % Refractive index of substrate
c0 = 299792458; % in (m/s)
n2_ini = [25,20];
n2 = zeros(N,1);

for i=1:N
    f = @(z) ThinFilmTheoTrans(n1,z,n3,d,delta_d,w(i)/c0) - H(i);
    c = @(x) complex(x(1),x(2));
    g = @(x) abs(f(c(x)));
    n2(i) = c(fminsearch(g,n2_ini));
end

figure(3);
sgtitle('Complex Refractive Index of Thin Film');
subplot(2,1,1);
plot(fTHz,real(n2),'Linewidth',0.9);
xlabel('Frequency(THz)');
ylabel('n');
title('Real Part');
axis([0,1.5,-inf,inf]);
subplot(2,1,2);
plot(fTHz,imag(n2),'Linewidth',0.9);
xlabel('Frequency(THz)');
ylabel('\kappa');
title('Imag Part');
axis([0,1.5,-inf,inf]);
%% Residue of Error Function Contour Plot at 0.5THz 
err = zeros(1000);
p = find(fTHz>=0.5,1);
k0_err = w(p)/c0;
n2_err = linspace(0.05*real(n2(p)),1.95*real(n2(p)),1000);
k2_err = linspace(-1*imag(n2(p)),3*imag(n2(p)),1000);
for ii=1:1000
    for jj=1:1000
        y = n2_err(ii) + 1i*k2_err(jj);
        err(ii,jj) = abs(ThinFilmTheoTrans(n1,y,n3,d,delta_d,k0_err) - H(p));
    end
end

figure(4);
contourf(n2_err,k2_err,err,50);
title('Error Function Contour at 500 GHz');
xlabel('Variation of Real Part');
ylabel('Variation of Imag Part');
colormap (jet(51))
colorbar
%% Write out to txt File
ReliableRange = and(fTHz>=0.2,fTHz<=1.7);
fTHz = fTHz(ReliableRange);
n2 = n2(ReliableRange);
O(1,:) = fTHz;
O(2,:) = real(n2);
O(3,:) = imag(n2);
formatSpec = '%g %g %g\n';
fprintf(fileID,formatSpec,O);
fclose('all');
%% Function body (Do not modify)
function Ttheo = ThinFilmTheoTrans(n1,n2,n3,d,deltad,k0)

    numer = 2*n2*(n1+n3)*exp(1i*(n2-1)*k0*d);
    denom = ((n2+n3)*(n2+n1)-(n2-n3)*(n2-n1)*exp(1i*2*n2*k0*d));
    thickdiff = exp(1i*(n3-1)*k0*deltad);
   
    Ttheo = numer./denom./thickdiff;
end    