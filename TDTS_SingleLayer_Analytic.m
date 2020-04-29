%%
% Analytical Extraction of Optical Constant in TDTS
% Single layered, neglecting multi-reflection
% Phaseshift corrected on 2019/12/09
% If the refractive index is negative, modify '-' to '+' in line 76
% Revised on 2020/04/10
%% Analysis Parameters
clear all;
Dref = fscanf(fopen('Purged-1-1024-10um-cut.txt','r'),'%f %f',[2,inf]);
Dsam = fscanf(fopen('Reference-SIGaAs-1-1024-10um-cut.txt','r'),'%f %f',[2,inf]);
fileID = fopen('RIX-SIGaAs.txt','w');  % Output file for complex RIX
d = 600e-6; % in m
%% Frequency Domain Analysis
t = Dref(1,:)*1e-12;
N = length(t);
dt = t(2)-t(1);
fs = 1/dt;
f = linspace(-0.5*fs,0.5*fs,N);
% f = f(f>0);
w = f*2*pi;
fTHz = f/1e+12;
Eref = Dref(2,:);
Esam = Dsam(2,:);
Fref = fftshift(fft(Eref));
Fsam = fftshift(fft(Esam));
%% Spline Interpolation to Smoothen Spectrums
% fTHz2 = fTHz + 0.5*(fTHz(2)-fTHz(1));
% Fref_db = db(abs(Fref).^2);
% Fsam_db = db(abs(Fsam).^2);
% Fref_db_spline = spline(fTHz,Fref_db,fTHz2);
% Fsam_db_spline = spline(fTHz,Fsam_db,fTHz2);
% figure(5);
% plot(fTHz,Fref_db,fTHz2,Fref_db_spline);
% xlabel('Frequnecy(THz)');
% title('Spectrums');
% axis([0,6.5,-inf,inf]);
%% Extract Modulus and Phase (unwrapped)
H = Fsam./Fref;
A = abs(H);
phaseshift = unwrap(angle(H));
%% Extrapolation to correct phase
ReliableRange = and(fTHz>0.2,fTHz<1.2);
p1 = find(fTHz>=0.5,1);
p2 = find(fTHz>=1.0,1);
m = (phaseshift(p2)-phaseshift(p1))/(fTHz(p2)-fTHz(p1));
b = phaseshift(p1) - m*fTHz(p1);
for i=1:N
    phaseshift(i) = phaseshift(i) - b;
end    
%% Plots for Verification
figure(1);
plot(t*1e+12,Eref,t*1e+12,Esam,'Linewidth',0.9);
legend('Reference Waveform','Sampled Waveform');
xlabel('Time Delay(ps)');
ylabel('Amplitude');
grid on;
title('Time Domain Plot');
figure(2);
subplot(2,1,1);
plot(fTHz,db(abs(Fref).^2),fTHz,db(abs(Fsam).^2),'Linewidth',0.9);
xlabel('Frequency(THz)');
ylabel('dB');
grid on;
title('Reference and Sample Spectrum');
legend('Reference','Sample');
axis([0,7,-inf,inf]);
subplot(2,1,2);
plot(fTHz,phaseshift,'Linewidth',0.9);
xlabel('Frequency');
ylabel('Phase');
grid on;
axis([-0.1,4,-inf,inf]);
title('Phase Shift');
%% Analytical Expression of Complex Refractive Index
c0 = 299792458; % in m/s
n = 1-c0./w.*(phaseshift/d);
alpha = log((4*n)./A./(1+n).^2)*(2/d); % in m^-1
alpha = alpha/100;
k = c0*alpha./(4*pi*f);
%% Plots of Complex Refractive Index
figure(3);
sgtitle('Complex Refractive Index');
subplot(2,1,1);
plot(fTHz,n,'Linewidth',0.8);
xlabel('Frequency(THz)');
ylabel('n');
title('Real Refractive Index');
axis([0.2,2.2,3,4]);
subplot(2,1,2);
plot(fTHz,k,'Linewidth',0.8);
xlabel('Frequency(THz)');
ylabel('\kappa');
title('Imaginary Refractive Index');
axis([0.3,2,-2e-3,2e-3]);
%% Write out to txt File
fTHz = fTHz(ReliableRange);
f = f(ReliableRange);
n = n(ReliableRange);
k = k(ReliableRange);
O(1,:) = fTHz;
O(2,:) = n;
O(3,:) = k;
formatSpec = '%g %g %g\n';
fprintf(fileID,formatSpec,O);
fclose(fileID);