%%
% To Process terahertz waveform and smoothen the tails
% Version 2.0 Support pre-cutting also
%%
[filename,path] = uigetfile('*.txt','Select a THz Waveform with Echoes');
D = fscanf(fopen(filename,'r'),'%f %f',[2,inf]);
t = D(1,:);
E = D(2,:);
N = length(t);

figure(1);
subplot(2,1,1);
plot(t,E,'Linewidth',0.9);
xlabel('Time');
ylabel('Amplitude');
title('Waveform before Processing');
subplot(2,1,2);
plot(E,'Linewidth',0.9);
grid on
xlabel('Point #');
ylabel('Amplitude');

Prompt1 = 'Enter the start of the THz peak ';
p1 = input(Prompt1);
Prompt2 = 'Enter the end of the THz peak ';
p2 = input(Prompt2);

for k=1:N
    if (k<=p1 | k>=p2)
        E(k) = 0;
    end    
end

figure(2);
subplot(2,1,1);
plot(t,E,'Linewidth',0.9);
xlabel('Time');
ylabel('Amplitude');
title('Waveform after Processing');
subplot(2,1,2);
plot(E,'Linewidth',0.9);
xlabel('Point #');
ylabel('Amplitude');

D2 = zeros(2,N);
for i=1:N
    D2(1,i) = t(i);
    D2(2,i) = E(i);
end
%%
fs = 1/(t(2)-t(1))*1e+12;
f = linspace(-0.5*fs,0.5*fs,N);
F = fftshift(fft(E));
fTHz = f/1e+12;

figure(3);
plot(fTHz,db(abs(F).^2),'Linewidth',0.9);
xlabel('Frequency(THz)');
ylabel('Spectral Power(dB)');
title('Spectrum after Cut');
axis([0,5,-inf,inf]);
%%
filename = strsplit(filename,'.');
filename = string(append(filename(1),'-cut.',filename(2)));
fileID = fopen(filename,'w');
formatSpec = '%g %g\n';
fprintf(fileID,formatSpec,D2);
fclose(fileID);