%%
% Provide Drude model Fitting with n,k data
% By comparing complex conductivity with classical EM theory
% Date : 2020/03/03
% Check epi_inf at line 12
% Initial guess at line 47
% Revised : 2020/04/10
%% Data I/O and Basic Treatment
D = fscanf(fopen('ThinFilm-RIX-RTA-1e14.txt','r'),'%f %f %f',[3,inf]);
global epi0 
epi0 = 8.85e-12;
epi_inf = 10.89;
fTHz = D(1,:);
n = D(2,:);
k = D(3,:);
Range = and(fTHz>=0.2,fTHz<=1.2);
fTHz = fTHz(Range);
n = n(Range);
k = k(Range);
f = fTHz*1e+12;
w = 2*pi*f;

figure(1);
title('Complex Refractive Index');
yyaxis left
plot(fTHz,n,'Linewidth',0.9);
xlabel('Frequency(THz)');
ylabel('Real Refractive Index');
yyaxis right
plot(fTHz,k,'Linewidth',0.9);
ylabel('Imag Refractive Index');
%% Complex Conductivity from Classical EM Wave Theory
ReCond_EM = 2*n.*k.*w*epi0;
ImCond_EM = epi0*w.*(epi_inf-n.^2+k.^2);
ComplexCond_EM = complex(ReCond_EM,ImCond_EM);

figure(2);
title('Complex Conductivity');
yyaxis left
plot(fTHz,ReCond_EM,'Linewidth',0.9);
xlabel('Frequency(THz)');
ylabel('Real Conductivity');
yyaxis right
plot(fTHz,ImCond_EM,'Linewidth',0.9);
ylabel('Imag Conductivity');
%% Drude Model Fitting
iniguess = [3e+14,10e-15];
% Initial Guess for plasma freq 
g = @(x) sum(abs(CondDrude(epi0,w,x(1),x(2)) - ComplexCond_EM));
[result,residue] = fminsearch(g,iniguess);
wp = result(1);
tau = result(2);
fprintf('Plasma Frequency = %g (rad*THz)\n',wp/1e+12);
fprintf('Collision Time = %g (fs)\n',tau*1e+15);
%% Comparison between Experimental Result and Drude Model Fitting
ReCondFit = real(CondDrude(epi0,w,wp,tau));
ImCondFit = imag(CondDrude(epi0,w,wp,tau));
CondFit = complex(ReCondFit,ImCondFit);

figure(3);
sgtitle('Comparsion of Experiment and Drude Fitting');
subplot(2,1,1);
plot(fTHz,ReCond_EM,'Linewidth',0.9);
xlabel('Frequency(THz)');
ylabel('Real Conductivity');
hold on
plot(fTHz,ReCondFit,'o');
hold off
legend('Experimental','Drude Fitting');
subplot(2,1,2);
plot(fTHz,ImCond_EM,'Linewidth',0.9);
xlabel('Frequency(THz)');
ylabel('Imag Conductivity');
hold on
plot(fTHz,ImCondFit,'o');
hold off
legend('Experimental','Drude Fitting');
%% Error Contour
% wp_err = linspace(0.5*wp,1.5*wp,1000);
% tau_err = linspace(0.5*tau,1.5*tau,1000);
% err = zeros(1000);
% 
% for a = 1:1000
%     for b = 1:1000
%         err(a,b) = abs(g([wp_err(a),tau_err(b)]));
%     end
% end
% 
% figure(4);
% contourf(wp_err,tau_err,err,50);
% title('Residue Contour of Drude Fitting');
% xlabel('Variation of Plasma Frequency');
% ylabel('Variation of Collision Time');
% colormap (jet(51))
% colorbar
%% Derivation of Carrier Concentration, Mobility and DC Conductivity
q = 1.6e-19;
m0 = 9.1e-31;
effm = 0.063*m0;
Ne = epi0*wp^2*effm/q^2;
mobility = q*tau/effm;
DCCond = epi0*wp^2*tau;
fprintf('Carrier Concentration = %g (cm^-3)\n',Ne/1e+6);
fprintf('Mobility = %g (cm^2/V/s)\n',mobility*1e+4);
fprintf('DC Conductivity = %g (S/cm)\n',DCCond/100);
%% Generation of Comparison Plots of Complex RIX
epi = epi_inf + 1i*CondFit./w/epi0;
n_pre = real(sqrt(epi));
k_pre = imag(sqrt(epi));

figure(4);
sgtitle('Comparsion of Refractive Index');
subplot(2,1,1);
plot(fTHz,n,'o');
xlabel('Frequency(THz)');
ylabel('n');
title('Real Refractive Index');
hold on
plot(fTHz,n_pre,'Linewidth',0.9);
legend('Experimental','Predicted');
hold off
subplot(2,1,2);
plot(fTHz,k,'o');
xlabel('Frequency(THz)');
ylabel('k');
title('Imag Refractive Index');
hold on
plot(fTHz,k_pre,'Linewidth',0.9);
legend('Experimental','Predicted');
hold off
fclose('all');
%% Function Body (Do Not Modify)
function CompCond = CondDrude(epi0,w,wp,tau)
    numer = epi0*wp^2*tau;
    denom = 1-1i*w*tau;
    CompCond = numer./denom;
end    