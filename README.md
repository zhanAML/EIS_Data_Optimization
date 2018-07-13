# EIS_Data_Optimization
Nonlinear least square method especially Gauss-Newton method.


clc;
close all;
clear;


%       Real      Imag    Frequency .
Out = [1.21E-02	-4.06E-03	1.01E+04
      1.21E-02	-2.73E-03	8.02E+03
      1.23E-02	-1.55E-03	6.33E+03
      1.24E-02	-5.24E-04	5.02E+03
      1.27E-02	3.99E-04	3.98E+03
      1.30E-02	1.24E-03	3.17E+03
      1.33E-02	2.02E-03	2.53E+03
      1.38E-02	2.87E-03	1.98E+03
      1.42E-02	3.63E-03	1.58E+03
      1.48E-02	4.35E-03	1.27E+03
      1.54E-02	5.19E-03	9.98E+02
      1.61E-02	6.00E-03	7.97E+02
      1.70E-02	6.91E-03	6.28E+02
      1.79E-02	7.80E-03	5.06E+02
      1.90E-02	8.84E-03	3.98E+02
      2.02E-02	9.91E-03	3.16E+02
      2.16E-02	1.10E-02	2.52E+02
      2.31E-02	1.23E-02	1.99E+02
      2.48E-02	1.36E-02	1.58E+02
      2.67E-02	1.50E-02	1.26E+02
      2.86E-02	1.64E-02	1.00E+02
      3.08E-02	1.80E-02	7.90E+01
      3.28E-02	1.96E-02	6.33E+01
      3.52E-02	2.17E-02	4.99E+01
      3.75E-02	2.40E-02	3.97E+01
      4.01E-02	2.69E-02	3.17E+01
      4.30E-02	3.03E-02	2.49E+01
      4.73E-02	3.42E-02	1.99E+01
      5.05E-02	3.91E-02	1.58E+01
      5.52E-02	4.52E-02	1.24E+01
      6.06E-02	5.15E-02	9.93E+00
      6.58E-02	5.90E-02	7.95E+00
      7.33E-02	6.76E-02	6.32E+00
      8.21E-02	7.74E-02	5.01E+00
      9.29E-02	8.86E-02	3.95E+00
      1.05E-01	9.98E-02	3.16E+00
      1.20E-01	1.11E-01	2.50E+00
      1.38E-01	1.24E-01	2.00E+00
      1.59E-01	1.36E-01	1.59E+00
      1.83E-01	1.47E-01	1.27E+00
      2.11E-01	1.55E-01	9.99E-01
      2.41E-01	1.61E-01	7.92E-01
      2.70E-01	1.63E-01	6.33E-01
      2.99E-01	1.61E-01	5.04E-01
      3.27E-01	1.56E-01	4.01E-01
      3.53E-01	1.49E-01	3.17E-01
      3.76E-01	1.41E-01	2.52E-01
      3.96E-01	1.33E-01	2.00E-01
      4.14E-01	1.26E-01	1.59E-01
      4.30E-01	1.19E-01	1.26E-01
      4.44E-01	1.14E-01	1.00E-01
      4.58E-01	1.11E-01	7.95E-02
      4.70E-01	1.09E-01	6.32E-02
      4.83E-01	1.09E-01	5.01E-02
      4.97E-01	1.11E-01	3.98E-02
      5.11E-01	1.14E-01	3.16E-02
      5.25E-01	1.19E-01	2.51E-02
      5.42E-01	1.25E-01	2.00E-02
      5.59E-01	1.32E-01	1.59E-02
      5.77E-01	1.41E-01	1.26E-02
      5.97E-01	1.51E-01	1.00E-02]


Zreal_50mA_Temp_30    =    Out(:,1);
Zimag_50mA_Temp_30    =    Out(:,2);
w = Out(:,3);

syms Re Rl Rt Qdl Ql alphadl alphal Wb;
Beta = [Re Rl Rt Qdl Ql alphadl alphal Wb]; 
Betak = [0.01256 0.06 0.43 0.32 0.4 0.86 0.644 40];
for k = 1:1:10 %iteration
  
    Zdl = 1./((w*1i).^Beta(6)*Beta(4)); % Zdl = 1./((w*1i).^alphadl*Qdl);
    Zl = 1./((w*1i).^Beta(7)*Beta(5));  % Zl = 1./((w*1i).^alphal*Ql);
    Zwb = coth(sqrt(1i.*w)*Beta(8))./(sqrt(1i.*w)*Beta(8)); % Zwb = coth(sqrt(1i.*w)*Wb)./(sqrt(1i.*w)*Wb);

    Z0 = Beta(2)+((Beta(3)+Zwb).*Zdl)./(Beta(3)+Zwb+Zdl); % Z0 = Rl+((Rt+Zwb).*Zdl)./(Rt+Zwb+Zdl);
    Z = Beta(1)+(Zl.*Z0)./(Zl+Z0); % Z = Re+(Zl.*Z0)./(Zl+Z0);
    
    Zreal = real(Z);
    Zimag = -imag(Z);
    
    ri = (Zreal_50mA_Temp_30 - Zreal).^2+(Zimag_50mA_Temp_30 - Zimag).^2;
    Ji = jacobian(ri,Beta); %15s
    
    Js = subs(Ji,Beta,Betak); %30s
    rs = subs(ri,Beta,Betak);
    
    J = eval(Js);   %50s
    r = eval(rs);   
    % inv(A)*B = A\B; B*inv(A) = B/A. This transform will increase the accuracy of calaculation
    % Beta = Beta - inv(J'*J)*J'*r
    Betak = Betak - (J'*J\J'*r)'
    
end

RMSE = sqrt(sum(r))

Zdl = 1./((w*1i).^Betak(6)*Betak(4)); % Zdl = 1./((w*1i).^alphadl*Qdl);
Zl = 1./((w*1i).^Betak(7)*Betak(5));  % Zl = 1./((w*1i).^alphal*Ql);
Zwb = coth(sqrt(1i.*w)*Betak(8))./(sqrt(1i.*w)*Betak(8)); % Zwb = coth(sqrt(1i.*w)*Wb)./(sqrt(1i.*w)*Wb);

Z0 = Betak(2)+((Betak(3)+Zwb).*Zdl)./(Betak(3)+Zwb+Zdl); % Z0 = Rl+((Rt+Zwb).*Zdl)./(Rt+Zwb+Zdl);
Z = Betak(1)+(Zl.*Z0)./(Zl+Z0); % Z = Re+(Zl.*Z0)./(Zl+Z0);

Zreal = real(Z);
Zimag = -imag(Z);

figure
plot(Zreal,Zimag,'*-',Zreal_50mA_Temp_30,Zimag_50mA_Temp_30,'o')
xlabel('Zreal(ohms)')
ylabel('Zimag(ohms)')
legend('Equivalent Circuit Model','Experimental Data');
title('Modeling of A123 Cell A002 at SOC=100%, -30°C with Galvanostatic EIS')
box on
grid on
