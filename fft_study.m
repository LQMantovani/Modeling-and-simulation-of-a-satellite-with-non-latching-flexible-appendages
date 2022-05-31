%Code to evaluate the frequency of vibration trough Fast Fourrier Tranform

% The FFT performed here is based on the example 
% provided in: https://www.mathworks.com/help/matlab/ref/fft.html

FR = Bodies.System.FR; %Sampling frequency
clear fft_matrice f_fft P1 pxx F 

for k = 1:3
    
    sample = 1/FR; %Sample period
    i0 = ceil(Bodies.System.LFFT*FR);
    t_ord = (1+i0):sample:(max(t)-sample); %Defines the time range

    Vec_ORD = sortrows([t, X(:,Bodies.B(1).pos_omega(k))]);     %Sorts the simulation output
    x_ord = interp1(Vec_ORD(1),Vec_ORD(2), t_ord, 'spline');    %Creates a vector with equally spaced time samples
    t_ord(end) = []; x_ord(end)=[];
    
    %FFT
    l_vec = length(t_ord);
    ft = fft(x_ord-mean(x_ord));        %Performs the FFT
    P2 = abs(ft/l_vec); 
    P1(:,k) = P2(1:floor(l_vec/2)+1);
    P1(2:end-1,k) = 2*P1(2:end-1,k);    %Creates a one-side FFT
    f_fft(:,k) = FR*(0:l_vec/2)/l_vec;  %Defines the frequency vector
    
    %Power Spectral Density
    [pxx_t, F_t] = pwelch(x_ord-mean(x_ord),[],[],[],FR);
    pxx(:,k) = pxx_t; F(:,k) = F_t;
end

%Generate plot
figure(nf);
subplot(231); stem(f_fft(:,1),P1(:,1)); grid minor; xlabel('Hz (p^1)'); ylabel('P1(Hz)'); xlim([0 50]);
subplot(234); stem(F(:,1),pxx(:,1)); grid minor; xlabel('Hz (p^1)'); ylabel('PSD'); xlim([0 50]);
subplot(232); stem(f_fft(:,2),P1(:,2)); grid minor; xlabel('Hz (q^1)'); ylabel('P1(Hz)'); xlim([0 50]);
subplot(235); stem(F(:,2),pxx(:,2)); grid minor; xlabel('Hz (q^1)'); ylabel('PSD'); xlim([0 50]);
subplot(233); stem(f_fft(:,3),P1(:,3)); grid minor; xlabel('Hz (r^1)'); ylabel('P1(Hz)'); xlim([0 50]);
subplot(236); stem(F(:,3),pxx(:,3)); grid minor; xlabel('Hz (r^1)'); ylabel('PSD'); xlim([0 50]);
sgtitle([Bodies.B(1).Name, ' vibration information']);
nf = nf + 1;

clear P1 f_fft pxx F


%Perform the FFT for the booms relative angular position and relative
%angular velocity
for j = 2:length(nflex)
    
    sample = 1/FR;                      %Sample period
    i0 = ceil(Bodies.System.LFFT*FR);
    vec_samples = 1+i0:FR*max(t);       %Defines the time interval
    
    t_ord = sample*vec_samples;
    %Creates an equally time spaced vector for the relative angular velocity
    v_ord = interp1(t,domega_v(:,Bodies.B(j).Main_rotation,j),sample*vec_samples);
    %Creates an equally time spaced vector for the relative angular position
    x_ord = interp1(t,d_Ang(:,Bodies.B(j).Main_rotation,j),sample*vec_samples);
    
    t_ord(end) = []; x_ord(end) = []; v_ord(end) = [];
    l_vec = length(t_ord);
    
    %FFT velocity
    ft = fft(v_ord-mean(v_ord));
    P2 = abs(ft/l_vec); P1 = P2(1:floor(l_vec/2)+1); P1_v=P1;
    P1_v(2:end-1) = 2*P1(2:end-1);
    
    %FFFT position
    ft = fft(x_ord-mean(x_ord));
    P2 = abs(ft/l_vec); P1 = P2(1:floor(l_vec/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    %Frequency for both FFTs are the same
    f_fft = FR*(0:l_vec/2)/l_vec;
    
    %Power Spectral Density
    [pxx_V, F_V] = pwelch(v_ord-mean(v_ord),[],[],[],FR);
    [pxx, F] = pwelch(x_ord-mean(x_ord),[],[],[],FR);
    
    figure(nf); 
    subplot(221); stem(f_fft,P1); grid minor; xlabel('Hz'); ylabel('P1(Hz) - Position'); xlim([0 50]);
    subplot(222); stem(f_fft,P1_v); grid minor; xlabel('Hz'); ylabel('P1(Hz) - Velocity'); xlim([0 50]);
    subplot(223); stem(F,pxx); grid minor; xlabel('Hz'); ylabel('PSD - Position'); xlim([0 50]);
    subplot(224); stem(F_V,pxx_V); grid minor; xlabel('Hz'); ylabel('PSD - Velocity'); xlim([0 50]);
    sgtitle([Bodies.B(j).Name, ' FFT - (Position, Velocity)']);
    nf = nf + 1;
end