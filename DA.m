%% basic housekeeping
close all; 
clear all; 
clc;
%% Loading EEG data
fs = 256                    % fs — Sampling frequency, 
load('EEG.mat')       % loading .mat file containing EEGsig 
N =length(EEGsig)           % returns the length of EEGsig
t =[0:length(EEGsig)-1]/fs; %sampling rate; time= N/frequency

%% plotting original signal in time domain and frequency domain
figure(1);
subplot(2,1,1) 
plot(t,EEGsig); %time domain  
title('EEG Signal'); xlabel('Time in sec'); ylabel('Voltage (uV)'); 
z= fftshift(EEGsig);
nfft=length(z);
freq = (0:1/nfft:1-1/nfft)*fs;  %frequency-domain
subplot(2,1,2)
plot(freq,abs(fftshift(z))) %done to observe the range of frequencies available in the signal
title('Frequency domain'); xlabel('Frequency'); ylabel('Magnitude');

%% filtering insignificant low frequency and high frequency components to reduce noise
EEG = highpass(EEGsig,1,fs); %cutoff freq is 1Hz
EEG= lowpass(EEG,50,fs) %cut off freq is 50 Hz 

figure(2);
plot(t,EEG);
xticks(0:0.2:10)
xlabel('Time in sec'); ylabel('Voltage (uV)');  title('EEG signal after filtering'); %EEG waveform

%% Sampling and Quantisation 
max=max(EEG)
min=min(EEG)
n1=16;%NO OF BITS PER SAMPLE
L=2^n1;
del=(max-min)/L %Value of each level. Since we have
partition=min:del:max % definition of decision lines/partition
codebook=min-(del/2):del:max+(del/2); 
%codebook tells the quantizer which common value to assign to inputs that fall into each range of the partition
%he lower range/upper range of the codebook must be equally distributed in both sides
[index,quants] = quantiz(EEGsig,partition,codebook);  
%The function returns quants, which contains the scalar quantization of EEGsig and depends on
%the quantization levels and prescribed values in the codebook
%and returns index i.e the quantization levels of EEGsig by using the scalar quantization partition specified in input partition i.e does not depend upon codebook value.

l1=length(index); % 1 to n index to 0 to n-1 index since in matlab arrays start from index 1. 
for i=1:l1
if (index(i)~=0)
index(i)=index(i)-1;
end
end
l2=length(quants);
for i=1:l2 % mapping the extremas to the required range (done based on codebook range)
if(quants(i)==min-(del/2))
quants(i)=min+(del/2);
end
if(quants(i)==max+(del/2))
quants(i)=max-(del/2);

end
end
figure(3)
plot(t,EEG);
hold on
scatter(t,quants)
xticks(0:0.2:10)
xlabel('Time in sec'); ylabel('Voltage (uV)');  title('Sampled EEG signal'); 
figure(4)
stairs(t,quants)
xticks(0:0.2:10)
xlabel('Time in sec'); ylabel('Voltage (uV)'); title('Quantized EEG Signal');
%% Quantisation error
figure(5)
plot(t,quants-EEG) %difference in amplitude between filtered signal and quantised signal
xticks(0:0.2:10)
xlabel('Time in sec'); ylabel('Voltage (uV)'); title('Quantization error');

%% Reconstruction
% Decimal to binary conversion
code=de2bi(index,'left-msb'); %decimal to binary conversion with left most bit being the msb
k=1;
for i=1:l1 % to convert column vector to row vector
for j=1:n1
coded(k)=code(i,j); % to convert column vector to row vector
j=j+1;
k=k+1;
end
i=i+1;
end
code1=reshape(coded,n1,(length(coded)/n1)); %converting the array to a single row array.
index1=bi2de(code1,'left-msb');  %converting from binary to decimal with left most bit being the msb
resignal=del*index+min+(del/2); %to recover the modulated signal
figure(6);
hold on;
plot(t,resignal);
xticks(0:0.2:10)
title('Reconstructed Signal');
xlabel('Time in sec');
ylabel('Voltage (uV)');
figure(7);
hold on;
plot(t,quants-resignal);
xticks(0:0.2:10)
title('Reconstruction Error'); xlabel('Time in sec'); ylabel('Voltage (uV)');
%% Alpha wave 
%alpha band: 8-13 Hz
alpha = highpass(EEG,8,fs); %cutoff freq is 8Hz
alpha= lowpass(alpha,13,fs) %cut off freq is 13 Hz  
figure(8);
hold on
subplot(2,1,1)
plot(t,alpha);
xticks(0:0.2:10)
xlabel('Time in sec'); ylabel('Voltage (uV)');; 
title('Alpha wave from original signal'); 
alpha1 = highpass(resignal,8,fs); %cutoff freq is 8Hz
alpha= lowpass(alpha1,13,fs) %cut off freq is 13 Hz  
subplot(2,1,2)
plot(t,alpha1);
xticks(0:0.2:10)
xlabel('Time in sec'); ylabel('Voltage (uV)');
title('Alpha wave from reconstructed signal'); 
hold off;
%% Delta wave 
%delta band: 0-3Hz
delta= lowpass(EEG,3,fs) ;
delta1= lowpass(resignal,3,fs) %cut off freq is 3 Hz 
figure(9);
hold on
subplot(2,1,1)
plot(t,delta);
xticks(0:0.2:10)
xlabel('Time in sec'); ylabel('Voltage (uV)');; 
title('Delta wave from original filtered signal'); 
 
subplot(2,1,2)
plot(t,delta1);
xticks(0:0.2:10)
xlabel('Time in sec'); ylabel('Voltage (uV)'); 
title('Delta wave from reconstructed signal'); 
hold off

%% Theta wave 
%theta band- 3.5-7.5Hz
Theta = highpass(EEG,3.5,fs); %cutoff freq is 1Hz
Theta= lowpass(Theta,7.5,fs) %cut off freq is 50 Hz
theta1 = highpass(resignal,8,fs); %cutoff freq is 8Hz
theta1= lowpass(theta1,13,fs) %cut off freq is 13 Hz  
figure(10);
hold on
subplot(2,1,1)
plot(t,Theta);
xticks(0:0.2:10)
xlabel('Time in sec'); ylabel('Voltage (uV)');; 
title('Theta wave from original filtered signal'); 
 
subplot(2,1,2)
plot(t,theta1);
xticks(0:0.2:10)
xlabel('Time in sec'); ylabel('Voltage (uV)'); 
title('Theta wave from reconstructed signal'); 
hold off

%% Beta wave 
%beta band > 13Hz
beta= highpass(EEG,13,fs)
beta1= highpass(resignal,13,fs) %cut off freq is 13 Hz 
figure(11);
hold on
subplot(2,1,1)
plot(t,beta);
xticks(0:0.2:10)
xlabel('Time in sec'); ylabel('Voltage (uV)');; 
title('Beta wave from original filtered signal'); 
 
subplot(2,1,2)
plot(t,beta1);
xticks(0:0.2:10)
xlabel('Time in sec'); ylabel('Voltage (uV)'); 
title('Beta wave from reconstructed signal'); 
hold off