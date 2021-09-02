clear,clc 
m=8;                                      % sensors ��Ԫ
n=2;                                      % sources �ź�Դ
theta=[-20 60];                            % in angle �Ƕ�
d=1/2;                                    % 1/2 lambada ��Ԫ���
N=1000;                                    % samples  ����
L=100;                                    % resolution in [-90' 90']  ����
Meann=0;                                  % mean of noise 
varn=1;                                   % variance of noise
SNR=10;                                   % signal-to-noise ratio  �����
INR=10;                                   % interference-to-noise ratio �����
rvar1=sqrt(varn) * 10^(SNR/20);           % variance of signal �ź���ɢ��
rvar2=sqrt(varn) * 10^(INR/20);           % variance of interference

% generate the source signals
s=[rvar1*exp(j*2*pi*50*0.001*[0:N-1])               %ѵ������
   rvar2*exp(j*2*pi*(100*0.001*[0:N-1]+rand))];
% generate the A matrix
A=exp(-j*2*pi*d*[0:m-1].'*sin(theta*pi/180));            %������������
% generate the noise component
e=sqrt(varn/2)*(randn(m,N)+j*randn(m,N));
% generate the ULA data
Y=A*s+e;

% initialize weight matrix and associated parameters for LMS predictor
de =s(1, :);
mu=1e-3;
w = zeros(m, 1);

for k = 1:N
    % predict next sample and error
    y(k) = w'*Y(:, k);
    e(k)  = de(k) - y(k);                          %û�������������
    % adapt weight matrix and step size
    w = w + mu * Y(:,k)*conj(e(k));
end 

% beamforming using the LMS method
beam=zeros(1,L);
for i = 1 : L
   a=exp(-j*2*pi*d*[0:m-1].'*sin(-pi/2 + pi*(i-1)/L));
   beam(i)=20*log10(abs(w'*a));
end
 
% plotting command followed
figure
angle=-90:180/L:(90-180/L);
plot(angle,beam);
xlabel('�Ƕ�/��');
ylabel('������Ӧ/dB');


