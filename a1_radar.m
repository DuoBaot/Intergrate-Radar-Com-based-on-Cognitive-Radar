%% 功能：雷达端只发送2类正负的调频，再解出信息，单目标
% 时间 2025 10 06 
% 已经证明，这种方法抗干扰好，主副瓣比高

clc;clear;close all

%% 参数设置
c=3e8;
para=struct();
para.B=1.2e7; % 线性调频信号的带宽（用复指数信号，因此也是调制后带宽）
para.fc=2.4e9; % 载频
para.Tc=40e-6; %单个脉冲持续时间，占空比1
para.fs=2.4e7; % 采样率
para.k_group=[1;-1]*para.B/para.Tc;
para.f_group=[-para.B/2;para.B/2]+para.fc;
para.text="电子科技大学";        
para.emit_codes=Qradcom.transform_words_into_codes(para);
para.bit_per_symbol=1;
para.num_chirp=length(para.emit_codes)/para.bit_per_symbol;
delta_R=c/(2*para.B);%距离分辨率,按带宽小的来算，如果按带宽大的，设定，那矮脉冲就测不全
delta_v=c/(2*para.fc*para.Tc*para.num_chirp);%速度分辨率
delta_f=1/para.Tc;
R_max=min(para.fs*para.Tc*c/(4*para.B),para.Tc*c/2);
v_max=c/(4*para.Tc*para.fc);
para.real_r=7*delta_R; % 目标距离，列向量
para.real_v=3*delta_v;%目标速度，列向量
para.fs_com=2e8; % 通信端采样率
para.num_samp_com=round(para.fs_com*para.Tc);
para.distance=0000; % 通信距离设为0，实际上通信距离对信号相位没有影响

% 创建实例
qradcom=Qradcom(para);

%% 混频信号生成,k斜率正负
[signal_IF, k_sequence, f_sequence] =qradcom.create_IF_signal();

%% 加入高斯白噪声
% SNR=-35;
SNR=-1;
[signal_IF_noise,noise] =qradcom.add_noise(signal_IF,SNR);

%% 每个脉冲的测距  历时 1.447798 秒。
[R_measured_each_chirp,fft1D,IF_matrix]= qradcom.measure_R(signal_IF_noise, k_sequence,f_sequence);
% qradcom.draw_fft1D(fft1D) %144X800

% 绘制某一个脉冲内的脉冲压缩结果
% 需要先补零再脉冲压缩
% 1006有问题，点数800-->2048周期从10000m-->4000m
%是不是时间变↑，频谱被压缩，但是观察时间为什么会变长？
% num_fft=2048*1;
% fft1D_zeros=qradcom.special_fft_for_range_zeros(IF_matrix.', k_sequence,f_sequence,num_fft);
% xk=fft1D_zeros(5,:);
% qradcom.draw_oneChirpCompr(xk)
% qradcom.draw_fft1D(fft1D_zeros)

%% 总的测距测速 历时 0.061572 秒。

[v_measured,R_measured,fft2D]= qradcom.measure_V(fft1D );
% qradcom.draw_fft2D(fft2D,delta_v)

%% 结果显示
fprintf("设定的距离=%f m ，速度=%f m/s \n",para.real_r,para.real_v);
fprintf("信号持续总时间 %.2f ms\n",para.Tc*para.num_chirp*1e6);
fprintf("实际测得的距离=%f m ，速度=%f m/s \n",R_measured,median(v_measured));


% 混频信号，时域，频域
qradcom.draw_t_and_f(signal_IF,para.fs,"原始混频信号");
qradcom.draw_stft(signal_IF);
qradcom.draw_t_and_f(signal_IF_noise,para.fs,"加噪混频信号");
% 绘制混频信号矩阵三维图
% qradcom.draw_signal_matrix(IF_matrix,"原始混频信号");
qradcom.draw_signal_matrix(IF_matrix,"加噪混频信号");
% 绘制脉压后1DFFT矩阵
qradcom.draw_fft1D(fft1D)

% 2dfft
qradcom.draw_fft2D(fft2D,delta_v)


