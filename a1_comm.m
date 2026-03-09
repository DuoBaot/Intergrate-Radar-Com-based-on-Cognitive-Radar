%% 功能：通信端产生接收信号，降采样并搞到发送的信息
% 时间 2025 10 06 

clc;clear;close all

%% 参数设置
c=3e8;
para=struct();
para.B=1.2e7; % 线性调频信号的带宽（用复指数信号，因此也是调制后带宽）
para.fc=2.4e9; % 载频
para.Tc=40e-6; %单个脉冲持续时间，占空比1
para.fs=2e7; % 采样率
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
para.distance=0; % 通信距离

% 创建实例
qradcom=Qradcom(para);
%%混频信号生成,k斜率正负
[~, k_sequence, f_sequence] =qradcom.create_IF_signal();

%% 构造接收信号
signal_reci=qradcom.create_recieve_signal_downsample(k_sequence,f_sequence);
%% 加入高斯白噪声
SNR=-5;
[signal_reci_noise,noise] =qradcom.add_noise(signal_reci,SNR);

%% 解调（STFT）
num_parts=8;
num_window=para.num_samp_com/num_parts;
num_fft=para.num_samp_com;
num_overlap = 0;        % 重叠点数
binary_decode=qradcom.b1_fx_decode_using_STFT(signal_reci_noise,num_parts,num_window,num_fft,num_overlap);
qradcom.draw_stft_decode(signal_reci_noise,num_overlap,num_fft,num_window);

%% 显示结果



