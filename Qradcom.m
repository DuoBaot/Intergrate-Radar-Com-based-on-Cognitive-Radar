%-----------------速度采用fc---------

classdef Qradcom <handle
    properties
        B % 线性调频信号的带宽（用复指数信号，因此也是调制后带宽）
        fc % 载频
        Tc %单个脉冲持续时间，占空比1
        num_samp %采样点
        num_chirp %chirp数
        real_r % 目标距离，列向量
        real_v %目标速度，列向量
        fs % 采样率
        Ts % 采样周期
        fs_com % 通信端采样率
        num_samp_com %通信段一个脉冲的采样点数
        text
        emit_codes
        k_group
        f_group
        t_vector
        bit_per_symbol
        distance % 通信距离
    end

    properties (Constant)
    c = 3e8;      % 光速 [m/s]
    end

    methods
        function obj=Qradcom(para)
            obj.B=para.B; % 线性调频信号的带宽（用复指数信号，因此也是调制后带宽）
            obj.fc=para.fc; % 载频
            obj.Tc=para.Tc; %单个脉冲持续时间，占空比1
            obj.num_samp=round(para.Tc*para.fs); %采样点
            obj.fs=para.fs; % 采样率
            obj.Ts=1/para.fs;
            obj.fs_com=para.fs_com;
            obj.num_samp_com=para.num_samp_com;
            obj.t_vector=(0:1/obj.fs:obj.Tc-1/obj.fs)';            
            obj.text=para.text;
            obj.emit_codes=para.emit_codes;
            obj.num_chirp=para.num_chirp;   %chirp数
            obj.bit_per_symbol=para.bit_per_symbol;
            obj.real_r=para.real_r; % 目标距离，列向量
            obj.real_v=para.real_v; %目标速度，列向量          
            obj.k_group=para.k_group;
            obj.f_group=para.f_group;
            obj.distance=para.distance;
        end

        function [k_sequence,f_sequence] = k_f_sequence(obj)
            %功能：根据编码，获得有序的调制序列
            k_sequence=zeros(obj.num_chirp,1);
            f_sequence=zeros(obj.num_chirp,1);            
            for ii=1:obj.num_chirp
                %斜率
                switch obj.emit_codes(ii)
                    case '1', k_sequence(ii) = obj.k_group(1); f_sequence(ii) = obj.f_group(1);
                    case '0', k_sequence(ii) = obj.k_group(2); f_sequence(ii) = obj.f_group(2);
                end
            end
        end
        %% ---------------------------雷达部分------------------------

        function [IF_signal, k_sequence, f_sequence] = create_IF_signal(obj)
            % 功能：根据2进制编码生成混频信号        
            % 输出：IF_signal（混频信号）,斜率序列，起始频率序列
            
            tau = 2*obj.real_r/obj.c;

            % 相位调制（处理第3位编码）% 幅度调制（处理第4位编码）
            [k_sequence, f_sequence] = obj.k_f_sequence();
            % 频率选择与信号生成
            IF_signal=zeros(obj.num_chirp*obj.num_samp,1);
            for jj = 1:obj.num_chirp
                k_now=k_sequence(jj);
                f_now = f_sequence(jj);
                for h=1:length(obj.real_r)
                    IF_static =  exp(1i*(2*pi*f_now*tau(h) + 2*pi*k_now*tau*obj.t_vector)) ;             
                    tar_h= IF_static *exp(1i*(2*pi*obj.fc*2*obj.real_v(h)*obj.Tc*(jj-1)/obj.c));
                    IF_signal((jj-1)*obj.num_samp+1:jj*obj.num_samp) =...
                        IF_signal((jj-1)*obj.num_samp+1:jj*obj.num_samp)+tar_h;                  
                end
            end      
           

        end

        function [r_measured,fft1D,IF_matrix]=measure_R(obj,signal_IF,k_sequence,f_sequence)
            % 功能：用另类FFT对混频信号测距
            % 输入：
            %      signal_IF 列混频信号,chirp*samp,1
            %      k_sequence 列斜率序列 chirp,1
            % 输出：
            %       r_measured,测得的距离序列 chirp,1
            %       fft1D, 1dfft矩阵 chirp,samp
            %       IF_matrix 混频信号矩阵 chirp,samp

            % 先得到信号矩阵
            IF_matrix=reshape(signal_IF,obj.num_samp,obj.num_chirp).'; % chirp*samp

            % 对矩阵的每一列作另类FFT
            fft1D = obj.special_fft_for_range(IF_matrix.',k_sequence,f_sequence); % (NsxM)

            % 对每一列fftshift再转置
            fft1D= fftshift(fft1D,1);% 每一列的上下交换
            fft1D=fft1D.';  %(MxNs)
            % 找出每个脉冲的距离
            if  mod(obj.num_samp,2)==0
                mid=obj.num_samp/2;
            else
                mid=(obj.num_samp-1)/2;
            end 
            % 保证索引是1~n-1   --shift-->  1~mid,mid+1~n-1
            % 对于距离 0~d_max-delta_d,-d_max~-delta_d -->-d_max~-delta_d,0~d_max-delta_d
            % 0 频分量在mid+1的索引处
            [~,index_1D]=max(abs(fft1D(:,:)),[],2);
            f_1D_max=(index_1D-mid-1)/(obj.num_samp/2)*(obj.fs/2);
            r_measured = obj.c * f_1D_max ./ (2 *abs(k_sequence));

        end

        function fft1D = special_fft_for_range(obj, IF_matrix, k_sequence,f_sequence)
            % 功能：对IF_matrix（矩阵）的每一列进行另类FFT，为测得距离
            % 输入： 
            %   IF_matrix  : 矩阵 (采样点×脉冲)，每列为一个信号
            %   k_sequence : 频率序列 (M×1)，与信号长度匹配
            % 输出：
            %   fft1D     : 未经shitf的FFT结果 (MxNs)，

            %% 输入维度检查
            [num_rows, num_cols] = size(IF_matrix);

            % 检查IF_matrix维度是否匹配预期
            if num_rows ~= obj.num_samp || num_cols ~= obj.num_chirp
                error('IF_matrix维度错误: 期望 [%d, %d], 实际 [%d, %d]', ...
                      obj.num_samp, obj.num_chirp, num_rows, num_cols);
            end

            % 检查k_sequence维度
            if length(k_sequence) ~= obj.num_chirp
                error('k_sequence长度错误: 期望 %d, 实际 %d', ...
                      obj.num_chirp, length(k_sequence));
            end

            %% 预计算参数
            Ns = obj.num_samp;
            k_vector = 0:(Ns-1);
            k_seq = k_sequence / abs(k_sequence(1));  % +-1序列归一化
            scale_factor = 1 / (obj.Ts * Ns * abs(k_sequence(1)));

            % 预分配结果矩阵
            fft1D = zeros(Ns, obj.num_chirp);

            %% 向量化计算 - 减少循环
            % 构造时间索引矩阵 (N×1)
            m_vec = (0:Ns-1)'; % 离散的时间向量

            for chirp = 1:obj.num_chirp
                k_now = k_seq(chirp);  % 当前chirp的k值,归一化后的
                f_now =f_sequence(chirp);
                % 向量化计算FFT矩阵 (一次计算整个矩阵)
                % 第一个相位项: -1j*2*pi*(m-1)*k_now*k/N 频率抵消
                phase1 = -1j * 2 * pi * (m_vec * k_now) * (k_vector / Ns);

                % 第二个相位项: -1j*m*k/(Ts*N*|k_sequence(1)|) 相位抵消
                phase2 = -1j * 2* pi *f_now* k_vector * scale_factor;      
                % 合并相位项
                total_phase = phase1 + phase2;

                % 计算FFT向量
                vector_fft = exp(total_phase); %num_sam*num_sam

                % 矩阵乘法计算FFT (一次计算所有频率点)
                fft1D(:, chirp) = vector_fft.' * IF_matrix(:, chirp); 
            end

        end
        
        function [v_measured,R_measured,fft2D]= measure_V(obj,fft1D)
            %功能：对某一帧的fft1d矩阵进行另类2dfft测速
            %输入：
            %   fft1D:已做FFTshift的FFT1D结果 (MxNs)，
            %返回：
            %   v_measured:峰值对应的速度,
            %   R_measured 峰值对应的速距离
            %   2dfft，chirp数*样本点数
                  
            %% 另类fft得到fft2D,
            %找中点索引，
            if mod(obj.num_chirp,2)==0 %脉冲个数偶数 k:-N/2+1~N/2
                mid_v=obj.num_chirp/2;
            else % 脉冲个数奇数k:-（N-1）/2~（N-1）/2
                mid_v=(obj.num_chirp+1)/2;
            end
            if  mod(obj.num_samp,2)==0
                mid_r=obj.num_samp/2;
            else
                mid_r=(obj.num_samp-1)/2;
            end 
            %% 对fft1D每一列进行一般的fft            
            fft2D=fft(fft1D); %每一行是一个脉冲，按列FFT
            fft2D= fftshift(fft2D,1);% 每一列的上下交换

            %% 测速
            % 直接找原矩阵 fft2D 中最大值的位置
            [~, linear_idx] = max(abs(fft2D), [], "all");  % 线性索引
            [index_v, index_r] = ind2sub(size(fft2D), linear_idx);  % 转换为行、列索引
            f_1D_max=(index_r-mid_r-1)/(obj.num_samp/2)*(obj.fs/2);
            R_measured = obj.c * f_1D_max ./ (2 *abs(obj.k_group(1)));
            v_measured=(index_v-mid_v-1)*obj.c/(obj.num_chirp*obj.fc*2*obj.Tc);

            % [~,index_v]=max(all_col);% 最大的行索引
            % v_measured=(index_v-mid_index-1)*obj.c/(obj.num_chirp*obj.fc*2*obj.Tc);
            %这里之所以-mid_index-1是因为，找到的最大值索引的范围永远是1~N，
            % 而事实上速度可以是负的 因此向左平移一半
        end

        function fft2D=special_fft_for_velocity(obj,fft1D,f_sequence)
            % 功能：对fft1Dn这样一个矩阵的每一列，进行另类fft
            % 输入：
            %   fft1D ：脉冲，采样点
            %   f_sequence 脉冲，1
            % 输出：
            %   fft2D：未做过FFTshift的fft2D: chirp，sampl

            %% 检验维度,fft1D,f_sequence

            M=obj.num_chirp;
            %构造fft向量，第一个维度是去相位抵消的匹配信号，第二个维度是k的取值数量          
            vector_fft=zeros(M,M);   

            k_vector=0:M-1;
            f_sequ=f_sequence/abs(f_sequence);
            %% 构造fft向量
            for k_index=1:M %一共M个fft向量
                for n=1:M %每个向量构造用到n循环
                    f_now=f_sequ(n);
                    vector_fft(n,k_index)=exp(-1j*2*pi*(n-1)*f_now*k_vector(k_index)/M)  ;
                end
            end
            %% 另类fft
            fft2D=NaN(M,obj.num_samp);

            for ns=1:obj.num_samp   %先对每个采样点对应的列循环
                for k_index=1:obj.num_chirp %再对该采样点列的所有脉冲循环
                    fft2D(k_index,ns)=fft1D(:,ns).'*vector_fft(:,k_index);
                end
            end
        end

        % 补零FFT目前存在问题
        function fft1D_zeros=special_fft_for_range_zeros(obj,IF_matrix, k_sequence,f_sequence,num_fft)
            % 功能：对IF_matrix（矩阵）的每一列进行另类FFT，补零搞多点数
            % 输入： 
            %   IF_matrix  : 矩阵 (采样点×脉冲)，每列为一个信号
            %   k_sequence : 斜率序列 (M×1)，与信号长度匹配
            %   f_sequence : 频率序列 (M×1)，与信号长度匹配
            %   num_fft    : 目标 FFT 点数（≥ 原始行数）
            % 输出：
            %   fft1D_zeros     : 未经shitf的FFT结果 (MxNs)，

            %% 输入维度检查
            [num_rows, num_cols] = size(IF_matrix);

            % 检查IF_matrix维度是否匹配预期
            if num_rows ~= obj.num_samp || num_cols ~= obj.num_chirp
                error('IF_matrix维度错误: 期望 [%d, %d], 实际 [%d, %d]', ...
                      obj.num_samp, obj.num_chirp, num_rows, num_cols);
            end

            % 检查k_sequence维度
            if length(k_sequence) ~= obj.num_chirp
                error('k_sequence长度错误: 期望 %d, 实际 %d', ...
                      obj.num_chirp, length(k_sequence));
            end
            %% 对混频信号的每一列补零
            %补零参数
            if num_fft < obj.num_samp
                warning("fft点数小于原点数，按原来点数")
                num_fft =obj.num_samp;   % 不允许缩减
            end
            % 仅补零
            IF_matrix_zeros = [IF_matrix; zeros(num_fft - size(IF_matrix,1), size(IF_matrix,2))];

            %% 预计算参数
            Ns = obj.num_samp;
            k_vector = 0:(num_fft-1);
            k_seq = k_sequence / abs(k_sequence(1));  % +-1序列归一化
            scale_factor = 1 / (obj.Ts * Ns * abs(k_sequence(1)));

            % 预分配结果矩阵
            fft1D_zeros = zeros(num_fft, obj.num_chirp);

            %% 向量化计算 - 减少循环
            % 构造时间索引矩阵 (N×1)
            m_vec = (0:num_fft-1)'; % 离散的时间向量

            for chirp = 1:obj.num_chirp
                k_now = k_seq(chirp);  % 当前chirp的k值,归一化后的
                f_now =f_sequence(chirp);
                % 向量化计算FFT矩阵 (一次计算整个矩阵)
                % 第一个相位项: -1j*2*pi*(m-1)*k_now*k/N 频率抵消
                phase1 = -1j * 2 * pi * (m_vec * k_now) * (k_vector / Ns);

                % 第二个相位项: -1j*m*k/(Ts*N*|k_sequence(1)|) 相位抵消
                phase2 = -1j * 2* pi *f_now* k_vector * scale_factor;      
                % 合并相位项
                total_phase = phase1 + phase2;

                % 计算FFT向量
                vector_fft = exp(total_phase); %num_sam*num_sam

                % 矩阵乘法计算FFT (一次计算所有频率点)
                fft1D_zeros(:, chirp) = vector_fft.' * IF_matrix_zeros(:, chirp); 
            end
            fft1D_zeros=fft1D_zeros.';
        end
        
        %% ---------------------------通信部分------------------------
        function signal_reci=create_recieve_signal_downsample(obj, k_sequence,f_sequence)
            % 功能：构造接收端已经降采样的接收信号
            % 输入：
            % 输出：
            signal_reci=zeros(obj.num_chirp*obj.num_samp_com,1);
            dely=obj.distance/obj.c;
            t_com=(0:1/obj.fs_com:obj.Tc-1/obj.fs_com)';    
            for ii =1:obj.num_chirp
                f_now=f_sequence(ii);
                k_now=k_sequence(ii);
                signal_reci((ii-1)*obj.num_samp_com+1:ii*obj.num_samp_com)=...
                    exp(   1j*(2*pi*f_now*(t_com-dely)+pi*k_now*(t_com-dely).^2)  );
            end

        end

        function binary_decode=b1_fx_decode_using_STFT(obj,signal,num_parts,num_window,num_fft,num_overlap)
            % 功能：短时傅里叶变换解码
            
            window = hamming(num_window); % 窗函数
   
            stft_matrix=stft(signal,obj.fs_com,'Window',window,'OverlapLength',num_overlap,'FFTLength',num_fft);
            
            binary_decode=blanks(length(signal)/obj.num_samp_com);%创建空白字符数组
            index=zeros(1,num_parts);
            k_door=0.2;%阈值
            for i=1:num_parts:size(stft_matrix,2)
                for jj=1:num_parts
                [~,index(jj)]= max(abs(stft_matrix(:,(i-1)+jj) ));
                end
                this_Symbol_k=diff(index)*obj.fs_com^2/( size(stft_matrix,1) )/num_window;
                %diff(index)只是单纯的对点数微分，要将横纵轴变为真实时间、频率，即求出横坐标相邻2个点的时间？纵坐标相邻2个点的频率？
                %横坐标相邻2个点的时间：1/fs_eff*num_samples/（num_samples/num_window）,一个码元1/fs_eff*num_samples的时间，而一个码元（窗中）又分成10个点
                % num_s纵坐标相邻2个点的频率，(fs_eff)/size(stft_matrix,1)    
                this_Symbol_k=mean(this_Symbol_k(3:num_parts-2));%此码元的斜率
                l=round(i/num_parts)+1;
                if  abs((this_Symbol_k-obj.k_group(1))/obj.k_group(1))<=k_door
                    binary_decode(l)='1';
                elseif abs((this_Symbol_k-obj.k_group(2))/obj.k_group(2))<=k_door
                    binary_decode(l)='0';
                else
                    binary_decode(l)=  NaN;
                    warning("*****发现未知斜率！！*****");
                    disp("此时i="+i+"，对应码元个数："+(i-1)/num_parts+"此时的斜率："+this_Symbol_k);
                end
            end
        end

        %% ---------------------------画图------------------------
        function draw_xk(obj,x_k,delta_v)
            % 绘制三维网格图
            samp_vector = round(linspace(0, obj.num_samp, size(x_k, 2))); % 采样向量，
            max_idx = size(x_k, 1) / 2;
            v_vector = (-max_idx : max_idx-1) * delta_v; % 速度向量，           
            [samp_grid,v_grid]=meshgrid(samp_vector ,v_vector);            
            % 绘制三维网格图
            figure;
            mesh(samp_grid,v_grid,abs(x_k)); % 创建三维网格图
            xlabel('sample'); % 设置X轴标签
            ylabel('Velocity(m/s) '); % 设置Y轴标签
            zlabel('Amplitude'); % 设置Z轴标签
            title('3D Mesh Plot of X[k]'); % 设置图表标题
            
        end

        function draw_signal_matrix(obj, IF_matrix_prec,mytitle)
            % 功能：绘制信号矩阵（取模值）的三维图
            % 输入：
            %   IF_matrix_prec - 信号矩阵，行是脉冲数，列是采样点数，值为复数
            % 
            % 绘制说明：
            %   x轴：时间（对应采样点）
            %   y轴：脉冲序号（对应慢时间）
            %   z轴：信号幅度（取模值）
            %   颜色：表示相位信息（取角度值）
            if nargin<3
                mytitle=['编码：',obj.text, ' - 信号（实部）矩阵三维图'];
            end
        
        
            % 计算实部和相位
            real_matrix = real(IF_matrix_prec);    % 幅度矩阵
            phase_matrix = angle(IF_matrix_prec); % 相位矩阵（弧度）
            
            % 创建坐标轴
            [num_pulses, num_samples] = size(IF_matrix_prec);
            time_axis = (0:num_samples-1)/obj.fs;      % 时间轴（秒）
            pulse_axis = 1:num_pulses;                 % 脉冲序号
            
            % 创建网格
            [T, P] = meshgrid(time_axis, pulse_axis);
            
            % 创建图形
            figure;
            set(gcf, 'Position', [100, 100, 800, 600]); % 设置图形大小
            
            % 绘制三维曲面
            surf(T, P, real_matrix, phase_matrix, ...
                 'EdgeColor', 'none', 'FaceAlpha', 0.8);
            
            % 设置视角
            view(45, 30);
            
            % 设置标签和标题
            xlabel('时间 (s)');
            ylabel('脉冲序号');
            zlabel('信号实部幅度');

            % title(['编码：',obj.text, ' - 信号（实部）矩阵三维图']);
            title(mytitle)
            
            % 添加颜色条
            cc = colorbar;
            cc.Label.String = '相位 (rad)';
            
            % 设置坐标轴范围
            xlim([0, obj.Tc]);
            ylim([1, num_pulses]);
            
            % 美化图形
            grid on;
            colormap(jet);
            shading interp;
            
            % 添加额外信息
            % annotation('textbox', [0.15, 0.85, 0.2, 0.1], ...
            %            'String', sprintf('带宽: %.1f MHz\n载频: %.1f MHz', ...
            %                             obj.B/1e6, obj.fc/1e6), ...
            %            'FitBoxToText', 'on', ...
            %            'BackgroundColor', 'white');
        end

        function draw_fft1D(obj,fft1D)
            %功能：绘制1DFFT矩阵三维图，同时再绘制频率轴变为距离的图（脉压处理）
            % 绘制三维网格图
            % 但是由于斜率k不一样，不能直接用距离表示，用频率
            y_vector=linspace(0, obj.num_chirp, size(fft1D, 1)); %y轴，脉冲数  
 
            %% 脉压处理 

            figure;
            % R=f*c/(2k)
            %本应-d_max~+d_max但是由于共2n+1个点，实际应该-d_max~+d_max-1个点
            %故先扩充一个点-->2n+1，再舍弃最后一个点
            x_vector_pul = linspace(-obj.fs*obj.c/(4*abs(obj.k_group(1))), obj.fs*obj.c/(4*abs(obj.k_group(1))), size(fft1D, 2)+1);  %x轴，距离唯，
            x_vector_pul=x_vector_pul(1:end-1);
            [x_grid_pul,y_grid_pul]=meshgrid(x_vector_pul ,y_vector);
            mesh(x_grid_pul,y_grid_pul,20*log10(abs(fft1D+eps))); % 创建三维网格图
            xlabel('距离(m)'); % 设置X轴标签
            ylabel('脉冲数 '); % 设置Y轴标签
            zlabel('幅度(dB)'); % 振幅20倍
            title('脉冲压缩处理'); % 设置图表标题

        end

        function draw_fft2D(obj,fft2D,delta_v)
            %功能：绘制2DFFT矩阵三维图，同时再绘制频率轴变为距离的图（脉压处理）
            % 绘制三维网格图

            figure;
            %本应-d_max~+d_max但是由于共2n+1个点，实际应该-d_max~+d_max-1个点
            %故先扩充一个点-->2n+1，再舍弃最后一个点
            x_vector_pul = linspace(-obj.fs*obj.c/(4*abs(obj.k_group(1))), obj.fs*obj.c/(4*abs(obj.k_group(1))), size(fft2D, 2)+1);  %x轴，距离唯，
            if mod(obj.num_chirp,2)==0 %脉冲个数偶数 k:-N/2+1~N/2
                mid_index=obj.num_chirp/2;
            else % 脉冲个数奇数k:-（N-1）/2~（N-1）/2
                mid_index=(obj.num_chirp+1)/2;
            end
            v_vector = (-mid_index : mid_index-1) * delta_v; % 速度向量， 
            x_vector_pul=x_vector_pul(1:end-1);
            [x_grid_pul,y_grid_pul]=meshgrid(x_vector_pul, v_vector);
            mesh(x_grid_pul,y_grid_pul,20*log10(abs(fft2D+eps))); % 创建三维网格图
            xlabel('距离(m)'); % 设置X轴标签
            ylabel('速度(m/s) '); % 设置Y轴标签
            zlabel('幅度(dB)'); % 设置Z轴标签
            title('2D-FFT矩阵'); % 设置图表标题

            end

        function draw_stft(obj,signal)
            %功能：绘制信号的STFT
            num_samples=length(signal)/obj.num_chirp;
            num_parts=16;%将一个码元分成几份，来算斜率，此数应可被num_samples整除!，影响横轴点数
            num_window=obj.num_samp/num_parts; %窗函数中的点数
            num_fft = 8*num_samples; % FFT点数
            figure
            stft(signal,obj.fs,'Window',hamming(num_window),'OverlapLength',0,'FFTLength',num_fft);
        end
        
        function draw_oneChirpCompr(obj,xk)
            % 功能：绘制某一个脉冲内的脉冲压缩结果，频域图
            x_vector_pul = linspace(-obj.fs/2, obj.fs/2, size(xk,2)+1);  %x轴，距离唯，
            x_vector_pul=x_vector_pul(1:end-1);
            % plot(x_vector_pul,20*log10(abs(xk)));
            xk=xk/max(xk);
            x_dB=20*log10(abs(xk));
            stem(x_vector_pul, x_dB, 'filled'); grid on;
            xlabel("频率（Hz）")
            ylabel("幅度（dB）")
            title("归一化后某脉冲的脉冲压缩结果")
        end
       
        function draw_stft_decode(obj,signal,num_overlap,num_fft,num_window)    
            % 功能：画对解调信号的短时傅里叶变换
            window = hamming(num_window); % 窗函数
            figure
            stft(signal,obj.fs_com,'Window',window,'OverlapLength',num_overlap,'FFTLength',num_fft);
        end
        %跳频,待修改！！！0925
        function signal_trans=d1_fx_create_trans_signal(obj)
            % 功能：构造接收信号，为了避免频率范围的差别，接收信号是原始的
            %调幅、调相项是没有包含在内的
            %时间： 2025 07 29
            
            signal_trans=zeros(1,obj.num_chirp*obj.num_samp);
            t=0:1/obj.fs:obj.Tc-1/obj.fs;

            tau = 2*target.real_R/obj.c;
               
            for ii=1:obj.num_chirp
                %以接收信号到达时间为参考，发射信号应该是t+tau
                phi=2*pi*f_sequence_hop(ii)*(t+tau)+pi*obj.k_sequence(ii)*(t+tau).^2;
                tmp=exp(1j*phi);
                signal_trans((ii-1)*num_sample+1:ii*num_sample)=tmp;
            
            end
   
        end

    end

    methods(Static)
        function emit_codes = transform_words_into_codes(obj)
            % 功能：将汉字变为二进制编码
            utf8Bytes = unicode2native(obj.text, 'UTF-8');
            %一个汉字3个3位数，一个数字占2个比特
            % 将每个字节转换为二进制字符串
            emit_codes = blanks(24*length(obj.text));
            for i = 1:length(utf8Bytes)
                binaryByte = dec2bin(utf8Bytes(i), 8); % 确保每个字节是 8 位
                emit_codes((i-1)*8+1:i*8) = binaryByte;
            end
            
        end
    
        function draw_t_and_f(signal, fs, mytitle)
            % 功能：在同一图形窗口中用上下子图显示信号的时域和频域图
            % 时间：2025-08-01
            % 输入：
            %   signal - 输入信号
            %   fs     - 采样频率
            %   mytitle - 图形标题（可选）
            
            if nargin < 3
                mytitle = "信号";
            end
            
            % 创建图形窗口（设置适当的大小）
            fig = figure('Position', [100, 100, 800, 500]); % 宽度800，高度500
            
            % ---------------- 时域图（上子图） ----------------
            subplot(2, 1, 1); % 2行1列，第1个子图
            t = (0:length(signal)-1)/fs;          % 时间轴
            plot(t*1e3, real(signal), 'LineWidth', 1);
            xlabel('时间 (ms)');
            ylabel('幅度');
            title(mytitle+' 实部时域波形');
            grid on;
            box on;
            
            % ---------------- 频域图（下子图） ----------------
            subplot(2, 1, 2); % 2行1列，第2个子图
            N = length(signal);
            f = (-N/2:N/2-1)*(fs/N);            % 双边频率轴      
            S = fftshift(fft(signal));
            plot(f/1e6, 20*log10(abs(S)+eps)); % MHz单位，dB
            xlabel('频率 (MHz)');
            ylabel('幅度 (dB)');
            title(mytitle+' 频谱');
            grid on;
            
            % 调整子图间距（可选）
            % h = gcf; % 获取当前图形句柄
            % set(h, 'Position', [100 100 800 500]); % 可以在这里调整图形大小
        end
     
        function [signal_IF_noise,noise]=add_noise(signal_IF,SNR)
            % 功能：给信号按照信噪比加入高斯白噪声
            % 输入：
            %   signal_IF,原始无噪的列信号
            %   SNR 信噪比dB SNR = 10*log10(信号功率 / 噪声功率)
            % 输出:
            %   signal_IF_noise, 加噪后的信号
            %   noise  纯噪声
            if isinf(SNR)
                signal_IF_noise=signal_IF;
                noise=zeros(size(signal_IF));
            else
                % 使用 awgn 函数添加噪声
                signal_IF_noise = awgn(signal_IF, SNR, 'measured');
                
                % 计算噪声
                noise = signal_IF_noise - signal_IF;
            end
        end


    end
end
