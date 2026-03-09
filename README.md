# Integrate-Radar-com-based-on-Cognitive-Radar - 认知雷达通信一体化技术与试验验证

## 项目说明
- 针对通感一体波形雷达感知性能恶化的问题，改进动目标检测（MTD）算法，在-35 dB的极端信噪比下实现无漏检，将单次检测耗时从约1分钟优化至10秒以内。为提升抗窄带压制干扰能力，探索脉间跳频技术的集成。

## 改进MTD算法原理

### 基于映射调制的通感一体
#### 信号构建原理
雷达通信系统采用正负调频的方法，脉冲起始频率：
$$
f_{0,n}=\left\{
\begin{matrix}
-\frac{B}{2}&,\text{发送码元“0”}\\
 \frac{B}{2}&, \text{发送码元“1”}
\end{matrix}
\right.
$$

脉冲的调频斜率：
$$
 k_{n}=\left\{
\begin{matrix}
\frac{B}{T_\text{c} }&,\text{发送码元“0”}\\
- \frac{B}{T_\text{c}}&, \text{发送码元“1”}
\end{matrix}
\right.
$$

定义门信号：
$$
\Pi (t)=\left\{\begin{matrix}
 1 & , 0 <t\le T_\text{c}  \\
 0 & , \text{else}
\end{matrix}\right.
$$
  
定义第一个脉冲内的参考发射信号为：
$$
  s_{n}( t ;f_{0,n},k_{n}) =
 A e^{j\left(2\pi f_{0,n}  t + \pi k_{n} t^{2} \right)}\cdot \Pi (t)
$$
其中$A$是信号的幅度，$f_{0,n}$是该脉冲的起始频率，$k_{n}$是该脉冲的调频斜率。

雷达发射信号为：
$$
  \begin{split} 
s_\text{tra}( t ;f_{0,n},k_{n}) &=\sum_{n=1}^{N} s_{n}( t -(n-1)T_\text{c};f_{0,n},k_{n})\\
 &=\sum_{n=1}^{N} A e^{j\left\{2\pi f_{0,n}  [t-(n-1)T_\text{c}] + \pi k_{n}  [t-(n-1)T_\text{c}]^{2} \right\}}\cdot \Pi (t-(n-1)T_\text{c})
\end{split}
$$

假设一个目标与雷达的距离为$R$，假设物体静止。雷达的接收信号：
$$
  \begin{split} 
s_\text{rec}^\text{rad}( t ;f_{0,n},k_{n}) &=\sum_{n=1}^{N} s_{n}( t -(n-1)T_\text{c}-\tau;f_{0,n},k_{n})\\
 &=\sum_{n=1}^{N} A e^{j\left\{2\pi f_{0,n}  [t-(n-1)T_\text{c}-\tau] + \pi k_{n}  [t-(n-1)T_\text{c}-\tau]^{2} \right\}}\cdot \Pi (t-(n-1)T_\text{c}-\tau)
\end{split}
$$
其中，时延$\tau=\frac{2R}{c}$。

如果目标相对于雷达，以速度$v$靠近，雷达的接收信号会增加一项多普勒频率。第$n$个脉冲相较于第$n-1$个脉冲而言，因速度的微小位移引起的时延$\Delta \tau=\frac{2vT_\text{c}}{c}$，由此造成的相位偏差近似建模为$\Delta \phi =2\pi f_{c}\Delta \tau $，以第一个脉冲为基准，第$n$个脉冲与第一个脉冲的相位偏差为：
$$
\Delta \phi_n =\frac{4\pi vT_\text{c} f_{c}}{c} \cdot (n-1)
$$
因此回波信号应修正为：
$$
  \begin{split} 
s_\text{rec}^\text{rad}( t ;f_{0,n},k_{n}) =&\sum_{n=1}^{N} s_{n}( t -(n-1)T_\text{c}-\tau;f_{0,n},k_{n})e^{j\frac{4 \pi(n-1) v T_\text{c}f_{c}}{c}}\\
 =&\sum_{n=1}^{N} A e^{j\left\{2\pi f_{0,n}  [t-(n-1)T_\text{c}-\tau] + \pi k_{n}  [t-(n-1)T_\text{c}-\tau]^{2} \right\}}e^{j\frac{4 \pi(n-1) v T_\text{c}f_{c}}{c}}\cdot  \\ 
 & \Pi (t-(n-1)T_\text{c}-\tau)
\end{split}
$$

混频信号由发送与接收信号相乘再经过低通滤波器，可以表示为：
$$
\begin{aligned} 
s_{\text{IF}}\left ( t;f_{0,n} ,k_{n}\right ) &=\sum_{n=1}^{N} A\cdot e^{j\left \{2\pi f_{0,n}\tau +2\pi k_{n}\tau [t-(n-1)T_\text{c}]-\pi k_{n} \tau^{2}\right \}}
e^{j\frac{4 \pi(n-1) v T_\text{c}f_{c}}{c}}  
\\ &\cdot \Pi (t-(n-1)T_\text{c})
\end{aligned}
$$
其中，$\pi k_{n} \tau^{2}$这一项很小，可以忽略，因此上式可以简化为：
$$
\begin{aligned} 
s_{\text{IF}}\left ( t;f_{0,n} ,k_{n}\right ) &=\sum_{n=1}^{N} A\cdot e^{j\left \{2\pi f_{0,n}\tau +2\pi k_{n}\tau [t-(n-1)T_\text{c}]\right \}}
e^{j\frac{4 \pi(n-1) v T_\text{c}f_{c}}{c}}  
\\ &\cdot \Pi (t-(n-1)T_\text{c})
\end{aligned}
$$

以周期$T_\text{s}$采样，将信号离散化处理，将每个脉冲内的采样点分别排列为一行，构成共$N_\text{c}$行、$N_\text{s}$列的矩阵$\mathbf{S}$，其中第$n$行$m$列对应的值为：
$$
\mathbf{S}_{n,m} =\ A\cdot e^{j\left [2\pi f_{0,n}\tau +2\pi k_{n}\tau T_\text{s}(m-1) \right ]}
e^{j\frac{4 \pi(n-1) v T_\text{c}f_{c}}{c}}  
$$

#### 雷达子系统
##### 测距原理
在单个脉冲周期内，由于$T_{c}$极小，故可以近似认为物体在一个脉冲周期内静止，则发射波经过时延$\tau$后，形成回波并被雷达接收。有下式成立：
$$
R=\frac{c\tau }{2} ,0< \tau < T_{c}
$$
雷达中的混频器将发射波与反射波信号相乘，再低通滤波，形成混频信号。从频域上看，混频信号的频率相当于发射波频率减去反射波频率。则混频信号的频率$f_\text{IF}$满足：
$$
f_\text{IF} =k _{n}\tau
$$
于是问题的关键在于通过频率计测得混频信号的频率$f_\text{IF}$。每个脉冲内，混频信号的频率是不变的，我们当然可以直接做FFT，根据峰值出现的位置判断频率。但是，不同于常规FMCW雷达，若将$n$个脉冲的所有$N_\text{s}$个采样点构成信号$x_{n}\left [ m \right ] $：
$$
x_{n}\left [ m \right ] = A\cdot e^{j\left [2\pi f_{0,n}\tau +2\pi k_{n}\tau T_\text{s}(m-1) \right ]} 
 e^{j\frac{4 \pi v T_\text{c}f_{c}(n-1)}{c}}
$$
对其直接做FFT：
$$
  \begin{split} 
x\left [ k \right ] &=\sum_{m=1}^{N_\text{s} } x_{n}\left [ m \right ] \cdot 
e^{-j\frac{2\pi k\left ( m-1 \right ) }{N_\text{s}} } \\
&= Ae^{j\frac{4 \pi v T_\text{c}f_{c}(n-1)}{c}}\cdot 
\sum_{m=1}^{N_\text{s} }e^{j\left [2\pi f_{0,n}\tau +2\pi k_{n}\tau T_\text{s}(m -1)\right ]}  
\cdot e^{-j\frac{2 \pi(m-1) k}{N_\text{s} }} 
\end{split}
$$
这样不能抵消脉冲间特定的随机相位$e^{j2\pi f_{0,n}\tau }$，也会因频率项$e^{j2\pi k_{n}\tau T_\text{s}m }$的频率有正有负，导致1D-FFT后峰值出现的位置不在一条直线上，进而无法实现能量的累积，对后续测速造成困扰。

本文提出一种改进MTD算法，在测距阶段，对一维信号$x\left [ m \right ] $做改进FFT：
$$
  \begin{split} 
x\left [ k \right ] &=\sum_{m=1}^{N_\text{s} } x\left [ m \right ] \cdot 
e^{-j\frac{2\pi k\left ( m-1 \right ) }{N_\text{s}}\cdot \frac{k_n}{\left | k_n \right | }}  
 e^{-j\frac{2\pi f_{0,n} k}{T_\text{s} N_\text{s}\left | k_n \right |}}  \\
&= Ae^{j\frac{4 \pi v T_\text{c}f_{c}(n-1)}{c}}  \cdot \sum_{m=1}^{N_\text{s} }
 e^{j2\pi k_{n}\tau T_\text{s}(m-1) } e^{-j\frac{2\pi k\left ( m-1 \right ) }{N_\text{s}}\cdot \frac{k_n}{\left | k_n \right | }} 
\cdot e^{j2\pi f_{0,n}\tau}  e^{-j\frac{2\pi f_{0,n} k}{T_\text{s} N_\text{s}\left | k_n \right |}}
\end{split}
$$
其中， $e^{-j\frac{2\pi k\left ( m-1 \right ) }{N_\text{s}}\cdot \frac{k_n}{\left | k_n \right | }} $项中$\frac{k_n}{\left | k_n \right | }$是为了抵消频率的同时使峰值出现在距离维的同一点处，即$k=\tau T_\text{s}N_\text{s}\left | k_n \right | $点处，此时 $  e^{-j\frac{2\pi f_{0,n} k}{T_\text{s} N_\text{s}\left | k_n \right |}}=e^{-j2\pi f_{0,n}\tau} $，恰好在这一点处将脉冲的初始相位抵消。

对混频信号的每一个脉冲作同样的处理可得到矩阵$\tilde{S}_{nk}$。如果用三维张量 $\mathbf{F} \in \mathbb{C}^{N_s \times N_s \times N_c}$ 与矩阵 $\mathbf{S} \in \mathbb{C}^{N_c \times N_s}$ 的模态乘积运算来描述整个过程，我们可以使用爱因斯坦求和约定更简洁地表示算法：
$$
\tilde{S}_{nk} = S_{nm} F_{mkn}
$$
其中：
- $\tilde{S}_{nk}$ 为输出矩阵 $\tilde{\mathbf{S}}$ 的元素
- $S_{nm}=A\cdot e^{j\left [2\pi f_{0,n}\tau +2\pi k_{n}\tau T_\text{s}(m-1) \right ]} 
 e^{j\frac{4 \pi v T_\text{c}f_{c}(n-1)}{c}}$ 为输入矩阵 $\mathbf{S}$ 的元素
- $F_{mkn}=e^{-j\frac{2\pi k\left ( m-1 \right ) }{N_\text{s}}\cdot \frac{k_n}{\left | k_n \right | }}  
 e^{-j\frac{2\pi f_{0,n} k}{T_\text{s} N_\text{s}\left | k_n \right |}}$ 为张量 $\mathbf{F}$ 的元素
- 重复下标 $m$ 表示对其求和

该运算等效于：
$$
\tilde{\mathbf{S}} = \mathbf{S} \times_2 \mathbf{F}
$$
其中，$\times_2$ 表示沿第二个模态的张量-矩阵乘积。

##### 测速原理
根据上述运算得到1D-FFT矩阵$\tilde{\mathbf{S}}$，其维度为$N_c \times N_s$，其中$N_c$为脉冲数，$N_s$为距离单元数。此时，假设在第$m_r$行是目标距离对应的索引，那么理想情况下这一行的与速度无关的相位就已经抵消，于是：
$$
\tilde{\mathbf{S}}_{n,m_r} =A'e^{j\frac{4 \pi v T_\text{c}f_{c}(n-1)}{c}}  ,\quad n = 1,2,\ldots,N_c
$$
其中$A'$是其幅度。

将第$m_r$行的抽取出来作为慢时间信号序列，并令：
$$
x_{m_r}[n] =\tilde{\mathbf{S}}_{n,m_r} =B e^{j\frac{4 \pi v T_\text{c}f_{c}(n-1)}{c}}  
$$
若希望获得速度，可对信号$x_{m_r}[n] $直接做FFT，即对矩阵$\tilde{\mathbf{S}}$的每一列作FFT：
$$
\tilde{\tilde{\mathbf{S}}} = \mathcal{F} \{\tilde{\mathbf{S} } \}
$$
其中，$\mathcal{F} \{\cdot \}$运算符表示对矩阵的每一列作FFT。

#### 通信子系统
接收端接收的信号为：
$$
    s_\text{rec}^\text{com}=\alpha s_\text{tra}( t ;f_{0,n},k_{n})
$$
其中，$\alpha$是信号衰减因子。

先对信号解调并降采样：
$$
s_{demodulate}(t,f_{0,n},k_n)=s_\text{rec}^\text{com}(t,f_{0,n},k_n)\cdot e^{-j2\pi f_ct}
$$
其中，$f_{c}$表示信号的载波频率，即中心频率。

之后，为获得当前发射波信号的斜率，采用短时傅里叶变换，选用合适的窗口长度与FFT点数，计算每个码元的斜率$\widehat{k_n}$，之后设置门限$k_{\tau}$，与设定的斜率$k_n$比较，如果下式成立：
$$
  \left| \frac{\widehat{k_n}-k_n}{k_n} \right |  \leq k_{\tau}  
$$
则认为该码元的斜率为$k_n$，由此解码。

## 致谢
- 该大创部分文档撰写、代码编写由项目团队共同完成，时间原因未在此项目中开源，特此感谢项目组全体成员的协作；
- 感谢电子科技大学崔国龙老师的指导与支持；
- 来源：电子科技大学大学生创新创业训练计划（2024/12-2025/10）。

## 免责声明
- 本仓库为电子科技大学大学生创新创业训练计划归档版本，核心的改进MTD算法为独立开发成果；
- 项目代码仅供学术交流使用，禁止未经授权的商业应用，商业使用需联系作者获得授权。