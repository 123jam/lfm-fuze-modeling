# LFM无线电引信数学建模基础

## 1. LFM信号数学模型

### 1.1 基本定义
线性调频（Linear Frequency Modulation, LFM）信号的复包络表示：

$$s(t) = \begin{cases}
A \cdot \exp\left(j2\pi\left(f_0 t + \frac{K}{2}t^2\right)\right) & 0 \leq t \leq T \\
0 & \text{otherwise}
\end{cases}$$

**参数说明：**
- $A$：信号幅度
- $f_0$：载波频率（Hz）
- $K = B/T$：调频斜率（Hz/s），其中 $B$ 为带宽，$T$ 为脉冲宽度
- $T$：脉冲持续时间（s）

**关键特性：**
$$\text{瞬时频率：} f_i(t) = f_0 + Kt$$

---

## 2. 目标回波建模

### 2.1 单个目标回波模型

假设目标距离为 $R$，径向速度为 $v$，则接收的目标回波为：

$$r(t) = \alpha \cdot s(t - \tau) \cdot \exp(j2\pi\Delta f \cdot t)$$

**其中：**
- $\tau = \frac{2R}{c}$：往返延迟（$c$ 为光速）
- $\alpha$：衰减系数（雷达方程确定）
- $\Delta f = \frac{2v f_0}{c}$：多普勒频移（一阶近似）

**展开形式：**
$$r(t) = \alpha \cdot \exp\left(j2\pi\left(f_0(t-\tau) + \frac{K}{2}(t-\tau)^2 + \Delta f \cdot t\right)\right)$$

$$= \alpha \cdot \exp(j\phi_0) \cdot \exp\left(j2\pi\left((f_0+\Delta f)t + \frac{K}{2}t^2 - \left(f_0\tau + K\tau t + \frac{K}{2}\tau^2\right)\right)\right)$$

**其中初始相位：**
$$\phi_0 = 2\pi f_0 \tau + \pi K \tau^2$$

### 2.2 衰减系数（雷达方程）

$$\alpha = \frac{\lambda^2 \sigma}{(4\pi)^3 R^4} \cdot G^2$$

- $\lambda = c/f_0$：波长
- $\sigma$：目标散射截面积（m²）
- $G$：天线增益
- $R$：目标距离

为简化起见，通常采用：
$$\alpha = \frac{1}{R^2} \text{ 或 } \alpha = \frac{1}{R^{2.5}}$$

---

## 3. 多普勒效应与距离-速度耦合

### 3.1 接收信号的时间-频率表示

在距离-速度同时变化的情况下，接收信号可表示为：

$$r(t) = \alpha \cdot \exp\left(j2\pi\left(f_d t + f_0(t-\tau) + \frac{K}{2}(t-\tau)^2\right)\right)$$

其中 $f_d = \Delta f = \frac{2vf_0}{c}$

### 3.2 距离-多普勒参数关系

- **距离分辨率：** $\Delta R = \frac{c}{2B}$
- **多普勒/速度分辨率：** $\Delta v = \frac{c}{2f_0 T}$
- **距离-速度耦合因子：** 同一个目标的 $(R, v)$ 会在不同的时频点产生响应

---

## 4. 多脉冲回波模型

在实际引信中，通常发射多个脉冲形成脉冲列。设发射周期为 $T_r$（重复周期），共 $N_p$ 个脉冲，则：

### 4.1 发射脉冲列

$$s_{\text{tx}}(t) = \sum_{n=0}^{N_p-1} s(t - nT_r) \cdot w_n$$

**其中：**
- $w_n$：第 $n$ 个脉冲的幅度加权（通常 $w_n=1$）

### 4.2 接收脉冲列

$$r_{\text{rx}}(t) = \sum_{n=0}^{N_p-1} \alpha_n \cdot s(t - nT_r - \tau) \cdot \exp(j2\pi\Delta f (nT_r + t))$$

**简化形式（假设短脉冲、衰减缓变）：**
$$r_{\text{rx}}(t) \approx \sum_{n=0}^{N_p-1} \alpha \cdot s(t - nT_r - \tau) \cdot \exp(j2\pi(\Delta f \cdot t + n\Delta f \cdot T_r))$$

---

## 5. 带噪接收信号模型

### 5.1 高斯白噪声

$$n(t) \sim \mathcal{CN}(0, \sigma_n^2)$$

其功率谱密度：$P_n = \sigma_n^2$

### 5.2 信噪比定义

$$\text{SNR} = \frac{E[|r(t)|^2]}{E[|n(t)|^2]} = \frac{|\alpha|^2}{2\sigma_n^2}$$

**以dB表示：**
$$\text{SNR(dB)} = 10\log_{10}(\text{SNR})$$

### 5.3 完整接收信号

$$y(t) = r(t) + n(t)$$

---

## 6. 匹配滤波与脉冲压缩

### 6.1 匹配滤波器设计

最优检测器的冲激响应为发射信号的反褶复共轭：

$$h(t) = s^*(-t) = A^* \cdot \exp\left(-j2\pi\left(-f_0 t + \frac{K}{2}t^2\right)\right)$$

或因果形式（延迟 $T$）：
$$h(t) = s^*(T-t), \quad 0 \leq t \leq T$$

### 6.2 脉冲压缩输出

$$y_{\text{comp}}(t) = \int_{-\infty}^{\infty} y(\tau) h^*(t-\tau) d\tau = y(t) * h(t)$$

### 6.3 压缩后的目标回波（理想情况）

设接收到单个延迟为 $\tau$ 的目标回波，经匹配滤波后，在 $t = \tau$ 处产生尖峰（距离脉冲）：

$$y_{\text{comp}}(\tau) = \sum_{n=0}^{N_p-1} \alpha \cdot \text{sinc}\left(\pi B (t - \tau - nT_r)\right) \cdot \exp(j2\pi\Delta f(t + nT_r))$$

**主瓣宽度（3dB）：** $\Delta t_{3dB} = \frac{1.2}{B}$

**对应距离分辨率：** $\Delta R = \frac{c}{2B}$

---

## 7. 距离-多普勒处理（FFT处理）

### 7.1 脉冲间相���处理

对多个脉冲的压缩输出进行FFT（沿脉冲维度）：

$$Y(m, k) = \sum_{n=0}^{N_p-1} y_{\text{comp}}(m, n) \cdot \exp\left(-j2\pi \frac{kn}{N_p}\right)$$

**其中：**
- $m$：距离单元索引
- $k$：多普勒单元索引（对应速度）
- $N_p$：脉冲数

### 7.2 距离-多普勒图（Range-Doppler Map）

$$Z(m, k) = |Y(m, k)|^2 \quad \text{(功率图)}$$

或

$$Z_{\text{dB}}(m, k) = 10\log_{10}\frac{|Y(m, k)|^2}{|Y_{\max}|^2} \quad \text{(归一化dB图)}$$

---

## 8. 模糊函数（Ambiguity Function）

### 8.1 定义

$$\chi(\tau, \nu) = \left|\int_{-\infty}^{\infty} s(t) s^*(t-\tau) \exp(j2\pi\nu t) dt\right|^2$$

**参数含义：**
- $\tau$：时间延迟（对应距离偏差）
- $\nu$：频率偏移（对应速度偏差）
- $\chi(\tau, \nu)$：相关函数的平方（功率）

### 8.2 LFM信号的模糊函数特性

对于LFM信号 $s(t) = \exp(j2\pi(f_0 t + \frac{K}{2}t^2))$，其模糊函数具有如下特点：

**（1）主瓣特性：**
- 距离维主瓣宽度：$\Delta \tau_{3dB} \approx \frac{1.2}{B}$
- 多普勒维主瓣宽度：$\Delta \nu_{3dB} \approx \frac{1.2}{T}$
- 主瓣积分体积守恒：$\int\int \chi(\tau, \nu) d\tau d\nu = T^2$（Parseval定理）

**（2）旁瓣特性：**
- 最高旁瓣电平（PSL）：$\text{PSL} \approx -13\text{ dB}$（对于矩形窗LFM）
- 旁瓣衰减率：-6 dB/octave（距离维）

**（3）Chirp Rate与分辨率关系：**
$$\text{距离分辨率} = \frac{c}{2B}$$
$$\text{速度分辨率} = \frac{c}{2f_0 T}$$

---

## 9. 三类欺骗干扰数学模型

### 9.1 DRFM（数字射频存储）干扰

**原理：** 对接收信号进行存储、延迟后重放

$$i_{\text{DRFM}}(t) = \beta \cdot y(t - t_d) \cdot \exp(j\phi_{\text{offset}})$$

**其中：**
- $\beta$：干扰幅度衰减因子
- $t_d$：重放延迟（可改变目标距离）
- $\phi_{\text{offset}}$：相位偏移（可为0，保持原相位）

**特点：** 完全保留信号调制结构，难以识别

### 9.2 RGPO（距离门脉冲振荡）干扰

**原理：** 产生多个距离门内的虚假回波脉冲

$$i_{\text{RGPO}}(t) = \sum_{k=0}^{N_g-1} A_k \cdot \text{rect}\left(\frac{t - t_{g,k}}{\Delta t_g}\right) \cdot \exp(j\phi_k)$$

**参数：**
- $N_g$：干扰门数（通常5~10个）
- $t_{g,k}$：第$k$个干扰门时刻
- $\Delta t_g$：干扰门宽度（通常等于脉冲压缩后的距离单元宽度）
- $\phi_k$：第$k$个门的相位（可随机变化）

**频域形式：**
$$I_{\text{RGPO}}(f) = \sum_{k=0}^{N_g-1} A_k \cdot \text{sinc}(f\Delta t_g) \cdot \exp(-j2\pi f t_{g,k})$$

### 9.3 VGPO（速度门脉冲振荡）干扰

**原理：** 产生多个多普勒频移分量，填充整个速度维

$$i_{\text{VGPO}}(t) = \sum_{m=0}^{N_v-1} A_m \cdot s(t) \cdot \exp(j2\pi f_m t)$$

**参数：**
- $N_v$：速度门数
- $f_m = f_{\text{center}} + m \cdot \Delta f_v$：多普勒频率
- $\Delta f_v = 1/T_r$：速度门间隔（由重复周期决定）

**等效表示：**
$$i_{\text{VGPO}}(t) = A(t) \cdot s(t) \cdot \text{comb}(f_m, \Delta f_v)$$

其中 $\text{comb}$ 为频率梳形函数

---

## 10. 干扰信噪比（JNR）定义

### 10.1 干扰功率

$$P_j = E[|i(t)|^2]$$

### 10.2 干扰信噪比

$$\text{JNR} = \frac{P_j}{P_s + P_n} = \frac{P_j}{P_s(1 + 1/\text{SNR})}$$

通常近似：
$$\text{JNR(dB)} = 10\log_{10}(P_j / P_s)$$

### 10.3 综合信干噪比（SINR）

$$\text{SINR} = \frac{P_s}{P_j + P_n}$$

---

## 11. 早炸机制分析

### 11.1 引信判决算法框架

```
输入：距离-多普勒图 Z(m,k)，距离-速度轴
过程：
  1. 在Z中搜索峰值
  2. 提取峰值对应的距离 R_det，速度 v_det
  3. 判决：
     IF R_det ≤ R_target - R_safe THEN
        引信起爆（判断为"到达引爆距离"）
        → 如果此时还未到达真实引爆距离，则早炸
     END IF
输出：起爆/不起爆 决策
```

### 11.2 早炸威胁评估

定义**早炸余度**：
$$\text{Margin} = R_{\text{detected,min}} - (R_{\text{target}} - R_{\text{safe}})$$

**含义：**
- $\text{Margin} > 0$：检测到的虚假目标在安全距离内 → **存在早炸风险**
- $\text{Margin} < 0$：检测正常 → **安全**

---

## 12. 符号总结表

| 符号 | 含义 | 单位 |
|------|------|------|
| $s(t)$ | 发射信号 | - |
| $r(t)$ | 目标回波 | - |
| $n(t)$ | 噪声 | - |
| $i(t)$ | 干扰信号 | - |
| $A, \alpha$ | 幅度/衰减系数 | - |
| $f_0$ | 载频 | Hz |
| $B$ | 带宽 | Hz |
| $K$ | 调频斜率 | Hz/s |
| $T$ | 脉冲宽度 | s |
| $T_r$ | 脉冲重复周期 | s |
| $R$ | 目标距离 | m |
| $v$ | 目标速度 | m/s |
| $c$ | 光速 | m/s |
| $\tau$ | 往返延迟 | s |
| $\Delta f$ | 多普勒频移 | Hz |
| $\sigma_n^2$ | 噪声功率 | - |
| $\chi(\tau, \nu)$ | 模糊函数 | - |

---

## 参考

- Skolnik, M. I. (2008). *Radar Handbook* (3rd ed.).
- Richards, M. A. (2014). *Fundamentals of Radar Signal Processing*.
- Mahafza, B. R. (2000). *Radar Systems Analysis and Design using MATLAB*.
