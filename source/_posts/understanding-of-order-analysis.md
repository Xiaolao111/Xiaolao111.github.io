---
title: 对阶次分析的一些理解
date: 2024-08-30 14:48:42
author: xiaolao
math: true
index_img: \img\blog_order\Figure_3.png
tags:
- FFT
- Python
---

本文是我在项目过程中，对振动加速度数据的阶次分析的一些理解

# 频谱分析
## 原理
在不考虑数据滤波的情况下，频谱分析为时域数据作FFT的结果，其中纵轴为振动幅度（常见单位有g，mm/s^2^等），横轴为频率，其中，能够采集到的最大频率为F~max~，频谱分辨率为F~n~，采样率为F~s~，采样点数为N
$$
F_{max}=\frac{Fs}{2},     F_n=\frac{Fs}{N}
$$
示例：
```
Fs = 1024
Fb = 10
x = np.arange(0, 1024)
t = x / Fb
data_len = 1024
y = np.sin(2 * np.pi * x / Fs * Fb)
fft_y = np.fft.fft(y)
fft_y = (np.abs(fft_y) / data_len)[range(int(data_len / 2))] * 2
fft_y[0] = fft_y[0] / 2
fig, ax = plt.subplots(2, 1)
line_0 = ax[0].plot(t, y)
ax[0].set_xlabel('time [s]')
line_1 = ax[1].plot(fft_y)
ax[1].set_xlabel('Frequency [Hz]')

mplcursors.cursor(line_0, hover=False, multiple=True).connect("add", lambda sel: sel.annotation.set_text(
    f"x={sel.target[0]:.2f}, y={sel.target[1]:.2f}"))
mplcursors.cursor(line_1, hover=False, multiple=True).connect("add", lambda sel: sel.annotation.set_text(
    f"x={sel.target[0]:.2f}, y={sel.target[1]:.2f}"))
plt.show()
```
![FFT演示](\img\blog_order\Figure_1.png)

# 阶次分析
阶次分析和频谱分析有完全不同的应用场景，阶次分析的好处是，一个转动系统，转速往往一直在变化，如果做频谱分析，我们会发现频谱会随转速变化，但是阶次分析与转速无关。频率分析适合系统中有固定频率存在的情况，而阶次分析适合系统中各部件互相耦合，转速一致变化的情况

## 基于时域采样的阶次分析
假设基准频率固定为*F*Hz，转速为*S*rpm，易得：
$$
F=\frac{S}{60}
$$
传统的频谱分析也是基于时域采样，所以，频谱到阶次谱只需要：
```
#  fft_y不变
fft_x = fft_x / F
```
这样，能够采集到的最大阶次为
$$
order_{max} =\frac{ F_{s}}{2* F}
$$
阶次分辨率为
$$
order_{n} =\frac{Fs}{N*F}
$$
假设我们在对一个传动系统的振动做阶次分析，以输入轴转速为基准转速，转速大约在300rpm到1500rpm，传感器的采集频率为5.12kHz，采集时间为0.5s，那么，在300rpm时，我们能计算的阶次范围为0-5120阶，分辨率为0.6阶，在1500rpm时，阶次范围为0-1024阶，分辨率为0.08阶
可以看出，当转速变化时，最大阶次和分辨率也会跟着变化，这样其实不利于分析和比较，而且时域采样的一个重要缺点是，无法处理频率随时间变化的非平稳信号，否则会有严重的拖尾效应，所以更多的阶次分析，会采用角度域采样而非时域采样

非平稳信号的拖尾效应示例：
```
Fb = 10
x = np.arange(0, 1024)
t = x / Fb
data_len = 1024
order_x = np.linspace(0, int(data_len / 2 / Fb), int(data_len / 2))
instantaneous_frequency = Fb + np.linspace(0, 20, data_len )
y = np.sin(2 * np.pi * np.cumsum(instantaneous_frequency) / len(t))
fft_y = np.fft.fft(y)
fft_y = (np.abs(fft_y) / data_len)[range(int(data_len / 2))] * 2
fft_y[0] = fft_y[0] / 2
fig, ax = plt.subplots(3, 1)
line_0 = ax[0].plot(t, y)
ax[0].set_xlabel('Time [s]')
line_1 = ax[1].plot(fft_y)
ax[1].set_xlabel('Frequency [Hz]')
line_2 = ax[2].plot(order_x, fft_y)
ax[2].set_xlabel('Order')
mplcursors.cursor(line_0, hover=False, multiple=True).connect("add", lambda sel: sel.annotation.set_text(
    f"x={sel.target[0]:.2f}, y={sel.target[1]:.2f}"))
mplcursors.cursor(line_1, hover=False, multiple=True).connect("add", lambda sel: sel.annotation.set_text(
    f"x={sel.target[0]:.2f}, y={sel.target[1]:.2f}"))
mplcursors.cursor(line_2, hover=False, multiple=True).connect("add", lambda sel: sel.annotation.set_text(
    f"x={sel.target[0]:.2f}, y={sel.target[1]:.2f}"))
plt.show()
```
![非平稳信号的分析示例](\img\blog_order\Figure_2.png)

## 基于角度域采样的阶次分析
角度域采样的原理是，选择旋转部件，该部件旋转一周，等间隔采样N个点，此时N为角度域采样的采样频率，对角度域信号做FFT，可得到基于旋转部件转速的阶次谱

### 采样时间内转速无变化
假设旋转部件的转动频率为F不变，转过一圈采样N个点，每采N×n做一次FFT，由于转速不变，可知时域采样率为F * N，FFT后的频谱，点数为N×n/2，最大频率为F*N/2，我们把频谱的横轴除以旋转部件的频率，得到阶次谱，阶次谱的最大阶次为Order~max~，谱线分辨率为Order~n~，可知：
$$
order_{max} =\frac{ N}{2} 
$$
$$
order_{n} =\frac{order_{max} }{\frac{N×n}{2}}  =\frac{1}{n}
$$
可以看到，相比于时域采样，角度域采样之后的阶次谱，最大阶次和阶次分辨率都与转速无关了，我们可以通过增大采样率N来扩大阶次谱范围，通过增大采样次数n来减小分辨率，提升精度
![角度域采样举例](\img\blog_order\Figure_3.png)

### 采样时间内转速变化
我的理解：阶次分析时需要把FFT后的频率轴转为阶次轴，也就是除以转频F，那么这里转频需要是一个固定的量，目前看到的处理方法基本为求采样时间内转速的平均值，作为该段时间的转速。如果在采样时间内转速变化很大，那么会导致求取阶次时误差也比较大

对于角度域采样，有*硬件*和*软件*两种处理方式，硬件处理是每转过恒定角度时触发采样，得到的原始数据就是角度域数据（但是不知道这种方法是否方便切换采样率）；软件处理指先采集时域信号，上位机端结合转速信号对时域信号进行基于角度域的重采样，这种做法好处是采样率很好控制，坏处是软件端的工作量会比较大
