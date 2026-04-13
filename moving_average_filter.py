# moving_average_filter.py
import numpy as np

class MovingAverageFilter:
    """
    滑动平均滤波器（循环缓冲区实现，O(1) 时间复杂度）

    公式：
        y[k] = (1 / n) * sum_{i=0}^{n-1} x[k-i]

    参数：
        window_size : 窗口长度 n，越大噪声抑制越强但响应越慢
    """

    def __init__(self, window_size):
        if window_size <= 0:
            raise ValueError("window_size must be a positive integer")
        self.window_size = window_size
        self.buffer = [0.0] * window_size
        self.index = 0          # 当前写入位置
        self.count = 0          # 已累计样本数（用于前 window_size 个样本的平均）
        self.sum = 0.0          # 窗口内元素之和

    def update(self, value):
        """
        输入新采样值，返回滤波后的值

        参数：
            value : 当前时刻的采样值

        返回：
            filtered_value : 滑动平均后的输出
        """
        # 减去将被覆盖的旧值
        old_value = self.buffer[self.index]
        self.sum -= old_value

        # 存入新值
        self.buffer[self.index] = value
        self.sum += value

        # 移动索引
        self.index = (self.index + 1) % self.window_size

        # 计算平均值（在窗口未满时用实际样本数）
        if self.count < self.window_size:
            self.count += 1
            return self.sum / self.count
        else:
            return self.sum / self.window_size

    def reset(self):
        """重置滤波器状态"""
        self.buffer = [0.0] * self.window_size
        self.index = 0
        self.count = 0
        self.sum = 0.0

    def get_buffer(self):
        """返回当前缓冲区内容的副本（按时间顺序，最旧的在前）"""
        if self.count == 0:
            return []
        ordered = []
        start = self.index if self.count == self.window_size else 0
        for i in range(self.count):
            idx = (start + i) % self.window_size
            ordered.append(self.buffer[idx])
        return ordered


# ========== 针对三轴加速度计的预配置类 ==========
class AccelerometerMAF:
    """
    三轴加速度计滑动平均滤波器

    对每个轴的加速度数据分别进行滑动平均滤波，适用于抑制扑翼飞行器的周期性振动。
    """

    def __init__(self, window_size):
        self.filter_x = MovingAverageFilter(window_size)
        self.filter_y = MovingAverageFilter(window_size)
        self.filter_z = MovingAverageFilter(window_size)

    def update(self, ax, ay, az):
        """
        输入三轴原始加速度值，返回滤波后的三轴加速度

        参数：
            ax, ay, az : 原始加速度计测量值（m/s² 或归一化值）

        返回：
            (ax_f, ay_f, az_f) : 滤波后的三轴加速度
        """
        return (self.filter_x.update(ax),
                self.filter_y.update(ay),
                self.filter_z.update(az))

    def reset(self):
        self.filter_x.reset()
        self.filter_y.reset()
        self.filter_z.reset()


# ========== 示例：仿真验证 ==========
if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # 生成带噪声的正弦信号（模拟加速度计受振动干扰）
    dt = 0.001           # 1 kHz 采样
    t = np.arange(0, 2.0, dt)
    freq = 50.0          # 振动基频 50 Hz
    amplitude = 1.0

    # 真实信号（缓慢变化的姿态相关加速度）
    true_signal = 0.5 * np.sin(2 * np.pi * 0.5 * t)

    # 高频振动噪声
    vibration = 2.0 * np.sin(2 * np.pi * freq * t)

    # 随机噪声
    noise = np.random.normal(0, 0.2, len(t))

    # 原始测量信号
    raw_signal = true_signal + vibration + noise

    # 应用滑动平均滤波（窗口大小 n=50，对应 50ms）
    maf = MovingAverageFilter(window_size=50)
    filtered_signal = np.array([maf.update(x) for x in raw_signal])

    # 绘图对比
    plt.figure(figsize=(12, 5))

    plt.subplot(2, 1, 1)
    plt.plot(t, true_signal, 'k-', linewidth=2, label='True signal')
    plt.plot(t, raw_signal, 'b-', alpha=0.4, linewidth=0.5, label='Raw measurement')
    plt.plot(t, filtered_signal, 'r-', linewidth=1.5, label='Moving average filtered')
    plt.ylabel('Acceleration')
    plt.title(f'Moving Average Filter (window_size = {maf.window_size})')
    plt.legend()
    plt.grid(True, linestyle=':')

    # 局部放大（前 0.2 秒）
    plt.subplot(2, 1, 2)
    mask = t < 0.2
    plt.plot(t[mask], true_signal[mask], 'k-', linewidth=2, label='True')
    plt.plot(t[mask], raw_signal[mask], 'b-', alpha=0.4, label='Raw')
    plt.plot(t[mask], filtered_signal[mask], 'r-', linewidth=1.5, label='Filtered')
    plt.xlabel('Time (s)')
    plt.ylabel('Acceleration')
    plt.legend()
    plt.grid(True, linestyle=':')
    plt.title('Zoomed view (first 0.2 s)')

    plt.tight_layout()
    plt.show()