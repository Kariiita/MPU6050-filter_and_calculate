# kalman_filter.py
import numpy as np

class KalmanFilter:
    """
    离散线性卡尔曼滤波器

    状态空间模型：
        x[k] = A * x[k-1] + B * u[k-1] + w[k-1]      (过程方程)
        z[k] = H * x[k] + v[k]                        (测量方程)

    其中：
        w ~ N(0, Q)   过程噪声
        v ~ N(0, R)   测量噪声

    使用方法：
        kf = KalmanFilter(A, B, H, Q, R, x_init, P_init)
        kf.predict(u)          # 先验估计
        kf.update(z)           # 后验估计（测量更新）
        x_est = kf.x           # 当前状态估计
    """

    def __init__(self, A, B, H, Q, R, x_init=None, P_init=None):
        """
        初始化卡尔曼滤波器

        参数：
            A : ndarray (n, n)   状态转移矩阵
            B : ndarray (n, m)   控制输入矩阵（若无控制输入则设为 None）
            H : ndarray (p, n)   观测矩阵
            Q : ndarray (n, n)   过程噪声协方差矩阵
            R : ndarray (p, p)   测量噪声协方差矩阵
            x_init : ndarray (n, 1)  初始状态估计（默认零向量）
            P_init : ndarray (n, n)  初始估计误差协方差（默认单位阵）
        """
        self.A = np.array(A, dtype=float)
        self.B = np.array(B, dtype=float) if B is not None else None
        self.H = np.array(H, dtype=float)
        self.Q = np.array(Q, dtype=float)
        self.R = np.array(R, dtype=float)

        self.n = self.A.shape[0]          # 状态维数
        self.p = self.H.shape[0]          # 测量维数

        # 初始状态
        if x_init is None:
            self.x = np.zeros((self.n, 1))
        else:
            self.x = np.array(x_init, dtype=float).reshape(-1, 1)

        # 初始误差协方差
        if P_init is None:
            self.P = np.eye(self.n)
        else:
            self.P = np.array(P_init, dtype=float)

    def predict(self, u=None):
        """
        预测步骤（时间更新）

        参数：
            u : ndarray (m, 1)  控制输入（若无控制则设为 None）

        返回：
            x_prior : 先验状态估计
            P_prior : 先验误差协方差
        """
        # 状态预测
        self.x = self.A @ self.x
        if u is not None and self.B is not None:
            u = np.array(u, dtype=float).reshape(-1, 1)
            self.x += self.B @ u

        # 协方差预测
        self.P = self.A @ self.P @ self.A.T + self.Q

        return self.x.copy(), self.P.copy()

    def update(self, z):
        """
        更新步骤（测量更新）

        参数：
            z : ndarray (p, 1)  测量向量

        返回：
            x_post : 后验状态估计
            P_post : 后验误差协方差
        """
        z = np.array(z, dtype=float).reshape(-1, 1)

        # 卡尔曼增益
        S = self.H @ self.P @ self.H.T + self.R
        K = self.P @ self.H.T @ np.linalg.inv(S)

        # 状态修正
        y = z - self.H @ self.x          # 测量残差
        self.x = self.x + K @ y

        # 协方差修正
        I = np.eye(self.n)
        self.P = (I - K @ self.H) @ self.P

        return self.x.copy(), self.P.copy()

    def filter(self, z, u=None):
        """
        执行一步完整的预测+更新

        参数：
            z : 测量向量
            u : 控制输入（可选）

        返回：
            后验状态估计
        """
        self.predict(u)
        self.update(z)
        return self.x

    def get_state(self):
        """返回当前状态估计值（一维数组）"""
        return self.x.flatten()

    def get_covariance(self):
        """返回当前误差协方差矩阵"""
        return self.P.copy()


# ========== 针对陀螺仪角速度估计的预配置类 ==========
class GyroKalmanFilter(KalmanFilter):
    """
    用于陀螺仪角速度估计与零偏补偿的卡尔曼滤波器

    状态向量 x = [omega, bias]^T
        omega : 真实角速度 (rad/s)
        bias  : 陀螺仪零偏 (rad/s)

    过程模型：
        omega[k] = omega[k-1] - bias[k-1] * dt + w_omega   (假设零偏缓慢变化)
        bias[k]  = bias[k-1] + w_bias

    测量模型：
        z[k] = omega[k] + bias[k] + v

    参数：
        dt : 采样周期 (s)
        q_omega : 角速度过程噪声方差
        q_bias  : 零偏过程噪声方差
        r_meas  : 测量噪声方差
    """

    def __init__(self, dt, q_omega=1e-4, q_bias=1e-8, r_meas=1e-3,
                 omega_init=0.0, bias_init=0.0):
        # 状态转移矩阵 A
        A = np.array([[1.0, -dt],
                      [0.0,  1.0]])

        # 控制矩阵 B（无外部控制输入）
        B = None

        # 观测矩阵 H
        H = np.array([[1.0, 1.0]])

        # 过程噪声协方差 Q
        Q = np.diag([q_omega, q_bias])

        # 测量噪声协方差 R
        R = np.array([[r_meas]])

        # 初始状态
        x_init = np.array([[omega_init], [bias_init]])

        # 初始误差协方差（可适当调大以加速收敛）
        P_init = np.eye(2) * 10.0

        super().__init__(A, B, H, Q, R, x_init, P_init)
        self.dt = dt

    def estimate_omega(self, z_meas):
        """
        输入原始陀螺仪测量值（rad/s），返回滤波后的角速度估计值

        参数：
            z_meas : 陀螺仪测量值 (rad/s)

        返回：
            omega_est : 滤波后的角速度 (rad/s)
            bias_est  : 估计的零偏 (rad/s)
        """
        z = np.array([[z_meas]])
        self.filter(z)
        omega = self.x[0, 0]
        bias  = self.x[1, 0]
        return omega, bias


# ========== 示例：如何使用 ==========
if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # 模拟参数
    dt = 0.01
    T = 10.0
    steps = int(T / dt)
    t = np.linspace(0, T, steps)

    # 真实角速度（正弦变化）
    omega_true = 0.5 * np.sin(2 * np.pi * 0.2 * t)

    # 陀螺仪零偏（缓慢漂移）
    bias_true = 0.1 + 0.05 * np.sin(2 * np.pi * 0.01 * t)

    # 测量噪声标准差
    meas_noise_std = 0.05
    # 陀螺仪测量值 = 真实角速度 + 零偏 + 噪声
    gyro_meas = omega_true + bias_true + np.random.normal(0, meas_noise_std, steps)

    # 创建滤波器
    kf = GyroKalmanFilter(dt, q_omega=1e-4, q_bias=1e-8, r_meas=meas_noise_std**2)

    # 存储估计结果
    omega_est = np.zeros(steps)
    bias_est = np.zeros(steps)

    for i in range(steps):
        omega_est[i], bias_est[i] = kf.estimate_omega(gyro_meas[i])

    # 绘图
    plt.figure(figsize=(12, 5))

    plt.subplot(2, 1, 1)
    plt.plot(t, omega_true, 'k-', label='True omega')
    plt.plot(t, gyro_meas, 'b-', alpha=0.5, label='Raw measurement')
    plt.plot(t, omega_est, 'r-', label='KF estimated omega')
    plt.ylabel('Angular velocity (rad/s)')
    plt.legend()
    plt.grid(True, ls=':')

    plt.subplot(2, 1, 2)
    plt.plot(t, bias_true, 'k-', label='True bias')
    plt.plot(t, bias_est, 'r-', label='KF estimated bias')
    plt.xlabel('Time (s)')
    plt.ylabel('Bias (rad/s)')
    plt.legend()
    plt.grid(True, ls=':')

    plt.tight_layout()
    plt.show()