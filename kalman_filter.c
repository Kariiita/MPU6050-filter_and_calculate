/**
  * @file    kalman_filter.c
  * @brief   一维卡尔曼滤波器实现
  */
#include "kalman_filter.h"

/**
 * @brief   初始化卡尔曼滤波器
 * @param   kfp: 滤波器结构体指针
 * @param   Q:   过程噪声协方差，一般取较小值，如0.001~0.01[reference:0]
 * @param   R:   测量噪声协方差，一般取较大值，如0.1~1，取决于传感器噪声水平[reference:1]
 */
void KalmanFilter_Init(KalmanFilter_t *kfp, float Q, float R)
{
    kfp->Q = Q;
    kfp->R = R;
    kfp->P = 1.0f;      // 初始估计误差协方差
    kfp->K = 0.0f;
    kfp->X = 0.0f;      // 初始状态估计值
}

/**
 * @brief   卡尔曼滤波器更新函数
 * @param   kfp:   滤波器结构体指针
 * @param   input: 本次的观测值（如加速度计计算的角度）
 * @retval  滤波后的最优估计值
 */
float KalmanFilter_Update(KalmanFilter_t *kfp, float input)
{
    // 1. 预测：先验估计误差协方差
    kfp->P = kfp->P + kfp->Q;

    // 2. 更新：计算卡尔曼增益
    kfp->K = kfp->P / (kfp->P + kfp->R);

    // 3. 更新：状态估计（最优估计值）
    kfp->X = kfp->X + kfp->K * (input - kfp->X);

    // 4. 更新：后验估计误差协方差
    kfp->P = (1 - kfp->K) * kfp->P;

    return kfp->X;
}
