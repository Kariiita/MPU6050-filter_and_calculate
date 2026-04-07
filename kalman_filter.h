/**
  * @file    kalman_filter.h
  * @brief   一维卡尔曼滤波器头文件
  */
#ifndef _KALMAN_FILTER_H_
#define _KALMAN_FILTER_H_

#include "stm32f10x.h"

/* 卡尔曼滤波器结构体 */
typedef struct
{
    float Q;        // 过程噪声协方差
    float R;        // 测量噪声协方差
    float P;        // 估计误差协方差
    float K;        // 卡尔曼增益
    float X;        // 滤波后的值
} KalmanFilter_t;

/* 函数声明 */
void KalmanFilter_Init(KalmanFilter_t *kfp, float Q, float R);
float KalmanFilter_Update(KalmanFilter_t *kfp, float input);

#endif /* _KALMAN_FILTER_H_ */
