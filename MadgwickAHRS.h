//=============================================================================================
// MadgwickAHRS.h
//=============================================================================================
//
// Implementation of Madgwick's IMU and AHRS algorithms in C.
// See: http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/
//
// This is a C translation of the original C++ implementation.
//=============================================================================================

#ifndef MADGWICKAHRS_H
#define MADGWICKAHRS_H

#include <math.h>

//--------------------------------------------------------------------------------------------
// 结构体定义：保存滤波器的所有状态变量
//--------------------------------------------------------------------------------------------
typedef struct {
    float beta;             // 算法增益，控制加速度计/磁力计修正强度
    float q0, q1, q2, q3;   // 四元数，表示传感器相对于参考系的姿态
    float invSampleFreq;    // 采样周期的倒数 (1 / 采样频率)
    float roll;             // 滚转角（弧度）
    float pitch;            // 俯仰角（弧度）
    float yaw;              // 偏航角（弧度）
    char anglesComputed;    // 标记欧拉角是否需要重新计算
} Madgwick;

//--------------------------------------------------------------------------------------------
// 函数声明
//--------------------------------------------------------------------------------------------

// 初始化滤波器，设置采样频率（Hz）
void Madgwick_init(Madgwick *filter, float sampleFrequency);

// 9轴AHRS更新（陀螺仪 + 加速度计 + 磁力计）
// 输入：gx, gy, gz 角速度（度/秒），ax, ay, az 加速度（任意单位，会被归一化），mx, my, mz 磁力计（任意单位）
void Madgwick_update(Madgwick *filter, float gx, float gy, float gz,
                     float ax, float ay, float az,
                     float mx, float my, float mz);

// 6轴IMU更新（陀螺仪 + 加速度计，无磁力计）
void Madgwick_updateIMU(Madgwick *filter, float gx, float gy, float gz,
                        float ax, float ay, float az);

// 获取欧拉角（度）
float Madgwick_getRoll(Madgwick *filter);
float Madgwick_getPitch(Madgwick *filter);
float Madgwick_getYaw(Madgwick *filter);

// 获取欧拉角（弧度）
float Madgwick_getRollRadians(Madgwick *filter);
float Madgwick_getPitchRadians(Madgwick *filter);
float Madgwick_getYawRadians(Madgwick *filter);

#endif // MADGWICKAHRS_H
