//=============================================================================================
// MadgwickAHRS.c
//=============================================================================================
//
// Implementation of Madgwick's IMU and AHRS algorithms in C.
// See: http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/
//
// This is a C translation of the original C++ implementation.
//=============================================================================================

#include "MadgwickAHRS.h"
#include <math.h>

// 快速平方根倒数（Quake III 算法，用于加速归一化）
static float invSqrt(float x) {
    float halfx = 0.5f * x;
    float y = x;
    long i = *(long*)&y;          // 将浮点数位模式解释为整数
    i = 0x5f3759df - (i >> 1);    // 魔术数近似
    y = *(float*)&i;              // 转换回浮点
    y = y * (1.5f - (halfx * y * y)); // 牛顿迭代一次
    y = y * (1.5f - (halfx * y * y)); // 再迭代一次提高精度
    return y;
}

// 将四元数转换为欧拉角（弧度），并存储到结构体中
static void computeAngles(Madgwick *filter) {
    filter->roll = atan2f(filter->q0 * filter->q1 + filter->q2 * filter->q3,
                          0.5f - filter->q1 * filter->q1 - filter->q2 * filter->q2);
    filter->pitch = asinf(-2.0f * (filter->q1 * filter->q3 - filter->q0 * filter->q2));
    filter->yaw = atan2f(filter->q1 * filter->q2 + filter->q0 * filter->q3,
                         0.5f - filter->q2 * filter->q2 - filter->q3 * filter->q3);
    filter->anglesComputed = 1;
}

//--------------------------------------------------------------------------------------------
// 初始化滤波器
//--------------------------------------------------------------------------------------------
void Madgwick_init(Madgwick *filter, float sampleFrequency) {
    filter->beta = 0.1f;                     // 默认beta = 0.1
    filter->q0 = 1.0f;                       // 初始四元数为单位四元数
    filter->q1 = 0.0f;
    filter->q2 = 0.0f;
    filter->q3 = 0.0f;
    filter->invSampleFreq = 1.0f / sampleFrequency;  // 采样周期倒数
    filter->anglesComputed = 0;
}

//--------------------------------------------------------------------------------------------
// 9轴AHRS算法更新（陀螺仪+加速度计+磁力计）
//--------------------------------------------------------------------------------------------
void Madgwick_update(Madgwick *filter,
                     float gx, float gy, float gz,
                     float ax, float ay, float az,
                     float mx, float my, float mz) {
    float recipNorm;
    float s0, s1, s2, s3;
    float qDot1, qDot2, qDot3, qDot4;
    float hx, hy;
    float _2q0mx, _2q0my, _2q0mz, _2q1mx, _2bx, _2bz, _4bx, _4bz;
    float _2q0, _2q1, _2q2, _2q3, _2q0q2, _2q2q3;
    float q0q0, q0q1, q0q2, q0q3, q1q1, q1q2, q1q3, q2q2, q2q3, q3q3;

    // 如果磁力计测量值全为0（无效），则退化为6轴IMU算法，避免NaN
    if ((mx == 0.0f) && (my == 0.0f) && (mz == 0.0f)) {
        Madgwick_updateIMU(filter, gx, gy, gz, ax, ay, az);
        return;
    }

    // 将陀螺仪单位从度/秒转换为弧度/秒（算法内部要求弧度制）
    gx *= 0.0174533f;
    gy *= 0.0174533f;
    gz *= 0.0174533f;

    // 由陀螺仪计算四元数的变化率（预测部分）
    qDot1 = 0.5f * (-filter->q1 * gx - filter->q2 * gy - filter->q3 * gz);
    qDot2 = 0.5f * ( filter->q0 * gx + filter->q2 * gz - filter->q3 * gy);
    qDot3 = 0.5f * ( filter->q0 * gy - filter->q1 * gz + filter->q3 * gx);
    qDot4 = 0.5f * ( filter->q0 * gz + filter->q1 * gy - filter->q2 * gx);

    // 加速度计有效时进行修正（避免除零）
    if (!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

        // 归一化加速度计测量值
        recipNorm = invSqrt(ax * ax + ay * ay + az * az);
        ax *= recipNorm;
        ay *= recipNorm;
        az *= recipNorm;

        // 归一化磁力计测量值
        recipNorm = invSqrt(mx * mx + my * my + mz * mz);
        mx *= recipNorm;
        my *= recipNorm;
        mz *= recipNorm;

        // 辅助变量，避免重复计算
        _2q0mx = 2.0f * filter->q0 * mx;
        _2q0my = 2.0f * filter->q0 * my;
        _2q0mz = 2.0f * filter->q0 * mz;
        _2q1mx = 2.0f * filter->q1 * mx;
        _2q0 = 2.0f * filter->q0;
        _2q1 = 2.0f * filter->q1;
        _2q2 = 2.0f * filter->q2;
        _2q3 = 2.0f * filter->q3;
        _2q0q2 = 2.0f * filter->q0 * filter->q2;
        _2q2q3 = 2.0f * filter->q2 * filter->q3;
        q0q0 = filter->q0 * filter->q0;
        q0q1 = filter->q0 * filter->q1;
        q0q2 = filter->q0 * filter->q2;
        q0q3 = filter->q0 * filter->q3;
        q1q1 = filter->q1 * filter->q1;
        q1q2 = filter->q1 * filter->q2;
        q1q3 = filter->q1 * filter->q3;
        q2q2 = filter->q2 * filter->q2;
        q2q3 = filter->q2 * filter->q3;
        q3q3 = filter->q3 * filter->q3;

        // 计算地球磁场的参考方向（将磁力计测量值从机体坐标系转换到地球坐标系）
        hx = mx * q0q0 - _2q0my * filter->q3 + _2q0mz * filter->q2 + mx * q1q1 + _2q1 * my * filter->q2 + _2q1 * mz * filter->q3 - mx * q2q2 - mx * q3q3;
        hy = _2q0mx * filter->q3 + my * q0q0 - _2q0mz * filter->q1 + _2q1mx * filter->q2 - my * q1q1 + my * q2q2 + _2q2 * mz * filter->q3 - my * q3q3;
        _2bx = sqrtf(hx * hx + hy * hy);
        _2bz = -_2q0mx * filter->q2 + _2q0my * filter->q1 + mz * q0q0 + _2q1mx * filter->q3 - mz * q1q1 + _2q2 * my * filter->q3 - mz * q2q2 + mz * q3q3;
        _4bx = 2.0f * _2bx;
        _4bz = 2.0f * _2bz;

        // 梯度下降算法的误差函数对四元数的偏导数（步进方向）
        s0 = -_2q2 * (2.0f * q1q3 - _2q0q2 - ax) + _2q1 * (2.0f * q0q1 + _2q2q3 - ay) - _2bz * filter->q2 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * filter->q3 + _2bz * filter->q1) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * filter->q2 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
        s1 = _2q3 * (2.0f * q1q3 - _2q0q2 - ax) + _2q0 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * filter->q1 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + _2bz * filter->q3 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * filter->q2 + _2bz * filter->q0) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * filter->q3 - _4bz * filter->q1) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
        s2 = -_2q0 * (2.0f * q1q3 - _2q0q2 - ax) + _2q3 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * filter->q2 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + (-_4bx * filter->q2 - _2bz * filter->q0) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * filter->q1 + _2bz * filter->q3) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * filter->q0 - _4bz * filter->q2) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
        s3 = _2q1 * (2.0f * q1q3 - _2q0q2 - ax) + _2q2 * (2.0f * q0q1 + _2q2q3 - ay) + (-_4bx * filter->q3 + _2bz * filter->q1) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * filter->q0 + _2bz * filter->q2) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * filter->q1 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
        // 归一化步长
        recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3);
        s0 *= recipNorm;
        s1 *= recipNorm;
        s2 *= recipNorm;
        s3 *= recipNorm;

        // 将修正项应用到四元数变化率（反馈）
        qDot1 -= filter->beta * s0;
        qDot2 -= filter->beta * s1;
        qDot3 -= filter->beta * s2;
        qDot4 -= filter->beta * s3;
    }

    // 对四元数变化率进行积分，得到新的四元数
    filter->q0 += qDot1 * filter->invSampleFreq;
    filter->q1 += qDot2 * filter->invSampleFreq;
    filter->q2 += qDot3 * filter->invSampleFreq;
    filter->q3 += qDot4 * filter->invSampleFreq;

    // 归一化四元数，确保模长为1
    recipNorm = invSqrt(filter->q0 * filter->q0 + filter->q1 * filter->q1 +
                        filter->q2 * filter->q2 + filter->q3 * filter->q3);
    filter->q0 *= recipNorm;
    filter->q1 *= recipNorm;
    filter->q2 *= recipNorm;
    filter->q3 *= recipNorm;
    filter->anglesComputed = 0;   // 欧拉角需要重新计算
}

//--------------------------------------------------------------------------------------------
// 6轴IMU算法更新（陀螺仪+加速度计，无磁力计）
//--------------------------------------------------------------------------------------------
void Madgwick_updateIMU(Madgwick *filter,
                        float gx, float gy, float gz,
                        float ax, float ay, float az) {
    float recipNorm;
    float s0, s1, s2, s3;
    float qDot1, qDot2, qDot3, qDot4;
    float _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2, _8q1, _8q2;
    float q0q0, q1q1, q2q2, q3q3;

    // 陀螺仪转换为弧度/秒
    gx *= 0.0174533f;
    gy *= 0.0174533f;
    gz *= 0.0174533f;

    // 陀螺仪积分预测四元数变化率
    qDot1 = 0.5f * (-filter->q1 * gx - filter->q2 * gy - filter->q3 * gz);
    qDot2 = 0.5f * ( filter->q0 * gx + filter->q2 * gz - filter->q3 * gy);
    qDot3 = 0.5f * ( filter->q0 * gy - filter->q1 * gz + filter->q3 * gx);
    qDot4 = 0.5f * ( filter->q0 * gz + filter->q1 * gy - filter->q2 * gx);

    // 加速度计有效时进行修正
    if (!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

        // 归一化加速度计
        recipNorm = invSqrt(ax * ax + ay * ay + az * az);
        ax *= recipNorm;
        ay *= recipNorm;
        az *= recipNorm;

        // 辅助变量
        _2q0 = 2.0f * filter->q0;
        _2q1 = 2.0f * filter->q1;
        _2q2 = 2.0f * filter->q2;
        _2q3 = 2.0f * filter->q3;
        _4q0 = 4.0f * filter->q0;
        _4q1 = 4.0f * filter->q1;
        _4q2 = 4.0f * filter->q2;
        _8q1 = 8.0f * filter->q1;
        _8q2 = 8.0f * filter->q2;
        q0q0 = filter->q0 * filter->q0;
        q1q1 = filter->q1 * filter->q1;
        q2q2 = filter->q2 * filter->q2;
        q3q3 = filter->q3 * filter->q3;

        // 计算误差函数梯度（仅利用重力参考）
        s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
        s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * filter->q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
        s2 = 4.0f * q0q0 * filter->q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
        s3 = 4.0f * q1q1 * filter->q3 - _2q1 * ax + 4.0f * q2q2 * filter->q3 - _2q2 * ay;
        // 归一化步长
        recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3);
        s0 *= recipNorm;
        s1 *= recipNorm;
        s2 *= recipNorm;
        s3 *= recipNorm;

        // 反馈修正
        qDot1 -= filter->beta * s0;
        qDot2 -= filter->beta * s1;
        qDot3 -= filter->beta * s2;
        qDot4 -= filter->beta * s3;
    }

    // 积分更新四元数
    filter->q0 += qDot1 * filter->invSampleFreq;
    filter->q1 += qDot2 * filter->invSampleFreq;
    filter->q2 += qDot3 * filter->invSampleFreq;
    filter->q3 += qDot4 * filter->invSampleFreq;

    // 归一化
    recipNorm = invSqrt(filter->q0 * filter->q0 + filter->q1 * filter->q1 +
                        filter->q2 * filter->q2 + filter->q3 * filter->q3);
    filter->q0 *= recipNorm;
    filter->q1 *= recipNorm;
    filter->q2 *= recipNorm;
    filter->q3 *= recipNorm;
    filter->anglesComputed = 0;
}

//--------------------------------------------------------------------------------------------
// 获取欧拉角（度）
//--------------------------------------------------------------------------------------------
float Madgwick_getRoll(Madgwick *filter) {
    if (!filter->anglesComputed) computeAngles(filter);
    return filter->roll * 57.29578f;   // 弧度转度
}

float Madgwick_getPitch(Madgwick *filter) {
    if (!filter->anglesComputed) computeAngles(filter);
    return filter->pitch * 57.29578f;
}

float Madgwick_getYaw(Madgwick *filter) {
    if (!filter->anglesComputed) computeAngles(filter);
    return filter->yaw * 57.29578f;
}

//--------------------------------------------------------------------------------------------
// 获取欧拉角（弧度）
//--------------------------------------------------------------------------------------------
float Madgwick_getRollRadians(Madgwick *filter) {
    if (!filter->anglesComputed) computeAngles(filter);
    return filter->roll;
}

float Madgwick_getPitchRadians(Madgwick *filter) {
    if (!filter->anglesComputed) computeAngles(filter);
    return filter->pitch;
}

float Madgwick_getYawRadians(Madgwick *filter) {
    if (!filter->anglesComputed) computeAngles(filter);
    return filter->yaw;
}
