#ifndef ATTITUDE_ESTIMATION_H
#define ATTITUDE_ESTIMATION_H

#include <stdint.h>

// 姿态角结构体（单位：度）
typedef struct {
    float roll;   // 滚转角
    float pitch;  // 俯仰角
    float yaw;    // 偏航角
} EulerAngle;



// 初始化姿态解算（包含陀螺仪零偏校准）
void attitude_init(void);

// 姿态解算更新函数（需要以固定周期调用）
// ax, ay, az: 加速度计原始值（单位：g）
// gx, gy, gz: 陀螺仪原始值（单位：度/秒）
// dt: 两次调用之间的时间间隔（秒）
// angle: 输出欧拉角
void attitude_update(float ax, float ay, float az,
                     float gx, float gy, float gz,
                     float dt, EulerAngle *angle);

#endif
					 