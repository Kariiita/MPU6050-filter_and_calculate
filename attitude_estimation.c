#include "attitude_estimation.h"
#include "OLED.h"
#include "MPU6050.h"
#include "Delay.h"
#include <math.h>

// 滤波器参数（可调）
#define Kp 0.5f   // 比例增益，控制加速度计修正速度
#define Ki 0.05f  // 积分增益，消除陀螺仪零偏漂移

#define M_PI 3.1415926f

// 全局变量
static float integralFBx = 0.0f;
static float integralFBy = 0.0f;
static float integralFBz = 0.0f;

// 陀螺仪零偏（单位：度/秒）
float gyro_bias[3] = {0};

// 四元数（全局姿态表示）
static float q0 = 1.0f, q1 = 0.0f, q2 = 0.0f, q3 = 0.0f;

// 辅助函数：将四元数转换为欧拉角
static void quaternion_to_euler(float q0, float q1, float q2, float q3, EulerAngle *angle)
{
    // 欧拉角计算（Z-Y-X 顺序，对应偏航-俯仰-滚转）
    float roll  = atan2f(2.0f * (q0 * q1 + q2 * q3), 1.0f - 2.0f * (q1 * q1 + q2 * q2));
    float pitch = asinf(2.0f * (q0 * q2 - q3 * q1));
    float yaw   = atan2f(2.0f * (q0 * q3 + q1 * q2), 1.0f - 2.0f * (q2 * q2 + q3 * q3));

    // 转换为度
    angle->roll  = roll  * (180.0f / M_PI);
    angle->pitch = pitch * (180.0f / M_PI);
    angle->yaw   = yaw   * (180.0f / M_PI);
}

// 陀螺仪零偏校准：静止状态下采集 N 组数据取平均
void attitude_calibrate_gyro(void)
{
    // 假设你已经实现了读取陀螺仪原始数据的函数，这里直接调用
    // 实际使用时，请在主函数中调用此函数前采集数据
    // 示例：采集 500 组数据，求平均作为零偏
    const int samples = 1000;
    float sum_gx = 0, sum_gy = 0, sum_gz = 0;

    for (int i = 0; i < samples; i++) {
        int16_t gx, gy, gz;
        // 读取原始数据（单位：度/秒）
		int16_t ax, ay, az;
		MPU6050_GetData(&ax, &ay, &az, &gx, &gy, &gz);	// 硬件读取

        sum_gx += gx / 16.384f + 4;
        sum_gy += gy / 16.384f;
        sum_gz += gz / 16.384f;

        // 适当延时，保证数据采集间隔均匀
        Delay_ms(1);  // 假设有延时函数
    }

    gyro_bias[0] = sum_gx / samples;
    gyro_bias[1] = sum_gy / samples;
    gyro_bias[2] = sum_gz / samples;
}

void attitude_init(void)
{
    // 重置四元数为初始姿态（水平、朝向默认）
    q0 = 1.0f;
    q1 = 0.0f;
    q2 = 0.0f;
    q3 = 0.0f;

    // 重置积分项
    integralFBx = 0.0f;
    integralFBy = 0.0f;
    integralFBz = 0.0f;

    // 执行陀螺仪零偏校准
    attitude_calibrate_gyro();
	
//	OLED_ShowSignedNum(2, 8, (int)(gyro_bias[0] * 10), 5);
//	OLED_ShowSignedNum(3, 8, (int)(gyro_bias[1] * 10), 5);
//	OLED_ShowSignedNum(4, 8, (int)(gyro_bias[2] * 10), 5);
	
}

void attitude_update(float ax, float ay, float az,
                     float gx, float gy, float gz,
                     float dt, EulerAngle *angle)
{
    // 去除陀螺仪零偏
    gx -= gyro_bias[0];
    gy -= gyro_bias[1];
    gz -= gyro_bias[2];

    // 加速度计归一化（得到重力方向单位向量）
    float norm = sqrtf(ax * ax + ay * ay + az * az);
    if (norm < 0.001f) return;  // 防止除零
    ax /= norm;
    ay /= norm;
    az /= norm;

    // 从四元数估算当前重力方向向量（即机体坐标系下的重力分量）
    float vx = 2.0f * (q1 * q3 - q0 * q2);
    float vy = 2.0f * (q0 * q1 + q2 * q3);
    float vz = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;

    // 计算误差：实际加速度计测量值 与 估算重力向量 的叉积
    float ex = (ay * vz - az * vy);
    float ey = (az * vx - ax * vz);
    float ez = (ax * vy - ay * vx);

    // 对误差进行积分
    integralFBx += Ki * ex * dt;
    integralFBy += Ki * ey * dt;
    integralFBz += Ki * ez * dt;

    // 使用比例+积分项修正陀螺仪角速度
    gx += Kp * ex + integralFBx;
    gy += Kp * ey + integralFBy;
    gz += Kp * ez + integralFBz;

    // 将角速度（度/秒）转换为弧度/秒
    gx = gx * (M_PI / 180.0f);
    gy = gy * (M_PI / 180.0f);
    gz = gz * (M_PI / 180.0f);

    // 四元数微分方程（一阶龙格-库塔）
    float q0_dot = -0.5f * (q1 * gx + q2 * gy + q3 * gz);
    float q1_dot =  0.5f * (q0 * gx + q2 * gz - q3 * gy);
    float q2_dot =  0.5f * (q0 * gy - q1 * gz + q3 * gx);
    float q3_dot =  0.5f * (q0 * gz + q1 * gy - q2 * gx);

    // 更新四元数
    q0 += q0_dot * dt;
    q1 += q1_dot * dt;
    q2 += q2_dot * dt;
    q3 += q3_dot * dt;

    // 四元数归一化
    norm = sqrtf(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    if (norm < 0.001f) return;
    q0 /= norm;
    q1 /= norm;
    q2 /= norm;
    q3 /= norm;

    // 转换为欧拉角
    quaternion_to_euler(q0, q1, q2, q3, angle);
}
