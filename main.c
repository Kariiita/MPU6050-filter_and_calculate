#include "stm32f10x.h"                  // Device header
#include "Delay.h"
#include "OLED.h"
#include "MPU6050.h"
#include "attitude_estimation.h"
#include <math.h>

#include "MadgwickAHRS.h"
#include "TIM3.h"

#include "kalman_filter.h"

#define PI 3.1415926
#define FILTER_BUFFERSIZE 8


extern volatile uint8_t tim3_flag;

int16_t AX, AY, AZ, GX, GY, GZ;			//定义用于存放各个数据的变量
int16_t filterdata[6] = {0};

EulerAngle angle = {0};


int16_t buffer[6][FILTER_BUFFERSIZE] = {0};
int16_t sum[6] = {0};
int i = 0;
void moving_average_filter(int16_t *datas, int16_t *filterdata)
{
	for (int j=0; j<6; j++)
	{
		sum[j] -= buffer[j][i];
		buffer[j][i] = datas[j];
		sum[j] += datas[j];
		filterdata[j] = sum[j] / FILTER_BUFFERSIZE;
	}
	i = (i + 1) % FILTER_BUFFERSIZE;
}

// gzy
void compute_euler_angle(int16_t ax, int16_t ay, int16_t az, int16_t wz, float dt, EulerAngle *angle)
{
	double pitch = atan2((double) -ax, sqrt(ay*ay+az*az))*180/PI;
	double roll  = atan2(ax, sqrt(ay*ay+az*az))*180/PI;
	angle->pitch = pitch;
	angle->roll  = roll;
	angle->yaw   += wz*dt;
}




Madgwick fltr;

KalmanFilter_t kalmanX, kalmanY, kalmanZ;  // 定义卡尔曼滤波器


float dt = 0.01f;

int main(void)
{
	/*模块初始化*/
	OLED_Init();		//OLED初始化
	
	OLED_ShowString(1, 1, "1");		//显示静态字符串
	OLED_ShowString(1, 8, "2");		//显示静态字符串
	
	MPU6050_Init();		//MPU6050初始化
	attitude_init();	//互补滤波初始化
	
	// 初始化Madgwick滤波器
	Madgwick_init(&fltr, 100.0f);   // 100Hz 采样率
	TIM3_Init(10);
	
	// 初始化卡尔曼滤波器
	KalmanFilter_Init(&kalmanX, 0.01f, 0.1f);
    KalmanFilter_Init(&kalmanY, 0.01f, 0.1f);
	KalmanFilter_Init(&kalmanZ, 0.01f, 0.1f);
	
	
	
	while (1)
	{
		if (TIM3_GetFlag())
		{
			TIM3_ClearFlag();
			
			MPU6050_GetData(&AX, &AY, &AZ, &GX, &GY, &GZ);		//获取MPU6050的数据
			
			
			AX = AX/2048.0f*10;
			AY = AY/2048.0f*10;
			AZ = AZ/2048.0f*10;	// 加速度 单位m/s^2
			GX /= 16.384f;
			GY /= 16.384f;
			GZ /= 16.384f;	// 角速度 单位°/s
			
			GX += 4;	// 硬件误差修正
			
			int16_t datas[6] = {AX, AY, AZ, GX, GY, GZ};
			moving_average_filter(datas, filterdata);
			
			float kalman_gx = KalmanFilter_Update(&kalmanX, GX);
			float kalman_gy = KalmanFilter_Update(&kalmanY, GY);
			float kalman_gz = KalmanFilter_Update(&kalmanZ, GZ);
			
			
//			OLED_ShowSignedNum(2, 1, filterdata[3], 5);					//OLED显示数据
//			OLED_ShowSignedNum(3, 1, filterdata[4], 5);
//			OLED_ShowSignedNum(4, 1, filterdata[5], 5);
//			OLED_ShowSignedNum(2, 8, kalman_gx, 5);
//			OLED_ShowSignedNum(3, 8, kalman_gy, 5);
//			OLED_ShowSignedNum(4, 8, kalman_gz, 5);
			
			attitude_update(filterdata[0], filterdata[1], filterdata[2], kalman_gx, kalman_gy, kalman_gz, dt, &angle);
			//compute_euler_angle(AX, AY, AZ, GZ, dt, &angle);
			Madgwick_updateIMU(&fltr, GX, GY, GZ, AX, AY, AZ);
			
			
			OLED_ShowSignedNum(2, 1, angle.roll, 5);
			OLED_ShowSignedNum(3, 1, angle.pitch, 5);
			OLED_ShowSignedNum(4, 1, angle.yaw, 5);
			
			
			float m_roll = Madgwick_getRoll(&fltr);
			float m_pitch = Madgwick_getPitch(&fltr);
			float m_yaw = Madgwick_getYaw(&fltr);
			
			OLED_ShowSignedNum(2, 8, m_roll, 5);
			OLED_ShowSignedNum(3, 8, m_pitch, 5);
			OLED_ShowSignedNum(4, 8, m_yaw, 5);
		
		
		}
	}
}

