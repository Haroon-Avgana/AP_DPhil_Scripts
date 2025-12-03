//例程2：边界模式
//可切换为:
//    无极-边界
//    棘轮边界

// DengFOC V0.2
// 灯哥开源，遵循GNU协议，转载请著名版权！
// GNU开源协议（GNU General Public License, GPL）是一种自由软件许可协议，保障用户能够自由地使用、研究、分享和修改软件。
// 该协议的主要特点是，要求任何修改或衍生的作品必须以相同的方式公开发布，即必须开源。此外，该协议也要求在使用或分发软件时，必须保留版权信息和许可协议。GNU开源协议是自由软件基金会（FSF）制定和维护的一种协议，常用于GNU计划的软件和其他自由软件中。
// 仅在DengFOC官方硬件上测试过，欢迎硬件购买/支持作者，淘宝搜索店铺：灯哥开源
// 你的支持将是接下来做视频和持续开源的经费，灯哥在这里先谢谢大家了

#include "DengFOC.h"

#define _2PI 6.28318530718f
/********************************************/
//0为无极-边界
//1为棘轮-边界
#define MODE 0

#define encoder1 4
#define encoder2 15

/********************************************/
int Sensor_DIR = 1;  // 传感器方向，若电机运动不正常，将此值取反
int Motor_PP = 7;    // 电机极对数
void Detents_Mode(void);
/********************************************/
//参数设置
float kp_pos1 = 1;
float kp_pos2 = 1;
float anlge_range = 0;     //边界范围90°
float anlge_trigger = 10;  //边界范围90°

float zero_angle = 0;  //零点位置
// zero_electric_angle
/********************************************/
void setup() {
  Serial.begin(115200);
  pinMode(12, OUTPUT);
  digitalWrite(12, HIGH);  // 使能，一定要放在校准电机前

  DFOC_Vbus(12.6);  // 设定驱动器供电电压
  DFOC_alignSensor(Motor_PP, Sensor_DIR);
  zero_angle = DFOC_M0_Angle();  //将当前位置设置为零点
  zero_angle = 3.14;             //将当前位置设置为零点

  pinMode(encoder1, OUTPUT);
  pinMode(encoder2, OUTPUT);
  digitalWrite(encoder1, HIGH);  // 初始化 encoder1 为high电平
  digitalWrite(encoder2, HIGH);  // 初始化 encoder2 为high电平
}
unsigned long lastTime = 0;   // 上次切换的时间
unsigned long interval = 5;  // 高低电平切换的时间（50ms）
bool isLow = true;            // 当前引脚B的状态


void loop() {
  runFOC();
  float angle = DFOC_M0_Angle();
  float angle_range2 = anlge_range / 2;
  float limit_range = angle_range2 / 360.0f * _2PI;
  float L1 = zero_angle + limit_range;
  float L2 = zero_angle - limit_range;
  float angle_trigger=anlge_trigger;
  float limit_range2 = angle_trigger / 360.0f * _2PI;
  float L3 = zero_angle + limit_range2;
  float L4 = zero_angle - limit_range2;
  //angle_range 范围内
  if (angle > L2 && angle < L1) {
    if (MODE == 0)
      DFOC_M0_setTorque(0);
    else if (MODE == 1)
      Detents_Mode();
  } else {  //angle_range 范围外
    if (angle < L2)
      DFOC_M0_setTorque(kp_pos1 * (L2 - angle));
    if (angle > L1)
      DFOC_M0_setTorque(kp_pos1 * (L1 - angle));
    if (angle > L3) {
      digitalWrite(encoder2, HIGH);
      unsigned long currentTime = millis();
      // 如果距离上次切换超过50毫秒，则切换引脚B的电平
      if (currentTime - lastTime >= interval) {
        if (isLow) {
          digitalWrite(encoder1, LOW);  // 设置引脚B为低电平
        } else {
          digitalWrite(encoder1, HIGH);  // 设置引脚B为高电平
        }
        lastTime = currentTime;  // 更新最后一次切换的时间
        isLow = !isLow;          // 切换状态
      }
    } 
    if (angle < L4) {
      digitalWrite(encoder2, LOW);
      unsigned long currentTime = millis();
      // 如果距离上次切换超过50毫秒，则切换引脚B的电平
      if (currentTime - lastTime >= interval) {
        if (isLow) {
          digitalWrite(encoder1, LOW);  // 设置引脚B为低电平
        } else {
          digitalWrite(encoder1, HIGH);  // 设置引脚B为高电平
        }
        lastTime = currentTime;  // 更新最后一次切换的时间
        isLow = !isLow;          // 切换状态
      }
    } 
  }
  // Serial.printf("%d,%d,%f\n", digitalRead(encoder1), digitalRead(encoder2),angle);

}



void Detents_Mode(void) {
  float attractor_distance = 15 * 3.1415926 / 180.0;
  float target = round(DFOC_M0_Angle() / attractor_distance) * attractor_distance;
  DFOC_M0_setTorque(kp_pos2 * (target - DFOC_M0_Angle()));
}
