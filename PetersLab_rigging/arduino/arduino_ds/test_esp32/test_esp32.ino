char receivedChars[32];    // 用于存储接收到的字符数据
int BonsaiValveTime = 0;   // 存储转换后的整数值
bool newData = false;      // 标志数据是否新接收

void GetBonsaiInput() {
  static byte ndx = 0;
  char endMarker = '\r';   // 结束标记是回车符
  char rc;

  // 检查串口是否有可用数据
  while (Serial2.available() > 0) {
    rc = Serial2.read();  // 读取一个字符

    // 如果没有遇到结束标记，则继续存储数据
    if (rc != endMarker) {
      if (ndx < sizeof(receivedChars) - 1) {
        receivedChars[ndx] = rc;  // 存储字符
        ndx++;                    // 更新索引
      } else {
        // 如果超出缓冲区，直接丢弃
        ndx = 0;
        Serial.println("Error: Buffer Overflow");
      }
    } else {
      receivedChars[ndx] = '\0';  // 添加字符串终止符
      ndx = 0;                   // 重置索引

      // 转换字符串为整数
      BonsaiValveTime = atoi(receivedChars);

      // 设置新数据标志
      newData = true;
    }
  }
}

void setup() {
  Serial.begin(115200);  // 调试串口
  Serial2.begin(9600);   // 与 Arduino 的串口通信
}

void loop() {
  GetBonsaiInput();  // 调用函数处理串口输入

  if (newData) {
    // 打印接收到的数据
    Serial.print("Received BonsaiValveTime: ");
    Serial.println(BonsaiValveTime);

    // 根据接收到的值执行不同的操作
    if (BonsaiValveTime == 1) {
      Serial.println("Action: Received 1, performing action...");
    } else if (BonsaiValveTime == -1) {
      Serial.println("Action: Received -1, performing another action...");
    } else {
      Serial.println("Action: Received an unsupported value.");
    }

    // 重置标志
    newData = false;
  }

  delay(10);  // 防止 WDT 触发重启
}
