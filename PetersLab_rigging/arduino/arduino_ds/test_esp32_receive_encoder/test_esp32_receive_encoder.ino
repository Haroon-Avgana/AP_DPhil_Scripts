
/*
  1) Read rotary encoder
    On an interrupt pin: for each wheel click, send +1 for CW and -1 for CCW
    Model: Kubler 05.2400.1122.1024

  2) Reward valve: on recieving positive integer from Bonsai, open REWARD valve for that length time

  3) Dummy valve: on recieving negative integer from Bonsai, open DUMMY valve for that length time

  (Based on code from Tomaso Muzzu - UCL - 25 May 2017)
  
*/

#define encoder0PinA 4    // sensor A of rotary encoder
#define encoder0PinB 15    // sensor B of rotary encoder


// variables for rotary encoder
volatile signed int encoder0Pos = 0;  // variable for counting ticks of rotary encoder


void setup() {

  pinMode(encoder0PinA, INPUT);  // rotary encoder sensor A
  pinMode(encoder0PinB, INPUT);  // rotary encoder sensor B


  // interrupts for rotary encoder
  attachInterrupt(digitalPinToInterrupt(encoder0PinA), doEncoderA, FALLING);

  Serial.begin(115200);
  Serial.setTimeout(5);


  delay(500);
}
int angle=0;
int count = 0;

void loop() {

  count++;
  if (count > 500) {
    count = 0;
    //Serial.printf("%f\n", DFOC_M0_Current());
Serial.println(angle); 
  }
  delay(1);
}

// Interrupt on A low to high transition
void doEncoderA() {
  if (digitalRead(encoder0PinB) == LOW) {
    encoder0Pos = 1;  
    angle= angle+encoder0Pos;

  } else {
    encoder0Pos = -1;
      angle= angle+encoder0Pos;

  }
  // Serial.println(encoder0Pos);  //
  // Serial.print("\n");
}


