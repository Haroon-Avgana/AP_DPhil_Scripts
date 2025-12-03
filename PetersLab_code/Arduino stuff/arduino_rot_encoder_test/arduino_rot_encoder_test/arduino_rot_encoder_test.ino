  
 /*
  Based on code from Tomaso Muzzu - UCL - 25 May 2017
  Script to read rotary encoder with quadrature encoding of position. Model Kubler 05.2400.1122.1024 (READ)

*/

#define encoder0PinA 2        // sensor A of rotary encoder
#define encoder0PinB 13        // sensor B of rotary encoder

// variables for rotary encoder
volatile signed int encoder0Pos = 0;    // variable for counting ticks of rotary encoder

void setup() {

  pinMode(encoder0PinA, INPUT);   // rotary encoder sensor A
  pinMode(encoder0PinB, INPUT);   // rotary encoder sensor B

  // interrupts for rotary encoder
  attachInterrupt(digitalPinToInterrupt(encoder0PinA), doEncoderA, RISING);

  Serial.begin (250000);
  Serial.setTimeout(5);

  delay(500);
}

void loop() {

}

// Interrupt on A low to high transition
void doEncoderA() {
    if (digitalRead(encoder0PinB)==HIGH) {
      encoder0Pos = 1;
    }
    else {
      encoder0Pos = - 1;
    }
    Serial.print(encoder0Pos);//
    Serial.print("\n");
}


// EOF
