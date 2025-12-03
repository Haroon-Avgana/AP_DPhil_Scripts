  
 /*
  Based on code from Tomaso Muzzu - UCL - 25 May 2017
  Script to read rotary encoder with quadrature encoding of position. 
  Model: Kubler 05.2400.1122.1024
*/

#define encoder0PinA 3        // sensor A of rotary encoder
#define encoder0PinB 7        // sensor B of rotary encoder
#define SValvePin 2           // digital pin controlling the solenoid valve


// variables for rotary encoder
volatile signed int encoder0Pos = 0;    // variable for counting ticks of rotary encoder

// variables for pintch valve of reward system
const byte numChars = 6;
char receivedChars[numChars];   // an array to store the received data
boolean newData = false;
int TimeON = 0;             // new for this version
int TempVar = 0;
boolean TimerFinished = false;
uint32_t StartTime = 0;      // variable to store temporary timestamps of previous iteration of the while loop


void setup() {

  pinMode(encoder0PinA, INPUT);   // rotary encoder sensor A
  pinMode(encoder0PinB, INPUT);   // rotary encoder sensor B

  pinMode(SValvePin, OUTPUT);     // solenoid valve

  // interrupts for rotary encoder
  attachInterrupt(digitalPinToInterrupt(encoder0PinA), doEncoderA, FALLING);

  Serial.begin (250000);
  Serial.setTimeout(5);

  delay(500);
}

void loop() {

  GetBonsaiInput();
  ActivatePV();

  delay(1);

}

// Interrupt on A low to high transition
void doEncoderA() {
    if (digitalRead(encoder0PinB)==LOW) {
      encoder0Pos = 1;
    }
    else {
      encoder0Pos = - 1;
    }
    Serial.print(encoder0Pos);//
    Serial.print("\n");
}


/////////////////////////////////////////////////////////////////////////////
// Read inputs from Bonsai. It is an integer specifing the amount of time in ms the pintch valve stays open
void GetBonsaiInput() { // part of code taken from http://forum.arduino.cc/index.php?topic=396450.0
  static byte ndx = 0;
  char endMarker = '\r';
  char rc;

  if (Serial.available() > 0) {
    rc = Serial.read();

    if (rc != endMarker) {
      receivedChars[ndx] = rc;
      ndx++;
      if (ndx >= numChars) {
        ndx = numChars - 1;
      }
    }
    else {
      receivedChars[ndx] = '\0'; // terminate the string
      ndx = 0;
      newData = true;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
// Timer for the output signal to the pintch valve to stay high
void ActivatePV() {
  if (newData == true) {
    TimeON = 0;             // zero previous time
    TimeON = atoi(receivedChars);   // convert array of chars to integers
    if (TempVar == 0) {
      StartTime = millis();
      TempVar = 1;
    }
    newData = false;
  }
  if (TempVar == 1) {
    // start checking the time and keep the valve open as long as you wish irrespective of what happens to Serial.read()
    if ((millis() - StartTime) <= (uint32_t)TimeON) {
      digitalWrite(SValvePin, HIGH); // open valve
      TimerFinished = false;
      TempVar = 1;
    }
    else {
      digitalWrite(SValvePin, LOW); // close valve
      TimerFinished = true;
      TempVar = 0;
    }
  }
}

// EOF
