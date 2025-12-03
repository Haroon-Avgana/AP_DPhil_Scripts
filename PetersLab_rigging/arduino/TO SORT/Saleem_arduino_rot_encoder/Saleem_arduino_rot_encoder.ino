  
 /*
  Tomaso Muzzu - UCL - 25 May 2017
  Script to communicate with the following devices from Matlab:
  - rotary encoder with quadrature encoding of position. Model Kubler 05.2400.1122.1024 (READ)
  - pintch valve for water reward. Model NResearch 225P011-21 (WRITE)
  - lick detector based on IR beam breaking circuit. Model OP550 and IR26-21C-L110-TR8 (READ)
*/

#include <Event.h>
#include <Timer.h>

#define encoder0PinA 11        // sensor A of rotary encoder
#define encoder0PinB 12        // sensor B of rotary encoder
#define LickPin 5             // digital pin of lick detector
#define SValvePin 1           // digital pin controlling the solenoid valve

// variables for rotary encoder
volatile unsigned int encoder0Pos = 0;    // variable for counting ticks of rotary encoder
unsigned int tmp_Pos = 1;                 // variable for counting ticks of rotary encoder
boolean A_set;
boolean B_set;
// variables for lick counter
volatile unsigned int LickCount = 0;      // variable for counting licks
unsigned int tmp_LickCount = 0;           // temporary variable for counting licks
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
  pinMode(LickPin, INPUT);        // lick detectot
  pinMode(SValvePin, OUTPUT);     // solenoid valve

  // interrupts for rotary encoder
  attachInterrupt(digitalPinToInterrupt(encoder0PinA), doEncoderA, CHANGE);
  attachInterrupt(digitalPinToInterrupt(encoder0PinB), doEncoderB, CHANGE);
  // interrupt for lick detector
  attachInterrupt(digitalPinToInterrupt(LickPin), Lick_Counter, FALLING);

  Serial.begin (250000);
  Serial.setTimeout(5);

  delay(500);
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void loop() {
  //Check for change in position and send to serial buffer
  if (tmp_Pos != encoder0Pos || (tmp_LickCount != LickCount)) {
    Serial.print(encoder0Pos);//
    Serial.print("\t");
    Serial.print(LickCount);//
    Serial.print("\n");
    tmp_Pos = encoder0Pos;
    tmp_LickCount = LickCount;
  }
  else {
    Serial.print(tmp_Pos);//
    Serial.print("\t");
    Serial.print(tmp_LickCount);//
    Serial.print("\n");
    
  }

  GetMatlabInput();
  ActivatePV();

  delay(1);
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// Read inputs from Matlab. It is an integer specifing the amount of time in ms the pintch valve stays open
void GetMatlabInput() { // part of code taken from http://forum.arduino.cc/index.php?topic=396450.0
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
      digitalWrite(SValvePin, LOW); // open valve
      TimerFinished = true;
      TempVar = 0;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
// Interrupt on A changing state
void doEncoderA() {
  // Low to High transition?
  if (digitalRead(encoder0PinA) == HIGH) {
    A_set = true;
    if (!B_set) {
      encoder0Pos = encoder0Pos + 1;
    }
  }
  // High-to-low transition?
  if (digitalRead(encoder0PinA) == LOW) {
    A_set = false;
  }
}
// Interrupt on B changing state
void doEncoderB() {
  // Low-to-high transition?
  if (digitalRead(encoder0PinB) == HIGH) {
    B_set = true;
    if (!A_set) {
      encoder0Pos = encoder0Pos - 1;
    }
  }
  // High-to-low transition?
  if (digitalRead(encoder0PinB) == LOW) {
    B_set = false;
  }
}

/////////////////////////////////////////////////////////////////////////////
// Interrupt for when IR beam breaking circuit goes down
void Lick_Counter() {
  // High-to-low transition?
  if (digitalRead(LickPin) == LOW) {
    LickCount = LickCount + 1;
  }
}



// EOF
