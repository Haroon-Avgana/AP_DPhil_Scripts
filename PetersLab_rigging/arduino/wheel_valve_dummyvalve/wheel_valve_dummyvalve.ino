
 /*
  1) Read rotary encoder
    On an interrupt pin: for each wheel click, send +1 for CW and -1 for CCW
    Model: Kubler 05.2400.1122.1024

  2) Reward valve: on recieving positive integer from Bonsai, open REWARD valve for that length time

  3) Dummy valve: on recieving negative integer from Bonsai, open DUMMY valve for that length time

  (Based on code from Tomaso Muzzu - UCL - 25 May 2017)
  
*/

#define encoder0PinA 3        // sensor A of rotary encoder
#define encoder0PinB 7        // sensor B of rotary encoder
#define RewardValvePin 2           // digital pin for reward valve
#define DummyValvePin 12           // digital pin for dummy valve
#define SpeakerMotor 4           // digital pin for reward valve


// variables for rotary encoder
volatile signed int encoder0Pos = 0;    // variable for counting ticks of rotary encoder

// variables for pinch valve of reward system
const byte numChars = 6;
char receivedChars[numChars];   // an array to store the received data
int BonsaiValveTime = 0; 
int ValveTimeOn = 0; 

// (reward valve)         
boolean rewardNewData = false;
int ValveClockRunning = 0;
boolean TimerFinished = false;
uint32_t StartTime = 0;      // variable to store temporary timestamps of previous iteration of the while loop

// (dummy valve)            
boolean dummyNewData = false;
int DummyValveClockRunning = 0;
boolean DummyTimerFinished = false;
uint32_t DummyStartTime = 0;      // variable to store temporary timestamps of previous iteration of the while loop


void setup() {

  pinMode(encoder0PinA, INPUT);   // rotary encoder sensor A
  pinMode(encoder0PinB, INPUT);   // rotary encoder sensor B

  pinMode(RewardValvePin, OUTPUT);     // reward valve
  pinMode(DummyValvePin, OUTPUT);     // dummy valve
  pinMode(SpeakerMotor, OUTPUT);     // dummy valve

  // interrupts for rotary encoder
  attachInterrupt(digitalPinToInterrupt(encoder0PinA), doEncoderA, FALLING);

  Serial.begin (250000);
  Serial.setTimeout(5);

  delay(500);
}

void loop() {

  GetBonsaiInput();
  OpenRewardValve();
  OpenDummyValve();

  delay(1);

}

// Interrupt on A low to high transition
void doEncoderA() {
    if (digitalRead(encoder0PinB)==LOW) {
      encoder0Pos = 1;
    }
    else {
      encoder0Pos = -1;
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

      BonsaiValveTime = atoi(receivedChars);

      if (BonsaiValveTime > 0) {
        rewardNewData = true;
        ValveTimeOn = BonsaiValveTime;
      }
      else if (BonsaiValveTime < 0) {
        dummyNewData = true;
        ValveTimeOn = -BonsaiValveTime;
      }
     
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
// REWARD VALVE TIMER
void OpenRewardValve() {
  if (rewardNewData == true) {
    if (ValveClockRunning == 0) {
      StartTime = millis();
      ValveClockRunning = 1;
    }
    rewardNewData = false;
  }
  if (ValveClockRunning == 1) {
    // start checking the time and keep the valve open as long as you wish irrespective of what happens to Serial.read()
    if ((millis() - StartTime) <= (uint32_t)ValveTimeOn) {
      digitalWrite(RewardValvePin, HIGH); // open valve
      TimerFinished = false;
      ValveClockRunning = 1;
    }
    else {
      digitalWrite(RewardValvePin, LOW); // close valve
      TimerFinished = true;
      ValveClockRunning = 0;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
// DUMMY VALVE TIMER
void OpenDummyValve() {
  if (dummyNewData == true) {
    if (DummyValveClockRunning == 0) {
      DummyStartTime = millis();
      DummyValveClockRunning = 1;
    }
    dummyNewData = false;
  }
  if (DummyValveClockRunning == 1) {
    // start checking the time and keep the valve open as long as you wish irrespective of what happens to Serial.read()
    if ((millis() - DummyStartTime) <= (uint32_t)ValveTimeOn) {
      digitalWrite(DummyValvePin, HIGH); // open valve
      DummyTimerFinished = false;
      DummyValveClockRunning = 1;
    }
    else {
      digitalWrite(DummyValvePin, LOW); // close valve
      DummyTimerFinished = true;
      DummyValveClockRunning = 0;
    }
  }
}

// EOF
