const int cameraOutPin = 29;

const int cameraInPin = 30;    // pin for "all lines exposing"
const int blueOutPin = 24;       // pin for blue's Gate1
const int purpleOutPin = 25;       // pin for purple's Gate1


int flipflopState = 0;        
int lastPCOstate = 0;
int currentPCOstate = 0;

void setup() {

// setup 70Hz square wave
analogWriteFrequency(cameraOutPin, 70); // change out pin to 70Hz
pinMode(cameraOutPin,OUTPUT);

pinMode(cameraInPin, INPUT);
pinMode(blueOutPin, OUTPUT);
pinMode(purpleOutPin, OUTPUT);
Serial.begin(9600);
flipflopState = 0;

analogWriteResolution(12);

}

void loop() {

analogWrite(cameraOutPin,127);


// analogWrite(A22, 0xfff);
// delay(2000); 
// analogWrite(A22, 0x7ff);
// delay(2000); 

currentPCOstate = digitalRead(cameraInPin);

if ((currentPCOstate==LOW) & (lastPCOstate==HIGH)) { 
  flipflopState = (flipflopState+1) % 2; 

  if (flipflopState==0) {
    digitalWrite(blueOutPin, HIGH);
    digitalWrite(purpleOutPin, LOW);
  } else {
    digitalWrite(blueOutPin, LOW);
    digitalWrite(purpleOutPin, HIGH);
  }
}

lastPCOstate = currentPCOstate;

}