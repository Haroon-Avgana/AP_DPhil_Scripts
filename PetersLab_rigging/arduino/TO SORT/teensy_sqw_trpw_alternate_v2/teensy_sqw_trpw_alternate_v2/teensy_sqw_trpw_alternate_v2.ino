#include "sq_wave.h"
#include "trp_wave.h"


volatile int i = 0;
volatile int j = 0;

const int goPin = 0; // pin for experiment

const int cameraOutPin = 29;

const int cameraInPin = 30;       // pin for "all lines exposing"
const int blueOutPin = 24;        // pin for blue's Gate1
const int violetOutPin = 25;      // pin for purple's Gate1
const uint8_t rampLightPin = A21;  // pin for ramping light

int flipflopState = 0;        
int lastPCOstate = 0;
int currentPCOstate = 0;

void setup() {

  pinMode(goPin, INPUT);
  pinMode(cameraInPin, INPUT);
  pinMode(blueOutPin, OUTPUT);
  pinMode(violetOutPin, OUTPUT);
  Serial.begin(9600);
  flipflopState = 0;

  analogWriteResolution(12);

  pinMode(rampLightPin,OUTPUT);
  analogWrite(rampLightPin,0);

  pinMode(cameraOutPin,OUTPUT);
  analogWrite(cameraOutPin,0);
}

void loop() {

  if (digitalRead(goPin)==HIGH) {
    digitalWrite(cameraOutPin, waveformsTable_sq_wave[i]);  // write the selected waveform on DAC
    i++;
    if (i==512)
      i=0;

    analogWrite(rampLightPin, waveformsTable_trp_wave[j]);  // write the selected waveform on DAC
    j++;
    if (j==512)
      j=0;

    delayMicroseconds(27.9); // to slow it down so it's at 70Hz


  } else {
    // If GO pin is low, write all pins to low
    digitalWrite(blueOutPin, LOW);
    digitalWrite(violetOutPin, LOW);
    digitalWrite(cameraOutPin, LOW);
    analogWrite(rampLightPin,0);
    flipflopState = 0;
  }

// Flip lights between blue and violet whenever cameraInPin goes low
  currentPCOstate = digitalRead(cameraInPin);

  if ((currentPCOstate==LOW) & (lastPCOstate==HIGH)) { 
    flipflopState = (flipflopState+1) % 2; 
    if (flipflopState==0) {
      digitalWrite(blueOutPin, HIGH);
      digitalWrite(violetOutPin, LOW);
    } else {
      digitalWrite(blueOutPin, LOW);
      digitalWrite(violetOutPin, HIGH);
    }
  }

  lastPCOstate = currentPCOstate;

}