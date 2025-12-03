#include "sq_wave.h"
#include "trp_wave.h"


volatile int i = 0;
volatile int j = 0;

const int cameraOutPin = 29;

void setup() {
  analogWriteResolution(12);
  
  pinMode(A22,OUTPUT);
  analogWrite(A22,0);

  pinMode(cameraOutPin,OUTPUT);
  analogWrite(cameraOutPin,0);

}


void loop() {

  analogWrite(cameraOutPin, waveformsTable_sq_wave[i]);  // write the selected waveform on DAC
  i++;
  if (i==512)
    i=0;

  analogWrite(A22, waveformsTable_trp_wave[j]);  // write the selected waveform on DAC
  j++;
  if (j==512)
    j=0;

  delayMicroseconds(27.9); // to slow it down so it's at 70Hz
}
