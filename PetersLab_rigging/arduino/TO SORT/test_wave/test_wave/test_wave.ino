#include <Audio.h>


AudioSynthWaveform       waveform;
AudioOutputI2S           i2s1;           //xy=360,98
AudioOutputAnalogStereo  DAC1;          //xy=372,173
AudioConnection          patchCord(waveform, DAC1);



int current_waveform=0;

// extern const int16_t myWaveform[256];  // defined in myWaveform.ino

void setup() {
  Serial.begin(9600);
  // Confirgure both to use "myWaveform" for WAVEFORM_ARBITRARY
  // waveform.arbitraryWaveform(myWaveform, 172.0);

  // configure waveform for 40 Hz and maximum amplitude
  waveform.frequency(440);
  waveform.amplitude(1.0);

  current_waveform = WAVEFORM_SQUARE;
  waveform.begin(current_waveform);

}

// USEFUL: File - Examples - Audio

void loop() {
  // put your main code here, to run repeatedly:
  AudioNoInterrupts();
  waveform.begin(current_waveform);
  AudioInterrupts(); 
}
