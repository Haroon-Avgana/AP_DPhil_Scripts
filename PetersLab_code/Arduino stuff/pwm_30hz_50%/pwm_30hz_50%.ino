void setup() {

TCCR2B = TCCR2B & 0b11111000 | 0x07; // for PWM frequency of 30.51 Hz

pinMode(11,OUTPUT);

}

void loop() {

analogWrite(11,127);

}