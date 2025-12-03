#include <DMAChannel.h>
#include "pdb.h"


DMAChannel dma(false);

static volatile uint16_t sinetable[] = {
   2047,    2147,    2248,    2348,    2447,    2545,    2642,    2737,
   2831,    2923,    3012,    3100,    3185,    3267,    3346,    3422,
   3495,    3564,    3630,    3692,    3750,    3804,    3853,    3898,
   3939,    3975,    4007,    4034,    4056,    4073,    4085,    4093,
   4095,    4093,    4085,    4073,    4056,    4034,    4007,    3975,
   3939,    3898,    3853,    3804,    3750,    3692,    3630,    3564,
   3495,    3422,    3346,    3267,    3185,    3100,    3012,    2923,
   2831,    2737,    2642,    2545,    2447,    2348,    2248,    2147,
   2047,    1948,    1847,    1747,    1648,    1550,    1453,    1358,
   1264,    1172,    1083,     995,     910,     828,     749,     673,
    600,     531,     465,     403,     345,     291,     242,     197,
    156,     120,      88,      61,      39,      22,      10,       2,
      0,       2,      10,      22,      39,      61,      88,     120,
    156,     197,     242,     291,     345,     403,     465,     531,
    600,     673,     749,     828,     910,     995,    1083,    1172,
   1264,    1358,    1453,    1550,    1648,    1747,    1847,    1948,
};

void setup() {
  dma.begin(true); // allocate the DMA channel first
  
  SIM_SCGC2 |= SIM_SCGC2_DAC0; // enable DAC clock
  DAC0_C0 = DAC_C0_DACEN | DAC_C0_DACRFS; // enable the DAC module, 3.3V reference
  // slowly ramp up to DC voltage, approx 1/4 second
  for (int16_t i=0; i<2048; i+=8) {
    *(int16_t *)&(DAC0_DAT0L) = i;
    delay(1);
  }
  
  // set the programmable delay block to trigger DMA requests
  SIM_SCGC6 |= SIM_SCGC6_PDB; // enable PDB clock
  PDB0_IDLY = 0; // interrupt delay register
  PDB0_MOD = PDB_PERIOD; // modulus register, sets period
  PDB0_SC = PDB_CONFIG | PDB_SC_LDOK; // load registers from buffers
  PDB0_SC = PDB_CONFIG | PDB_SC_SWTRIG; // reset and restart
  PDB0_CH0C1 = 0x0101; // channel n control register?
  
  dma.sourceBuffer(sinetable, sizeof(sinetable));
  dma.destination(*(volatile uint16_t *)&(DAC0_DAT0L));
  dma.triggerAtHardwareEvent(DMAMUX_SOURCE_PDB);
  dma.enable();

}


void loop() {
}
