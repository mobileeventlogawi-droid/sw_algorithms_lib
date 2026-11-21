/***********************************************************************//**
 * @file		microphone_analysis.h
 * @brief		IGC Environmental Noise Level (ENL) and band levels
 *			from the microphone sample stream
 * @license 		This project is released under the GNU Public License GPL-3.0

    <Larus Flight Sensor Firmware>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 **************************************************************************/
#ifndef MICROPHONE_ANALYSIS_H_
#define MICROPHONE_ANALYSIS_H_

#include <stdint.h>
#include <math.h>

// Signal chain, record layout and the IGC requirements behind the numbers:
// sw_sensor documentation/microphone_ENL.md

#define MIC_SAMPLE_RATE_HZ		8000.0f
#define MIC_RECORD_SAMPLES		2000	// one record every 250 ms
#define MIC_SETTLING_SAMPLES		4000	// microphone start-up and filter transients

#define ENL_CENTER_FREQUENCY_HZ		150.0f	// IGC: peak sensitivity 100..200 Hz
#define ENL_FILTER_Q			1.2f	// per stage, two identical stages
#define ENL_WINDOW_RECORDS		4	// ENL = RMS over the last 1 s
#define ENL_BASELINE			10	// IGC: "010" = serviceable, no noise
#define ENL_MAXIMUM			999
#define ENL_UNSERVICEABLE		0	// IGC: "000" = system not working
// Uncalibrated default: ENL 900 at an in-band level of -26 dBFS (RMS
// 0.03544), i.e. 94 dB SPL with the ICS-43434's sensitivity.
#define ENL_DEFAULT_GAIN		(890.0f / 0.03544f)

#define MIC_OCTAVE_BANDS		7	// 31.5 Hz .. 2 kHz
#define MIC_OCTAVE_Q			1.414f	// one octave bandwidth

#define MIC_STATUS_NO_SIGNAL		0x0001	// constant samples: microphone missing or stuck
#define MIC_STATUS_LOST_BLOCKS_SHIFT	16	// DMA blocks not processed in time (set by driver)

//! LRSX record MICROPHONE_DATA, 11 words
typedef struct
{
  uint32_t ENL;				//!< 0..999, 1 s window
  uint32_t status;			//!< MIC_STATUS_xxx
  float ENL_band_dBFS;			//!< ENL filter output, 250 ms window
  float broadband_dBFS;			//!< DC removed, 250 ms window
  float octave_dBFS[MIC_OCTAVE_BANDS];	//!< 31.5 .. 2000 Hz, 250 ms window
} microphone_record_t;

//! RBJ band-pass (0 dB peak gain), direct form I: robust in float for
//! the low f0/fs ratios used here
class mic_biquad_t
{
public:
  void setup (float f0, float q, float fs)
  {
    float w0 = 2.0f * (float) M_PI * f0 / fs;
    float alpha = sinf (w0) / (2.0f * q);
    float a0 = 1.0f + alpha;
    b0 = alpha / a0;
    b2 = -alpha / a0;
    a1 = -2.0f * cosf (w0) / a0;
    a2 = (1.0f - alpha) / a0;
    x1 = x2 = y1 = y2 = 0.0f;
  }
  float filter (float x)
  {
    float y = b0 * x + b2 * x2 - a1 * y1 - a2 * y2; // b1 == 0
    x2 = x1;
    x1 = x;
    y2 = y1;
    y1 = y;
    return y;
  }
private:
  float b0, b2, a1, a2;
  float x1, x2, y1, y2;
};

class microphone_analysis_t
{
public:
  microphone_analysis_t (float ENL_center_frequency = ENL_CENTER_FREQUENCY_HZ,
			 float ENL_gain = ENL_DEFAULT_GAIN)
  : ENL_gain (ENL_gain)
  {
    ENL_stage[0].setup (ENL_center_frequency, ENL_FILTER_Q, MIC_SAMPLE_RATE_HZ);
    ENL_stage[1].setup (ENL_center_frequency, ENL_FILTER_Q, MIC_SAMPLE_RATE_HZ);
    for (unsigned i = 0; i < MIC_OCTAVE_BANDS; ++i)
      octave[i].setup (octave_center_frequency (i), MIC_OCTAVE_Q, MIC_SAMPLE_RATE_HZ);
    reset ();
  }

  static float octave_center_frequency (unsigned band)
  {
    return 31.25f * (float) (1 << band); // nominal 31.5, 63, 125 ... 2000 Hz
  }

  //! restart after a gap in the sample stream
  void reset (void)
  {
    settling = MIC_SETTLING_SAMPLES;
    dc_x1 = dc_y1 = 0.0f;
    clear_accumulators ();
    filled_records = 0;
    next_ENL_slot = 0;
    record.ENL = ENL_UNSERVICEABLE;
    record.status = 0;
  }

  //! @param raw left-justified 32-bit sample, full scale = 2^31
  //! @return true when a new record is available
  bool update (int32_t raw)
  {
    float x = (float) raw * (1.0f / 2147483648.0f);

    // DC blocker, corner ~10 Hz
    float hp = x - dc_x1 + 0.992f * dc_y1;
    dc_x1 = x;
    dc_y1 = hp;

    float enl = ENL_stage[1].filter (ENL_stage[0].filter (hp));
    float band[MIC_OCTAVE_BANDS];
    for (unsigned i = 0; i < MIC_OCTAVE_BANDS; ++i)
      band[i] = octave[i].filter (hp);

    if (settling > 0)
      {
	--settling;
	return false;
      }

    if (raw < raw_min)
      raw_min = raw;
    if (raw > raw_max)
      raw_max = raw;
    sum_broadband += hp * hp;
    sum_ENL += enl * enl;
    for (unsigned i = 0; i < MIC_OCTAVE_BANDS; ++i)
      sum_octave[i] += band[i] * band[i];

    if (++count < MIC_RECORD_SAMPLES)
      return false;

    finish_record ();
    return true;
  }

  const microphone_record_t & get_record (void) const	{ return record; }
  uint16_t get_ENL (void) const				{ return (uint16_t) record.ENL; }

private:
  static float dBFS (float sum)
  {
    // +3 dB: RMS relative to a full-scale sine, as usual for dBFS
    return 10.0f * log10f (2.0f * sum / (float) MIC_RECORD_SAMPLES + 1e-20f);
  }

  void clear_accumulators (void)
  {
    count = 0;
    raw_min = INT32_MAX;
    raw_max = INT32_MIN;
    sum_broadband = sum_ENL = 0.0f;
    for (unsigned i = 0; i < MIC_OCTAVE_BANDS; ++i)
      sum_octave[i] = 0.0f;
  }

  void finish_record (void)
  {
    ENL_mean_square[next_ENL_slot] = sum_ENL / (float) MIC_RECORD_SAMPLES;
    next_ENL_slot = (next_ENL_slot + 1) % ENL_WINDOW_RECORDS;
    if (filled_records < ENL_WINDOW_RECORDS)
      ++filled_records;

    float mean_square = 0.0f;
    for (unsigned i = 0; i < filled_records; ++i)
      mean_square += ENL_mean_square[i];
    mean_square /= (float) filled_records;

    // a working microphone always shows self-noise
    bool no_signal = (raw_min == raw_max);

    if (no_signal)
      record.ENL = ENL_UNSERVICEABLE;
    else
      {
	float value = (float) ENL_BASELINE + ENL_gain * sqrtf (mean_square);
	record.ENL = value >= (float) ENL_MAXIMUM ? ENL_MAXIMUM : (uint32_t) value;
      }
    record.status = no_signal ? MIC_STATUS_NO_SIGNAL : 0;
    record.ENL_band_dBFS = dBFS (sum_ENL);
    record.broadband_dBFS = dBFS (sum_broadband);
    for (unsigned i = 0; i < MIC_OCTAVE_BANDS; ++i)
      record.octave_dBFS[i] = dBFS (sum_octave[i]);

    clear_accumulators ();
  }

  mic_biquad_t ENL_stage[2];
  mic_biquad_t octave[MIC_OCTAVE_BANDS];
  float ENL_gain;
  unsigned settling;
  float dc_x1, dc_y1;
  unsigned count;
  int32_t raw_min, raw_max;
  float sum_broadband, sum_ENL, sum_octave[MIC_OCTAVE_BANDS];
  float ENL_mean_square[ENL_WINDOW_RECORDS];
  unsigned filled_records;
  unsigned next_ENL_slot;
  microphone_record_t record;
};

#endif /* MICROPHONE_ANALYSIS_H_ */
