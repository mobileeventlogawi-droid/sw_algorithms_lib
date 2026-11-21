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

// third-octave bands, base-2 centers 1000 Hz * 2^(k/3): 24.8 Hz .. 2520 Hz
// (nominal 25 .. 2500 Hz), the highest band edge stays below fs/2
#define MIC_BANDS			21
#define MIC_BAND_1KHZ			16	// index of the 1 kHz band
#define MIC_BAND_SECTIONS		3	// 6th order Butterworth band-pass

#define MIC_STATUS_NO_SIGNAL		0x0001	// constant samples: microphone missing or stuck
#define MIC_STATUS_LOST_BLOCKS_SHIFT	16	// DMA blocks not processed in time (set by driver)

//! LRSX record MICROPHONE_DATA, 27 words. Levels in dBFS: 0 dBFS = sine
//! with full-scale peak.
typedef struct
{
  uint32_t record_count;		//!< +1 per record, shows gaps
  uint32_t ENL;				//!< 0..999, 1 s window
  uint32_t status;			//!< MIC_STATUS_xxx
  float ENL_band_dBFS;			//!< ENL filter output, 250 ms RMS
  float broadband_dBFS;			//!< DC removed, 250 ms RMS
  float peak_dBFS;			//!< largest raw sample in 250 ms, 0 = clipping
  float band_dBFS[MIC_BANDS];		//!< third-octave bands, 250 ms RMS
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

//! Third-octave band-pass: 3rd order Butterworth low-pass prototype,
//! low-pass to band-pass transform, bilinear transform (band edges
//! prewarped). Every section has the numerator 1 - z^-2; the overall gain
//! is applied to the mean square only, saving a multiplication per sample.
//! Cascaded direct form I sections share their state between neighbours.
class mic_third_octave_t
{
public:
  void setup (float center, float fs)
  {
    struct cplx
    {
      float re, im;
      cplx operator+ (cplx b) const { return { re + b.re, im + b.im }; }
      cplx operator- (cplx b) const { return { re - b.re, im - b.im }; }
      cplx operator* (cplx b) const { return { re * b.re - im * b.im, re * b.im + im * b.re }; }
      cplx operator/ (cplx b) const
      {
	float d = b.re * b.re + b.im * b.im;
	return { (re * b.re + im * b.im) / d, (im * b.re - re * b.im) / d };
      }
      cplx root (void) const
      {
	float r = sqrtf (hypotf (re, im)), phi = atan2f (im, re) / 2.0f;
	return { r * cosf (phi), r * sinf (phi) };
      }
    };

    // edges of the third-octave band, prewarped for the bilinear transform
    float w1 = tanf ((float) M_PI * center / powf (2.0f, 1.0f / 6.0f) / fs);
    float w2 = tanf ((float) M_PI * center * powf (2.0f, 1.0f / 6.0f) / fs);
    float w0 = sqrtf (w1 * w2);
    float bw = w2 - w1;

    // upper half-plane band-pass poles: two from the prototype pole at
    // 120 degrees, one from the real prototype pole
    cplx p_complex = { cosf (2.0f * (float) M_PI / 3.0f) * bw / 2.0f, sinf (2.0f * (float) M_PI / 3.0f) * bw / 2.0f };
    cplx p_real = { -bw / 2.0f, 0.0f };
    cplx w0_squared = { w0 * w0, 0.0f };
    cplx d = (p_complex * p_complex - w0_squared).root ();
    cplx poles[MIC_BAND_SECTIONS] = { p_complex + d, p_complex - d, p_real + (p_real * p_real - w0_squared).root () };

    // section response at the center frequency to normalize the gain
    float omega0 = 2.0f * atanf (w0);
    cplx z_inv = { cosf (omega0), -sinf (omega0) };
    cplx z_inv2 = z_inv * z_inv;
    float gain = 1.0f;
    cplx one = { 1.0f, 0.0f };

    for (unsigned k = 0; k < MIC_BAND_SECTIONS; ++k)
      {
	cplx s = poles[k];
	if (s.im < 0.0f)
	  s.im = -s.im; // the conjugate belongs to the same section
	cplx z = (one + s) / (one - s);
	float c1 = -2.0f * z.re;
	float c2 = z.re * z.re + z.im * z.im;
	a1[k] = c1;
	a2[k] = c2;
	cplx denominator = one + cplx { c1, 0.0f } * z_inv + cplx { c2, 0.0f } * z_inv2;
	cplx response = (one - z_inv2) / denominator;
	gain /= hypotf (response.re, response.im);
      }
    gain_squared = gain * gain;

    for (unsigned k = 0; k <= MIC_BAND_SECTIONS; ++k)
      state[k][0] = state[k][1] = 0.0f;
  }

  //! @return unscaled output, multiply its mean square by gain_squared
  float filter (float x)
  {
    float v = x;
    for (unsigned k = 0; k < MIC_BAND_SECTIONS; ++k)
      {
	float y = v - state[k][1] - a1[k] * state[k + 1][0] - a2[k] * state[k + 1][1];
	state[k][1] = state[k][0];
	state[k][0] = v;
	v = y;
      }
    state[MIC_BAND_SECTIONS][1] = state[MIC_BAND_SECTIONS][0];
    state[MIC_BAND_SECTIONS][0] = v;
    return v;
  }

  float gain_squared;

private:
  float a1[MIC_BAND_SECTIONS], a2[MIC_BAND_SECTIONS];
  float state[MIC_BAND_SECTIONS + 1][2]; // node k: last two values
};

class microphone_analysis_t
{
public:
  microphone_analysis_t (float ENL_center_frequency = ENL_CENTER_FREQUENCY_HZ,
			 float ENL_gain = ENL_DEFAULT_GAIN)
  : ENL_gain (ENL_gain),
    record_count (0)
  {
    ENL_stage[0].setup (ENL_center_frequency, ENL_FILTER_Q, MIC_SAMPLE_RATE_HZ);
    ENL_stage[1].setup (ENL_center_frequency, ENL_FILTER_Q, MIC_SAMPLE_RATE_HZ);
    for (unsigned i = 0; i < MIC_BANDS; ++i)
      band[i].setup (band_center_frequency (i), MIC_SAMPLE_RATE_HZ);
    reset ();
  }

  static float band_center_frequency (unsigned band)
  {
    return 1000.0f * powf (2.0f, ((float) band - (float) MIC_BAND_1KHZ) / 3.0f);
  }

  //! restart after a gap in the sample stream, record_count keeps counting
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

    // DC blocker, corner ~2.5 Hz: keeps the 25 Hz band unattenuated
    float hp = x - dc_x1 + 0.998f * dc_y1;
    dc_x1 = x;
    dc_y1 = hp;

    float enl = ENL_stage[1].filter (ENL_stage[0].filter (hp));
    float band_output[MIC_BANDS];
    for (unsigned i = 0; i < MIC_BANDS; ++i)
      band_output[i] = band[i].filter (hp);

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
    for (unsigned i = 0; i < MIC_BANDS; ++i)
      sum_band[i] += band_output[i] * band_output[i];

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
    for (unsigned i = 0; i < MIC_BANDS; ++i)
      sum_band[i] = 0.0f;
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
    record.record_count = ++record_count;
    record.status = no_signal ? MIC_STATUS_NO_SIGNAL : 0;
    record.ENL_band_dBFS = dBFS (sum_ENL);
    record.broadband_dBFS = dBFS (sum_broadband);

    // |INT32_MIN| does not fit into int32_t
    float peak = fmaxf (-(float) raw_min, (float) raw_max) * (1.0f / 2147483648.0f);
    record.peak_dBFS = 20.0f * log10f (peak + 1e-10f);

    for (unsigned i = 0; i < MIC_BANDS; ++i)
      record.band_dBFS[i] = dBFS (sum_band[i] * band[i].gain_squared);

    clear_accumulators ();
  }

  mic_biquad_t ENL_stage[2];
  mic_third_octave_t band[MIC_BANDS];
  float ENL_gain;
  uint32_t record_count;
  unsigned settling;
  float dc_x1, dc_y1;
  unsigned count;
  int32_t raw_min, raw_max;
  float sum_broadband, sum_ENL, sum_band[MIC_BANDS];
  float ENL_mean_square[ENL_WINDOW_RECORDS];
  unsigned filled_records;
  unsigned next_ENL_slot;
  microphone_record_t record;
};

#endif /* MICROPHONE_ANALYSIS_H_ */
