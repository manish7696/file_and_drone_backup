/* -*- c++ -*- */
/*
 * Copyright 2020 Free Software Foundation, Inc.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/fft/fft.h>
#include <gnuradio/fft/window.h>
#include <gnuradio/io_signature.h>

#include "fft_burst_tagger_impl.h"

#include <volk/volk.h>

#include <inttypes.h>
#include <stdio.h>
#include <chrono>
#include <math.h>

namespace gr {
namespace iridium {

fft_burst_tagger::sptr fft_burst_tagger::make(double center_frequency,
                                              int fft_size,
                                              int sample_rate,
                                              int burst_pre_len,
                                              int burst_post_len,
                                              int burst_width,
                                              int max_bursts,
                                              int max_burst_len,
                                              float threshold,
                                              int history_size,
                                              bool offline,
                                              bool debug)
{
    return gnuradio::get_initial_sptr(new fft_burst_tagger_impl(center_frequency,
                                                                fft_size,
                                                                sample_rate,
                                                                burst_pre_len,
                                                                burst_post_len,
                                                                burst_width,
                                                                max_bursts,
                                                                max_burst_len,
                                                                threshold,
                                                                history_size,
                                                                offline,
                                                                debug));
}


/*
 * The private constructor
 */
fft_burst_tagger_impl::fft_burst_tagger_impl(double center_frequency,
                                             int fft_size,
                                             int sample_rate,
                                             int burst_pre_len,
                                             int burst_post_len,
                                             int burst_width,
                                             int max_bursts,
                                             int max_burst_len,
                                             float threshold,
                                             int history_size,
                                             bool offline,
                                             bool debug)
    : gr::sync_block("fft_burst_tagger",
                     gr::io_signature::make(1, 1, sizeof(gr_complex)),
                     gr::io_signature::make(1, 1, sizeof(gr_complex))),
      d_center_frequency(center_frequency),
      d_sample_rate(sample_rate),
      d_fft_size(fft_size),
      d_burst_pre_len(burst_pre_len),
      d_max_burst_len(max_burst_len),
      d_burst_id(0),
      d_sample_count(0),
      d_n_tagged_bursts(0),
      d_squelch_count(0),
      d_fft(NULL),
      d_history_size(history_size),
      d_peaks(std::vector<peak>()),
      d_bursts(std::vector<burst>()),
      d_history_primed(false),
      d_history_index(0),
      d_burst_post_len(burst_post_len),
      d_debug(debug),
      d_burst_debug_file(NULL),
      d_last_rx_time_offset(0),
      d_last_rx_time_timestamp(0),
      d_offline(offline)
{
    const int nthreads = 1; //........
    d_fft = new fft::fft_complex_fwd(d_fft_size, nthreads); //........

    set_output_multiple(d_fft_size); //.........

    s_cycle=0;//..........
    s_file1 = fopen("/home/ktk/Documents/abc_fft16384_noise_290mhz_band_hackrf_20msps.bin", "wb");//............

    // We need to keep d_burst_pre_len samples
    // in the buffer to be able to tag a burst at it's start.
    // Set the history to this + 1, so we always have
    // this amount of samples available at the start of
    // our input buffer.
    set_history(d_burst_pre_len + 1);//..........

    d_window_f = (float*)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());//.......************
    std::vector<float> window =
        fft::window::build(fft::window::WIN_BLACKMAN, d_fft_size, 0);
    memcpy(d_window_f, &window[0], sizeof(float) * d_fft_size); //..............*******************

    // To get better SNR and noise floor estimates we apply the scaling factor and
    // Equivalent noise bandwidth of our window.
    // See
    // https://www.sjsu.edu/people/burford.furman/docs/me120/FFT_tutorial_NI.pdf
    // page 15 and
    // https://www.ap.com/blog/fft-spectrum-and-spectral-densities-same-data-different-scaling/

    // Scaling factor*****************************************
    for (int j = 0; j < d_fft_size; j++) {
        d_window_f[j] /= 0.42;
    } //******************************************

    // Equivalent noise bandwidth
    d_window_enbw = 1.72;

    L=d_fft_size;
    frequency=95e6;

    h15=(float*) malloc (sizeof(float)*20); //*********************
      for (int ii = 0; ii <= 20-1; ++ii) {
        if(ii<=9) h15[ii]=1/20.0;
        else h15[ii]=-1/20.0;
      }
      y15=(float*) malloc (sizeof(float)*L);
      y15_abs=(float*) malloc (sizeof(float)*L);
      del_meanVec_log = (float*) malloc (sizeof(float)*L);
      temp11 = (float*) malloc (sizeof(float)*L);

      pPos = (float*) malloc (sizeof(float)*L);
      pPosMax=(float*) malloc (sizeof(float)*L);
      pNev=(float*) malloc (sizeof(float)*L);
      pNevMin = (float*) malloc (sizeof(float)*L);
      posVec = (float*) malloc (sizeof(float)*L);
      posLoc = (int*) malloc (sizeof(int)*L);
      nevVec = (float*) malloc (sizeof(float)*L);
      nevLoc = (int*) malloc (sizeof(int)*L);
      nevPosIdx = (int*) malloc (sizeof(int)*L);
      for (int ii = 0; ii <= L-1; ++ii) {
      	pPos[ii]=pPosMax[ii]=pNev[ii]=pNevMin[ii]=posVec[ii]=posLoc[ii]=nevVec[ii]=nevLoc[ii]=nevPosIdx[ii]=0;}

      chLtedge = (int*) malloc (sizeof(int)*L);
      chRtedge = (int*) malloc (sizeof(int)*L);
      chCen = (int*) malloc (sizeof(int)*L);
      chBw = (float*) malloc (sizeof(float)*L);
      chMag = (float*) malloc (sizeof(float)*L);
      for (int ii = 0; ii <= L-1; ++ii) {chLtedge[ii]=chRtedge[ii]=chCen[ii]=chBw[ii]=chMag[ii]=0;}
      tempBuf = (gr_complex*) malloc (sizeof(gr_complex)*L);


    d_baseline_history_f = (float*)volk_malloc(
        sizeof(float) * d_fft_size * d_history_size, volk_get_alignment());
    d_baseline_sum_f =
        (float*)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());
    d_magnitude_shifted_f =
        (float*)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());
    d_relative_magnitude_f =
        (float*)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());
    d_burst_mask_f =
        (float*)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());
    d_ones_f = (float*)volk_malloc(sizeof(float) * d_fft_size, volk_get_alignment());



    memset(d_baseline_history_f, 0, sizeof(float) * d_fft_size * d_history_size);
    memset(d_baseline_sum_f, 0, sizeof(float) * d_fft_size);
    memset(d_magnitude_shifted_f, 0, sizeof(float) * d_fft_size); //..
    memset(d_relative_magnitude_f, 0, sizeof(float) * d_fft_size); //......

    for (int i = 0; i < d_fft_size; i++) {
        d_ones_f[i] = 1.0;
        d_burst_mask_f[i] = 1.0;
    }
    meanVec = (float*) malloc (sizeof(float)*L);     //....
	for (int ii = 0; ii <= L-1; ++ii) meanVec[ii]=0;     //............
    meanVec_log = (float*) malloc (sizeof(float)*L);     //....
    meanVec_log_s = (float*) malloc (sizeof(float)*L);



    // Divide by the ENBW as the calculation of d_relative_magnitude_f does
    // not take the ENBW of the FFT into account.
    d_threshold = pow(10, threshold / 10) / d_history_size / d_window_enbw;

    d_peaks.reserve(d_fft_size);

    if (max_bursts) {
        d_max_bursts = max_bursts;
    } else {
        // Consider the signal to be invalid if more than 80%
        // of all channels are in use.
        d_max_bursts = (sample_rate / burst_width) * 0.8;
    }

    // Area to ignore around an already found signal in FFT bins
    // Internal representation is in FFT bins
    d_burst_width = burst_width / (sample_rate / fft_size);

}

/*
 * Our virtual destructor.
 */
fft_burst_tagger_impl::~fft_burst_tagger_impl() //...........
{
    delete d_fft;
    volk_free(d_window_f);
    volk_free(d_baseline_history_f);
    volk_free(d_baseline_sum_f);
    volk_free(d_relative_magnitude_f);
    if (d_burst_debug_file) {
        fclose(d_burst_debug_file);
    }
      free (y15);
      free (y15_abs);
      free (meanVec_log);
      free (meanVec_log_s);
      free (del_meanVec_log);
      free (temp11);
      free (tempBuf);
      free (pPos);
      free (pPosMax);
      free (pNev);
      free (pNevMin);
      free (posVec);
      free (posLoc);
      free (nevVec);
      free (nevLoc);
      free (nevPosIdx);
      free (chLtedge);
      free (chRtedge);
      free (chCen);
      free (chBw);
      free (chMag);

    fclose(s_file1);
}

// *********************************************************************************************



bool fft_burst_tagger_impl::update_filters_pre(void)
{
    if (!d_history_primed) {
// return back in initial 512 cycles....
        return false;
    }
// divide for 513 onwards cycles.....

    volk_32f_x2_divide_32f(
        d_relative_magnitude_f, d_magnitude_shifted_f, d_baseline_sum_f, d_fft_size);
    return true;
}

#define HIST(i) (d_baseline_history_f + i * d_fft_size)

void fft_burst_tagger_impl::update_filters_post(bool force)
{
    // We only update the average if there is no burst going on at the moment
    if (d_bursts.size() == 0 || force) {
        volk_32f_x2_subtract_32f(
            d_baseline_sum_f, d_baseline_sum_f, HIST(d_history_index), d_fft_size);
        volk_32f_x2_add_32f(
            d_baseline_sum_f, d_baseline_sum_f, d_magnitude_shifted_f, d_fft_size);
        memcpy(HIST(d_history_index), d_magnitude_shifted_f, sizeof(float) * d_fft_size);

        d_history_index++;

        if (d_history_index == d_history_size) {
            d_history_primed = true;
            d_history_index = 0;
        }
    }
}

void fft_burst_tagger_impl::update_bursts(void)
{
    auto b = std::begin(d_bursts);
    while (b != std::end(d_bursts)) {
        if (d_relative_magnitude_f[b->center_bin - 1] > d_threshold ||
            d_relative_magnitude_f[b->center_bin] > d_threshold ||
            d_relative_magnitude_f[b->center_bin + 1] > d_threshold) {
            b->last_active = d_index;
        }
        ++b;
    }
}

void fft_burst_tagger_impl::delete_gone_bursts(void)
{
    bool update_noise_floor = false;
    auto b = std::begin(d_bursts);

    while (b != std::end(d_bursts)) {
        bool long_burst = false;

        if (d_max_burst_len) {
            if (b->last_active - b->start > d_max_burst_len) {
                update_noise_floor = true;
                long_burst = true;
            }
        }

        if ((b->last_active + d_burst_post_len) <= d_index || long_burst) {
            // printf("Deleting gone burst %" PRIu64 " (start=%" PRIu64 ", d_index=%"
            // PRIu64 ")\n", b->id, b->start, d_index);
            b->stop = d_index;
            d_gone_bursts.push_back(*b);
            b = d_bursts.erase(b);
        } else {
            ++b;
        }
    }

    if (update_noise_floor) {
        update_filters_post(true);
    }
}

void fft_burst_tagger_impl::create_new_bursts(void)
{
    for (peak p : d_peaks) {
        if (d_burst_mask_f[p.bin]) {
            burst b;
            b.id = d_burst_id;
            b.center_bin = p.bin;

            // Allow downstream blocks to split this burst
            // and assign sub ids
            d_burst_id += 10;

            // Normalize the relative magnitude
            // relative_magnitude relates to the uncorrected (regarding ENBW) noise floor.
            // We apply the ENBW here to have a more accurate SNR estimate
            b.magnitude =
                10 * log10(p.relative_magnitude * d_history_size * d_window_enbw);
            // The burst might have started one FFT earlier
            b.start = d_index - d_burst_pre_len;
            b.last_active = b.start;
            // Keep noise level around (dbFS/Hz)
            // Need to divide by the fft size twice as d_baseline_sum_f is a square of the
            // FFT's magnitude Apply the ENBW again to get an accurate estimate
            b.noise = 10 * log10(d_baseline_sum_f[b.center_bin] / d_history_size /
                                 (d_fft_size * d_fft_size) / d_window_enbw /
                                 (d_sample_rate / d_fft_size));

            d_bursts.push_back(b);
            d_new_bursts.push_back(b);
            mask_burst(b);
        }
    }
    if (d_max_bursts > 0 && d_bursts.size() > d_max_bursts) {
        fprintf(
            stderr, "Detector in burst squelch at %f\n", d_index / float(d_sample_rate));
        d_new_bursts.clear();
        for (burst b : d_bursts) {
            if (b.start != d_index - d_burst_pre_len) {
                b.stop = d_index;
                d_gone_bursts.push_back(b);
            }
        }
        d_bursts.clear();
        update_burst_mask();

        d_squelch_count += 3;
        // If we get too many squelches our noise floor might be significantly off
        // and the normal adjustments are not good enough.
        if (d_squelch_count >= 10) {
            fprintf(stderr, "Resetting noise estimate\n");
            d_history_index = 0;
            d_history_primed = 0;
            memset(d_baseline_history_f, 0, sizeof(float) * d_fft_size * d_history_size);
            memset(d_baseline_sum_f, 0, sizeof(float) * d_fft_size);
            d_squelch_count = 0;
        }
    } else {
        if (d_squelch_count) {
            d_squelch_count--;
        }
    }
}

void fft_burst_tagger_impl::mask_burst(burst& b)
{
    int clear_start = std::max(b.center_bin - d_burst_width / 2, 0);
    int clear_stop = std::min(b.center_bin + d_burst_width / 2, d_fft_size - 1);
    memset(
        d_burst_mask_f + clear_start, 0, (clear_stop - clear_start + 1) * sizeof(float));
}

void fft_burst_tagger_impl::update_burst_mask(void)
{
    memcpy(d_burst_mask_f, d_ones_f, sizeof(float) * d_fft_size);
    for (burst b : d_bursts) {
        mask_burst(b);
    }
}

void fft_burst_tagger_impl::remove_peaks_around_bursts(void)
{
    volk_32f_x2_multiply_32f(
        d_relative_magnitude_f, d_relative_magnitude_f, d_burst_mask_f, d_fft_size);
}

void fft_burst_tagger_impl::extract_peaks(void)
{

    d_peaks.clear();

    for (int bin = d_burst_width / 2; bin < (d_fft_size - d_burst_width / 2); bin++) {
        if (d_relative_magnitude_f[bin] > d_threshold) {
            peak p;
            p.bin = bin;
            p.relative_magnitude = d_relative_magnitude_f[bin];
            d_peaks.push_back(p);
            // printf("ts %" PRIu64 " bin %d val %f\n", d_index, p.bin,
            // p.relative_magnitude);
        }
    }

    struct {
        bool operator()(peak a, peak b)
        {
            return a.relative_magnitude > b.relative_magnitude;
        }
    } mag_gt;

    std::sort(d_peaks.begin(), d_peaks.end(), mag_gt);
}

void fft_burst_tagger_impl::save_peaks_to_debug_file(char* filename)
{
    FILE* file = fopen(filename, "a");
    for (peak p : d_peaks) {
        fprintf(file, "%" PRIu64 ",%d,x\n", d_index, p.bin);
        // float f_rel = (p.bin - d_fft_size / 2) / float(d_fft_size);
        // fprintf(file, "%f,%f,x\n", d_index/4e6, f_rel * 4e6 + 1624800000);
    }
    fclose(file);
}

void fft_burst_tagger_impl::tag_new_bursts(void)
{
    for (burst b : d_new_bursts) {
        // printf("new burst %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", nitems_read(0),
        // b.start, nitems_read(0) - b.start);
        pmt::pmt_t key = pmt::string_to_symbol("new_burst");
        float relative_frequency = (b.center_bin - d_fft_size / 2) / float(d_fft_size);


        const uint64_t offset = b.start - d_last_rx_time_offset;

        // Convert the offset in samples into a timestamp in nanoseconds.
        // offset * 1000000000 can become larger than 2**64, so we need to scale it
        // together with the sample rate. This means the sample rate must be a multiple of
        // 100000
        const uint64_t timestamp =
            d_last_rx_time_timestamp +
            (offset * (1000000000ULL / 100000)) / (d_sample_rate / 100000);

        pmt::pmt_t value = pmt::make_dict();
        value = pmt::dict_add(value, pmt::mp("id"), pmt::from_uint64(b.id));
        value = pmt::dict_add(
            value, pmt::mp("relative_frequency"), pmt::from_float(relative_frequency));
//        value = pmt::dict_add(
//            value, pmt::mp("center_frequency"), pmt::from_double(d_center_frequency));
        value = pmt::dict_add(value, pmt::mp("magnitude"), pmt::from_float(b.magnitude));
        value =
//            pmt::dict_add(value, pmt::mp("sample_rate"), pmt::from_float(d_sample_rate));
        value = pmt::dict_add(value, pmt::mp("timestamp"), pmt::from_uint64(timestamp));
        value = pmt::dict_add(value, pmt::mp("noise"), pmt::from_float(b.noise));

        // Our output is lagging by d_burst_pre_len samples.
        // Compensate by moving the tag into the past
        // printf("Tagging new burst %" PRIu64 " on sample %" PRIu64 " (nitems_read(0)=%"
        // PRIu64 ")\n", b.id, b.start + d_burst_pre_len, nitems_read(0));
        add_item_tag(0, b.start + d_burst_pre_len, key, value);
    }
    d_new_bursts.clear();
}

void fft_burst_tagger_impl::tag_gone_bursts(int noutput_items)
{
    auto b = std::begin(d_gone_bursts);

    while (b != std::end(d_gone_bursts)) {
        uint64_t output_index = b->stop + d_burst_pre_len;

        if (nitems_read(0) <= output_index &&
            output_index < nitems_read(0) + noutput_items) {
            pmt::pmt_t key = pmt::string_to_symbol("gone_burst");
            pmt::pmt_t value = pmt::make_dict();
            value = pmt::dict_add(value, pmt::mp("id"), pmt::from_uint64(b->id));
            // printf("Tagging gone burst %" PRIu64 " on sample %" PRIu64 "
            // (nitems_read(0)=%" PRIu64 ", noutput_items=%u)\n", b->id, output_index,
            // nitems_read(0), noutput_items);
            add_item_tag(0, output_index, key, value);
            d_n_tagged_bursts++;

            b = d_gone_bursts.erase(b);
        } else {
            ++b;
        }
    }
}

uint64_t fft_burst_tagger_impl::get_n_tagged_bursts() { return d_n_tagged_bursts; }

uint64_t fft_burst_tagger_impl::get_sample_count() { return d_sample_count; }

int fft_burst_tagger_impl::work(int noutput_items,                   //*******************************************************************************
                                gr_vector_const_void_star& input_items,
                                gr_vector_void_star& output_items)
{
    d_sample_count += noutput_items;
    int chanCount=-1;

    // We keep d_burst_pre_len additional samples in front of our data
    const gr_complex* in = (const gr_complex*)input_items[0] + 0;
    gr_complex* out = (gr_complex*)output_items[0];

    assert(noutput_items % d_fft_size == 0);

    std::vector<tag_t> rx_freq_tags;
    get_tags_in_window(rx_freq_tags, 0, 0, noutput_items, pmt::mp("rx_freq"));
    if (!rx_freq_tags.empty()) {            // std::cout<< cCount<<" ";
    cCount=0;
    for (int ii = 0; ii <= L-1; ++ii) meanVec[ii]=0;
    std::sort(rx_freq_tags.begin(), rx_freq_tags.end(), tag_t::offset_compare);
    const auto& rx_freq_tag = rx_freq_tags.back();
    const pmt::pmt_t& value_freq = rx_freq_tag.value;
    frequency = pmt::to_double(value_freq); // std::cout<<frequency<<"\n";
    }
    else cCount++;

	//  std::cout<<frequency<<"  "<<noutput_items<<"\n";
	
    for (int i = 0; i < noutput_items; i += d_fft_size) {
        d_index = nitems_read(0) + i;

        volk_32fc_32f_multiply_32fc(d_fft->get_inbuf(), &in[i], d_window_f, d_fft_size);
        d_fft->execute();
        volk_32fc_magnitude_squared_32f(&d_magnitude_shifted_f[0],
                                        &d_fft->get_outbuf()[d_fft_size / 2],
                                        d_fft_size / 2);
        volk_32fc_magnitude_squared_32f(&d_magnitude_shifted_f[d_fft_size / 2],
                                        &d_fft->get_outbuf()[0],
                                        d_fft_size / 2);
       s_cycle=d_index/d_fft_size;

//      if(s_cycle>0 && s_cycle<5001)
// fwrite(d_magnitude_shifted_f,sizeof(float),d_fft_size,s_file1);

      for (int ii = 0; ii <= L-1; ++ii) {
        meanVec[ii] = (999/1000.0)*meanVec[ii] + (1/1000.0)*d_magnitude_shifted_f[ii];
        meanVec_log[ii] = log10(meanVec[ii]);
        meanVec_log_s[ii]=0;
        temp11[ii]=meanVec_log[ii];
        if(ii>=20) del_meanVec_log[ii]=meanVec_log[ii-20];
        else del_meanVec_log[ii]=1;
    }
    // if(cCount>0 && (cCount==49)){
      if(cCount>0 && (cCount==2)) {
        // std::cout<<frequency<<"\n";
        for(int ii=0; ii<=L-1; ++ii) {
        out[ii].real(1); out[ii].imag(0);
        if(ii==L-1) out[ii].real(frequency);
        tempBuf[ii].real(meanVec_log[ii]); tempBuf[ii].imag(0);
        }
      }
        // else if(cCount>0 && (cCount==50)) {
      else if(cCount>0 && (cCount==3)) {

      float smoothwidth=30.0;
      int w=round(smoothwidth);
      int halfw=round(w/2.0);

      float sumPoints=0,tempSum=0,localSum=0,minValue=100;
      int k;

      for (int ii = 0; ii <= w-1; ++ii) sumPoints+=temp11[ii];
      for (k = 0; k <= L-1-w; ++k) {
           meanVec_log_s[k+halfw-1]=sumPoints;
           sumPoints = sumPoints - temp11[k] + temp11[k+w];
          }
      for (int kk = L-1-w+1; kk <= L-1; ++kk) tempSum+=temp11[kk];
      meanVec_log_s[k+halfw-1]=tempSum;
      for (int ii = 0; ii <= L-1; ++ii){
          meanVec_log_s[ii]/=w;
          if(meanVec_log_s[ii]<minValue) minValue=meanVec_log_s[ii]; }

      for (int ii=0; ii<= L-1; ++ii) { meanVec_log_s[ii]+= std::abs(minValue);
        localSum+=meanVec_log_s[ii];
         }
      localSum/=float(L);


      float y15_abs_mean=0, fSum=0;
      for (int ii = 0; ii <= L-1; ++ii) {
          for (int ij = 0; ij <= 20-1; ++ij) {
              if((ii-ij)>=0) fSum+=meanVec_log_s[ii-ij] * h15[ij];
          }
          y15[ii]=fSum;
          y15_abs[ii]=std::abs(fSum);
          y15_abs_mean+=y15_abs[ii];
          fSum=0;
      }
      y15_abs_mean/=float(L);

      float s_threshold=4*y15_abs_mean, maxx=0, minn=0;
      int maxIdx=-1, minIdx=-1, state=0, prevState=0, posCount=-1, nevCount=-1;

      for(int j1=0; j1<=L-1; ++j1){ out[j1]=tempBuf[j1]; tempBuf[j1]=0;
         if (y15[j1]>s_threshold)
            { state=1;
              pPos[j1]= y15[j1];
              if (pPos[j1]>maxx){
            	maxx=pPos[j1]; maxIdx=j1; }
            }
   	 else if (y15[j1]<-1*s_threshold)
   	    { state=-1;
              pNev[j1]= y15[j1];
              if (pNev[j1]<minn)
            	{minn=pNev[j1]; minIdx=j1;}
            }
    	 else state=0;
    	 if (prevState==1 && state<=0) {
              pPosMax[maxIdx]=maxx;
       	      posCount=posCount+1; posVec[posCount]=maxx; posLoc[posCount]=maxIdx;
              maxx=0;maxIdx=-1;
            }
         if (prevState==-1 && state>=0) {
              pNevMin[minIdx]=minn;
              nevCount=nevCount+1;
              nevVec[nevCount]=minn; nevLoc[nevCount]=minIdx; nevPosIdx[nevCount]=posCount;
              minn=0;minIdx=-1;
              }
         prevState=state;
      }
      int tmpPosIdx=-1, posIdxLimit=0, writeIdx; 
      float L1=L;
      for (int jj=0; jj<=nevCount-1; jj++) {
    	tmpPosIdx=nevPosIdx[jj];
    	while (tmpPosIdx>posIdxLimit) {
        	float t1=posVec[tmpPosIdx]; float t2=nevVec[jj];
           
        	if (std::abs(t1+t2)<0.6*std::max(t1,-1*t2)) {
           	   chanCount=chanCount+1;chLtedge[chanCount]=posLoc[tmpPosIdx];chRtedge[chanCount]=nevLoc[jj];
            	   chCen[chanCount]=std::floor((chLtedge[chanCount]+chRtedge[chanCount])/2.0);
                   chMag[chanCount]=del_meanVec_log[chCen[chanCount]];
                   chBw[chanCount]=chRtedge[chanCount]-chLtedge[chanCount]+1;
                   writeIdx=chanCount*4;
                   tempBuf[writeIdx+0].real(chanCount);
                   tempBuf[writeIdx+1].real((frequency + ((chCen[chanCount]-(L1/2))*(d_sample_rate/L1)))/1000000.0);
                   tempBuf[writeIdx+2].real((10*chMag[chanCount]));
                   tempBuf[writeIdx+3].real((chBw[chanCount]*(d_sample_rate/L1))/1000.0);//	std::cout<< tempBuf[writeIdx + 1] << "\n";
                   posIdxLimit=nevPosIdx[jj];
                   break;
        	}
        	tmpPosIdx=tmpPosIdx-1;
    }
} tmp_chanCount=chanCount+1;

    }
    // else if(cCount>0 && (cCount==51)) {
      else if(cCount>0 && (cCount==4)) {
      
      	      for(int j2=0; j2<=L-1; ++j2) {	
      	      	if(j2==0) { out[j2]=tmp_chanCount;  
      	      	//  std::cout<<tmp_chanCount<<"  "<<frequency<<"  "<<tempBuf[j2]<<"  "<<tempBuf[j2+1]<<"  "<<tempBuf[j2+2]<<"  "<<tempBuf[j2+3]<<"  "<<tempBuf[j2+4]<<"  "<<tempBuf[j2+5]<<"  "<<tempBuf[j2+6]<<"  "<<tempBuf[j2+7]<<"\n";
      	      	}
      	      	else out[j2]=tempBuf[j2-1];      	      
      	      }
      }
      
      else {
      	      for(int j2=0; j2<=L-1; ++j2) out[j2]=0;
      }


    }
    

//    memcpy(out, in - d_burst_pre_len, sizeof(gr_complex) * noutput_items);

    // Tell runtime system how many output items we produced.
    return noutput_items;
}

} /* namespace iridium */
} /* namespace gr */
