/*******************************************************************************
  file name: CollectOneEvent.cxx
  author: Zhe Yang
  created: 01/25/2019
  last modified: 02/05/2019

  description:
  Collect new 32-bit data and manage signal length and trigger length

  remark:
  Not used in current DecodeRawData.cxx code
  
*******************************************************************************/

#include <iostream>
#include <math.h>

const bool GOTNEWEVENT = true;
const bool NOTGOTNEWEVENT = false;
static unsigned int processed_data_count = 0;

void CollectOneEvent(unsigned int *trigger_cache_length, 
                    unsigned int trigger_cache[128][5],
                    unsigned int *signal_cache_length,
                    unsigned int signal_cache[128][5],
                    unsigned int *new_data_line) {
  bool new_trigger_flag = (new_data_line[0] == 4 || new_data_line[0] == 5) 
                          && (new_data_line[1] == 1);
  bool new_signal_flag = (new_data_line[0] == 4 || new_data_line[0] == 5) 
                         && (new_data_line[1] != 1);

  if (new_trigger_flag) {
    trigger_cache[*trigger_cache_length][0] = new_data_line[0];
    trigger_cache[*trigger_cache_length][1] = new_data_line[1];
    trigger_cache[*trigger_cache_length][2] = new_data_line[2];
    trigger_cache[*trigger_cache_length][3] = new_data_line[3];
    trigger_cache[*trigger_cache_length][4] = new_data_line[4];
    (*trigger_cache_length)++;
  }
  if (new_signal_flag) {
    signal_cache[*signal_cache_length][0] = new_data_line[0];
    signal_cache[*signal_cache_length][1] = new_data_line[1];
    signal_cache[*signal_cache_length][2] = new_data_line[2];
    signal_cache[*signal_cache_length][3] = new_data_line[3];
    signal_cache[*signal_cache_length][4] = new_data_line[4];
    (*signal_cache_length)++;
  }

  if (*trigger_cache_length == 0 && *signal_cache_length == 0) {
    processed_data_count = 0;
  } else {
    processed_data_count++;
  }

}