/*******************************************************************************
  file name: GetEventTrigger.cxx
  author: Zhe Yang
  created: 01/25/2019
  last modified: 02/05/2019

  description:
  Get best trigger signal from collected trigger data

  remark:
  Not used in current DecodeRawData.cxx code
*******************************************************************************/

#include <iostream>
#include <math.h>

int GetEventTrigger (unsigned int trigger_cache_length,
                      unsigned int trigger_cache[128][5], 
                      unsigned int signal_mean_coarse,
                      unsigned int *trigger_data_line) {
  int selected_trigger_id = 128;
  long current_coarse_difference = 2048;
  if (trigger_cache_length == 0) {
    return 1;
  }

  long temp_coarse = 0;
  long temp_mean_coarse = signal_mean_coarse;
  for (int trigger_id = 0; trigger_id < trigger_cache_length; trigger_id++) {
    temp_coarse = trigger_cache[trigger_id][3];
    if (abs(temp_coarse - temp_mean_coarse) < current_coarse_difference) {
      current_coarse_difference = abs(temp_coarse - temp_mean_coarse);
      selected_trigger_id = trigger_id;
    }
    cout << "current_coarse_difference= " << current_coarse_difference;
    cout << "  selected_trigger_id=" << selected_trigger_id << endl;
  }
  
  trigger_data_line[0] = trigger_cache[selected_trigger_id][0];
  trigger_data_line[1] = trigger_cache[selected_trigger_id][1];
  trigger_data_line[2] = trigger_cache[selected_trigger_id][2];
  trigger_data_line[3] = trigger_cache[selected_trigger_id][3];
  trigger_data_line[4] = trigger_cache[selected_trigger_id][4];

  return 0;
}