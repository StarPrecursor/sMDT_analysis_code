#include <iostream>

#include "CheckEvent.cxx"
#include "GetEventTrigger.cxx"

using namespace std;

int testmain() {
  unsigned int signal_cache_length;
  unsigned int signal_cache[128][5];
  unsigned int mean_coarse;

  unsigned int trigger_cache_length;
  unsigned int trigger_cache[128][5];
  unsigned int trigger_data_line[5];

  signal_cache_length = 3;
  signal_cache[0][3] = 10;
  signal_cache[1][3] = 10;
  signal_cache[2][3] = 10;

  trigger_cache_length = 2;
  trigger_cache[0][3] = 11;
  trigger_cache[1][3] = 12;
  trigger_cache[2][3] = 9;

  bool result;

  cout << "testing Check Event... " << endl;

  result = CheckEvent(signal_cache_length, signal_cache, &mean_coarse);

  cout << "result= " << result << endl;

  cout << "testing GetEventTrigger... mean coarse=" << mean_coarse << endl;

  int selected_trigger_id = 128;
  int current_coarse_difference = 2048;
  if (trigger_cache_length == 0) {
    return 1;
  }

  int temp_coarse = 0;
  int temp_mean_coarse = mean_coarse;
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
  
  cout << "trigger_data_line coarse: " << trigger_data_line[3] << endl;

  return 0;
}