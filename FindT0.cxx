/*******************************************************************************
  file name: FindT0.cxx
  author: Zhe Yang
  created: 02/25/2019
  last modified: 02/25/2019

  description:
  -Find T0, Tmax, and other parameters for the ntuple generated by
   DecodeRawData.cxx

  remark:
  -Use the T0fit class wrote by Edward Diehl
*******************************************************************************/

#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>

#include "DecodeRawData.h"
#include "FitTDCSlope.h"
//#include "T0Fit.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

//string s_input_filename = "run00187664_20190222.dat";
string s_input_filename = "run00187691_20190301.dat";

using namespace std;

int TdcFit(TH1F *input_hist, Double_t fit_output_value[2][5],
           Double_t fit_initial_value_t0[5], Double_t fit_initial_value_tmax[5],Bool_t plot_flag);

void FindT0() {
  // get histogram name
  char directory_name[256];
  char hist_name[256];
  char temp_name[256];
  char input_directory_name[256];
  char input_path_name[256];
  char output_directory_name[256];
  char output_path_name[256];
  Double_t output_data[2][5];
  TCanvas *plot_canvas = new TCanvas();

  sprintf(input_path_name, s_input_filename.c_str());
  strcat(input_path_name, ".dir/");
  strcpy(input_directory_name, input_path_name);
  strcat(input_path_name, s_input_filename.c_str());
  strcat(input_path_name, ".out.root");
  cout << "input path name: " << input_path_name << endl;
  TFile *p_input_rootfile = new TFile(input_path_name);
  TH1F *p_mezzanine_drift_time_hist[MAX_TDC_QUANTITY];
  TH1F *p_channel_drift_time_hist[MAX_TDC_QUANTITY][MAX_TDC_CHANNEL_QUANTITY];

  // fit and plot
  ////////////chdir(input_directory_name);
  strcpy(output_directory_name, input_directory_name);
  strcat(output_directory_name, "T0Fit_output");
  if (mkdir(output_directory_name, 0777) == -1) {
    cerr << strerror(errno) << endl;
  }
  for(Int_t tdc_id = 0; tdc_id < MAX_TDC_QUANTITY; tdc_id++) {
    sprintf(directory_name, "TDC_%02d_of_%02d_Time_Spectrum", tdc_id, MAX_TDC_QUANTITY);
    //sprintf(temp_name, "tdc_%d_tdc_time_spectrum", tdc_id);
    sprintf(temp_name, "tdc_%d_tdc_time_spectrum_corrected", tdc_id);
    strcpy(hist_name, directory_name);
    strcat(hist_name, "/");
    strcat(hist_name, temp_name);
    p_input_rootfile->GetObject(hist_name, p_mezzanine_drift_time_hist[tdc_id]);

    char output_subdirectory_name[256];
    sprintf(temp_name, "TDC_%02d_of_%02d_Time_Spectrum", tdc_id,
            MAX_TDC_QUANTITY);
    strcpy(output_subdirectory_name, output_directory_name);
    strcat(output_subdirectory_name, "/");
    strcat(output_subdirectory_name, temp_name);
    if (mkdir(output_subdirectory_name, 0777) == -1) {
      cerr << strerror(errno) << endl;
    }

    strcpy(output_path_name, output_subdirectory_name);
    strcat(output_path_name, "/");
    //sprintf(temp_name, "tdc_%d_tdc_time_spectrum", tdc_id);
    sprintf(temp_name, "tdc_%d_tdc_time_spectrum_corrected", tdc_id);
    strcat(output_path_name, temp_name);
    strcat(output_path_name, ".png");

    if (tdc_id >= 8 && tdc_id <= 11) { // currently we only use MEZZ#8~11
      plot_canvas->cd();
      p_mezzanine_drift_time_hist[tdc_id]->Draw();
      TdcFit(p_mezzanine_drift_time_hist[tdc_id], output_data, nullptr, nullptr,
             false);
      plot_canvas->SaveAs(output_path_name);
    }

    for (Int_t ch_id = 0; ch_id < MAX_TDC_CHANNEL_QUANTITY; ch_id++) {

      //sprintf(temp_name, "tdc_%d_channel_%d_tdc_time_spectrum", tdc_id, ch_id);
      sprintf(temp_name, "tdc_%d_channel_%d_tdc_time_spectrum_corrected", tdc_id, ch_id);
      strcpy(hist_name, directory_name);
      strcat(hist_name, "/");
      strcat(hist_name, temp_name);
      cout << "hist name:" << hist_name << endl;
      //getchar();

      p_input_rootfile->GetObject(hist_name,
                                  p_channel_drift_time_hist[tdc_id][ch_id]);

      if (tdc_id >= 8 && tdc_id <= 11) { // currently we only use MEZZ#8~11
        plot_canvas->cd();
        p_channel_drift_time_hist[tdc_id][ch_id]->Draw();
        TdcFit(p_channel_drift_time_hist[tdc_id][ch_id], output_data, nullptr,
              nullptr, false);

        strcpy(output_path_name, output_subdirectory_name);
        strcat(output_path_name, "/");
        //sprintf(temp_name, "tdc_%d_channel_%d_tdc_time_spectrum", tdc_id, ch_id);
        sprintf(temp_name, "tdc_%d_channel_%d_tdc_time_spectrum_corrected", tdc_id, ch_id);
        strcat(output_path_name, temp_name);
        strcat(output_path_name, ".png");
        plot_canvas->SaveAs(output_path_name);
      }
    }
  }

  cout << "finished" << endl;
}

Int_t TdcFit(TH1F *input_hist,
             Double_t fit_output_value[2][5],
             Double_t fit_initial_value_t0[5] = nullptr,
             Double_t fit_initial_value_tmax[5] = nullptr,
             Bool_t plot_flag = false) {
  // t0 fitting parameters
  Double_t t0_amplitude, t0, t0_slope, t0_top_slope, t0_noise, t0_fit_range_min, t0_fit_range_max;

  // histogram features
  Double_t max_bin_value;

  // Protect for null histograms and low statistics time spectra
  if (input_hist == nullptr || input_hist->GetEntries() < 3000.0) return 1;

  // default text box locations
  gStyle->SetOptStat(0);
  TPaveText *entry_pave_text = new TPaveText(0.65, 0.8, 0.9, 0.9, "NDC");
  TPaveText *t0_pave_text = new TPaveText(0.65, 0.7, 0.9, 0.8, "NDC");
  TPaveText *tmax_pave_text = new TPaveText(0.65, 0.6, 0.9, 0.7, "NDC");
  TPaveText *t_drift_time_pave_text = new TPaveText(0.65, 0.5, 0.9, 0.6, "NDC");

  // Perform T0 fit
  // initialize parameters for t0 fittind, use default if not specified
  if (fit_initial_value_t0 != nullptr) {
    t0_amplitude = fit_initial_value_t0[0];
    t0 = fit_initial_value_t0[1];
    t0_slope = fit_initial_value_t0[2];
    t0_top_slope = fit_initial_value_t0[3];
    t0_noise = fit_initial_value_t0[4];
  } else {
    max_bin_value = input_hist->GetMaximum();
    t0_amplitude = max_bin_value / 1.1;
    t0 = input_hist->GetBinCenter(input_hist->FindFirstBinAbove(0.45 * max_bin_value));
    t0_slope = 2.5;
    t0_top_slope = 0;
    t0_noise = t0_amplitude / 100.0;
  }
  t0_fit_range_min = t0 - 100.0;
  t0_fit_range_max = input_hist->GetBinCenter(input_hist->GetMaximumBin());
  // Define TF1 function for T0 fit and perform fitting
  TF1 *t0_fit_function =  new TF1("t0_fit_function", FermiDiracFunction,
                                t0_fit_range_min, t0_fit_range_max, 5);
  t0_fit_function->SetParameters(t0_amplitude, t0, t0_slope, t0_top_slope);
  t0_fit_function->SetParLimits(5, 0, t0_amplitude / 5.0); //Do not allow negative t0_noise
  input_hist->Fit("t0_fit_function", "R");

  if (fit_output_value != nullptr) {
    fit_output_value[0][0] = t0_fit_function->GetParameter(0);
    fit_output_value[0][1] = t0_fit_function->GetParameter(1);
    fit_output_value[0][2] = t0_fit_function->GetParameter(2);
    fit_output_value[0][3] = t0_fit_function->GetParameter(3);
    fit_output_value[0][4] = t0_fit_function->GetParameter(4);
  }

  t0_pave_text->AddText(Form("T0 = %.2lf #pm %.2lf ns",
                      t0_fit_function->GetParameter(1),
                      t0_fit_function->GetParError(1)));
  t0_pave_text->AddText(Form("Slope = %.2lf #pm %.2lf ns",
                      t0_fit_function->GetParameter(2),
                      t0_fit_function->GetParError(2)));
  t0_pave_text->SetTextColor(kRed);

  // Perform Tmax fit //
  // tmax fitting parameters
  Double_t tmax_amplitude, tmax, tmax_slope, tmax_top_slope, tmax_noise, tmax_fit_range_min, tmax_fit_range_max;
  // initialize parameters for tmax fittind, use default if not specified
  if (fit_initial_value_t0 != nullptr) {
    tmax_amplitude = fit_initial_value_tmax[0];
    tmax = fit_initial_value_tmax[1];
    tmax_slope = fit_initial_value_tmax[2];
    tmax_top_slope = fit_initial_value_tmax[3];
    tmax_noise = fit_initial_value_tmax[4];
  } else {
    tmax_amplitude = max_bin_value / 2.0;
    tmax = input_hist->GetBinCenter(input_hist->FindLastBinAbove(0.1 * max_bin_value));
    tmax_slope = -8.0;
    tmax_top_slope = 0;
    tmax_noise = t0_fit_function->GetParameter(2);
  }
  tmax_fit_range_min = tmax - 80.;
  tmax_fit_range_max = tmax + 200.;
  // Define TF1 function for Tmax fit
  TF1 *tmax_fit_function = new TF1("tmax_fit_function", FermiDiracFunction,
                                    tmax_fit_range_min, tmax_fit_range_max, 5);
  tmax_fit_function->SetParameters(tmax_amplitude, tmax, tmax_slope,
                                   tmax_top_slope);
  tmax_fit_function->SetParLimits(1, t0_fit_function->GetParameter(1),
                                  t0_fit_function->GetParameter(1) + 300);
  tmax_fit_function->SetLineColor(kGreen);
  input_hist->Fit("tmax_fit_function", "r+"); //the "+" is to add the function to the list
  tmax_fit_function->SetLineColor(kGreen);

  if (fit_output_value != nullptr) {
    fit_output_value[1][0] = tmax_fit_function->GetParameter(0);
    fit_output_value[1][1] = tmax_fit_function->GetParameter(1);
    fit_output_value[1][2] = tmax_fit_function->GetParameter(2);
    fit_output_value[1][3] = tmax_fit_function->GetParameter(3);
    fit_output_value[1][4] = tmax_fit_function->GetParameter(4);
  }

  tmax_pave_text->AddText(Form("Tmax %.1lf #pm %.1lf ns",
                      tmax_fit_function->GetParameter(1),
                      tmax_fit_function->GetParError(1)));
  tmax_pave_text->AddText(Form("Slope %.1lf #pm %.1lf /ns",
                      tmax_fit_function->GetParameter(2),tmax_fit_function->GetParError(2)));
  tmax_pave_text->SetTextColor(kGreen);

  // Compute max drift time of time spectra: Tmax-T0
  Double_t drift_time_max = tmax_fit_function->GetParameter(1) -
                            t0_fit_function->GetParameter(1);
  Double_t drift_time_max_err = sqrt(tmax_fit_function->GetParError(0) *
                                     tmax_fit_function->GetParError(0) +
                                     t0_fit_function->GetParError(0) *
                                     t0_fit_function->GetParError(0));
  t_drift_time_pave_text->AddText(Form("Max drift time"));
  t_drift_time_pave_text->AddText(Form("%.1lf #pm %.1lf ns", drift_time_max,
                                        drift_time_max_err));
  t_drift_time_pave_text->SetTextColor(kViolet);
  t_drift_time_pave_text->Draw();

  // draw TPaveText in order to avoid overlap between each TPaveText
  entry_pave_text->AddText(Form("Total entries: %.0f",
                                 input_hist->GetEntries()));
  entry_pave_text->SetTextColor(kBlack);
  entry_pave_text->Draw();
  t0_pave_text->Draw();
  tmax_pave_text->Draw();
  t_drift_time_pave_text->Draw();

  return 0;
}
// end TdcFit //////////////////////////////////////////////////////////////////