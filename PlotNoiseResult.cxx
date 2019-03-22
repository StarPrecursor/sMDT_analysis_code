#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TMultiGraph.h>

using namespace std;

Int_t PlotNoiseResult() {
  // open data file
  ifstream input_data_file;
  ifstream input_data_file2;
  //input_data_file.open("noise_result.txt");
  //input_data_file.open("noise_result2.txt");
  input_data_file.open("noise_result3.txt");
  input_data_file2.open("noise_result_m1.txt");

  // Plot noise distribution
  Double_t noise_rate;
  Double_t noise_rate2;
  TH1F *hist_noise_distribution08_m1 = new TH1F("Mezz08 noise method1", "Mezz08 noise method1", 24, 0, 24);
  TH1F *hist_noise_distribution08_m2 = new TH1F("Mezz08 noise method2", "Mezz08 noise method2", 24, 0, 24);
  TH1F *hist_noise_distribution09_m1 = new TH1F("Mezz09 noise method1", "Mezz09 noise method1", 24, 0, 24);
  TH1F *hist_noise_distribution09_m2 = new TH1F("Mezz09 noise method2", "Mezz09 noise method2", 24, 0, 24);
  TH1F *hist_noise_distribution10_m1 = new TH1F("Mezz10 noise method1", "Mezz10 noise method1", 24, 0, 24);
  TH1F *hist_noise_distribution10_m2 = new TH1F("Mezz10 noise method2", "Mezz10 noise method2", 24, 0, 24);
  TH1F *hist_noise_distribution11_m1 = new TH1F("Mezz11 noise method1", "Mezz11 noise method1", 24, 0, 24);
  TH1F *hist_noise_distribution11_m2 = new TH1F("Mezz11 noise method2", "Mezz11 noise method2", 24, 0, 24);
  for (Int_t tube_id = 0; tube_id != 24; tube_id++ ) {
    input_data_file >> noise_rate;
    input_data_file2 >> noise_rate2;
    hist_noise_distribution08_m1->Fill(hist_noise_distribution08_m1->GetBinCenter(tube_id + 1), noise_rate);
    hist_noise_distribution08_m2->Fill(hist_noise_distribution08_m2->GetBinCenter(tube_id + 1), noise_rate2);
  }
  for (Int_t tube_id = 0; tube_id != 24; tube_id++ ) {
    input_data_file >> noise_rate;
    input_data_file2 >> noise_rate2;
    hist_noise_distribution09_m1->Fill(hist_noise_distribution09_m1->GetBinCenter(tube_id + 1), noise_rate);
    hist_noise_distribution09_m2->Fill(hist_noise_distribution09_m2->GetBinCenter(tube_id + 1), noise_rate2);
  } for (Int_t tube_id = 0; tube_id != 24; tube_id++ ) {
    input_data_file >> noise_rate;
    input_data_file2 >> noise_rate2;
    hist_noise_distribution10_m1->Fill(hist_noise_distribution10_m1->GetBinCenter(tube_id + 1), noise_rate);
    hist_noise_distribution10_m2->Fill(hist_noise_distribution10_m2->GetBinCenter(tube_id + 1), noise_rate2);
  } for (Int_t tube_id = 0; tube_id != 24; tube_id++ ) {
    input_data_file >> noise_rate;
    input_data_file2 >> noise_rate2;
    hist_noise_distribution11_m1->Fill(hist_noise_distribution11_m1->GetBinCenter(tube_id + 1), noise_rate);
    hist_noise_distribution11_m2->Fill(hist_noise_distribution11_m2->GetBinCenter(tube_id + 1), noise_rate2);
  }
  TCanvas *canvas_mezz08 = new TCanvas("Mezz08 noise", "Mezz08 noise method1");
  hist_noise_distribution08_m1->SetTitle("Mezz08 noise level;tube id;rate/Hz");
  hist_noise_distribution08_m1->SetLineColor(kBlack);
  hist_noise_distribution08_m1->SetLineWidth(3);
  hist_noise_distribution08_m1->SetStats(0);
  hist_noise_distribution08_m1->Draw("hist");
  hist_noise_distribution08_m2->Draw("hist same");
  canvas_mezz08->SaveAs("Mezz08_noise.png");
  TCanvas *canvas_mezz09 = new TCanvas("Mezz09 noise", "Mezz09 noise");
  hist_noise_distribution09_m1->SetTitle("Mezz09 noise level;tube id;rate/Hz");
  hist_noise_distribution09_m1->SetLineColor(kBlack);
  hist_noise_distribution09_m1->SetLineWidth(3);
  hist_noise_distribution09_m1->SetStats(0);
  hist_noise_distribution09_m1->Draw("hist");
  hist_noise_distribution09_m2->Draw("hist same");
  canvas_mezz09->SaveAs("Mezz09_noise.png");  
  TCanvas *canvas_mezz10 = new TCanvas("Mezz10 noise", "Mezz10 noise method1");
  hist_noise_distribution10_m1->SetTitle("Mezz10 noise level;tube id;rate/Hz");
  hist_noise_distribution10_m1->SetLineColor(kBlack);
  hist_noise_distribution10_m1->SetLineWidth(3);
  hist_noise_distribution10_m1->SetStats(0);
  hist_noise_distribution10_m1->Draw("hist");
  hist_noise_distribution10_m2->Draw("hist same");
  canvas_mezz10->SaveAs("Mezz10_noise.png");
  TCanvas *canvas_mezz11 = new TCanvas("Mezz11 noise", "Mezz11 noise");
  hist_noise_distribution11_m1->SetTitle("Mezz11 noise level;tube id;rate/Hz");
  hist_noise_distribution11_m1->SetLineColor(kBlack);
  hist_noise_distribution11_m1->SetLineWidth(3);
  hist_noise_distribution11_m1->SetStats(0);
  hist_noise_distribution11_m1->Draw("hist");
  hist_noise_distribution11_m2->Draw("hist same");
  canvas_mezz11->SaveAs("Mezz11_noise.png");

  return 0;
}