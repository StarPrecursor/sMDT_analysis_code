/*******************************************************************************
  file name: FitTDCSlope.h
  author: Zhe Yang
  created: 02/25/2019
  last modified: 02/25/2019

  description:
  -Header for FitTDCSlope.cxx
*******************************************************************************/

#ifndef FITTDCSLOPE_H_
#define FITTDCSLOPE_H_

#include <string>
#include <math.h>

Double_t FermiDiracFunction(Double_t *x, Double_t *par) {
  return (par[0] + par[3] * x[0]) / (1 + exp((par[1] - x[0]) / par[2]));
}

#endif // FITTDCSLOPE_H_