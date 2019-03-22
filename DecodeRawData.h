/*******************************************************************************
  file name: DecodeRawData.h
  author: Zhe Yang
  created: 01/25/2019
  last modified: 01/31/2019

  description:
  -Header for DecodeRawData.cxx
*******************************************************************************/

#ifndef DECODERAWDATA_H_
#define DECODERAWDATA_H_

int DecodeRawData();

const Int_t MAX_TDC_QUANTITY = 18;
const Int_t MAX_TDC_CHANNEL_QUANTITY = 24;
const Int_t MAX_TUBE_LAYER = 8;
const Int_t MAX_TUBE_COLUMN = 54;
const Long_t MAX_WORD_QUANTITY = 2000000;

#endif // DECODERAWDATA_H_