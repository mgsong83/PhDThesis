
class CCOR {

private:
  TGraphAsymmErrors *mz62, *mz44,  *mz31,  *pttkt;
  TGraphErrors *pttjt;
  TGraphErrors *fff53, *fff500, *fff1000;
  double msz, lwd;
  
public:
  
  CCOR(double size=1.3){
    msz = size;
    lwd = 2.1;
    
    const int NC = 8;  //CCOR sqrt{s}=62.4 GeV
    double s = 62.4;
    double pt[NC] = { 3.39173, 4.42291, 5.44249, 6.43632, 7.44221, 8.44988, 9.458, 10.7077};
    double ept[NC] = {NC * 0};
    double zt[NC] = {0.785493, 0.817104, 0.837532, 0.845597, 0.863079, 0.873503, 0.882162, 0.905629};
    double eztl[NC] = {0.00529367, 0.00882279, 0.00352912, 0.00529367, 0.00294093, 0.00529367, 0.00705823, 0.00352912};
    double ezth[NC] = {0.00529367, 0.00999916, 0.00352912, 0.00647005, 0.00235274, 0.00470549, 0.00529367, 0.00294093};
    for (int i = 0; i < NC; i++) pt[i] /= s / 2;
    mz62 = new TGraphAsymmErrors(NC, pt, zt, ept, ept, eztl, ezth);
    mz62->SetMarkerStyle(24); mz62->SetMarkerSize(msz);

    const int NC = 8; //CCOR sqrt{s}=44.8 GeV
    double s = 44.8;
    double pt[NC] = { 3.39917, 4.44329, 5.42285, 6.42814, 7.42138, 8.4545, 9.41858, 10.6976};
    double zt[NC] = {0.813144, 0.850642, 0.858113, 0.877947, 0.888365, 0.912329, 0.923911, 0.945625};
    double eztl[NC] = {0.0082346, 0.00235274, 0.00352912, 0.00294093, 0.00235274, 0.00294093, 0.00529367, 0.00764642};
    double ezth[NC] = {0.00705823, 0.00294093, 0.00176456, 0.00352912, 0.00294093, 0.00235274, 0.00647005, 0.00999916};
    for (int i = 0; i < NC; i++) pt[i] /= s / 2;
    mz44 = new TGraphAsymmErrors(NC, pt, zt, ept, ept, eztl, ezth);
    mz44->SetMarkerStyle(25); mz44->SetMarkerSize(msz);

    const int NC = 5;   //CCOR sqrt{s}=31 GeV
    double s = 29.7641;
    double pt[NC] = { 3.39069, 4.43615, 5.45305, 6.44287, 7.44102};
    double zt[NC] = { 0.84667, 0.878875, 0.909891, 0.933837, 0.924845};
    double eztl[NC] = {0.00588186, 0.00588186, 0.00588186, 0.00529367, 0.0117637};
    double ezth[NC] = {0.00705823, 0.0082346, 0.00588186, 0.00588186, 0.0141165};
    for (int i = 0; i < NC; i++) pt[i] /= s / 2;
    mz31 = new TGraphAsymmErrors(NC, pt, zt, ept, ept, eztl, ezth);
    mz31->SetMarkerStyle(26); mz31->SetMarkerSize(msz);

    const int NC=4;
    double ptt[NC] = {   3.413,  5.43453,  7.42318,  9.41464};
    double eptt[NC] = {0};
    double jt[NC] = {0.386444, 0.408286, 0.394021, 0.395808};
    double ejt[NC] = {0.0548461, 0.014046, 0.00668855, 0.0234099};
    for(int ii=0;ii<NC;ii++){
      jt[ii]  *= sqrt(pi);
      ejt[ii] *= sqrt(pi);
    }
    pttjt = new TGraphErrors(NC, ptt, jt, eptt, ejt);
    pttjt->SetMarkerStyle(26); pttjt->SetMarkerSize(msz);

    const int CNC = 4;
    double CCORpt [CNC] = {3.29904, 5.30672, 7.27197, 9.29834};
    double CCORkt [CNC] = {0.539998, 0.651385, 0.689406, 0.800769};
    double elCCORkt[CNC] = {0.0303784, 0.0177207, 0.0151892, 0.0354415};
    double ehCCORkt[CNC] = {0.0253154, 0.0202523, 0.0177207, 0.0405046};
    double empt_l [CNC] = {CNC * 0};
    double mz_assoc3[CNC] = {0.6, 0.72, 0.754, 0.77 };
    for (int ik = 0; ik < CNC; ik++){
      CCORkt[ik] *= sqrt(pi);
      elCCORkt[ik] *= sqrt(pi);
    }
    pttkt = new TGraphAsymmErrors(CNC, CCORpt, CCORkt, empt_l, empt_l, elCCORkt, ehCCORkt);
    pttkt->SetMarkerStyle(26); pttkt->SetMarkerSize(msz);

    //Feynman==================================
    const int NC = 4;   //FFFPRD18, 3320 (1978)
    double s = 53;
    double pt[] = {2, 4, 7, 9};
    double zt[] = {0.68, 0.75, 0.81, 0.84};
    double ezt[] = {0, 0, 0, 0};
    for (int i = 0; i < NC; i++) pt[i] /= s / 2;
    fff53 = new TGraphErrors(NC, pt, zt, ept, ezt);
    fff53->SetLineStyle(2); fff53->SetLineColor(4); fff53->SetLineWidth(lwd);
    
    const int NC = 2;   //FFFPRD18, 3320 (1978)
    double s = 500;
    double pt[] = {10, 30};
    double zt[] = {0.58, 0.69};
    double ezt[] = {0};
    for (int i = 0; i < NC; i++) pt[i] /= s / 2;
    fff500 = new TGraphErrors(NC, pt, zt, ept, ezt);
    fff500->SetLineStyle(2); fff500->SetLineColor(4);fff500->SetLineWidth(lwd);

    const int NC = 2;   //FFFPRD18, 3320 (1978)
    double s = 1000;
    double pt[NC] = {10, 30};
    double zt[NC] = {0.53, 0.61};
    double ezt[NC] = {0};
    for (int i = 0; i < NC; i++) pt[i] /= s / 2;
    fff1000 = new TGraphErrors(NC, pt, zt, ept, ezt);
    fff1000->SetLineStyle(2); fff1000->SetLineColor(4);fff1000->SetLineWidth(lwd);

  }
  
public:

  TGraphAsymmErrors *GetMz62(){return mz62;}
  TGraphAsymmErrors *GetMz44(){return mz44;}
  TGraphAsymmErrors *GetMz31(){return mz31;}
  TGraphErrors *GetJt()  {return pttjt;}
  TGraphAsymmErrors *GetKt()  {return pttkt;}

  TGraphErrors *GetFFF53()  {return fff53;}
  TGraphErrors *GetFFF500() {return fff500;}
  TGraphErrors *GetFFF1000(){return fff1000;}
  
};
