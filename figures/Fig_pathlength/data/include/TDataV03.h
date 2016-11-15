#include "Tjtkt.h"
#define MAXD 100;


class TDataV03 {

public:

  TDataV03(){
    jtkt = new Tjtkt();
    //cout<<"Init: "<<jtkt<<endl;
  }
  
public:
  
  Tjtkt *jtkt;
  TGraphErrors *grr[MAXD][MAXD][MAXD];  
  double ptlim[MAXD][2];
  double emptt[MAXD], empta[MAXD];
  double x, sys[3][100];
  
  enum enGrPT  { ptt1540, ptt1530, ptt24, ptt25, ptt1850, ptt1851, ptt1852, pttjt, pttjt36,
		  pta34, pta45, pta5580, pta510, pta710,
		  MCkt2_ptt, MCkt3_ptt};
  enum enGrWhat { tSigmaN, tSigmaF, tJt, tPout, tZKT, tKT, tZT, tXhq, tKTa3, tZTa3, tXhqa3, tKTa5, tZTa5, tXhqa5 };
  enum enGrSys { tReal, tSysL, tSysH };
  enGrPT grtyp;
  
  TGraphErrors *GetGrr(enGrPT grtyp, enGrWhat twhat, enGrSys type){
    return grr[grtyp][twhat][type];
  }
  TGraphErrors *GetGrr(enGrPT grtyp, enGrWhat twhat){
    return grr[grtyp][twhat][0];
  }
  void WhatGetPtLims(){ cout<<"????"<<endl;}
  double GetPTs(enGrPT grtyp, int i){ return ptlim[grtyp][i]; }
  
  //--------------------------------
  void GetPttDist() {
  //--------------------------------
    //piv03_jt
    const int NC =  7 ;
    double sigN  [NC]= {0.235212, 0.217575, 0.200074, 0.192325, 0.202032, 0.189764, 0.175386}; 
    double esigN [NC]= {0.0011917, 0.00188785, 0.00294503, 0.0044955, 0.00659598, 0.00684463, 0.0097398}; 
    double sigF  [NC]= {0.515988, 0.456308, 0.402666, 0.384932, 0.348294, 0.309224, 0.211282}; 
    double esigF [NC]= {0.00539758, 0.00737429, 0.00989722, 0.0134954, 0.016223, 0.0147617, 0.0147013}; 
    double pout  [NC]= { 1.17121, 1.0298, 0.912181, 0.886833, 0.806236, 0.717642, 0.508798}; 
    double epout [NC]= {0.0162147, 0.0198087, 0.0240964, 0.0337279, 0.038672, 0.0338021, 0.0356728}; 
    double mptt  [NC]= { 3.3613, 4.37786, 5.39766, 6.41083, 7.42298, 8.7485, 11.2319}; 
    double mpta  [NC]= { 2.13976, 2.19094, 2.23755, 2.27325, 2.29128, 2.32995, 2.40858}; 
    double SoverB  [NC]= { 2.29017, 3.44142, 4.93059, 6.45895, 9.82344, 11.3273, 18.6435}; 
    double bpt[NC+1] = {1.50, 4.00};
    grtyp = pttjt;
    ptlim[grtyp][0]=bpt[0]; ptlim[grtyp][1]=bpt[1]; 
    jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mptt, sigN, emptt, esigN);
    grr[grtyp][tSigmaN][tReal]->SetMarkerStyle(20);
    grr[grtyp][tJt][tReal]    = jtkt->GetJTy_ptt();
    grr[grtyp][tJt][tReal]->SetMarkerStyle(20);

    const int NC = 5 ;
    double sigN [NC] = {0.169453, 0.163354, 0.152024, 0.1428, 0.128654};
    double esigN [NC] = {0.00140279, 0.00215899, 0.00310616, 0.00373597, 0.00567196};
    double sigF [NC] = {0.469352, 0.366998, 0.35992, 0.347772, 0.209094};
    double esigF [NC] = {0.0071926, 0.00976757, 0.0121463, 0.0145803, 0.0109607};
    double pout [NC] = { 1.59801, 1.45702, 1.40399, 1.26933, 0.781135};
    double epout [NC] = {0.0286966, 0.0414733, 0.0484922, 0.0544153, 0.0414274};
    double mptt [NC] = { 3.3613, 4.37786, 5.39766, 6.68388, 9.31824};
    double mpta [NC] = { 3.82484, 3.87342, 3.88995, 3.97298, 4.03606};
    double bpt[NC + 1] = {3.00, 6.00};    
    grtyp = pttjt36;
    ptlim[grtyp][0]=bpt[0]; ptlim[grtyp][1]=bpt[1]; 
    jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mptt, sigN, emptt, esigN);
    grr[grtyp][tSigmaN][tReal]->SetMarkerStyle(20);
    grr[grtyp][tJt][tReal]    = jtkt->GetJTy_ptt();
    grr[grtyp][tJt][tReal]->SetMarkerStyle(20);
    
    //piv03_kt
    const int NC =  6 ;
    double sigN  [NC]= {0.235212, 0.217575, 0.200074, 0.192325, 0.199454, 0.172846}; 
    double esigN [NC]= {0.0011917, 0.00188785, 0.00294503, 0.0044955, 0.00516161, 0.00787827}; 
    double sigF  [NC]= {0.515988, 0.456308, 0.402666, 0.384932, 0.332543, 0.261762}; 
    double esigF [NC]= {0.00539758, 0.00737429, 0.00989722, 0.0134954, 0.0122275, 0.0137989}; 
    double pout  [NC]= { 1.17121, 1.0298, 0.912181, 0.886833, 0.771245, 0.621531}; 
    double epout [NC]= {0.0162147, 0.0198087, 0.0240964, 0.0337279, 0.029248, 0.0339336}; 
    double YN  [NC]= {0.0581068, 0.0827414, 0.110073, 0.144111, 0.178117, 0.247617}; 
    double eYN [NC]= {0.000298434, 0.000703738, 0.00158538, 0.00329689, 0.0049161, 0.0111272}; 
    double YF  [NC]= {0.0445665, 0.0648264, 0.0891338, 0.122702, 0.165568, 0.236338}; 
    double eYF [NC]= {0.000574708, 0.00113778, 0.00222146, 0.00436133, 0.00601838, 0.0130114}; 
    double mptt  [NC]= { 3.3613, 4.37786, 5.39766, 6.41083, 7.7178, 10.3154}; 
    double mpta  [NC]= { 2.13976, 2.19094, 2.23755, 2.27325, 2.30173, 2.37748}; 
    double mNpta  [NC]= { 2.16634, 2.21295, 2.26033, 2.30687, 2.32726, 2.37382}; 
    double mFpta  [NC]= { 2.09927, 2.15807, 2.20453, 2.22759, 2.26901, 2.38224}; 
    double SoverB  [NC]= { 2.29017, 3.44142, 4.93059, 6.45895, 10.3384, 13.7838}; 
    double bpt[NC+1] = {1.50, 4.00};
    grtyp = ptt1540;
    ptlim[grtyp][0]=bpt[0]; ptlim[grtyp][1]=bpt[1]; 
    jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mptt, sigN, emptt, esigN); grr[grtyp][tSigmaN]->SetMarkerStyle(20);
    grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mptt, sigF, emptt, esigF); grr[grtyp][tSigmaF]->SetMarkerStyle(21);
    grr[grtyp][tJt][tReal]     = jtkt->GetJTy_ptt(); grr[grtyp][tJt]->SetMarkerStyle(20);
    grr[grtyp][tPout][tReal] = new TGraphErrors(NC, mptt, pout, emptt, epout); grr[grtyp][tPout]->SetMarkerStyle(21);
    grr[grtyp][tZKT][tReal]  = jtkt->GetKTyz_ptt(); grr[grtyp][tZKT]->SetMarkerStyle(21);
    
    for(int i=0;i<NC;i++){
      grr[grtyp][tZKT][tReal]->GetPoint(i,x, sys[0][i]);
      sys[1][i] = sys[0][i]*0.05;
      sys[2][i] = sys[0][i]*0.05;
      emptt[i]=0.2;
    }
    grr[grtyp][tZKT][tSysL] = new TGraphErrors(NC, mptt, sys[0], emptt, sys[1]);
    grr[grtyp][tZKT][tSysL]->SetFillColor(41);

    const int NC = 5 ;
    double sigN [NC] = { 0.21129, 0.194303, 0.175706, 0.179157, 0.191894};
    double esigN [NC] = {0.00130761, 0.00211195, 0.00308274, 0.00378276, 0.00684469};
    double sigF [NC] = {0.487505, 0.423028, 0.388436, 0.369033, 0.270805};
    double esigF [NC] = {0.00645191, 0.00794819, 0.0113457, 0.0132413, 0.0131555};
    double pout [NC] = { 1.32657, 1.13326, 1.0547, 0.997665, 0.744404};
    double epout [NC]= {0.0216152, 0.0246219, 0.0334601, 0.0357628, 0.0362193}; 
    double YN [NC] = {0.0301311, 0.0453028, 0.0637024, 0.0923693, 0.127404 };
    double eYN [NC] = {0.000190391, 0.000489064, 0.00113445, 0.00208535, 0.00533067};
    double YF [NC] = {0.0208344, 0.0327481, 0.0474368, 0.07707, 0.136064};
    double eYF [NC] = {0.000304273, 0.000639679, 0.00140997, 0.00256304, 0.00659996};
    double mptt [NC] = { 3.3613, 4.37786, 5.39766, 6.68388, 9.12797};
    double mpta [NC] = { 2.61375, 2.64392, 2.69149, 2.70511, 2.74831};
    double mNpta [NC] = { 2.62942, 2.65484, 2.70593, 2.73088, 2.7421};
    double mFpta [NC] = { 2.58714, 2.62616, 2.66913, 2.66758, 2.75627};
    double SoverB [NC] = { 3.24186, 5.10616, 7.38133, 12.2124, 24.3161};
    double bpt[2] = {2.00, 4.00};
    grtyp = ptt24;
    ptlim[grtyp][0]=bpt[0]; ptlim[grtyp][1]=bpt[1]; 
    jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mptt, sigN, emptt, esigN); grr[grtyp][tSigmaN]->SetMarkerStyle(20);
    grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mptt, sigF, emptt, esigF); grr[grtyp][tSigmaF]->SetMarkerStyle(21);
    grr[grtyp][tJt][tReal]     = jtkt->GetJTy_ptt(); grr[grtyp][tJt]->SetMarkerStyle(20);
    grr[grtyp][tPout][tReal] = new TGraphErrors(NC, mptt, pout, emptt, epout); grr[grtyp][tPout]->SetMarkerStyle(21);
    grr[grtyp][tZKT][tReal]  = jtkt->GetKTyz_ptt(); grr[grtyp][tZKT]->SetMarkerStyle(21);
    
    for(int i=0;i<NC;i++){
      grr[grtyp][tZKT][tReal]->GetPoint(i,x, sys[0][i]);
      sys[1][i] = sys[0][i]*0.05;
      sys[2][i] = sys[0][i]*0.05;
      emptt[i]=0.2;
    }
    grr[grtyp][tZKT][tSysL] = new TGraphErrors(NC, mptt, sys[0], emptt, sys[1]);
    grr[grtyp][tZKT][tSysL]->SetFillColor(41);

    const int NC = 5 ;
    double sigN [NC] = {0.205404, 0.188575, 0.171199, 0.173097, 0.179451};
    double esigN [NC] = {0.00117361, 0.00183119, 0.00277372, 0.00348642, 0.00655649};
    double sigF [NC] = {0.844553, 0.499883, 0.388455, 0.355308, 0.292405};
    double esigF [NC] = {0.00812949, 0.00723097, 0.0106277, 0.0115472, 0.0144891};
    double pout [NC] = { 1.36988, 1.18932, 1.14229, 1.04456, 0.813267};
    double epout [NC] = {0.0209125, 0.0242663, 0.0326486, 0.0341407, 0.0405123};
    double YN [NC] = {0.0326686, 0.0500063, 0.0724826, 0.105217, 0.139015};
    double eYN [NC] = {0.000190477, 0.000500614, 0.00117138, 0.00220174, 0.00620163};
    double YF [NC] = {0.0224514, 0.0360154, 0.0537579, 0.0877465, 0.132193};
    double eYF [NC] = {0.000306431, 0.000660704, 0.00140119, 0.00245922, 0.00694972};
    double mptt [NC] = { 3.3613, 4.37786, 5.39766, 6.68388, 8.7485};
    double mpta [NC] = { 2.74542, 2.80071, 2.87061, 2.91797, 2.98446};
    double mNpta [NC] = { 2.76695, 2.82195, 2.89559, 2.94842, 2.98636};
    double mFpta [NC] = { 2.70842, 2.76549, 2.83109, 2.87303, 2.98201};
    double SoverB [NC] = { 33.6134, 12.6676, 7.94347, 13.2734, 25.3648};
    double bpt[2] = {2.00, 5.00};
    grtyp = ptt25;
    ptlim[grtyp][0]=bpt[0]; ptlim[grtyp][1]=bpt[1]; 
    jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mptt, sigN, emptt, esigN); grr[grtyp][tSigmaN]->SetMarkerStyle(20);
    grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mptt, sigF, emptt, esigF); grr[grtyp][tSigmaF]->SetMarkerStyle(21);
    grr[grtyp][tJt][tReal]     = jtkt->GetJTy_ptt(); grr[grtyp][tJt]->SetMarkerStyle(20);
    grr[grtyp][tPout][tReal] = new TGraphErrors(NC, mptt, pout, emptt, epout); grr[grtyp][tPout]->SetMarkerStyle(21);
    grr[grtyp][tZKT][tReal]  = jtkt->GetKTyz_ptt(); grr[grtyp][tZKT]->SetMarkerStyle(21);
    
    for(int i=0;i<NC;i++){
      grr[grtyp][tZKT][tReal]->GetPoint(i,x, sys[0][i]);
      sys[1][i] = sys[0][i]*0.05;
      sys[2][i] = sys[0][i]*0.05;
      emptt[i]=0.2;
    }
    grr[grtyp][tZKT][tSysL] = new TGraphErrors(NC, mptt, sys[0], emptt, sys[1]);
    grr[grtyp][tZKT][tSysL]->SetFillColor(41);

    const int NC = 6 ;
    double sigN [NC] = {0.214924, 0.195291, 0.177628, 0.185717, 0.183041, 0.162512};
    double esigN [NC] = {0.00110908, 0.00179689, 0.00271126, 0.00427842, 0.00627432, 0.00599166};
    double sigF [NC] = {0.489397, 0.434137, 0.398009, 0.361572, 0.328464, 0.279952};
    double esigF [NC] = {0.00554644, 0.00751139, 0.00999885, 0.0131291, 0.014668, 0.0160401};
    double pout [NC] = { 1.29789, 1.14803, 1.07321, 1.0022, 0.906725, 0.782218};
    double epout [NC] = {0.018266, 0.0222768, 0.0294446, 0.0380342, 0.0422164, 0.0445693};
    double YN [NC] = {0.0418667, 0.0628001, 0.0883947, 0.114525, 0.142776, 0.168431};
    double eYN [NC] = {0.000223973, 0.000577659, 0.00132515, 0.00261378, 0.00493203, 0.00711082};
    double YF [NC] = {0.0297372, 0.0467976, 0.0663763, 0.0994319, 0.125709, 0.141652};
    double eYF [NC] = {0.000372411, 0.000811083, 0.00166045, 0.00317097, 0.00578639, 0.00720541};
    double mptt [NC] = { 3.36129, 4.37784, 5.39766, 6.41086, 7.42285, 8.74846};
    double mpta [NC] = { 2.53947, 2.59908, 2.67533, 2.73355, 2.75575, 2.80777};
    double mNpta [NC] = { 2.56646, 2.62719, 2.70249, 2.76857, 2.77533, 2.8204};
    double mFpta [NC] = { 2.49507, 2.5543, 2.63324, 2.6825, 2.72735, 2.79196};
    double SoverB [NC] = { 3.06138, 4.88066, 7.35541, 14.503, 15.6307, 16.0835};
    double bpt[NC + 1] = {1.80, 5.00};
    grtyp = ptt1850;
    ptlim[grtyp][0]=bpt[0]; ptlim[grtyp][1]=bpt[1]; 
    jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mptt, sigN, emptt, esigN); grr[grtyp][tSigmaN]->SetMarkerStyle(20);
    grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mptt, sigF, emptt, esigF); grr[grtyp][tSigmaF]->SetMarkerStyle(21);
    grr[grtyp][tJt][tReal]     = jtkt->GetJTy_ptt(); grr[grtyp][tJt]->SetMarkerStyle(20);
    grr[grtyp][tPout][tReal] = new TGraphErrors(NC, mptt, pout, emptt, epout); grr[grtyp][tPout]->SetMarkerStyle(21);
    grr[grtyp][tZKT][tReal]  = jtkt->GetKTyz_ptt(); grr[grtyp][tZKT]->SetMarkerStyle(21);
    
    for(int i=0;i<NC;i++){
      grr[grtyp][tZKT][tReal]->GetPoint(i,x, sys[0][i]);
      sys[1][i] = sys[0][i]*0.05;
      sys[2][i] = sys[0][i]*0.05;
      emptt[i]=0.2;
    }
    grr[grtyp][tZKT][tSysL] = new TGraphErrors(NC, mptt, sys[0], emptt, sys[1]);
    grr[grtyp][tZKT][tSysL]->SetFillColor(41);

    const int NC = 6 ;
    double sigN [NC] = {0.214967, 0.195409, 0.177925, 0.184939, 0.183284, 0.166631};
    double esigN [NC] = {0.00110674, 0.00179394, 0.00271635, 0.00430891, 0.00623973, 0.00510273};
    double sigF [NC] = {0.489533, 0.434251, 0.397221, 0.360226, 0.324241, 0.240508};
    double esigF [NC] = {0.00553681, 0.00750027, 0.00990232, 0.013018, 0.0145185, 0.0135813};
    double pout [NC] = { 1.29838, 1.14834, 1.07024, 0.997914, 0.897093, 0.682666};
    double epout [NC] = {0.0182349, 0.0222621, 0.0291634, 0.0377463, 0.0415638, 0.0382978};
    double YN [NC] = {0.0418557, 0.0627524, 0.088317, 0.113802, 0.144431, 0.171918};
    double eYN [NC] = {0.000223499, 0.000576098, 0.00132388, 0.00260028, 0.00499207, 0.00605132};
    double YF [NC] = {0.0297628, 0.046769, 0.0666517, 0.0997213, 0.12523, 0.158343};
    double eYF [NC] = {0.000371816, 0.000811121, 0.0016531, 0.0031673, 0.00571722, 0.0066716};
    double mptt [NC] = { 3.3613, 4.37786, 5.39766, 6.41083, 7.42298, 9.31824};
    double mpta [NC] = { 2.5398, 2.59902, 2.67443, 2.73345, 2.75505, 2.82089};
    double mNpta [NC] = { 2.56693, 2.6269, 2.70091, 2.76893, 2.77489, 2.82143};
    double mFpta [NC] = { 2.49515, 2.55465, 2.63341, 2.6818, 2.72634, 2.82022};
    double SoverB [NC] = { 3.06584, 4.87483, 7.43273, 14.5883, 15.7605, 15.9191};
    double bpt[NC + 1] = {1.80, 5.00};
    grtyp = ptt1851;
    ptlim[grtyp][0]=bpt[0]; ptlim[grtyp][1]=bpt[1]; 
    jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mptt, sigN, emptt, esigN); grr[grtyp][tSigmaN]->SetMarkerStyle(20);
    grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mptt, sigF, emptt, esigF); grr[grtyp][tSigmaF]->SetMarkerStyle(21);
    grr[grtyp][tJt][tReal]     = jtkt->GetJTy_ptt(); grr[grtyp][tJt]->SetMarkerStyle(20);
    grr[grtyp][tPout][tReal] = new TGraphErrors(NC, mptt, pout, emptt, epout); grr[grtyp][tPout]->SetMarkerStyle(21);
    grr[grtyp][tZKT][tReal]  = jtkt->GetKTyz_ptt(); grr[grtyp][tZKT]->SetMarkerStyle(21);
    
    for(int i=0;i<NC;i++){
      grr[grtyp][tZKT][tReal]->GetPoint(i,x, sys[0][i]);
      sys[1][i] = sys[0][i]*0.05;
      sys[2][i] = sys[0][i]*0.05;
      emptt[i]=0.2;
    }
    grr[grtyp][tZKT][tSysL] = new TGraphErrors(NC, mptt, sys[0], emptt, sys[1]);
    grr[grtyp][tZKT][tSysL]->SetFillColor(41);

    const int NC = 6 ;
    double sigN [NC] = {0.215082, 0.196028, 0.179104, 0.188499, 0.186395, 0.178329};
    double esigN [NC] = {0.00112678, 0.00182231, 0.00277082, 0.00455375, 0.00659044, 0.00738062};
    double sigF [NC] = {0.490787, 0.434691, 0.397507, 0.356715, 0.334835, 0.25976};
    double esigF [NC] = {0.00578143, 0.00749725, 0.0099102, 0.0128286, 0.0152801, 0.0204667};
    double pout [NC] = { 1.30296, 1.14926, 1.07103, 0.987582, 0.924561, 0.736178};
    double epout [NC] = {0.0185704, 0.0227089, 0.0293327, 0.0376041, 0.0432348, 0.0551019};
    double YN [NC] = {0.0419453, 0.0630872, 0.0888439, 0.113814, 0.139815, 0.160708};
    double eYN [NC] = {0.000227927, 0.000588847, 0.00135243, 0.00266354, 0.00486525, 0.0085144};
    double YF [NC] = {0.0299164, 0.0470265, 0.0671515, 0.100756, 0.127701, 0.125128};
    double eYF [NC] = {0.000379772, 0.000826457, 0.00169261, 0.00321948, 0.0058248, 0.0078393};
    double mptt [NC] = { 3.3612, 4.37785, 5.39776, 6.41104, 7.42264, 8.43726};
    double mpta [NC] = { 2.54012, 2.59811, 2.67501, 2.7358, 2.75397, 2.8272};
    double mNpta [NC] = { 2.56731, 2.62582, 2.70089, 2.7711, 2.77155, 2.84534};
    double mFpta [NC] = { 2.49544, 2.55412, 2.6348, 2.6844, 2.72858, 2.80483};
    double SoverB [NC] = { 3.08752, 4.93624, 7.49493, 15.1579, 16.1433, 17.0543};
    double bpt[NC + 1] = {1.80, 5.00};
    grtyp = ptt1852;
    ptlim[grtyp][0]=bpt[0]; ptlim[grtyp][1]=bpt[1]; 
    jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mptt, sigN, emptt, esigN); grr[grtyp][tSigmaN]->SetMarkerStyle(20);
    grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mptt, sigF, emptt, esigF); grr[grtyp][tSigmaF]->SetMarkerStyle(21);
    grr[grtyp][tJt][tReal]     = jtkt->GetJTy_ptt(); grr[grtyp][tJt]->SetMarkerStyle(20);
    grr[grtyp][tPout][tReal] = new TGraphErrors(NC, mptt, pout, emptt, epout); grr[grtyp][tPout]->SetMarkerStyle(21);
    grr[grtyp][tZKT][tReal]  = jtkt->GetKTyz_ptt(); grr[grtyp][tZKT]->SetMarkerStyle(21);
    
    for(int i=0;i<NC;i++){
      grr[grtyp][tZKT][tReal]->GetPoint(i,x, sys[0][i]);
      sys[1][i] = sys[0][i]*0.05;
      sys[2][i] = sys[0][i]*0.05;
      emptt[i]=0.2;
    }
    grr[grtyp][tZKT][tSysL] = new TGraphErrors(NC, mptt, sys[0], emptt, sys[1]);
    grr[grtyp][tZKT][tSysL]->SetFillColor(41);


    const int NC = 6 ;
    double sigN [NC] = {0.245751, 0.229378, 0.212863, 0.204513, 0.216378, 0.209021};
    double esigN [NC] = {0.00135102, 0.00226807, 0.00357528, 0.00556264, 0.00819616, 0.00714859};
    double sigF [NC] = {0.533841, 0.471472, 0.428381, 0.409132, 0.366702, 0.324337};
    double esigF [NC] = {0.00608982, 0.0084974, 0.0121991, 0.0155371, 0.0175264, 0.0145955};
    double pout [NC] = { 1.14464, 0.988007, 0.888612, 0.853214, 0.765476, 0.674255};
    double epout [NC] = {0.0179899, 0.0216376, 0.0282725, 0.0356674, 0.0393086, 0.0310634};
    double YN [NC] = {0.0512984, 0.0720957, 0.0937819, 0.119826, 0.145356, 0.177219};
    double eYN [NC] = {0.000291396, 0.000711352, 0.00152143, 0.00304825, 0.00540891, 0.0066234};
    double YF [NC] = {0.0413653, 0.0583387, 0.0777228, 0.10871, 0.126987, 0.173625};
    double eYF [NC] = {0.000598348, 0.00115809, 0.00220069, 0.00418322, 0.00636195, 0.00768409};
    double mptt [NC] = { 3.36121, 4.37787, 5.39769, 6.41107, 7.42268, 9.12764};
    double mpta [NC] = { 1.99241, 2.02067, 2.03505, 2.04876, 2.06178, 2.08436};
    double mNpta [NC] = { 2.00694, 2.03226, 2.04688, 2.06518, 2.09038, 2.08941};
    double mFpta [NC] = { 1.9709, 2.00381, 2.01822, 2.02736, 2.02545, 2.07796};
    double SoverB [NC] = { 2.18721, 3.28157, 4.68157, 6.63627, 10.312, 19.6147};
    double bpt[NC + 1] = {1.50, 3.00};
    grtyp = ptt1530;
    ptlim[grtyp][0]=bpt[0]; ptlim[grtyp][1]=bpt[1]; 
    jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mptt, sigN, emptt, esigN); grr[grtyp][tSigmaN]->SetMarkerStyle(20);
    grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mptt, sigF, emptt, esigF); grr[grtyp][tSigmaF]->SetMarkerStyle(21);
    grr[grtyp][tJt][tReal]     = jtkt->GetJTy_ptt(); grr[grtyp][tJt]->SetMarkerStyle(20);
    grr[grtyp][tPout][tReal] = new TGraphErrors(NC, mptt, pout, emptt, epout); grr[grtyp][tPout]->SetMarkerStyle(21);
    grr[grtyp][tZKT][tReal]  = jtkt->GetKTyz_ptt(); grr[grtyp][tZKT]->SetMarkerStyle(21);
    
    for(int i=0;i<NC;i++){
      grr[grtyp][tZKT][tReal]->GetPoint(i,x, sys[0][i]);
      sys[1][i] = sys[0][i]*0.05;
      sys[2][i] = sys[0][i]*0.05;
      emptt[i]=0.2;
    }
    grr[grtyp][tZKT][tSysL] = new TGraphErrors(NC, mptt, sys[0], emptt, sys[1]);
    grr[grtyp][tZKT][tSysL]->SetFillColor(41);

    
  }

  //--------------------------------
  void GetPtaDist() {
  //--------------------------------
    //piv03_pta
    const int NC =  5 ;
    double sigN  [NC]= { 0.26983, 0.231358, 0.206467, 0.180784, 0.148719}; 
    double esigN [NC]= {0.00219424, 0.00220827, 0.00241551, 0.00188543, 0.00217478}; 
    double sigF  [NC]= {0.552870, 0.542573, 0.4648335, 0.423873, 0.383274}; 
    double esigF [NC]= {0.00891034, 0.0108586, 0.0111493, 0.0104516, 0.0100599}; 
    double pout  [NC]= { 1.03967, 1.29717, 1.28438, 1.46419, 1.77942}; 
    double epout [NC]= {0.0241539, 0.0367924, 0.0366588, 0.0378912, 0.0502167}; 
    double YN  [NC]= {0.0285678, 0.0149818, 0.00797014, 0.00714051, 0.00354311}; 
    double eYN [NC]= {0.000242226, 0.000145236, 9.49959e-05, 7.64715e-05, 5.47512e-05}; 
    double YF  [NC]= {0.024739, 0.0116362, 0.00534485, 0.00434568, 0.00209421}; 
    double eYF [NC]= {0.000538787, 0.000299971, 0.000138475, 0.000106704, 5.8308e-05}; 
    double mptt  [NC]= { 3.3613, 3.3613, 3.3613, 3.3613, 3.3613}; 
    double mpta  [NC]= { 1.71662, 2.21921, 2.7232, 3.4095, 4.68893}; 
    double SoverB  [NC]= { 2.75402, 2.64296, 3.66263, 5.33595, 9.17837}; 
    double bpt[NC+1] = {1.50, 2.00, 2.50, 3.00, 4.00, 6.00};
    grtyp = pta34;
    ptlim[grtyp][0]=3; ptlim[grtyp][1]=4; 
    jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mpta, sigN, empta, esigN);
    grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mpta, sigF, empta, esigF);
    grr[grtyp][tPout][tReal]    = new TGraphErrors(NC, mpta, pout, empta, epout);
    grr[grtyp][tJt][tReal]      = jtkt->GetJTy_pta();
    grr[grtyp][tZKT][tReal]     = jtkt->GetKTyz_pta();
    int marker = 21;
    grr[grtyp][tSigmaN][tReal]->SetMarkerStyle(marker);
    grr[grtyp][tSigmaF][tReal]->SetMarkerStyle(marker);
    grr[grtyp][tPout][tReal]  ->SetMarkerStyle(marker);
    grr[grtyp][tJt][tReal]    ->SetMarkerStyle(marker);
    grr[grtyp][tZKT][tReal]   ->SetMarkerStyle(marker);

    const int NC =  5 ;
    double sigN  [NC]= {0.254752, 0.214598, 0.188864, 0.169248, 0.152667}; 
    double esigN [NC]= {0.00378349, 0.00371848, 0.00383629, 0.00308092, 0.00367133}; 
    double sigF  [NC]= {0.517601, 0.46193, 0.432401, 0.40559, 0.333183}; 
    double esigF [NC]= {0.0144318, 0.0126335, 0.0146411, 0.0134332, 0.012747}; 
    double pout  [NC]= {0.945121, 1.197298, 1.20376, 1.39132, 1.53995}; 
    double epout [NC]= {0.034468, 0.0328399, 0.0469905, 0.0487018, 0.0616629}; 
    double YN  [NC]= {0.0384788, 0.0216408, 0.0120397, 0.011249, 0.00642419}; 
    double eYN [NC]= {0.000582038, 0.000369237, 0.000247783, 0.000205634, 0.000161075}; 
    double YF  [NC]= {0.0332892, 0.016163, 0.00948276, 0.00820146, 0.0043183}; 
    double eYF [NC]= {0.00113026, 0.000517991, 0.000334728, 0.000263832, 0.000170131}; 
    double mptt  [NC]= { 4.37786, 4.37786, 4.37786, 4.37786, 4.37786}; 
    double mpta  [NC]= { 1.72256, 2.22357, 2.72825, 3.42195, 4.68967}; 
    double mNpta  [NC]= { 1.72333, 2.22626, 2.72816, 3.422, 4.67644}; 
    double mFpta  [NC]= { 1.72152, 2.21935, 2.7284, 3.42187, 4.71547}; 
    double SoverB  [NC]= { 2.4799, 4.07931, 6.17897, 9.76308, 20.1161}; 
    double bpt[NC+1] = {1.50, 2.00, 2.50, 3.00, 4.00, 6.00};
    grtyp = pta45;
    ptlim[grtyp][0]=4; ptlim[grtyp][1]=5; 
    jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mpta, sigN, empta, esigN);
    grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mpta, sigF, empta, esigF);
    grr[grtyp][tPout][tReal]    = new TGraphErrors(NC, mpta, pout, empta, epout);
    grr[grtyp][tJt][tReal]      = jtkt->GetJTy_pta();
    grr[grtyp][tZKT][tReal]     = jtkt->GetKTyz_pta();
    int marker = 25;
    grr[grtyp][tSigmaN][tReal]->SetMarkerStyle(marker);
    grr[grtyp][tSigmaF][tReal]->SetMarkerStyle(marker);
    grr[grtyp][tPout][tReal]  ->SetMarkerStyle(marker);
    grr[grtyp][tJt][tReal]    ->SetMarkerStyle(marker);
    grr[grtyp][tZKT][tReal]   ->SetMarkerStyle(marker);

    const int NC =  5 ;
    double sigN  [NC]= { 0.22934, 0.187448, 0.176848, 0.154741, 0.134747}; 
    double esigN [NC]= {0.00575869, 0.00569099, 0.00610231, 0.00451015, 0.00395705}; 
    double sigF  [NC]= {0.417101, 0.423744, 0.31629, 0.298417, 0.339967}; 
    double esigF [NC]= {0.015356, 0.0183441, 0.0176440, 0.0116635, 0.0147683}; 
    double pout  [NC]= {0.724789, 0.934351, 0.865934, 1.0525, 1.60319}; 
    double epout [NC]= {0.0307428, 0.045972, 0.0487456, 0.0426619, 0.0765677}; 
    double YN  [NC]= {0.0585078, 0.0337564, 0.0207421, 0.0223372, 0.0164396}; 
    double eYN [NC]= {0.00147732, 0.000963678, 0.000719428, 0.000685176, 0.000588772}; 
    double YF  [NC]= {0.0522264, 0.0317654, 0.0163221, 0.0203796, 0.0115058}; 
    double eYF [NC]= {0.00208839, 0.00144154, 0.000857869, 0.000834862, 0.000580878}; 
    double mptt  [NC]= { 6.24682, 6.24682, 6.24682, 6.24682, 6.24682}; 
    double mpta  [NC]= { 1.71982, 2.22285, 2.72548, 3.43889, 4.72026}; 
    double mNpta  [NC]= { 1.72392, 2.229, 2.72745, 3.44497, 4.70157}; 
    double mFpta  [NC]= { 1.71466, 2.21471, 2.72227, 3.42939, 4.75205}; 
    double SoverB  [NC]= { 4.72399, 8.81738, 303.956, 28.0712, 28.3373}; 
    double bpt[NC+1] = {1.50, 2.00, 2.50, 3.00, 4.00, 6.00};
    grtyp = pta5580;
    ptlim[grtyp][0]=5.5; ptlim[grtyp][1]=8; 
    jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mpta, sigN, empta, esigN);
    grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mpta, sigF, empta, esigF);
    grr[grtyp][tPout][tReal]    = new TGraphErrors(NC, mpta, pout, empta, epout);
    grr[grtyp][tJt][tReal]      = jtkt->GetJTy_pta();
    grr[grtyp][tZKT][tReal]     = jtkt->GetKTyz_pta();
    int marker = 22;
    grr[grtyp][tSigmaN][tReal]->SetMarkerStyle(marker);
    grr[grtyp][tSigmaF][tReal]->SetMarkerStyle(marker);
    grr[grtyp][tPout][tReal]  ->SetMarkerStyle(marker);
    grr[grtyp][tJt][tReal]    ->SetMarkerStyle(marker);
    grr[grtyp][tZKT][tReal]   ->SetMarkerStyle(marker);

    const int NC =  5 ;
    double sigN  [NC]= {0.225402, 0.19559, 0.173119, 0.15193, 0.110517}; 
    double esigN [NC]= {0.010695, 0.00912837, 0.0104797, 0.00868875, 0.00547745}; 
    double sigF  [NC]= {0.349627, 0.364686, 0.313225, 0.275312, 0.227544}; 
    double esigF [NC]= {0.0238558, 0.0255262, 0.0205065, 0.0192433, 0.0137075}; 
    double pout  [NC]= {0.615771, 0.809745, 0.85063, 0.941835, 1.07443}; 
    double epout [NC]= {0.044723, 0.0585744, 0.0561736, 0.0677997, 0.066315}; 
    double YN  [NC]= {0.0718054, 0.0522536, 0.0257537, 0.0272917, 0.0285725}; 
    double eYN [NC]= {0.00349727, 0.00269477, 0.00169076, 0.00169394, 0.00180596}; 
    double YF  [NC]= {0.0668274, 0.047942, 0.0298742, 0.030101, 0.0245675}; 
    double eYF [NC]= {0.0044592, 0.0034413, 0.00225357, 0.00220832, 0.00199752}; 
    double mptt  [NC]= { 7.91816, 7.91816, 7.91816, 7.91816, 7.91816}; 
    double mpta  [NC]= { 1.72616, 2.22843, 2.72858, 3.43881, 4.7724}; 
    double SoverB  [NC]= { 6.24245, 11.8844, 25.0802, 24.0522, 54.2041}; 
    double bpt[NC+1] = {1.50, 2.00, 2.50, 3.00, 4.00, 6.00};
    grtyp = pta710;
    ptlim[grtyp][0]=7; ptlim[grtyp][1]=10; 
    jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mpta, sigN, empta, esigN);
    grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mpta, sigF, empta, esigF);
    grr[grtyp][tPout][tReal]    = new TGraphErrors(NC, mpta, pout, empta, epout);
    grr[grtyp][tJt][tReal]      = jtkt->GetJTy_pta();
    grr[grtyp][tZKT][tReal]     = jtkt->GetKTyz_pta();
    int marker = 26;
    grr[grtyp][tSigmaN][tReal]->SetMarkerStyle(marker);
    grr[grtyp][tSigmaF][tReal]->SetMarkerStyle(marker);
    grr[grtyp][tPout][tReal]  ->SetMarkerStyle(marker);
    grr[grtyp][tJt][tReal]    ->SetMarkerStyle(marker);
    grr[grtyp][tZKT][tReal]   ->SetMarkerStyle(marker);
  }


  //--------------------------------
  void GetMCDist() {
  //--------------------------------
    //rcf_kt3.0I0F0V1.0
    const int NC =  9 ;
    double sigN  [NC]= {0.279601, 0.256517, 0.240819, 0.230552, 0.217021, 0.204037, 0.186437, 0.192801, 0.175649}; 
    double esigN [NC]= {0.000661104, 0.000798326, 0.000999109, 0.00131703, 0.00139006, 0.00249416, 0.00382185, 0.00727508, 0.00829801}; 
    double sigF  [NC]= {0.620922, 0.58732, 0.554556, 0.53633, 0.506112, 0.479666, 0.429235, 0.384783, 0.311632}; 
    double esigF [NC]= {0.007952, 0.00979387, 0.0116786, 0.0128858, 0.0113308, 0.0150914, 0.0176987, 0.0237479, 0.0200451}; 
    double YN  [NC]= {0.854934, 0.886377, 0.873904, 0.818234, 0.737375, 0.648546, 0.580041, 0.638938, 0.828225}; 
    double eYN [NC]= {0.00227566, 0.00313039, 0.00408838, 0.00527788, 0.00526524, 0.00938566, 0.0139097, 0.0254754, 0.0424444}; 
    double YF  [NC]= { 0.14268, 0.139446, 0.141678, 0.14812, 0.171071, 0.223552, 0.299506, 0.34024, 0.534785}; 
    double eYF [NC]= {0.00201385, 0.00249648, 0.00310314, 0.00378087, 0.0039541, 0.00702162, 0.0120584, 0.0217349, 0.036637}; 
    double mptt  [NC]= { 2.22004, 2.72404, 3.22586, 3.72644, 4.41066, 5.4131, 6.40833, 7.41237, 8.69772}; 
    double mpta  [NC]= { 2.11777, 2.11323, 2.11407, 2.12122, 2.11927, 2.1404, 2.18127, 2.20048, 2.3429}; 
    double mNpta  [NC]= { 2.14105, 2.13754, 2.13913, 2.14669, 2.14609, 2.1632, 2.20573, 2.2231, 2.27119}; 
    double mFpta  [NC]= { 2.0408, 2.02732, 2.02631, 2.03567, 2.04327, 2.09009, 2.14071, 2.16322, 2.44562}; 
//     for(int j=0;j<NC;j++) {YN[j]*= 1.0; eYN[j]*= 1.0;}
//     for(int j=0;j<NC;j++) {YF[j]*= 2.0; eYF[j]*= 2.0;}
//     jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    enGrPT grtyp = MCkt3_ptt;
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mptt, sigN, emptt, esigN);
    grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mptt, sigF, emptt, esigF);

//     pytYTrigN = new TGraphErrors(NC-1, mptt, YN, emptt, eYN);
//     pytYTrigF = new TGraphErrors(NC-1, mptt, YF, emptt, eYF);
//     pytYTrigN->SetLineWidth(2);  pytYTrigF->SetLineWidth(2);
//     pytTrigSigN = new TGraphErrors(NC, mptt, sigN, emptt, esigN);
//     pytTrigSigF = new TGraphErrors(NC, mptt, sigF, emptt, esigF);
//     //   pytTrigSigN->SetFillColor(1);  pytTrigSigF->SetFillColor(10);
//     pytTrigSigN->SetLineWidth(2);   pytTrigSigF->SetLineWidth(2);
    
    
    //rcf_kt2.0I0F0V1_trigg.root
    const int NC =  9 ;
    double sigN  [NC]= {0.296123, 0.270335, 0.251138, 0.23536, 0.215885, 0.195755, 0.199738,               0.173866, 0.16174}; 
    double esigN [NC]= {0.00106432, 0.0013774, 0.00184228, 0.00244328, 0.00256413, 0.0039821, 0.0063094,   0.008948, 0.00876596}; 
    double sigF  [NC]= {0.545048, 0.528385, 0.481404, 0.4437, 0.404859, 0.349816, 0.326035,                0.271085, 0.249243}; 
    double esigF [NC]= {0.00550298, 0.00695167, 0.00827349, 0.00910347, 0.00783672, 0.0108368, 0.0132557,  0.013767, 0.0213012}; 
    double mptt  [NC]= { 2.21606, 2.71929, 3.22129, 3.72091, 4.38965, 5.39127, 6.41301, 7.42143,                     8.71778}; 
    double mpta  [NC]= { 2.01974, 2.04502, 2.06938, 2.11165, 2.16768, 2.27635, 2.36079, 2.44171,                     2.46449}; 
    double mNpta  [NC]= { 2.02231, 2.04413, 2.06689, 2.10677, 2.16423, 2.25867, 2.33759, 2.39894,                    2.38207}; 
    double mFpta  [NC]= { 2.01389, 2.04702, 2.07448, 2.12072, 2.17314, 2.30058, 2.39083, 2.48812,                    2.57237}; 
    grtyp = MCkt2_ptt;
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mptt, sigN, emptt, esigN);
    grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mptt, sigF, emptt, esigF);
//     pytTrigSigN3 = new TGraphErrors(NC-1, mptt, sigN, emptt, esigN);
//     pytTrigSigF3 = new TGraphErrors(NC-1, mptt, sigF, emptt, esigF);
//     pytTrigSigN3->SetLineWidth(2);  pytTrigSigF3->SetLineWidth(2);    
//     jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
//     grr[tSigmaPttN] = new TGraphErrors(NC, mptt, sigN, emptt, esigN);
//     grr[tSigmaPttF] = new TGraphErrors(NC, mptt, sigF, emptt, esigF);

    
    //pythia rcf_kt2.0I0F0V1_trigg.root
    const int NC =  9 ;
    double pyptt [NC] = { 2.21606,  2.71929,  3.22129,  3.72091,  4.38965,  5.39127,  6.41301,  7.42143,  8.71778};
    double pyzkt [NC] = { 1.10358,  1.32113,  1.37313,  1.41096,  1.51186,   1.5366,  1.66693,  1.55607,  1.60923};
    double epyptt[NC] = {     0.1,      0.1,      0.1,      0.1,      0.1,      0.1,      0.1,      0.1,      0.1};
    double epyzkt[NC] = {0.0216568, 0.0311547, 0.0381869, 0.0464232, 0.0438894, 0.0721763, 0.107043, 0.145824, 0.238137};
    grtyp = MCkt2_ptt;
    grr[grtyp][tPout][tReal] = new TGraphErrors(NC-1, pyptt, pyzkt, epyptt, epyzkt);
    
    //pythia rcf_kt3.0I0F0V1.0
    const int NC =  8 ;
    double pyptt [NC] = { 2.22004,  2.72404,  3.22586,  3.72644,  4.41066,   5.4131,  6.40833,  7.41237};
    double pyzkt [NC] = { 1.41818,    1.611,  1.76383,  1.96637,    2.138,  2.46192,   2.5257,  2.58226};
    double epyptt[NC] = {     0.1,      0.1,      0.1,      0.1,      0.1,      0.1,      0.1,      0.1};
    double epyzkt[NC] = {0.0358274,  0.04902, 0.0628342, 0.0790482, 0.0736502, 0.110399, 0.133004, 0.215087};
    grtyp = MCkt3_ptt;
    grr[grtyp][tPout][tReal] = new TGraphErrors(NC, pyptt, pyzkt, epyptt, epyzkt);
  }

  //--------------------------------
  void Get_03_PttDist() {
  //--------------------------------

    const int NC =  6 ; //partially recovered and from pi03_NewNew_14pta50.roo
    double sigN  [NC]= {0.236589, 0.216757, 0.2109972, 0.195693, 0.184941,          0.099999}; 
    double esigN [NC]= {0.00358158, 0.00342673, 0.00587283, 0.00803991, 0.00998691, 0.099999};
    double sigF  [NC]= {   0.521,    0.464, 0.431,    0.396,     0.35,              0.99999};
    double esigF [NC]= {0.0158789, 0.0141421, 0.0207893, 0.0220255, 0.0244331,      0.99999}; 
    double pout  [NC]= { 1.18757, 1.03306, 1.06972, 0.899801, 0.645105,             0.99999}; 
    double epout [NC]= {0.0490876, 0.0374822, 0.056701, 0.0563907, 0.0593627,       0.99999};
    double mptt  [NC]= { 3.22152, 3.89, 4.9044, 5.91016, 7.23529,                   8.34193}; 
    double mpta  [NC]= { 2.10403, 2.15981, 2.21359, 2.27439, 2.41041,               2.49999}; 
    double jt[NC]  = {   0.587,    0.577,      0.6,    0.596,    0.597,    0.632};
    double ejt[NC] = {   0.009,    0.009,    0.017,    0.026,    0.038,    0.085};
    double bpt[NC+1] = {1.40, 5.00};
    grtyp = pttjt;
    ptlim[grtyp][0]=bpt[0]; ptlim[grtyp][1]=bpt[1]; 
    jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mptt, sigN, emptt, esigN);
    grr[grtyp][tSigmaN][tReal]->SetMarkerStyle(20);
    grr[grtyp][tJt][tReal]    = new TGraphErrors(NC, mptt, jt, emptt, ejt); //jtkt->GetJTy_ptt();
    grr[grtyp][tJt][tReal]->SetMarkerStyle(20);

    //run03
    const int NC = 6;
    double mptt[NC] = { 3.22152, 3.89,    4.9044,  5.91016, 7.23529, 8.34193 };
    double xzkt[NC] = { 1.62953, 1.66402, 1.88903, 2.06215, 2.17462, 2.52895 };
    double exzkt[NC]= {0.0828566, 0.0773639, 0.132192, 0.193457, 0.24573, 0.617122 };
    double bpt[NC+1] = {1.40, 5.00};
    grtyp = ptt1540;
    ptlim[grtyp][0]=bpt[0]; ptlim[grtyp][1]=bpt[1]; 
    grr[grtyp][tZKT][tReal] = new TGraphErrors(NC, mptt, xzkt, emptt, exzkt);  //jtkt->GetKTyz_ptt();
    grr[grtyp][tZKT][tReal]->SetMarkerStyle(20);
    //grr[grtyp][tPout][tReal] = new TGraphErrors(NC, mptt, pout, emptt, epout); grr[grtyp][tPout]->SetMarkerStyle(21);
    
    for(int i=0;i<NC;i++){
      grr[grtyp][tZKT][tReal]->GetPoint(i,x, sys[0][i]);
      sys[1][i] = sys[0][i]*0.05;
      sys[2][i] = sys[0][i]*0.05;
      emptt[i]=0.2;
    }
    grr[grtyp][tZKT][tSysL] = new TGraphErrors(NC, mptt, sys[0], emptt, sys[1]);
    grr[grtyp][tZKT][tSysL]->SetFillColor(41);

    //=== FINAL kT ptt from ppfin.C and above
    // the seam grtyp and ptlims
    const int NC=6;
    double mptt[NC]={ 3.22152, 3.89, 4.9044, 5.91016, 7.23529, 8.34193, };
    double mpta[NC]={ 2.10403, 2.15981, 2.21359, 2.27439, 2.41, 2.40717, };
  
    double kT[NC]={   2.78856, 2.56872, 2.66249, 2.73924, 2.82573, 3.11495, };
    double ekTl[NC]={ 0.132664, 0.105996, 0.168991, 0.19505, 0.254544, 0.601929, };
    double ekTh[NC]={ 0.147192, 0.106377, 0.138255, 0.200021, 0.235688, 0.491619, };
  
    double mz[NC]={  0.512864, 0.56048, 0.605206, 0.635991, 0.660047, 0.682522, };
    double emzl[NC]={0.00409504, 0.00202162, 0.0013198, -7.15624e-05, -0.00128994, -0.00566139, };
    double emzh[NC]={0.00491338, 0.00216728, 0.00119696, -0.000269158, -0.00222622, -0.00908892, };
  
    double xhq[NC]={  0.877438, 0.868518, 0.85076, 0.848901, 0.859358, 0.841775, };
    double exhql[NC]={ 0.00679598, 0.00772673, 0.0148326, 0.0176072, 0.0211194, 0.0499965, };
    double exhqh[NC]={0.00715344, 0.00770412, 0.0124174, 0.0190238, 0.020743, 0.0479045, };
    //-------------------------------------------------------------
    double syskth[NC]= { 3.13574, 2.90125, 3.00391, 3.08819, 3.15449, 3.45073, };
    //double ekTl[NC]= { 0.161454, 0.112493, 0.169952, 0.227228, 0.284492, 0.671501, };
    double sysmzh[NC]= { 0.459596, 0.500076, 0.542144, 0.572449, 0.597002, 0.623431, };
    //double emzl[NC]= { 0.00310951, 0.0006694, -0.000625784, -0.00218767, -0.00279109, -0.00761746, };
    double sysxhqh[NC]={ 0.881252, 0.875918, 0.860641, 0.858783, 0.868835, 0.851811, };
    //double exhql[NC]={ 0.00823916, 0.00759532, 0.0131869, 0.0173032, 0.0191454, 0.0455805, };    
    double sysktl[NC]= { 2.45496, 2.24549, 2.31433, 2.41893, 2.50641, 2.78778, };
    //double ekTl[NC]= { 0.117885, 0.0968849, 0.13264, 0.180775, 0.220924, 0.519263, };
    double sysmzl[NC]= { 0.575572, 0.630781, 0.676804, 0.705684, 0.727081, 0.744532, };
    //double emzl[NC]= { 0.00533346, 0.00333942, 0.00261934, 0.00145364, -0.000533472, -0.00513088, };  
    double sysxhql[NC]={ 0.866767, 0.85208, 0.831968, 0.827071, 0.841193, 0.823419, };
    //double exhql[NC]={ 0.00601329, 0.00772176, 0.0136357, 0.0202763, 0.0240169, 0.0555096, };

    double zero[100]={0};
    grr[grtyp][tKT][tReal] = new TGraphErrors(NC, mptt, kT, zero, ekTl);
    grr[grtyp][tKT][tReal]->SetMarkerStyle(20);
    grr[grtyp][tKT][tSysL] = new TGraphErrors(NC, mptt, sysktl, zero, zero);
    grr[grtyp][tKT][tSysH] = new TGraphErrors(NC, mptt, syskth, zero, zero);
    
    grr[grtyp][tZT][tReal] = new TGraphErrors(NC, mptt, mz, zero, emzl);
    grr[grtyp][tZT][tReal]->SetMarkerStyle(20);
    grr[grtyp][tZT][tSysL] = new TGraphErrors(NC, mptt, sysmzl, zero, zero);
    grr[grtyp][tZT][tSysH] = new TGraphErrors(NC, mptt, sysmzh, zero, zero);

    grr[grtyp][tXhq][tReal] = new TGraphErrors(NC, mptt, xhq, zero, exhql);
    grr[grtyp][tXhq][tReal]->SetMarkerStyle(22);
    grr[grtyp][tXhq][tSysL] = new TGraphErrors(NC, mptt, sysxhql, zero, zero);
    grr[grtyp][tXhq][tSysH] = new TGraphErrors(NC, mptt, sysxhqh, zero, zero);
  }

  //--------------------------------
  void Get_03_PtaDist() {
  //--------------------------------
    
  // 3 < pTtrigg < 4 pi0  pi03_3t4_new.root 
    const int NC =  6 ;
    double sigN  [NC]= {0.259379, 0.227474, 0.190342, 0.178368, 0.166084, 0.127265}; 
    double esigN [NC]= {0.00525097, 0.00531542, 0.00551501, 0.00612342, 0.00460053, 0.00723908}; 
    double sigF  [NC]= {0.538939, 0.527376, 0.438261, 0.456289, 0.454404, 0.327827}; 
    double esigF [NC]= {0.0222664, 0.0250665, 0.0302659, 0.0306966, 0.0255876, 0.026706}; 
    double pout  [NC]= {0.996454, 1.24423, 1.22196, 1.4959, 1.79262, 1.67537}; 
    double epout [NC]= {0.0559085, 0.0785499, 0.0948472, 0.105082, 0.115185, 0.14148}; 
    double mptt  [NC]= { 3.38892, 3.38892, 3.38892, 3.38892, 3.38892, 3.38892}; 
    double mpta  [NC]= { 1.71639, 2.22036, 2.72949, 3.22803, 3.92835, 5.03744}; 
    grtyp = pta34;
    ptlim[grtyp][0]=3; ptlim[grtyp][1]=4; 
    jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mpta, sigN, empta, esigN);
    grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mpta, sigF, empta, esigF);
    grr[grtyp][tPout][tReal]    = new TGraphErrors(NC, mpta, pout, empta, epout);
    grr[grtyp][tJt][tReal]      = jtkt->GetJTy_pta();
    grr[grtyp][tZKT][tReal]     = jtkt->GetKTyz_pta();
    int marker = 24;
    grr[grtyp][tSigmaN][tReal]->SetMarkerStyle(marker);
    grr[grtyp][tSigmaF][tReal]->SetMarkerStyle(marker);
    grr[grtyp][tPout][tReal]  ->SetMarkerStyle(marker);
    grr[grtyp][tJt][tReal]    ->SetMarkerStyle(marker);
    grr[grtyp][tZKT][tReal]   ->SetMarkerStyle(marker);
    for(int i=0;i<NC;i++){
      grr[grtyp][tZKT][tReal]->GetPoint(i,x, sys[0][i]);
      sys[1][i] = sys[0][i]*0.05;
      sys[2][i] = sys[0][i]*0.05;
      empta[i]=0.1;
    }
    grr[grtyp][tZKT][tSysL] = new TGraphErrors(NC, mpta, sys[0], empta, sys[1]);
    grr[grtyp][tZKT][tSysL]->SetFillColor(41);

    // 5 < pTtrigg < 10  pi03_5t10_new2.root 
    const int NC =  5 ;
    double sigN  [NC]= {0.242694, 0.234205, 0.163402, 0.187071, 0.127472}; 
    double esigN [NC]= {0.0134617, 0.012483, 0.0105629, 0.0120335, 0.00603555}; 
    double pout  [NC]= {0.959736, 1.09961, 1.08775, 1.28503, 1.26768}; 
    double epout [NC]= {0.102137, 0.102808, 0.110133, 0.135535, 0.209666}; 
    double mptt  [NC]= { 5.86593, 5.86593, 5.86593, 5.86593, 5.86593}; 
    double mpta  [NC]= { 1.84732, 2.23604, 2.73176, 3.4433, 4.65253}; 
    grtyp = pta510;
    ptlim[grtyp][0]=5; ptlim[grtyp][1]=10; 
    jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
    grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mpta, sigN, empta, esigN);
    grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mpta, sigF, empta, esigF);
    grr[grtyp][tPout][tReal]    = new TGraphErrors(NC, mpta, pout, empta, epout);
    grr[grtyp][tJt][tReal]      = jtkt->GetJTy_pta();
    grr[grtyp][tZKT][tReal]     = jtkt->GetKTyz_pta();
    int marker = 20;
    grr[grtyp][tSigmaN][tReal]->SetMarkerStyle(marker);
    grr[grtyp][tSigmaF][tReal]->SetMarkerStyle(marker);
    grr[grtyp][tPout][tReal]  ->SetMarkerStyle(marker);
    grr[grtyp][tJt][tReal]    ->SetMarkerStyle(marker);
    grr[grtyp][tZKT][tReal]   ->SetMarkerStyle(marker);
    for(int i=0;i<NC;i++){
      grr[grtyp][tZKT][tReal]->GetPoint(i,x, sys[0][i]);
      sys[1][i] = sys[0][i]*0.05;
      sys[2][i] = sys[0][i]*0.05;
      empta[i]=0.1;
    }
    grr[grtyp][tZKT][tSysL] = new TGraphErrors(NC, mpta, sys[0], empta, sys[1]);
    grr[grtyp][tZKT][tSysL]->SetFillColor(30);

    // ===== kT =============
    const int NC=6;
    double mptt[NC]={3.38892, 3.38892, 3.38892, 3.38892, 3.38892, 3.38892, };
    double mpta[NC]={1.71639, 2.22036, 2.72949, 3.22803, 3.92835, 5.03744, };
    
    double kT[NC]={2.65971, 2.93875, 2.57151, 2.92909, 3.19017, 2.68119, };
    double ekTl[NC]={0.188477, 0.215353, 0.232712, 0.233007, 0.226005, 0.250578, };
    double ekTh[NC]={0.171709, 0.230893, 0.237499, 0.245125, 0.241247, 0.237131, };
    
    double mz[NC]={0.539397, 0.51772, 0.506838, 0.487994, 0.467395, 0.411453, };
    double emzl[NC]={0.00817732, 0.0057948, 0.000945396, -0.00165861, -0.00463482, -0.00580352, };
    double emzh[NC]={0.00783731, 0.00731948, 0.00140203, -0.000698139, -0.00378033, -0.00703773, };

    double xhq[NC]={0.810886, 0.876485, 0.951451, 0.987946, 1.03728, 1.06104, };
    double exhql[NC]={0.0146864, 0.0108862, 0.00638619, 0.00147691, -0.00440618, -0.00915709, };
    double exhqh[NC]={0.0123396, 0.0107569, 0.00664045, 0.00159581, -0.00498729, -0.00934807, };
    //-------------------------------------------------------------
    double syskth[NC]={2.96346, 3.30025, 2.90937, 3.29392, 3.5752, 3.06352, };
    //double ekTl[NC]={0.199783, 0.246828, 0.246051, 0.241654, 0.230637, 0.27826, };
    double sysmzh[NC]={0.48597, 0.464223, 0.448644, 0.433275, 0.415064, 0.359334, };
    //double emzl[NC]={0.00603293, 0.00393503, -0.00173301, -0.00384298, -0.00615228, -0.00683699, };
    double sysxhqh[NC]={0.81708, 0.880514, 0.954824, 0.988818, 1.03432, 1.05605, };
    //double exhql[NC]={0.0154719, 0.0124152, 0.00611894, 0.00142423, -0.00418864, -0.00828728, };
    double sysktl[NC]={2.33745, 2.58582, 2.24895, 2.58073, 2.82881, 2.34298, };
    //double ekTl[NC]={0.156152, 0.191639, 0.214306, 0.217885, 0.214581, 0.229078, };
    double sysmzl[NC]={0.601715, 0.581278, 0.576375, 0.554486, 0.532216, 0.477152, };
    //double emzl[NC]={0.00910618, 0.00779648, 0.0038341, 0.000803405, -0.00276841, -0.00519443, };
    double sysxhql[NC]={0.798852, 0.865699, 0.943797, 0.986, 1.04393, 1.07517, };
    //double exhql[NC]={0.0123178, 0.00970025, 0.00656033, 0.00150673, -0.00462773, -0.0108602, };    
    grtyp = pta34;
    ptlim[grtyp][0]=3; ptlim[grtyp][1]=4; 
    double zero[100]={0};
    grr[grtyp][tKTa3][tReal] = new TGraphErrors(NC, mpta, kT, zero, ekTl);
    grr[grtyp][tKTa3][tReal]->SetMarkerStyle(20);
    grr[grtyp][tKTa3][tSysL] = new TGraphErrors(NC, mpta, sysktl, zero, zero);
    grr[grtyp][tKTa3][tSysH] = new TGraphErrors(NC, mpta, syskth, zero, zero);
    
    grr[grtyp][tZTa3][tReal] = new TGraphErrors(NC, mpta, mz, zero, emzl);
    grr[grtyp][tZTa3][tReal]->SetMarkerStyle(20);
    grr[grtyp][tZTa3][tSysL] = new TGraphErrors(NC, mpta, sysmzl, zero, zero);
    grr[grtyp][tZTa3][tSysH] = new TGraphErrors(NC, mpta, sysmzh, zero, zero);

    grr[grtyp][tXhqa3][tReal] = new TGraphErrors(NC, mpta, xhq, zero, exhql);
    grr[grtyp][tXhqa3][tReal]->SetMarkerStyle(22);
    grr[grtyp][tXhqa3][tSysL] = new TGraphErrors(NC, mpta, sysxhql, zero, zero);
    grr[grtyp][tXhqa3][tSysH] = new TGraphErrors(NC, mpta, sysxhqh, zero, zero);

    const int NC=5;
    double mptt[NC]={5.86593, 5.86593, 5.86593, 5.86593, 5.86593, };
    double mpta[NC]={1.84732, 2.23604, 2.73176, 3.4433, 4.65253, };
    
    double kT[NC]={3.08755, 3.19366, 3.03564, 3.03952, 2.64219, };
    double ekTl[NC]={0.304571, 0.296834, 0.302184, 0.384397, 0.560654, };
    double ekTh[NC]={0.27555, 0.262244, 0.278299, 0.344797, 0.53019, };
    
    double mz[NC]={0.657809, 0.636684, 0.609251, 0.571949, 0.519472, };
    double emzl[NC]={0.00383682, 0.000456943, -0.00193182, -0.00196791, 0.0052766, };
    double emzh[NC]={0.00601082, 0.00209165, -0.00260897, -0.00592326, -0.00106922, };
    
    double xhq[NC]={0.751108, 0.798761, 0.86846, 0.91933, 0.976078, };
    double exhql[NC]={0.0389635, 0.0304537, 0.020735, 0.0153872, 0.00788006, };
    double exhqh[NC]={0.0348043, 0.0281079, 0.0211909, 0.0154559, 0.00780096, };
    //-------------------------------------------------------------
    double syskth[NC]={3.43313, 3.56196, 3.41428, 3.43384, 3.03504, };
    //double ekTl[NC]={0.333775, 0.30953, 0.328588, 0.404667, 0.645655, };
    double sysmzh[NC]={0.599962, 0.578271, 0.547568, 0.509082, 0.453992, };
    //double emzl[NC]={0.000990045, -0.00273702, -0.0056658, -0.00616929, 0.00254761, };
    double sysxhqh[NC]={0.76447, 0.811669, 0.878767, 0.926212, 0.977463, };
    //double exhql[NC]={0.0395657, 0.0291846, 0.0191513, 0.0131486, 0.00757706, };
    double sysktl[NC]={2.7234, 2.80461, 2.66108, 2.64428, 2.28795, };
    //double ekTl[NC]={0.269123, 0.254184, 0.273911, 0.345952, 0.498483, };
    double sysmzl[NC]={0.721616, 0.702246, 0.679252, 0.645156, 0.59821, };
    //double emzl[NC]={0.00558996, 0.00303193, 0.00148695, 0.00251244, 0.00866513, };
    double sysxhql[NC]={0.729355, 0.777068, 0.847998, 0.904459, 0.970422, };
    //double exhql[NC]={0.0387191, 0.0294914, 0.0228741, 0.0176065, 0.00912529, };
    grtyp = pta510;
    ptlim[grtyp][0]=5; ptlim[grtyp][1]=10; 
    double zero[100]={0};
    grr[grtyp][tKTa5][tReal] = new TGraphErrors(NC, mpta, kT, zero, ekTl);
    grr[grtyp][tKTa5][tReal]->SetMarkerStyle(24);
    grr[grtyp][tKTa5][tSysL] = new TGraphErrors(NC, mpta, sysktl, zero, zero);
    grr[grtyp][tKTa5][tSysH] = new TGraphErrors(NC, mpta, syskth, zero, zero);
    
    grr[grtyp][tZTa5][tReal] = new TGraphErrors(NC, mpta, mz, zero, emzl);
    grr[grtyp][tZTa5][tReal]->SetMarkerStyle(24);
    grr[grtyp][tZTa5][tSysL] = new TGraphErrors(NC, mpta, sysmzl, zero, zero);
    grr[grtyp][tZTa5][tSysH] = new TGraphErrors(NC, mpta, sysmzh, zero, zero);

    grr[grtyp][tXhqa5][tReal] = new TGraphErrors(NC, mpta, xhq, zero, exhql);
    grr[grtyp][tXhqa5][tReal]->SetMarkerStyle(26);
    grr[grtyp][tXhqa5][tSysL] = new TGraphErrors(NC, mpta, sysxhql, zero, zero);
    grr[grtyp][tXhqa5][tSysH] = new TGraphErrors(NC, mpta, sysxhqh, zero, zero);
  }

  
};

