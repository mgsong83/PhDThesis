#include "/home/rak/nanoDST/pqm2002/paper/TDataV03.h"
#include "/home/rak/nanoDST/pqm2002/include/CCOR.h"

void paper(){
      
  enum enGrPT  { ptt1540, ptt1530, ptt24, ptt25, ptt1850, ptt1851, ptt1852, pttjt, pttjt36,
		  pta34, pta45, pta5580, pta510, pta710,
		  MCkt2_ptt, MCkt3_ptt};
  enum enGrWhat { tSigmaN, tSigmaF, tJt, tPout, tZKT, tKT, tZT, tXhq, tKTa3, tZTa3, tXhqa3, tKTa5, tZTa5, tXhqa5 };
  enum enGrSys { tReal, tSysL, tSysH };
  
  double hpttl=2, hptth=10;
  double sc=1.1, sdx=300, sdy;
  
  TDataV03 *v3data = new TDataV03();
  ccor = new CCOR(1.1);

  //v3data->GetPttDist();
  //v3data->GetPtaDist();
  v3data->Get_03_PttDist();
  v3data->Get_03_PtaDist();
  v3data->GetMCDist();
  int what[13]={
    0, //all 
    0, //sigma_ptt 1
    0, //jt_ptt 2
    0, //pout_ptt 3
    0, //xzkt_ptt 4
    0, //sigma_pta 5
    0, //jt_pta 6
    0, //pout_pta 7
    0, //xzkt_pta 8
    0, //fin kT_ptt 9
    0, //fin zt_xhq_ptt 10
    1, //fin kt_pta 11
    1, //fin zt_xhq_pta 12
  };
  
  
  enGrPT grtypN = pttjt;    // pttjt36
  enGrPT grtypF = ptt1540;  // 1530 1540, 1850;

  enGrPT grtyp1 = pta34;
  enGrPT grtyp2 = pta510;  //45
  enGrPT grtyp3 = pta710;
  //enGrPT grtyp3 = pta5580;

  
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  gStyle->SetMarkerSize(1.2);
  char chname[100];

  //--------------------------------
  if(what[0] || what[1]){
    mc(1, sc);
    mpad->SetLogy(0);  mpad->SetGridx(0);
    gPad->SetLeftMargin(0.17);
    hfr = new TH2F("hfr"," ", 10, hpttl, hptth, 10, 0, 0.8);
    hset( *hfr, "p_{Tt}    [GeV/c]","#sigma_{N,A}   [rad]",
	  1.15,1.1, 0.06, 0.06, 0.01,0.001, 0.05,0.05, 505,505);
    hfr->Draw();
    gStyle->SetFuncWidth(4);

    TGraphErrors *pyn2, *pyf2, *pyn3, *pyf3;
    pyn2 = v3data->GetGrr(MCkt2_ptt, tSigmaN);
    pyf2 = v3data->GetGrr(MCkt2_ptt, tSigmaF);
    pyn3 = v3data->GetGrr(MCkt3_ptt, tSigmaN);
    pyf3 = v3data->GetGrr(MCkt3_ptt, tSigmaF);
  
    pyn2->SetLineStyle(2);  pyf2->SetLineStyle(2);
    pyn2->Draw("Ce3");pyn2->SetFillColor(41);pyn2->Draw("Ce2");
    pyf2->Draw("Ce3");pyf2->SetFillColor(41);pyf2->Draw("Ce2");
    pyn3->Draw("Ce3");pyn3->SetFillColor(43);pyn3->Draw("Ce2");
    pyf3->Draw("Ce3");pyf3->SetFillColor(43);pyf3->Draw("Ce2");
    
    v3data->GetGrr(grtypN, tSigmaN)->Draw("p");
    v3data->GetGrr(grtypF, tSigmaF)->Draw("p");
    
    sprintf(chname,"%2.1f<p_{Ta}<%2.1f GeV/c", v3data->GetPTs(grtypF,0), v3data->GetPTs(grtypF,1));
    leg = new TLegend(0.5,0.62,0.72,0.87,chname,"brNDC");
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.04);
    leg->AddEntry(v3data->GetGrr(grtypF, tSigmaN), "Near-side  #sigma_{N}, data","p");
    leg->AddEntry(v3data->GetGrr(grtypF, tSigmaF), "Away-side #sigma_{A}, data","p");
    leg->AddEntry(pyn3, "pythia k_{T}=3.0 GeV/c","l");
    leg->AddEntry(pyf2, "              =2.0 GeV/c","l");
    leg->Draw();
    pPrint("/home/JanRak/run5pp_jTkT_pi0h/run35/sigmaPTTrig","c1");
  }

  //--------------------------------
  if(what[0] || what[2]){
    mc(2, sc, sdx);
    mpad->SetLogy(0);  mpad->SetGridx(0);
    gPad->SetLeftMargin(0.17);
    hfr = new TH2F("hfr"," ", 10, hpttl, hptth, 10, 0, 0.9);
    hset( *hfr, "p_{Tt}    [GeV/c]","#sqrt{#LTj^{2}_{T}#GT}   [GeV/c]",
	  1.15,1.1, 0.06,0.06, 0.01,0.001, 0.05,0.05, 505,505);
    hfr->Draw();
    
    v3data->GetGrr(grtypN,tJt)->Fit("pol0"," "," ",3,10);
    v3data->GetGrr(grtypN,tJt)->GetFunction("pol0")->SetLineStyle(2);
    v3data->GetGrr(grtypN,tJt)->Draw("p");
    ccor->GetJt()->Draw("p");
    
    sprintf(chname,"#sqrt{s}=200 GeV (%2.1f<p_{Ta}<%2.1f GeV/c)", v3data->GetPTs(grtypN,0), v3data->GetPTs(grtypN,1));
    leg = new TLegend(0.2,0.17,0.6,0.33, NULL,"brNDC");
    leg->SetFillStyle(0); leg->SetTextSize(0.04);  leg->SetBorderSize(0); 
    leg->AddEntry(v3data->GetGrr(grtypN,tJt), chname,"p");
    leg->AddEntry(ccor->GetJt(), "#sqrt{s}=62.4 GeV CCOR","p");
    leg->Draw();
    pPrint("/home/JanRak/run5pp_jTkT_pi0h/run35/jt_ptt","c2");
  }
    
  //--------------------------------
  if(what[0] || what[3]){
    mc(3, sc, sdx);
    mpad->SetLogy(0);  mpad->SetGridx(0);
    gPad->SetLeftMargin(0.17);
    hfr = new TH2F("hfr"," ", 10, hpttl, hptth, 10, 0, 2.5);
    hset( *hfr, "p_{Tt}    [GeV/c]","#sqrt{#LTp^{2}_{out}#GT}   [GeV/c]",
	  1.15,1.1, 0.06,0.06, 0.01,0.001, 0.05,0.05, 505,505);
    hfr->Draw();
    
    v3data->GetGrr(grtypF,tPout)->Draw("p");
    //  v3data->GetGrr(ptt24,tPout)->Draw("p");
    
    sprintf(chname,"%2.1f<p_{Ta}<%2.1f GeV/c", v3data->GetPTs(grtypF,0), v3data->GetPTs(grtypF,1));
    leg = new TLegend(0.5,0.62,0.72,0.87,chname,"brNDC");
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.04);
    leg->AddEntry(v3data->GetGrr(grtypF, tPout), "Away-side","p");
    leg->Draw();
    pPrint("/home/JanRak/run5pp_jTkT_pi0h/run35/poutt","c3");
  }
  
  //--------------------------------
  if(what[0] || what[4]){
    mc(4, sc, sdx);
    mpad->SetLogy(0);  mpad->SetGridx(0);
    gPad->SetLeftMargin(0.17);
    hfr = new TH2F("hfr"," ", 10, hpttl, hptth, 10, 0, 3.2);
    hset( *hfr, "p_{Tt}    [GeV/c]","#hat{x}^{-1}_{h}#LTz_{t}#GT#sqrt{#LTk^{2}_{T}#GT}   [GeV/c]",
	  1.15,1.1, 0.06,0.06, 0.01,0.001, 0.05,0.05, 505,505);
    hfr->Draw();
    
    TGraphErrors *py2, *py3;
    py2 = v3data->GetGrr(MCkt2_ptt, tPout, tReal);
    py3 = v3data->GetGrr(MCkt3_ptt, tPout, tReal);
    py2->SetLineStyle(2);
    py2->Draw("Ce3");py2->SetFillColor(18);py2->Draw("Ce2"); //41
    py3->Draw("Ce3");py3->SetFillColor(18);py3->Draw("Ce2");  //43
    
    v3data->GetGrr(grtypF, tZKT, tSysL)->Draw("2");
    v3data->GetGrr(grtypF, tZKT, tReal)->Draw("pz");
    //v3data->GetGrr(grtypF, tZKT, tReal)->Print();
    ccor->GetKt()->SetMarkerSize(1.4);
    ccor->GetKt()->Draw("p");
        
    leg = new TLegend(0.3,0.18,0.6,0.3, NULL,"brNDC");
    leg->SetFillStyle(0); leg->SetTextSize(0.04);  leg->SetBorderSize(0); 
    sprintf(chname,"#sqrt{s}=200 GeV; %2.1f<p_{Ta}<%2.1f GeV/c", v3data->GetPTs(grtypF,0), v3data->GetPTs(grtypF,1));
    leg->AddEntry(v3data->GetGrr(grtypF, tZKT, tReal), chname,"p");
    leg->AddEntry(ccor->GetKt(), "#sqrt{k^{2}_{T}} at #sqrt{s}=62.4 GeV","p");
//     leg->AddEntry(py3,     "             #sqrt{k^{2}_{T}}=3 GeV/c; PYTHIA","l");
//     leg->AddEntry(py2,     "                   2 GeV/c; PYTHIA","l");
    leg->Draw();
    pPrint("/home/JanRak/run5pp_jTkT_pi0h/run35/zkt","c4");
  }



  sdy=250;
  //--------------------------------
  if(what[0] || what[5]){
    mc(5, sc, sdx, sdy, 5);
    mpad->SetLogy(0);  mpad->SetGridx(0);
    gPad->SetLeftMargin(0.17);
    hfr = new TH2F("hfr"," ", 10, 1, 6, 10, 0, 0.7);
    hset( *hfr, "p_{Ta}    [GeV/c]","#sigma_{N,A}   [rad]",
	  1.15,1.1, 0.06,0.06, 0.01,0.001, 0.05,0.05, 505,505);
    hfr->Draw();
        
    v3data->GetGrr(grtyp1, tSigmaN)->Draw("p");
    v3data->GetGrr(grtyp1, tSigmaF)->Draw("p");
    
    v3data->GetGrr(grtyp2, tSigmaN)->Draw("p");
    v3data->GetGrr(grtyp2, tSigmaF)->Draw("p");
    
    //   v3data->GetGrr(pta5580, tSigmaN)->Draw("p");
    //   v3data->GetGrr(pta5580, tSigmaF)->Draw("pl");
    
    v3data->GetGrr(grtyp3, tSigmaN)->Draw("p");
    v3data->GetGrr(grtyp3, tSigmaF)->Draw("p");
    
    //run03
    //   const int NC = 6; //3-4
    //   double pta[NC] = {1.71639, 2.22036, 2.72949, 3.22803, 3.92835, 5.03744 };
    //   double sf[NC] = {0.538939, 0.527376, 0.438261, 0.456289, 0.454404, 0.327827 };
    //   double esf[NC] = {0.0222664, 0.0250665, 0.0302659, 0.0306966, 0.0255876, 0.026706 };
    //   double epta[NC]={0};
    //   grr03f_34 = new TGraphErrors(NC, pta, sf, epta, esf);
    //   grr03f_34->Draw("l");
    //   double pta[NC] = {1.71639, 2.22036, 2.72949, 3.22803, 3.92835, 5.03744 };
    //   double sn[NC] = {0.259379, 0.227474, 0.190342, 0.178368, 0.166084, 0.127265 };
    //   double esn[NC] = {0.00525097, 0.00531542, 0.00551501, 0.00612342, 0.00460053, 0.00723908 };
    //   grr03n_34 = new TGraphErrors(NC, pta, sn, epta, esn);
    //   grr03n_34->Draw("l");
    
  
    leg = new TLegend(0.5,0.62,0.72,0.87,NULL,"brNDC");
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.04);

    enGrTyp  grtyp = grtyp1;
    sprintf(chname,"%2.1f<p_{Tt}<%2.1f GeV/c", v3data->GetPTs(grtyp,0), v3data->GetPTs(grtyp,1));
    leg->AddEntry(v3data->GetGrr(grtyp, tSigmaN), chname,"p");
    grtyp = grtyp2;
    sprintf(chname,"%2.1f<p_{Tt}<%2.1f GeV/c", v3data->GetPTs(grtyp,0), v3data->GetPTs(grtyp,1));
    leg->AddEntry(v3data->GetGrr(grtyp, tSigmaN), chname,"p");
    //   grtyp = pta5580;
    //   sprintf(chname,"%2.1f<p_{Tt}<%2.1f GeV/c", v3data->GetPTs(grtyp,0), v3data->GetPTs(grtyp,1));
    //   leg->AddEntry(v3data->GetGrr(grtyp, tSigmaN), chname,"p");
    grtyp = grtyp3;
    sprintf(chname,"%2.1f<p_{Tt}<%2.1f GeV/c", v3data->GetPTs(grtyp,0), v3data->GetPTs(grtyp,1));
    leg->AddEntry(v3data->GetGrr(grtyp, tSigmaN), chname,"p");
    leg->Draw();
    pPrint("/home/JanRak/run5pp_jTkT_pi0h/run35/sigma_pta","c5");
  }
  
  //--------------------------------
  if(what[0] || what[6]){
    mc(6, sc, sdx, sdy, 5);
    mpad->SetLogy(0);  mpad->SetGridx(0);
    gPad->SetLeftMargin(0.17);
    hfr = new TH2F("hfr"," ", 10, 1, 6, 10, 0, 0.9);
    hset( *hfr, "p_{Ta}    [GeV/c]","#sqrt{#LTj^{2}_{T}#GT}   [GeV/c]",
	  1.15,1.1, 0.06,0.06, 0.01,0.001, 0.05,0.05, 505,505);
    hfr->Draw();
    
    grtyp = grtyp1;
    v3data->GetGrr(grtyp, tJt)->Fit("pol0","","",3,5);  v3data->GetGrr(grtyp, tJt)->Draw("p");
    grtyp = grtyp2;
    v3data->GetGrr(grtyp, tJt)->Fit("pol0","","",3,5);  v3data->GetGrr(grtyp, tJt)->Draw("p");
    grtyp = grtyp3;
    v3data->GetGrr(grtyp, tJt)->Fit("pol0","","",3,5);  v3data->GetGrr(grtyp, tJt)->Draw("p");
    pPrint("/home/JanRak/run5pp_jTkT_pi0h/run35/jt_pta","c6");
  }


  //--------------------------------
  if(what[0] || what[7]){
    mc(7, sc, sdx, sdy, 5);
    mpad->SetLogy(0);  mpad->SetGridx(0);
    gPad->SetLeftMargin(0.17);
    hfr = new TH2F("hfr"," ", 10, 1, 6, 10, 0, 2.5);
    hset( *hfr, "p_{Ta}    [GeV/c]","#sqrt{#LTp^{2}_{out}#GT}   [GeV/c]",
	  1.15,1.1, 0.06,0.06, 0.01,0.001, 0.05,0.05, 505,505);
    hfr->Draw();
    
    grtyp = grtyp1;
    v3data->GetGrr(grtyp, tPout)->Draw("p");
    grtyp = grtyp2;
    v3data->GetGrr(grtyp, tPout)->Draw("p");
    grtyp = grtyp3 ;
    v3data->GetGrr(grtyp, tPout)->Draw("p");
    pPrint("/home/JanRak/run5pp_jTkT_pi0h/run35/pouta","c7");
  }
  
  //--------------------------------
  if(what[0] || what[8]){
    mc(8, sc, sdx, sdy, 5);
    mpad->SetLogy(0);  mpad->SetGridx(0);
    gPad->SetLeftMargin(0.17);
    hfr = new TH2F("hfr"," ", 10, 1, 6, 10, 0, 3.0);
    hset( *hfr, "p_{Ta}    [GeV/c]","#hat{x}^{-1}_{h}#LTz_{t}#GT#sqrt{#LTk^{2}_{T}#GT}   [GeV/c]",
	  1.0, 1.1, 0.07,0.06, 0.01,0.001, 0.05,0.05, 505,505);
    hfr->Draw();
    
    v3data->GetGrr(grtyp1,  tZKT, tSysL)->Draw("2");
    v3data->GetGrr(grtyp2,  tZKT, tSysL)->Draw("2");
    v3data->GetGrr(grtyp1,  tZKT)->Draw("pz");
    v3data->GetGrr(grtyp2,  tZKT)->Draw("pz");

    //v3data->GetGrr(grtyp2,  tZKT)->Print();
    // v3data->GetGrr(grtyp3, tZKT)->Draw("p");
    
    //printGrr(v3data->GetGrr(grtyp1,  tZKT));
    
    leg = new TLegend(0.2, 0.2, 0.6, 0.35,NULL,"brNDC");
    leg->SetFillStyle(0); leg->SetTextSize(0.045);  leg->SetBorderSize(0); 
    sprintf(chname,"%2.0f<p_{Ta}(#pi^{0})<%2.0f GeV/c", v3data->GetPTs(grtyp1,0), v3data->GetPTs(grtyp1,1));
    leg->AddEntry(v3data->GetGrr(grtyp1, tZKT), chname,"p");
    sprintf(chname,"%2.0f<p_{Ta}(#pi^{0})<%2.0f GeV/c", v3data->GetPTs(grtyp2,0), v3data->GetPTs(grtyp2,1));
    leg->AddEntry(v3data->GetGrr(grtyp2, tZKT), chname,"p");
    leg->Draw();
    pPrint("/home/JanRak/ppg029/draft/final/zkt_pta","c8");
    //pPrint("/home/JanRak/run5pp_jTkT_pi0h/run35/zkt_pta","c8");
  }
    
  sdy=450;
  //--------------------------------
  if(what[0] || what[9]){
    mc(9, sc, sdx, sdy, 9);
    mpad->SetLogy(0);  mpad->SetGridx(0);
    gPad->SetLeftMargin(0.17);
    hfr = new TH2F("hfr"," ", 10, 0, 10, 1, 0.0, 4);
    hset( *hfr, "p_{Tt}  (GeV/#font[72]{c})", "#sqrt{#LTk^{2}_{T}#GT}  (GeV/#font[72]{c}) ",
	  0.9,1.2, 0.06,0.06, 0.01,0.001, 0.05,0.05, 505,505);
    hfr->Draw();
    
    v3data->GetGrr(grtypF, tKT, tReal)->Draw("pz");
    v3data->GetGrr(grtypF, tKT, tSysH)->Draw("l");
    v3data->GetGrr(grtypF, tKT, tSysL)->Draw("l");
    //v3data->GetGrr(grtypF, tKT, tSysL)->Print();
    //f0 = new TF1("f0","pol0",3.0,9); f0->SetLineStyle(2);
    //v3data->GetGrr(grtypF, tKT, tSysL)->Fit("f0","R");

    ccor->GetKt()->Draw("p");

    leg = new TLegend(0.13,0.15,0.7,0.25,NULL,"brNDC");
    leg->SetFillStyle(0); leg->SetTextSize(0.037);  leg->SetBorderSize(0); 
    leg->AddEntry(ccor->GetKt(), "#sqrt{s}=63 GeV","p");
    leg->AddEntry(v3data->GetGrr(grtypF, tKT, tReal), "#sqrt{s}=200 GeV 1.4<p_{Ta}<5 GeV/c","p");
    leg->Draw();
    pPrint("/home/JanRak/ppg029/draft/final/kt_ptt","c9");
  }

  //--------------------------------
  if(what[0] || what[10]){
    mc(10, sc, sdx, sdy, 9);
    mpad->SetLogy(0);  mpad->SetGridx(0);
    gPad->SetLeftMargin(0.17);
    hfr = new TH2F("hfr"," ", 10, 0, 10, 1, 0.0, 1.0);
    hset( *hfr, "p_{Tt}  (GeV/#font[72]{c})", "#LTz_{t}#GT, #hat{x}_{h}",
	  1.1, 1.2, 0.06,0.06, 0.01,0.001, 0.05,0.05, 505,505);
    hfr->Draw();

    v3data->GetGrr(grtypF, tZT, tReal)->Draw("pz");
    v3data->GetGrr(grtypF, tZT, tSysH)->Draw("l");
    v3data->GetGrr(grtypF, tZT, tSysL)->Draw("l");

    v3data->GetGrr(grtypF, tXhq, tReal)->Draw("pz");
    v3data->GetGrr(grtypF, tXhq, tSysH)->Draw("l");
    v3data->GetGrr(grtypF, tXhq, tSysL)->Draw("l");
    //v3data->GetGrr(grtypF, tXhq, tReal)->Print();

    leg = new TLegend(0.31,0.2,0.79,0.35,NULL,"brNDC");
    leg->SetFillStyle(0); leg->SetTextSize(0.045);  leg->SetBorderSize(0); 
    sprintf(chname,"#LTz_{t}#GT %2.1f<p_{Ta}<%2.1f GeV/c",
	    v3data->GetPTs(grtypF, 0), v3data->GetPTs(grtypF, 1));
    leg->AddEntry(v3data->GetGrr(grtypF, tZT, tReal), chname,"p");
    sprintf(chname,"#hat{x}_{h}   %2.1f<p_{Ta}<%2.1f GeV/c",
	    v3data->GetPTs(grtypF, 0), v3data->GetPTs(grtypF, 1));
    leg->AddEntry(v3data->GetGrr(grtypF, tXhq, tReal), chname,"p");
    leg->Draw();
    pPrint("/home/JanRak/ppg029/draft/final/zt_xhq_ptt","c10");    
  }
  
  //--------------------------------
  if(what[0] || what[11]){
    mc(11, sc, sdx, sdy, 9);
    mpad->SetLogy(0);  mpad->SetGridx(0);
    gPad->SetLeftMargin(0.17);
    hfr = new TH2F("hfr"," ", 10, 0, 6, 1, 0.0, 4);
    hset( *hfr, "p_{Ta}  (GeV/#font[72]{c})", "#sqrt{#LTk^{2}_{T}#GT}  (GeV/#font[72]{c}) ",
	  0.9,1.2, 0.06,0.06, 0.01,0.001, 0.05,0.05, 505,505);
    hfr->Draw();
    
    v3data->GetGrr(grtyp1, tKTa3, tReal)->Draw("pz");
    //v3data->GetGrr(grtyp1, tKTa3, tReal)->Print();
    v3data->GetGrr(grtyp1, tKTa3, tSysH)->Draw("l");
    v3data->GetGrr(grtyp1, tKTa3, tSysL)->Draw("l");

    v3data->GetGrr(grtyp2, tKTa5, tReal)->Draw("pz");
    //v3data->GetGrr(grtyp2, tKTa5, tReal)->Print();
    v3data->GetGrr(grtyp2, tKTa5, tSysH)->SetLineStyle(2);
    v3data->GetGrr(grtyp2, tKTa5, tSysH)->Draw("l");
    v3data->GetGrr(grtyp2, tKTa5, tSysL)->SetLineStyle(2);
    v3data->GetGrr(grtyp2, tKTa5, tSysL)->Draw("l");
//     f0 = new TF1("f0","pol0",1.5,5); f0->SetLineStyle(2);
//     v3data->GetGrr(grtyp2, tKTa5, tReal)->Fit("f0","R");

    leg = new TLegend(0.3,0.2,0.7,0.4,NULL,"brNDC");
    leg->SetFillStyle(0); leg->SetTextSize(0.045);  leg->SetBorderSize(0); 
    sprintf(chname,"%2.0f<p_{Tt}(#pi^{0})<%2.0f GeV/c",
	    v3data->GetPTs(grtyp1,0), v3data->GetPTs(grtyp1,1));
    leg->AddEntry(v3data->GetGrr(grtyp1, tKTa3), chname,"p");
    sprintf(chname,"%2.0f<p_{Tt}(#pi^{0})<%2.0f GeV/c",
	    v3data->GetPTs(grtyp2,0), v3data->GetPTs(grtyp2,1));
    leg->AddEntry(v3data->GetGrr(grtyp2, tKTa5), chname,"p");
    leg->Draw();
    pPrint("/home/JanRak/ppg029/draft/final/kt_pta","c11");    
  }

  //--------------------------------
  if(what[0] || what[12]){
    mc(12, sc, sdx, sdy, 9);
    mpad->SetLogy(0);  mpad->SetGridx(0);
    gPad->SetLeftMargin(0.17);
    hfr = new TH2F("hfr"," ", 10, 0, 6, 1, 0.0, 1.2);
    hset( *hfr, "p_{Ta}  (GeV/#font[72]{c})", "#LTz_{t}#GT, #hat{x}_{h}",
	  1.1, 1.2, 0.06,0.06, 0.01,0.001, 0.05,0.05, 505,505);
    hfr->Draw();

    v3data->GetGrr(grtyp1, tZTa3, tReal)->Draw("pz");
    v3data->GetGrr(grtyp1, tZTa3, tSysH)->Draw("l");
    v3data->GetGrr(grtyp1, tZTa3, tSysL)->Draw("l");

    v3data->GetGrr(grtyp1, tXhqa3, tReal)->Draw("pz");
    v3data->GetGrr(grtyp1, tXhqa3, tSysH)->Draw("l");
    v3data->GetGrr(grtyp1, tXhqa3, tSysL)->Draw("l");

    v3data->GetGrr(grtyp2, tZTa5, tReal)->Draw("pz");
    v3data->GetGrr(grtyp2, tZTa5, tSysH)->SetLineStyle(2);
    v3data->GetGrr(grtyp2, tZTa5, tSysL)->SetLineStyle(2);
    v3data->GetGrr(grtyp2, tZTa5, tSysH)->Draw("l");
    v3data->GetGrr(grtyp2, tZTa5, tSysL)->Draw("l");

    v3data->GetGrr(grtyp2, tXhqa5, tReal)->Draw("pz");
    v3data->GetGrr(grtyp2, tXhqa5, tSysH)->SetLineStyle(2);
    v3data->GetGrr(grtyp2, tXhqa5, tSysL)->SetLineStyle(2);
    v3data->GetGrr(grtyp2, tXhqa5, tSysH)->Draw("l");
    v3data->GetGrr(grtyp2, tXhqa5, tSysL)->Draw("l");

    v3data->GetGrr(grtyp2, tXhqa5, tReal)->Print();

    
    leg = new TLegend(0.31,0.15,0.79,0.35,NULL,"brNDC");
    leg->SetFillStyle(0); leg->SetTextSize(0.035);  leg->SetBorderSize(0);
    
    sprintf(chname,"#LTz_{t}#GT  %2.1f<p_{Tt}<%2.1f GeV/c",
	    v3data->GetPTs(grtyp1, 0), v3data->GetPTs(grtyp1, 1));
    leg->AddEntry(v3data->GetGrr(grtyp1, tZTa3, tReal), chname,"p");
    leg->AddEntry(v3data->GetGrr(grtyp1, tXhqa3, tReal), "#hat{x}_{h}","p");

    sprintf(chname,"#LTz_{t}#GT %2.1f<p_{Tt}<%2.1f GeV/c",
	    v3data->GetPTs(grtyp2, 0), v3data->GetPTs(grtyp2, 1));
    leg->AddEntry(v3data->GetGrr(grtyp2, tZTa5, tReal), chname,"p");
    leg->AddEntry(v3data->GetGrr(grtyp2, tXhqa5, tReal), "#hat{x}_{h}","p");

    leg->Draw();
    pPrint("/home/JanRak/ppg029/draft/final/zt_xhq_pta","c12");    
  }

}
