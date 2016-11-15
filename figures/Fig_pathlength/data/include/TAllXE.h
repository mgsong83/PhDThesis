#define MAXD 100;

class TAllXE {

    public:
        TGraphErrors *xegrr[MAXD], *ptgrr[MAXD], *pi0cross;
        double ptt[MAXD], sc[MAXD], normAt[MAXD], norm[MAXD], ptNormAt[MAXD], ptNorm[MAXD];
        double pttl[MAXD], ptth[MAXD], SoverB;
        int ixegrr;
        double msz;
        double InclNormAt, InclNorm;

    private:
        TH1D     *hxEFNew[2][MAXD], *hxEF[MAXD], *hpTN[MAXD], *hpTF[2][MAXD], *hpTFnomix[MAXD];
        TProfile *pxEpTN[MAXD], *pxExEN[MAXD], *ppTpTN[MAXD], *pdphiN[MAXD];
        TProfile *pxEpTF[MAXD], *pxExEF[MAXD], *ppTpTF[MAXD], *pdphiF[MAXD];
        TH1D     *hTriggPtBin[MAXD], *hptt[MAXD], *hInclusive;
        TF1 *effPt;  

    public:

        TAllXE(){
            ixegrr=0;
            msz=1.2;
            effPt = new TF1("effPt","0.187*([0]/x+[1]+[2]*x+[3]*x*x)",0,6);
            effPt->SetParameters(4.017,23.52,-0.3135,-0.0108);
        }

        TGraphErrors *GetGrr(int i){   check(i); return xegrr[i]; }
        TGraphErrors *GetPtGrr(int i){ check(i); return ptgrr[i]; }
        double GetPtt(int i){ check(i); return ptt[i]; }
        double GetLPtt(int i){ check(i); return pttl[i]; }
        double GetHPtt(int i){ check(i); return ptth[i]; }
        double GetSc(int i){ check(i); return sc[i]; }
        double GetNormAt(int i) {check(i); return normAt[i];}
        double GetNorm(int i) {check(i); return norm[i];}
        double GetPtNormAt(int i) {check(i); return ptNormAt[i];}
        double GetPtNorm(int i) {check(i); return ptNorm[i];}

        TGraphErrors *GetInclGrr(){return pi0cross; }
        double GetInclNormAt() {return InclNormAt;}
        double GetInclNorm() {return InclNorm;}

        int GetNoGrr(){ return ixegrr;}

        void check(int i){
            if(i<ixegrr)
                return;
            else {
                cout<<"i = "<<i<<" no such xe graph"<<endl;
                exit(0);
            }
        }

        void PrintData(int ingrr=0){
            double x[500],y[500],ey[500];
            for(int ii=0; ii<xegrr[ingrr]->GetN();ii++){
                xegrr[ingrr]->GetPoint(ii,x[ii],y[ii]);
                ey[ii] = xegrr[ingrr]->GetErrorY(ii);
                //       cout<<setprecision(2)<<x[ii]<<" & "<<setprecision(2)<<y[ii]/sc[ingrr]
                // 	  <<" $\pm$ "<<setprecision(2)<<ey[ii]/sc[ingrr]<<" & "<<endl;
                printf("%1.2f & %1.1e $\pm$ %1.1e &\n",x[ii],y[ii]/sc[ingrr],ey[ii]/sc[ingrr]);
            }
        }


        void AddInclusExp(){
            // pi0 cross section =================================================
            const int ND =  23 ;
            double diff_data[ND], ediff_data[ND], ecross[ND];
            double pt[ND]    = {   1.215,    1.719,    2.223,    2.725,    3.228,     3.73,    4.231,    4.732,    5.234,    5.735,    6.237,    6.738,    7.238,    7.739,     8.24,     8.74,    9.241,    9.741,    10.88,     12.9,    14.91,    16.92,    18.93};
            double cross[ND] = {  0.3912,  0.06237,  0.01278,  0.00333, 0.000984, 0.0003348, 0.0001199, 4.915e-05, 2.161e-05, 1.049e-05, 5.265e-06, 2.752e-06, 1.57e-06, 9.109e-07, 5.487e-07, 3.514e-07, 2.142e-07, 1.387e-07, 5.796e-08, 1.482e-08, 3.891e-09, 1.185e-09, 6.309e-10};
            double estat[ND] = {0.0003148, 8.048e-05, 2.788e-05, 1.215e-05, 5.773e-06, 3.031e-06, 2.471e-07, 1.467e-07, 9.113e-08, 6.025e-08, 4.051e-08, 2.792e-08, 2.008e-08, 1.48e-08, 1.103e-08, 8.766e-09, 6.336e-09, 5.099e-09, 1.554e-09, 7.144e-10, 3.247e-10, 1.719e-10, 1.304e-10};
            double esys[ND]  = { 0.02444, 0.003502, 0.0007218, 0.0001923, 5.837e-05, 2.047e-05, 1.141e-05, 4.74e-06, 2.11e-06, 1.037e-06, 5.271e-07, 2.786e-07, 1.608e-07, 9.427e-08, 5.741e-08, 3.717e-08, 2.288e-08, 1.491e-08, 6.504e-09, 1.787e-09, 5.056e-10, 1.741e-10, 1.028e-10};
            //   const int ND = 17;
            //   double pt[ND]     = { 1.215, 1.719, 2.223, 2.726, 3.228, 3.73, 4.232, 4.733, 5.234, 5.735, 6.236, 6.737, 7.452, 8.457, 9.46, 10.861, 13.25};
            //   double cross[ND]  = { 0.3733, 0.06052, 0.01221, 0.003308, 0.0009978, 0.0003385, 0.0001187, 4.726e-05, 2.206e-05, 1.113e-05, 4.999e-06, 3.003e-06, 1.08e-06, 4.853e-07, 1.643e-07, 5.227e-08, 1.19e-08};
            //   double estat[ND] = { 0.006117, 0.001075, 0.0003011, 0.0001181, 5.652e-05, 2.456e-05, 2.906e-06, 1.99e-06, 1.096e-06, 5.005e-07, 3.172e-07, 2.304e-07, 9.53e-08, 5.832e-08, 3.166e-08, 1.169e-08, 4.91e-09};
            //   double esys[ND]  = { 0.02731, 0.004284, 0.0008644, 0.0002375, 7.278e-05, 2.591e-05, 9.888e-06, 4.001e-06, 1.929e-06, 1.028e-06, 4.759e-07, 2.95e-07, 1.092e-07, 5.253e-08, 1.813e-08, 6.137e-09, 1.894e-09};
            for(int i=0;i<ND;i++) ecross[i] = sqrt(estat[i]*estat[i]+esys[i]*esys[i]);
            double ept[ND]    = {0};
            pi0cross = new TGraphErrors(ND, pt, cross, ept, ecross );
            pi0cross->SetMarkerStyle(20);

            InclNormAt=8;     //was at 5, InclNorm = fpow->Eval(InclNormAt);  
            fpow = new TF1("fpow","[0]/pow(x,[1])",3,20);
            fpow->SetParameters(1.39346e+01, 8.09682e+00); //run05 fit
            InclNorm = fpow->Eval(InclNormAt);
            // ====================================================================
        }

        //------------------------------------------
        void AddFromFile( TString InFile ){
            //------------------------------------------

            const int NH=13;
            double insc[NH];
            for(int i=0;i<NH;i++) insc[i]=1e2*pow(10,-i);      
            double XENormsAt[NH] = {0.892512, 0.685266, 0.555796, 0.467958, 0.40415, 0.355583, 0.31785, 0.278366, 0.230191, 0.180739, 0.230191, 0.180739, 0.180739  };
            double XENorms[NH] = {0.750194, 0.131934, 0.0217233, 0.00325984, 0.000406407, 6.15334e-05, 6.61735e-06, 7.55714e-07, 1.06371e-07, 1.15507e-08 , 1.06371e-07, 1.15507e-08, 1.15507e-08 };
            double Pi0PtNormsAt[NH] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
            double Pi0PtNorms[NH] = {0.0487263, 0.0113764, 0.00172654, 0.000205437, 2.24061e-05, 3.0811e-06 , 3e-7, 3e-7, 3e-7, 3e-7, 3e-7, 3e-7, 3e-7};
            double inSoverB  [NH]= { 1.45587, 2.0981,  2.88227, 4.04818, 4.99095, 4.2165, 20, 20, 20, 20, 20, 20, 20}; 
            int marker[NH]= {29, 24, 21, 25, 22, 27, 28, 23, 24, 26, 20, 21, 21};

            cout<<"Reading from "<<InFile<<" ";
            fin = new TFile(InFile);
            fin->cd("JCard");
            const int NumCent = CentBinBorders->GetNoElements()-2; 
            const int NumPtt  = TriggPtBorders->GetNoElements()-2;
            const int NumPta  = AssocPtBorders->GetNoElements()-2;
            cout<<endl<<"+++++++++++++++++++++++++++++++"<<endl;
            cout<<InFile<<" cent=0-"<<NumCent<<"  #trig_pt=0-"<<NumPtt<<"  #assoc_pt=0-"<<NumPta<<endl<<endl;
            cout<<endl<<"+++++++++++++++++++++++++++++++"<<endl;
            
            int lowih=0;
            if(NumPtt-lowih > NH) {
                cout<<"you have a wrong file, NumPtt="<<NumPtt-lowih<<endl;
                exit(0);
            }

            TString StrName;
            for(int ih=lowih; ih<=NumPtt; ih++){
                StrName = "hxEF0"; StrName+=ih; hxEFNew[0][ih-lowih] = (TH1D *)fin->Get(StrName);
                //cout<<endl<<StrName<<" ";
                StrName = "hxEF1"; StrName+=ih; hxEFNew[1][ih-lowih] = (TH1D *)fin->Get(StrName);
                StrName = "hpTF0"; StrName+=ih; hpTF[0][ih-lowih] = (TH1D *)fin->Get(StrName);      
                StrName = "hpTF1"; StrName+=ih; hpTF[1][ih-lowih] = (TH1D *)fin->Get(StrName);      
                StrName = "pxEpTF"; StrName+=ih; pxEpTF[ih-lowih] = (TProfile *)fin->Get(StrName);
                StrName = "pxExEF"; StrName+=ih; pxExEF[ih-lowih] = (TProfile *)fin->Get(StrName);
                StrName = "ppTpTF"; StrName+=ih; ppTpTF[ih-lowih] = (TProfile *)fin->Get(StrName);
                StrName = "pdphiF"; StrName+=ih; pdphiF[ih-lowih] = (TProfile *)fin->Get(StrName);      
                StrName = "hTriggPtBin0"; StrName+=ih; hTriggPtBin[ih-lowih] = (TH1D *)fin->Get(StrName);
                //cout<<ih<<" "<< hxEFNew[0][ih-lowih]<<" "<<pdphiF[ih-lowih]<<endl;
            }
            //    for(int i=0;i<NH;i++){
            for(int i=0;i<=NumPtt-lowih;i++){
                sc[ixegrr]=insc[i];
                ptt[ixegrr]=hTriggPtBin[i]->GetMean();
                SoverB = inSoverB[i];
                //=============
                //xegrr[ixegrr] = getXEgrr_New(i,sc[ixegrr], 1.0, 4.8); //1. 4.8
                xegrr[ixegrr] = getXEgrr(i,sc[ixegrr], 0.5, 10.0); //1. 4.8
                //=============

                xegrr[ixegrr]->SetMarkerStyle(marker[i]);
                ptgrr[ixegrr] = getPTgrr_New(i,sc[ixegrr]);
                ptgrr[ixegrr]->SetMarkerStyle(marker[i]);
                normAt[ixegrr]= XENormsAt[i];
                norm[ixegrr]  = XENorms[i];
                ptNormAt[ixegrr]= Pi0PtNormsAt[i];
                ptNorm[ixegrr]  = Pi0PtNorms[i];
                pttl[ixegrr] = hTriggPtBin[i]->GetXaxis()->GetXmin();
                ptth[ixegrr] = hTriggPtBin[i]->GetXaxis()->GetXmax();
                ixegrr++;
            }
            cout<<ixegrr<<" histos in "<<endl;
        }


        //------------------------------------------
        TGraphErrors *getXEgrr(int ih, double sc=1, double pTa_l=2.0, double pTa_h=5.0){
            //------------------------------------------
            TH1D *hfxe     =  hxEFNew[0][ih];
            TH1D *hbxe     =  hxEFNew[1][ih];
            TProfile *pxept =  pxEpTF[ih];
            TProfile *pxexe =  pxExEF[ih];
            TProfile *pdphi =  pdphiF[ih];
            int NC=hfxe->GetNbinsX();
            double xE[100], exE[100]={0};
            double normYF[100], enormYF[100];
            double triggs = hTriggPtBin[ih]->GetEntries();
            double ptt    = hTriggPtBin[ih]->GetMean(), wpta, wxe, w;
            double bw  = hfxe->GetBinWidth(1);
            int good=0, rebin=1;
            if(ptt>5) rebin=1;

            double forSum = hfxe->GetEntries();
            double bckSum = 1;//hbxe->GetEntries();
            double forToMix = bckSum>0 ? forSum/bckSum/SoverB : 1;

            for(int ib=0; ib<NC+1-rebin; ib+=rebin){
                //      if(hfxe->GetBinContent(ib+1)>0){
                wpta=0;
                wxe =0;
                w=0;
                for(int ii=0; ii<rebin; ii++){
                    double bck  = 0; //hbxe->GetBinContent(ib+1+ii)*forToMix;
                    double forg = hfxe->GetBinContent(ib+1+ii);// - bck;
                    wpta += pxept->GetBinContent(ib+1+ii)*forg;
                    wxe  += pxexe->GetBinContent(ib+1+ii)*forg;
                    //wxe  += hfxe->GetBinCenter(ib+1+ii)*forg;
                    w    += forg;
                }
                if(w<=0) continue;
                wpta /= w;
                wxe /= w;
                if(wxe < pTa_l/ptt || wxe > pTa_h/ptt ) continue;
                //double bw  = hfxe->GetBinWidth(ib+1);
                double eff = 1; //effPt->Eval(wpta);
                double norm=sc*eff/triggs/bw/ptt/rebin;   // 1/ptt is Jacobian dN/dpTassoc = 1/ptt dN/dxE
                xE[good]     =  wxe; //pxexe->GetBinContent(ib+1)/cos(pdphi->GetBinError(ib+1));
                normYF[good] =  w*norm;  //hfxe->GetBinContent(ib+1)*norm*sc/mixed;
                enormYF[good]=  sqrt(w)*norm;  //hfxe->GetBinError(ib+1)*norm*sc/mixed;
                //cout<<good<<" "<<triggs<<" "<<bw<<" "<<norm<<" "<<xE[good]<<" "<<normYF[good]<<endl;
                //cout<<good<<" "<<ib<<" "<<xE[good]<<" "<<pta<<endl;
                //printf("%3.2f %3.2f & %3.1e $\\pm$ %3.1e \n", ptt, xE[good],normYF[good],enormYF[good]);
                good++;
            }
            TGraphErrors *gr = new TGraphErrors(good, xE, normYF, exE, enormYF);  
            gr->SetTitle(hfxe->GetTitle());
            return gr;
            }

            //------------------------------------------
            TGraphErrors  *getPTgrr(int ih, double sc=1){
                //------------------------------------------

                TH1D *hid = (TH1D*)hpTFnomix[ih]->Clone();
                TProfile *pptpt =  ppTpTF[ih];
                hid->Sumw2();
                double pta[100], epta[100]={0}, YF[100], eYF[100];

                double ptt=hTriggPtBin[ih]->GetMean(), wpt, w;
                int good=0, rebin=1;
                if(ptt>5) rebin=2;
                if(ptt>6) rebin=3;
                if(ptt>7) rebin=3;

                int NC=hid->GetNbinsX();
                double bw  = hid->GetBinWidth(1);
                double triggs = hTriggPtBin[ih]->GetEntries();
                int first=1;

                double forSum = hid->GetEntries();

                for(int ib=0; ib<NC+1-rebin; ib+=rebin){
                    wpt=0;
                    w=0;
                    for(int ii=0; ii<rebin; ii++){
                        double forg = hid->GetBinContent(ib+1+ii);
                        //wpt += pptpt->GetBinContent(ib+1+ii)*forg;
                        wpt += hid->GetBinCenter(ib+1+ii)*forg;
                        w   += forg;
                    }
                    if(w<=0 ) continue;
                    if(first){
                        first=0;
                        continue;
                    }
                    pta[good] = wpt/w;
                    if(pta[good]>5.5) continue;
                    if(ptt>8 && pta[good]>4.7) continue;
                    double norm=sc/triggs/bw*effPt->Eval(pta[good])/rebin;
                    YF[good]   = w*norm;
                    eYF[good]  = sqrt(w)*norm; //hid->GetBinError(ib+1)*norm;
                    good++;
                }

                TGraphErrors *g = new TGraphErrors(good, pta, YF, epta, eYF);
                g->SetTitle(hid->GetTitle());
                return g;
            }


            //------------------------------------------
            TGraphErrors  *getPTgrr_New(int ih, double sc=1){
                //------------------------------------------

                TH1D *hid = (TH1D*)hpTF[0][ih]->Clone();
                TH1D *hib = (TH1D*)hpTF[0][ih]->Clone();
                TProfile *pptpt =  ppTpTF[ih];
                hid->Sumw2();
                double pta[100], epta[100]={0}, YF[100], eYF[100];

                double ptt=hTriggPtBin[ih]->GetMean(), wpt, w;
                int good=0, rebin=1;
                if(ptt>5) rebin=2;
                if(ptt>6) rebin=3;
                if(ptt>7) rebin=3;

                int NC=hid->GetNbinsX();
                double bw  = hid->GetBinWidth(1);
                double triggs = hTriggPtBin[ih]->GetEntries();
                int first=1;

                double forSum = hid->GetEntries();
                double bckSum = 1;//hib->GetEntries();
                double forToMix = bckSum>0 ? forSum/bckSum/SoverB : 1;

                for(int ib=0; ib<NC+1-rebin; ib+=rebin){
                    wpt=0;
                    w=0;
                    for(int ii=0; ii<rebin; ii++){
                        double bck  = 0;//hib->GetBinContent(ib+1+ii)*forToMix;
                        double forg = hid->GetBinContent(ib+1+ii) - bck;
                        //wpt += pptpt->GetBinContent(ib+1+ii)*forg;
                        wpt += hib->GetBinCenter(ib+1+ii)*forg;
                        w   += forg;
                    }
                    if(w<=0 ) continue;
                    if(first){
                        first=0;
                        continue;
                    }
                    pta[good] = wpt/w;
                    if(pta[good]>5.5) continue;
                    if(ptt>8 && pta[good]>4.7) continue;
                    double norm=sc/triggs/bw*effPt->Eval(pta[good])/rebin;
                    YF[good]   = w*norm;
                    eYF[good]  = sqrt(w)*norm; //hid->GetBinError(ib+1)*norm;
                    good++;
                }

                TGraphErrors *g = new TGraphErrors(good, pta, YF, epta, eYF);
                g->SetTitle(hid->GetTitle());
                return g;
            }


            //------------------------------------------
            AddJia(){
                //------------------------------------------

                const int NCT=6;
                double inptt[NCT] = {5.43, 6.31, 7.28, 8.61, 10.65, 13.27};
                double inpttl[NCT] = {5,6,7,8,10,12};
                double inptth[NCT] = {6,7,8,10,12,16};
                double insc[NCT]  = {1,    1e-1, 1e-2, 1e-3, 1e-4,  1e-5};
                double PipmNormAt=0.1;
                double PipmNorms[NCT] = {8.35497, 0.838858, 0.0787992, 0.0117858, 0.00110482, 9.58872e-05};
                double zero[50] = {0};
                for(int i=0;i<NCT;i++){
                    ptt[ixegrr+i]    = inptt[i];
                    sc[ixegrr+i]     = insc[i];
                    normAt[ixegrr+i] = PipmNormAt;
                    norm[ixegrr+i]   = PipmNorms[i];
                    pttl[ixegrr+i] = inpttl[i];
                    ptth[ixegrr+i] = inptth[i];
                }

                const int NC = 33;
                double xE[NC]       = {   0.09,      0.11,    0.13,    0.15,    0.17,    0.19,       0.21,     0.23,    0.25,     0.27,      0.29,    0.31,    0.33,       0.35,     0.37,    0.39,      0.41,     0.43,     0.45,     0.47,     0.49,     0.51,     0.53,     0.55,     0.57,     0.59,    0.61,     0.63,     0.65,     0.67,      0.69,      0.71,    0.75};
                double xEyield[NC]  = {10.4629,   7.37825, 5.79812, 4.22275, 4.98068, 3.57347,    2.84904,  2.61984,  2.0275,  1.65985,   1.27239, 1.43563, 1.23164,   0.825873, 0.571647, 1.02383,  0.385127, 0.321005, 0.581726, 0.802587, 0.420161, 0.390575, 0.570991, 0.344496, 0.289838, 0.254824, 0.34605, 0.383122, 0.263067, 0.292372, 0.0917321, 0.0841873, 0.118473};
                double exEyield[NC] = {0.943269, 0.731912, 0.618526, 0.501761, 0.49362, 0.40975, 0.375422, 0.333672, 0.290257, 0.246909, 0.220229, 0.231007, 0.230915, 0.186823, 0.180878, 0.222021, 0.154512, 0.13543, 0.165204, 0.200602,  0.13884,  0.139991, 0.160031, 0.127793,  0.11485, 0.108744, 0.12613, 0.131419, 0.109294, 0.112017, 0.0672707, 0.0616134, 0.0702532};
                for(int i=0;i<NC;i++){
                    xEyield[i]  *= sc[ixegrr];
                    exEyield[i] *= sc[ixegrr];
                }
                xegrr[ixegrr] = new TGraphErrors(NC, xE, xEyield, zero, exEyield);
                xegrr[ixegrr]->SetMarkerSize(msz);
                xegrr[ixegrr++]->SetMarkerStyle(20);

                const int NC = 29;
                double xE[NC] = {0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.27, 0.29, 0.31, 0.33, 0.35, 0.37, 0.39, 0.41, 0.43, 0.45, 0.47, 0.49, 0.51, 0.53, 0.55, 0.57, 0.59, 0.61, 0.63};
                double xEyield[NC] = {14.6233, 10.223, 6.92365, 5.14871, 4.60475, 3.0561, 3.18087, 3.16883, 2.68896, 1.9568, 2.08802, 1.38273, 1.32056, 1.52548, 0.651595, 0.660068, 0.749598, 0.564136, 0.896963, 0.554129, 0.635135, 0.352991, 0.555512, 0.340856, 0.383412, 0.497981, 0.125165, 0.234853, 0.0912686};
                double exEyield[NC] = {1.15862, 0.906145, 0.687555, 0.575504, 0.515905, 0.417809, 0.399777, 0.393402, 0.328804, 0.28983, 0.344181, 0.290566, 0.284763, 0.303205, 0.209128, 0.201171, 0.205636, 0.179032, 0.217299, 0.174571, 0.184417, 0.141584, 0.175176, 0.133885, 0.138685, 0.160709, 0.0772643, 0.108247, 0.0676516};
                for(int i=0;i<NC;i++){
                    xEyield[i]  *= sc[ixegrr];
                    exEyield[i] *= sc[ixegrr];
                }
                xegrr[ixegrr] = new TGraphErrors(NC, xE, xEyield, zero, exEyield);
                xegrr[ixegrr]->SetMarkerSize(msz);
                xegrr[ixegrr++]->SetMarkerStyle(24);

                const int NC = 24;
                double xE[NC]       = {0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.27, 0.29, 0.31, 0.33, 0.35, 0.37, 0.39, 0.41, 0.43, 0.45, 0.47, 0.49, 0.51};
                double xEyield[NC]  = {9.47935, 13.7164, 10.8773, 8.18869, 5.08611, 5.19384, 2.46006, 2.81683, 1.7995, 2.67155, 2.04583, 1.11442, 1.25879, 0.867667, 1.2775, 0.607706, 1.08027, 1.19047, 0.824809, 0.627602, 0.327297, 0.25548, 0.358955, 0.181038};
                double exEyield[NC] = {1.72404, 1.55349, 1.1563, 0.942171, 0.741234, 0.703735, 0.503727, 0.502558, 0.378209, 0.471264, 0.437384, 0.343657, 0.363776, 0.296347, 0.372475, 0.25728, 0.330956, 0.342373, 0.306162, 0.244693, 0.17213, 0.151335, 0.184497, 0.134466};
                for(int i=0;i<NC;i++){
                    xEyield[i]  *= sc[ixegrr];
                    exEyield[i] *= sc[ixegrr];
                }
                xegrr[ixegrr] = new TGraphErrors(NC, xE, xEyield, zero, exEyield);
                xegrr[ixegrr]->SetMarkerSize(msz);
                xegrr[ixegrr++]->SetMarkerStyle(21);

                const int NC = 21;
                double xE[NC] = {0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.27, 0.29, 0.31, 0.33, 0.35, 0.37, 0.39, 0.41, 0.43, 0.45};
                double xEyield[NC] = {21.9236, 24.2565, 13.0568, 9.05137, 4.96246, 6.54269, 4.98254, 4.34761, 2.2251, 2.14975, 1.79706, 2.09865, 1.4095, 1.18981, 0.900859, 0.634165, 0.714002, 0.136205, 0.259855, 0.3995, 0.420841};
                double exEyield[NC] = {2.30704, 1.77138, 1.35268, 1.08396, 0.815885, 0.89722, 0.849842, 0.763927, 0.563995, 0.638988, 0.550295, 0.566441, 0.490492, 0.426342, 0.362251, 0.304166, 0.333205, 0.150663, 0.191904, 0.233285, 0.246208};
                for(int i=0;i<NC;i++){
                    xEyield[i]  *= sc[ixegrr];
                    exEyield[i] *= sc[ixegrr];
                }
                xegrr[ixegrr] = new TGraphErrors(NC, xE, xEyield, zero, exEyield);
                xegrr[ixegrr]->SetMarkerSize(msz);
                xegrr[ixegrr++]->SetMarkerStyle(25);

                const int NC = 15;
                double xE[NC]       = {   0.05,    0.07,    0.09,    0.11,    0.13,    0.15,    0.17,    0.19,    0.21,    0.23,    0.25,    0.27,    0.29,     0.31,      0.37};
                double xEyield[NC]  = { 31.913, 18.4999, 14.5577, 5.97833, 4.66367, 6.03208, 2.85665, 2.27552, 3.48028, 3.36806, 2.23274, 3.05916, 2.89312, 0.767766,  0.660111};
                double exEyield[NC] = {4.40949, 3.40248, 2.81009, 1.79487, 1.52814, 1.65853, 1.17074, 1.41102, 1.60123, 1.54323,  1.3099, 1.54398, 1.45368,  0.77651,  0.661979};
                for(int i=0;i<NC;i++){
                    xEyield[i]  *= sc[ixegrr];
                    exEyield[i] *= sc[ixegrr];
                }
                xegrr[ixegrr] = new TGraphErrors(NC, xE, xEyield, zero, exEyield);
                xegrr[ixegrr]->SetMarkerSize(msz);
                xegrr[ixegrr++]->SetMarkerStyle(29);

                const int NC = 11;
                double xE[NC]       = {0.05,       0.07,    0.09,    0.11, 0.13,       0.15, 0.17,       0.19,    0.21,    0.25,    0.27};
                double xEyield[NC]  = {31.8547, 7.12175, 9.16751, 11.3278, 5.63912, 4.01009, 6.41103, 3.61058, 2.38152,  2.6378, 1.36064};
                double exEyield[NC] = {6.60093, 4.00456,  3.4402, 3.73535, 2.91223, 2.36194, 3.02606, 2.18727,  1.7424, 1.88909, 1.38614};
                for(int i=0;i<NC;i++){
                    xEyield[i]  *= sc[ixegrr];
                    exEyield[i] *= sc[ixegrr];
                }
                xegrr[ixegrr] = new TGraphErrors(NC, xE, xEyield, zero, exEyield);
                xegrr[ixegrr]->SetMarkerSize(msz);
                xegrr[ixegrr++]->SetMarkerStyle(30);
            }


            //------------------------------------------
            AddMCFromFile(TString InFile="../ktsmearing/MC/rform_rcf_D6.root"){
                //------------------------------------------

                const int NH=6;
                //double insc[NH]  = {1, 3e-1, 1e-1, 3e-2, 1e-2,  3e-3};
                //double insc[NH]  = {1, 1, 1, 1, 1, 1};
                double insc[NH]  = {1, 1e-1, 1e-2, 1e-3, 1e-4,  1e-5};
                //     double XENormsAt[NH] = {0.883929, 0.682367, 0.555999, 0.469196, 0.405844, 0.347895 };
                //     double XENorms[NH] = {0.0409001, 0.00661064, 0.00112794, 0.000202491, 3.88212e-05, 8.5431e-06 };
                double Pi0PtNormsAt[NH] = {3, 3, 3, 3, 3, 3 };
                double Pi0PtNorms[NH] = {0.0487263, 0.0113764, 0.00172654, 0.000205437, 2.24061e-05, 3.0811e-06 };
                double inSoverB  [NH]= { 1.45587, 2.0981,  2.88227, 4.04818, 4.99095, 4.2165}; 
                int marker[NH]= {29, 24, 21, 25, 22, 27};

                cout<<"Reading from "<<InFile<<endl;
                fin = new TFile(InFile);
                int NumPtt   = NH;
                TString StrName;
                for(int ih=0; ih<NumPtt; ih++){
                    StrName = "hxe"; StrName+=ih; hxEFNew[0][ih] = (TH1D *)fin->Get(StrName);
                    StrName = "hptassoc"; StrName+=ih; hpTF[0][ih] = (TH1D *)fin->Get(StrName);      
                    StrName = "hpttrigg"; StrName+=ih; hptt[ih] = (TProfile *)fin->Get(StrName);
                }
                StrName = "hfragpt"; hInclusive = (TProfile *)fin->Get(StrName);

                cout<<" histos in "<<endl;

                for(int i=0;i<NumPtt;i++){
                    sc[ixegrr]=insc[i];
                    ptt[ixegrr]=hptt[i]->GetMean();
                    //=============
                    xegrr[ixegrr] = getXEgrr_MC(i,sc[ixegrr], 1.5, 5.0); //2. 4.8
                    //=============

                    xegrr[ixegrr]->SetMarkerStyle(marker[i]);
                    ptgrr[ixegrr] = getPtgrr_MC(i,sc[ixegrr]);
                    ptgrr[ixegrr]->SetMarkerStyle(marker[i]);
                    //       normAt[ixegrr]= XENormsAt[i];
                    //       norm[ixegrr]  = XENorms[i];
                    ptNormAt[ixegrr]= Pi0PtNormsAt[i];
                    ptNorm[ixegrr]  = Pi0PtNorms[i];
                    pttl[ixegrr] = hptt[i]->GetXaxis()->GetXmin();
                    ptth[ixegrr] = hptt[i]->GetXaxis()->GetXmax();
                    ixegrr++;
                }
            }

            void AddInclusMC(){
                const int ND =  50 ;
                int good=0;
                double pt[ND], ept[ND] = {0};
                double cross[ND],  ecross[ND];
                for(int i=0;i<ND;i++) {
                    double w   = hInclusive->GetBinContent(i+1);
                    double wpt = hInclusive->GetBinCenter(i+1);
                    if(w<=0 || wpt>20) continue;
                    pt[good]    = wpt;
                    cross[good] = hInclusive->GetBinContent(i+1)*1e-2;
                    ecross[good++] = hInclusive->GetBinError(i+1)*1e-2;
                }
                pi0cross = new TGraphErrors(good, pt, cross, ept, ecross );
                pi0cross->SetMarkerStyle(20);

                InclNormAt=10;     //was at 5, InclNorm = fpow->Eval(InclNormAt);  
                fpow = new TF1("fpow","[0]/pow(x,[1])",10, 20);
                fpow->SetParameters( 498.149,  8.9598 );
                pi0cross->Fit(fpow,"RQN");
                InclNorm = fpow->Eval(InclNormAt);
                // ====================================================================
            }



            //------------------------------------------
            TGraphErrors *getXEgrr_MC(int ih, double sc=1, double pTa_l=2.0, double pTa_h=5.0){
                //------------------------------------------

                TH1D *hfxe     =  hxEFNew[0][ih];
                double sum = hfxe->GetEntries();
                double invsum = sum>0 ? 1/sum : 1;
                hfxe->Sumw2();
                hfxe->Scale(invsum*sc*5e1);

                int NC=hfxe->GetNbinsX();
                double xE[100], exE[100]={0};
                double normYF[100], enormYF[100];
                double ptt    = hptt[ih]->GetMean(), wpta, wxe, w, ew;
                double bw  = hfxe->GetBinWidth(1);
                int good=0;

                for(int ib=0; ib<NC; ib++){
                    wxe  = hfxe->GetBinCenter(ib+1);
                    w    = hfxe->GetBinContent(ib+1);
                    ew   = hfxe->GetBinError(ib+1);
                    //cout<<ib<<" "<<good<<" "<<ptt<<" xe="<<wxe<<" xe1="<<pTa_l/ptt<<" "<<w<<" "<<ew<<" "<<endl;
                    if(wxe < pTa_l/ptt || wxe > pTa_h/ptt || w<=0 ) continue;
                    xE[good]     =  wxe; //pxexe->GetBinContent(ib+1)/cos(pdphi->GetBinError(ib+1));
                    normYF[good] =  w;
                    enormYF[good]=  ew;
                    good++;
                }
                TGraphErrors *gr = new TGraphErrors(good, xE, normYF, exE, enormYF);  
                gr->SetTitle(hfxe->GetTitle());
                gr->Fit("expo","Q","",  pTa_l/ptt,  pTa_h/ptt );
                normAt[ih]= (pTa_h/ptt - pTa_l/ptt)/2;
                norm[ih]  = gr->GetFunction("expo")->Eval(normAt[ih]);
                gr->GetFunction("expo")->SetLineStyle(3);
                return gr;
            }
            //------------------------------------------
            TGraphErrors *getPtgrr_MC(int ih, double sc=1, double pTa_l=2.0, double pTa_h=5.0){
                //------------------------------------------

                TH1D *hfxe     =  hpTF[0][ih];
                int NC=hfxe->GetNbinsX();
                double xE[100], exE[100]={0};
                double normYF[100], enormYF[100];
                double ptt    = hptt[ih]->GetMean(), wpta, wxe, w, ew;
                double bw  = hfxe->GetBinWidth(1);
                int good=0;

                for(int ib=0; ib<NC; ib++){
                    wxe  = hfxe->GetBinContent(ib+1);
                    w    = hfxe->GetBinContent(ib+1);
                    ew   = hfxe->GetBinError(ib+1);
                    if(wxe < pTa_l/ptt || wxe > pTa_h/ptt || w<=0 ) continue;
                    xE[good]     =  wxe; //pxexe->GetBinContent(ib+1)/cos(pdphi->GetBinError(ib+1));
                    normYF[good] =  w;
                    enormYF[good]=  ew;
                    good++;
                }
                TGraphErrors *gr = new TGraphErrors(good, xE, normYF, exE, enormYF);  
                gr->SetTitle(hfxe->GetTitle());
                return gr;
            }


        };

