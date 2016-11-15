#define MAXD 100;

class TDataTaxi {

    public:

        TDataTaxi(){
            jtkt = new Tjtkt();
            //cout<<"Init: "<<jtkt<<endl;
        }

    public:

        Tjtkt *jtkt;
        TGraphErrors *grr[MAXD][MAXD][MAXD];  

        double ptlim[MAXD][2];
        double emptt[MAXD], empta[MAXD];
        double x, sys[3][100];

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

            const int NM=20;
            double sigN[NM],esigN[NM],sigF[NM],esigF[NM],pout[NM],epout[NM];
            double YN[NM],eYN[NM], YF[NM],eYF[NM],mptt[NM],mpta[NM],SoverB[NM],bpt[NM];
        
            //--------------------------------
            const int NC =  11 ; // pp6_Npizero_Ahadron_jt.root
            double sigN  [NC]= {0.294686, 0.277295, 0.271135, 0.261978, 0.263525, 0.278508, 0.219882, 0.23246, 0.324618, 0.257874, 0.652381}; 
            double esigN [NC]= {0.00100039, 0.00162661, 0.00251275, 0.00390201, 0.00614865, 0.00996239, 0.0104128, 0.0102172, 0.0261076, 0.0364272, 0.251221}; 

            double pout  [NC]= {0.849089, 0.779225, 0.741395, 0.664183, 0.634216, 0.610628, 0.595014, 0.558851, 0.643295, 0.467767, 0.157895}; 
            double epout [NC]= {0.00794477, 0.0103738, 0.0142399, 0.016839, 0.0237591, 0.0296114, 0.0354045, 0.0360833, 0.0613235, 0.081254, 0.0240484}; 
            double YN  [NC]= { 6.08345, 6.21155, 6.50243, 6.54369, 7.34828, 6.9952, 7.33229, 7.42498, 8.33791, 7.17684, 17.3941}; 
            double eYN [NC]= {0.0233352, 0.039999, 0.0651345, 0.0984677, 0.166894, 0.230245, 0.300327, 0.326393, 0.586513, 0.836175, 6.59262}; 
            double YF  [NC]= { 6.31923, 6.90899, 7.33065, 7.08871, 7.33799, 7.41272, 6.99124, 7.40086, 8.47174, 7.05675, 8.70344}; 
            double eYF [NC]= {0.0562168, 0.0866206, 0.132024, 0.167921, 0.258258, 0.326227, 0.43836, 0.445659, 0.734538, 0.964574, 2.95701}; 
            double mptt  [NC]= { 3.40315, 4.39576, 5.41133, 6.42175, 7.43451, 8.43425, 9.45605, 10.804, 12.8201, 14.8722, 16.7782}; 
            double mpta  [NC]= { 1.35993, 1.37184, 1.38266, 1.39298, 1.39471, 1.40244, 1.39939, 1.4235, 1.41857, 1.4308, 1.38909}; 
           //double SoverB  [NC]= { 1.62658, 2.30681, 3.30907, 4.45112, 6.94908, 8.06288, 10.9932}; 
            double bpt[2] = {1, 2, 3, 5};
            grtyp = ptt1020_pi0h_jt;
            int mtyp=27;
            ptlim[grtyp][0]=bpt[0]; ptlim[grtyp][1]=bpt[1]; 
            jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
            grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mptt, sigN, emptt, esigN); 
            grr[grtyp][tSigmaN][tReal]->SetMarkerStyle(mtyp);
            //grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mptt, sigF, emptt, esigF); 
            //grr[grtyp][tSigmaF][tReal]->SetMarkerStyle(mtyp);
            grr[grtyp][tJt][tReal]     = jtkt->GetJTy_ptt(); grr[grtyp][tJt][tReal]->SetMarkerStyle(mtyp);
            grr[grtyp][tPout][tReal] = new TGraphErrors(NC, mptt, pout, emptt, epout); grr[grtyp][tPout][tReal]->SetMarkerStyle(21);
            grr[grtyp][tZKT][tReal]  = jtkt->GetKTyz_ptt(); grr[grtyp][tZKT][tReal]->SetMarkerStyle(21);

            //--------------------------------
            const int NC =  11 ; // pp6_Npizero_Ahadron_jt.root
            const int NC =  11 ;
            double sigN  [NC]= {0.209494, 0.195288, 0.186991, 0.184245, 0.183066, 0.173298, 0.177779, 0.155226, 0.192722, 0.165345, 0.0967664}; 
            double esigN [NC]= {0.000914526, 0.00146266, 0.0023699, 0.00329391, 0.00498641, 0.00808826, 0.0116429, 0.00927141, 0.039349, 0.022052, 0.0256287}; 

            double pout  [NC]= { 1.18092, 1.07993, 0.981963, 0.903389, 0.936277, 0.892644, 0.828675, 0.829067, 0.747175, 0.434274, 0.592808}; 
            double epout [NC]= {0.0135359, 0.0171395, 0.0217662, 0.0282011, 0.0431812, 0.0563819, 0.0630941, 0.0717956, 0.121337, 0.152662, 0.0938639}; 
            double YN  [NC]= { 1.68059, 1.92388, 2.23722, 2.43303, 2.67429, 2.62122, 2.69159, 3.05797, 2.32173, 5.89193, 4.00763}; 
            double eYN [NC]= {0.0076734, 0.0147407, 0.027488, 0.0458639, 0.0769124, 0.104468, 0.143035, 0.169916, 0.330626, 0.851483, 1.36955}; 
            double YF  [NC]= { 1.27187, 1.55634, 1.88198, 2.24312, 2.28459, 2.53114, 3.1272, 3.05121, 2.38969, 4.88471, 10.7598}; 
            double eYF [NC]= {0.0132336, 0.0227478, 0.0396394, 0.0625941, 0.0971057, 0.144247, 0.197663, 0.205564, 0.326368, 0.922953, 3.91221}; 
            double mptt  [NC]= { 3.40315, 4.39576, 5.41133, 6.42175, 7.43451, 8.43425, 9.45605, 10.804, 12.8201, 14.8722, 16.7782}; 
            double mpta  [NC]= { 2.39294, 2.40002, 2.40651, 2.41658, 2.41806, 2.42806, 2.43828, 2.43975, 2.46265, 2.41156, 2.47844}; 
            //double SoverB  [NC]= { 1.62658, 2.30681, 3.30907, 4.45112, 6.94908, 8.06288, 10.9932}; 
            double bpt[2] = {1, 2, 3, 5};
            grtyp = ptt2030_pi0h_jt;
            int mtyp=26;
            ptlim[grtyp][0]=bpt[1]; ptlim[grtyp][1]=bpt[2]; 
            jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
            grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mptt, sigN, emptt, esigN); 
            grr[grtyp][tSigmaN][tReal]->SetMarkerStyle(mtyp);
            //grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mptt, sigF, emptt, esigF); 
            //grr[grtyp][tSigmaF][tReal]->SetMarkerStyle(mtyp);
            grr[grtyp][tJt][tReal]     = jtkt->GetJTy_ptt(); grr[grtyp][tJt][tReal]->SetMarkerStyle(mtyp);
            grr[grtyp][tPout][tReal] = new TGraphErrors(NC, mptt, pout, emptt, epout); grr[grtyp][tPout][tReal]->SetMarkerStyle(21);
            grr[grtyp][tZKT][tReal]  = jtkt->GetKTyz_ptt(); grr[grtyp][tZKT][tReal]->SetMarkerStyle(21);

            const int NC =  11 ;
            double sigN  [NC]= {0.158005, 0.152481, 0.143909, 0.142707, 0.135133, 0.128145, 0.146626, 0.153624, 0.116832, 0.116643, 0.543867}; 
            double esigN [NC]= {0.00087998, 0.0012621, 0.00174107, 0.00334242, 0.00455377, 0.00647253, 0.00608359, 0.01238, 0.0142216, 0.0270471, 0.769534}; 

            double pout  [NC]= { 1.4871, 1.35168, 1.30074, 1.16849, 1.27428, 1.28042, 0.827645, 1.37539, 0.857094, 0.795705, 2.36123}; 
            double epout [NC]= {0.0187687, 0.0237107, 0.0303744, 0.0415734, 0.0475148, 0.0631577, 0.0610966, 0.17966, 0.0713772, 0.193899, 2.71181}; 
            double YN  [NC]= {0.793607, 1.01789, 1.34828, 1.70736, 1.53288, 1.6085, 1.81325, 1.77316, 2.26606, 3.30994, 5.32828}; 
            double eYN [NC]= {0.00456835, 0.0095261, 0.019603, 0.0376862, 0.0493003, 0.0731402, 0.112794, 0.122224, 0.27676, 0.662298, 7.5066}; 
            double YF  [NC]= {0.494091, 0.718896, 0.988493, 1.10677, 1.58979, 1.72372, 2.06385, 1.49577, 2.33508, 0.466849, 12.5764}; 
            double eYF [NC]= {0.00567582, 0.0115206, 0.0220399, 0.0358557, 0.0667329, 0.0988372, 0.168371, 0.137943, 0.349254, 0.154471, 22.3819}; 
            double mptt  [NC]= { 3.40315, 4.39576, 5.41133, 6.42175, 7.43451, 8.43425, 9.45605, 10.804, 12.8201, 14.8722, 16.7782}; 
            double mpta  [NC]= { 3.67493, 3.71612, 3.732, 3.75969, 3.74713, 3.77945, 3.77466, 3.80343, 3.80653, 3.77914, 3.7199}; 
            //double SoverB  [NC]= { 1.62658, 2.30681, 3.30907, 4.45112, 6.94908, 8.06288, 10.9932}; 
            double bpt[2] = {1, 2, 3, 5};
            grtyp = ptt3050_pi0h_jt;
            int mtyp=25;
            ptlim[grtyp][0]=bpt[2]; ptlim[grtyp][1]=bpt[3]; 
            jtkt->Compute(NC, mpta, mptt, sigN, esigN, pout, epout);
            grr[grtyp][tSigmaN][tReal] = new TGraphErrors(NC, mptt, sigN, emptt, esigN); 
            grr[grtyp][tSigmaN][tReal]->SetMarkerStyle(mtyp);
            //grr[grtyp][tSigmaF][tReal] = new TGraphErrors(NC, mptt, sigF, emptt, esigF); 
            //grr[grtyp][tSigmaF][tReal]->SetMarkerStyle(mtyp);
            grr[grtyp][tJt][tReal]     = jtkt->GetJTy_ptt(); grr[grtyp][tJt][tReal]->SetMarkerStyle(mtyp);
            grr[grtyp][tPout][tReal] = new TGraphErrors(NC, mptt, pout, emptt, epout); grr[grtyp][tPout][tReal]->SetMarkerStyle(21);

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
            grtyp = pptt1050_pi0;
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

        }


};

