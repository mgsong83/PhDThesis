class Tjtkt {

private:
  double jt[50], ejt[50];
  TGraphErrors *jgrr_pta, *kgrrz_pta, *jgrr_ptt, *kgrrz_ptt;
  double Con;
  
public:
  Tjtkt(){;}
  
  void Compute(const int NC, double mpta[],  double mptt[],  
	       double sigN[], double esigN[],
	       //double sigF[], double esigF[],
	       double pout[], double epout[] ){
    double jt[NC], ejt[NC], mpt_t, mpt_a;
    double empta[NC]={0}, zkt[NC], ezkt[NC];
    
    for(int i=0; i<NC; i++) {
      mpt_a = mpta[i];                                //associated pt
      mpt_t = mptt[i];
      empta[i]=0.1; //for painting purposes ?
      double fact = mpt_a*mpt_t/sqrt(mpt_a*mpt_a+mpt_t*mpt_t)*sqrt(2); //sqrt(2) because rms(jty)->rms(jt)
      jt[i]  = sigN[i]*fact;
      ejt[i] = esigN[i]*fact;
      //kt ===========================
      //double mpout=mpta[i]*sin(sigF[i]);//pout[i];
      double mpout=pout[i];
      double empout=epout[i];
      double xh = mpt_a/mpt_t;
      double jty = jt[i]/sqrt(2);
      double ejty= ejt[i]/sqrt(2);
      zkt[i] = sqrt(mpout*mpout-(1+xh*xh)*jty*jty)/xh;
      
      double nom = pow(empout*mpout,2)+pow(jty*ejty*(1+xh*xh),2);
      double denom = pow(xh,2)*(mpout*mpout-jty*jty*(1+xh*xh));
      ezkt[i] = sqrt(nom/denom);
//       cout<<"pta="<<mpt_a<<" ptt="<<mpt_t<<" xh="<<xh<<" jty="<<jty<<" +- "<<ejty<<endl;
//       cout<<"zkt="<<zkt[i]<<" +- "<<ezkt[i]<<" pout="<<pout[i]<<" +- "<<epout[i]<<" nom="<<nom<<" denom="<<denom<<endl<<endl;;
    }
    jgrr_pta  = new TGraphErrors(NC, mpta, jt,  empta, ejt);
    kgrrz_pta = new TGraphErrors(NC, mpta, zkt, empta, ezkt);
    jgrr_ptt  = new TGraphErrors(NC, mptt, jt,  empta, ejt);
    kgrrz_ptt = new TGraphErrors(NC, mptt, zkt, empta, ezkt);
 }
  TGraphErrors *GetJTy_pta() {return jgrr_pta;}
  TGraphErrors *GetKTyz_pta(){return kgrrz_pta;}
  TGraphErrors *GetJTy_ptt() {return jgrr_ptt;}
  TGraphErrors *GetKTyz_ptt(){return kgrrz_ptt;}
};
