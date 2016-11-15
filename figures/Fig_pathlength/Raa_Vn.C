#define Do_DEBUG
#ifdef Do_DEBUG
#define DEBUG(x) cout << "DEBUGE " << x << endl
#else
#define DEBUG(x)
#endif 

#include "TF1.h"
#include "TFile.h"
#include "data/include/JmxnTPad.h" // for multi canvas
#include "data/include/alicelogo.h"
int gMarker[11] = {20,24,21,25,22,26,27,28,29,30,31};
int OpenMarker[6] = {24, 25, 26, 27, 28, 30};
int gColor[11] = {1,2,4,kCyan+4,6,12,8,9,kViolet,kCyan+1,3};

const int NC = 6; // nunmber of cent  // just test now
const int NH = 1; // 0 :v2, 1: v3
int cent[NC+1] = {0,10,20,30,40,50,60};

// Load Ldep from Glabuer
TProfile *hLeccvsDphiRP[NH][NC];
//TGraphErrors gLvsDphi[NC];
TString Plane = "RP"; // or "RP"
double lL=0.5;
double hL=8.5;


// -------------------------------------------------------------------- 
void LoadGlabuer(){
		cout <<"Loading Glaber L with "<< Plane.Data() << endl;
//		TFile fin_Ldep("data/Glauber_PathLength23_PbPb_65mb.root");
//		TFile fin_Ldep("data/Glauber_PathLength23_PbPb_65mb_bnew.root");
//		TFile fin_Ldep("data/Glauber_PathLength23_PbPb_65mb_outside.root");
		TFile fin_Ldep("data/Glauber_PathLength23_PbPb_65mb_outside_mean.root");
		gROOT->cd();
		for(int ih=0; ih<NH; ih++){
				for(int ic=0;ic<NC;ic++ ){
						TString hname = Form("hL%deccvsDphi%s%02d",ih+2,Plane.Data(), ic);
						printf("hL%deccvsDphi%s%02d",ih,Plane.Data(), ic);
						hLeccvsDphiRP[ih][ic] = (TProfile*)fin_Ldep.Get(hname.Data())->Clone();
				}
		}
		cout << "load done" << endl;
} // get TGraph  Path-Length <-> Angle
// -------------------------------------------------------------------- 


// Main Code ----------------------------------------------------------
void Raa_Vn(){
	// OPEN input files
	 
	DEBUG("file loading");
	TFile *fin_Raa = new TFile("data/ALICE_data.root","read");
	TFile *fin_Vn = new TFile("data/cms_vn_hadrons.root","read");
	TFile *fout = new TFile("ALICE_RAA_Ldep.root","recreate");
	DEBUG("file load done");

	const int NEP = 6; // number of angle bin
	const int NPT= 12;
	//int IPT[NPT+1] =   {10,11,12,13,14,15,16,17,18,19,20,21,22};
	//int IPT[NPT+1] =   {20,21,22,23,24,25,26,27,28,29,30,31,32};
	//int IPT[NPT+1] =   {30,31,32,33,34,35,36,37,38,39,40,41,42};
	//int IPT[NPT+1] =   {40,41,42,43,44,45,46,47,48,49,50,51,52};
	//int IPT[NPT+1] =   {50,51,52,53,54,55,56,57,58,59,60,61,62};
//	int IPT[NPT+1] =   {37,38,39,40,41,42,43,44,45,46,47,48,49};
	int IPT[NPT+1] =   {28,29,30,31,32,33,34,35,36,37,38,39,40};
	double dphi[NEP] = {7.5, 22.5, 37.5, 52.5, 67.5, 82.5}; // angle?? (why?)
//	double dphi[NEP] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};

	

	for(int i=0; i<NEP; i++){
		dphi[i] = (TMath::Pi() / 180.) * dphi[i]; 
	} // change angle unit (degree -> radian) 
	DEBUG("angle changed");

	// read from data txt file in TGraphErrors or TH1D ??
	TH1D *hRAAInclusive[NC];
	TString sCENTs[NC];
	for(int ic=0; ic<NC; ic++){
		hRAAInclusive[ic] =(TH1D*)fin_Raa->Get(Form("hRAAPt_%02d",ic));
		sCENTs[ic] = Form("%d-%d %%",cent[ic], cent[ic+1]);
	} // load RAA histo/graphi
	DEBUG("RAA histo load done");

	TF1 *fVN[NH][NC];
	TH1D *hVN_err[NH][NC]; // will be used when fit is ok.
	for(int ih=0; ih<NH; ih++){
			for(int ic=0; ic<NC; ic++){
					fVN[ih][ic] = (TF1*)fin_Vn->Get(Form("fit_v%d_C%02d%02d_cmshadrons",ih+2,cent[ic], cent[ic+1]));
					hVN_err[ih][ic] = (TH1D*)fin_Vn->Get(Form("fit_v%d_C%02d%02d_error",ih+2,cent[ic], cent[ic+1]));
				
			} // load Vn graphi // get Vn from fit fucntion
	}
	DEBUG("read vn data done");
	TH1D *hRAA_VN[NH][NC][NEP];

	for(int ih=0; ih<NH; ih++){
			for(int ic=0; ic<NC; ic++){
					for(int iep=0; iep<NEP; iep++){
							hRAA_VN[ih][ic][iep] = GetRAAEPUnfold( hRAAInclusive[ic], fVN[ih][ic], dphi[iep], ih+2 );
							//hRAA_VN[ih][ic][iep] = GetRAAEPUnfold( hRAAInclusive[ic], hVN_err[ih][ic], dphi[iep] );
					}
			} // calculate Raa * vn(ratio) value.
	} // 
	DEBUG("RAA * vn calculated");


	// Making New TGraphErrors 

	
	TF1  *ffFlowVN[NH];
	for (int ih=0; ih<NH; ih++){	
		if(ih==0) ffFlowVN[ih] = new TF1("ffFlowV2",FlowFitFunctionV2,0.-0.02,TMath::Pi()/2+0.02,2);
		if(ih==1) ffFlowVN[ih] = new TF1("ffFlowV3",FlowFitFunctionV3,0.-0.02,TMath::Pi()/2+0.02,2);
		ffFlowV2->SetLineColor(9);ffFlowV2->SetLineStyle(2);
	}
	DEBUG("ffFlow ready");

	TGraphErrors *gRAA_Dphi[NH][NC][NPT];
	TGraphErrors *gRAA_L[NH][NC][NPT];
	TString sPT[NH][NC][NPT];

	DEBUG("TGraphs ready");

	for(int ih=0; ih<NH; ih++){
			for(int ic=0; ic<NC; ic++){
					for(int ipt=0; ipt<NPT; ipt++){
							gRAA_Dphi[ih][ic][ipt] = new TGraphErrors(NEP);
							gRAA_Dphi[ih][ic][ipt]->SetName(Form("gRAA_Dphi%02d%02d%02d", ih, ic, ipt));
							sPT[ih][ic][ipt] = Form("%.2f < p_{T} < %.2f", hRAA_VN[ih][ic][0]->GetBinLowEdge(IPT[ipt]), hRAA_VN[ih][ic][0]->GetBinLowEdge(IPT[ipt])+hRAA_VN[ih][ic][0]->GetBinWidth(IPT[ipt]    ));

							for(int iep=0; iep<NEP; iep++){
									double raa_ep = hRAA_VN[ih][ic][iep]->GetBinContent( IPT[ipt] );
									double raa_ep_err = hRAA_VN[ih][ic][iep]->GetBinError( IPT[ipt] ); 
									gRAA_Dphi[ih][ic][ipt]->SetPoint(iep, dphi[iep], raa_ep );
									gRAA_Dphi[ih][ic][ipt]->SetPointError(iep, 0, raa_ep_err);  
							}
					}
			}
	}
	DEBUG("TGraph filled");
	LoadGlabuer(); // read Path-length vs phi info.
	DEBUG("GLAUBER simulation loaded");

	for(int ih=0;ih<NH; ih++){
			for(int ic=0;ic<NC;ic++) {
					for(int ipt=0;ipt<NPT;ipt++) {
							TGraphErrors * gr = new TGraphErrors(gRAA_Dphi[ih][ic][ipt]->GetN());
							TGraphErrors * gr2 = gRAA_Dphi[ih][ic][ipt];
							gr->SetName(Form("gRAA_Psi%d_L%02d%02d",ih+2 ,ic, ipt));
							for( int ip=0;ip<gr2->GetN();ip++){ // ip is number of phi angle
									double _dphi = gr2->GetX()[ip];
									double dphi_bin = hLeccvsDphiRP[ih][ic]->GetXaxis()->FindBin( _dphi );
									double _l = hLeccvsDphiRP[ih][ic]->GetBinContent(dphi_bin);
									double errx = hLeccvsDphiRP[ih][ic]->GetBinWidth(dphi_bin); 
									double erry = hLeccvsDphiRP[ih][ic]->GetBinError(dphi_bin); 
									double _lerr = TMath::Sqrt(errx*errx + erry*erry); 
									gr->SetPoint(ip, _l, gr2->GetY()[ip]);
									gr->SetPointError(ip, _lerr, gr2->GetErrorY( ip ) );
//									if(ic==2 && ipt==0) cout << ic <<"\t"<<_dphi <<"\t" << _l << "\t" << _lerr << endl;
							}//ip
							gRAA_L[ih][ic][ipt]=gr;
					}//ipt
			}//ic
	} // ih
	//----------------------- 
	DEBUG("path length calculated");

	///////////////////////////////////////////////////////////////////////
	//// Raa as a function of dphi in pt bins
	///////////////////////////////////////////////////////////////////////
	int iC = 2; // for drawing
	int iH = 2; // for drawing
	int ih = iH-2; // ih=0 : v2, ih=1 : v3;  

	const char *toptitle = "#sqrt{s_{NN}}= 2.76 TeV Pb+Pb";
	TLatex *toptex = new TLatex(0.30,0.92,Form("%s, Centrality %d-%d %%",toptitle,cent[iC],cent[iC+1]));
	toptex->SetTextSize(0.03);
	TLatex *toptex2 = new TLatex(0.40,0.92,toptitle);
	toptex2->SetTextSize(0.03);

	const int nx=3;
	const int ny=3; // last ptt bin to be shown
	int PTbin1[9] = {3,4,5,6,7,8,9,10,11};

	double lrX[nx]={-0.2};      double urX[nx]={TMath::Pi()/2 + 0.2};
	double lraa = 0, hraa = 1.;
	double lrY[nx][ny];
	double urY[nx][ny];
	double scale[nx][ny];
	for (int ix=0;ix<nx;ix++){
		lrX[ix] = -0.22;
		urX[ix] = TMath::Pi()/2 + 0.22;
		for (int iy=0;iy<ny;iy++){
			lrY[ix][iy]= lraa;
			urY[ix][iy]= hraa;
		}
	}
	DEBUG("draw par setted");

	// Drowing Raa - dPhi 3x3 Canvas ------------------------------- 
	JmxnTPad  *pad1 =new JmxnTPad(nx, ny, lrX, urX, lrY[0], urY[0]);
	pad1->SetOptionsLabelAxisX(0.008,0.58, 42);
	pad1->SetOptionsLabelAxisY(0.008,0.58, 42);
	pad1->SetSizeOfLabelPadsAxisY(0.05,0.017);
	pad1->SetTitleX(Form("| #Phi_{hadron} - #psi_{%d} | [rad]",ih+2),1,0.04,42,1.4);
	pad1->SetTitleY("R_{AA}",1,0.04,42,1.6);
	TString pLeg[nx][ny];
	int count = 0;
	for(int iy=0;iy<ny;iy++) {
		for(int ix=0;ix<nx;ix++) {
			pLeg[ix][iy]=sPT[ih][iC][PTbin1[count]];
			pad1->SetCaption(ix, iy, pLeg[ix][iy]);  // puting caption to each pad
			pad1->SetOptionsCaption(20, 43, 0.9,4./5);  // puting caption to each pad
			count++;
		}
	}
	pad1->SetYCaption(8.8/10.);
	pad1->Draw();
	toptex->Draw();

	count = 0;
	for(int iy=0;iy<ny;iy++) {
		for(int ix=0;ix<nx;ix++) {
			TPad* p = pad1->GetPad(ix,iy); //include graph to the upper left pad
			p->cd();
			p->SetGridx(0);
			p->SetGridy(1);
			gRAA_Dphi[ih][iC][PTbin1[count]]->SetMarkerStyle(20);
			gRAA_Dphi[ih][iC][PTbin1[count]]->Draw("p,same");
			gRAA_Dphi[ih][iC][PTbin1[count]]->Fit(Form("ffFlowV%d",ih+2),"R");
			p->Update();
			count++;
		}
	}
	c_9997->SaveAs("fig3.eps");
	// Drowing Raa - dPhi 3x3 Canvas ------------------------------- done. 

	// Draw L
	const int nx2=3;
	const int ny2=3; // last ptt bin to be shown
	double lrY2[nx2][ny2];
	double urY2[nx2][ny2];
	double lrX2[nx2]; 
	double urX2[nx2];
	for (int ix=0;ix<nx2;ix++){
		lrX2[ix] = lL;
		urX2[ix] = hL;
		for (int iy=0;iy<ny2;iy++){
			lrY2[ix][iy]= lraa;
			urY2[ix][iy]= hraa;
		}
	} // initiationg done. 

	JmxnTPad  *pad2 =new JmxnTPad(nx2, ny2, lrX2, urX2, lrY2[0], urY2[0]);
	pad2->SetOptionsLabelAxisX(0.008,0.58, 42);
	pad2->SetOptionsLabelAxisY(0.008,0.58, 42);
	pad2->SetSizeOfLabelPadsAxisY(0.05,0.017);
	pad2->SetTitleX("L_{#epsilon} [fm]",1,0.04,42,1.4);
	pad2->SetTitleY("R_{AA}",1,0.04,42,1.6);
	TString pLeg2[nx2][ny2];
	int count = 0;
	for(int iy=0;iy<ny2;iy++) {
		for(int ix=0;ix<nx2;ix++) {
			pad2->SetCaption(ix, iy, sPT[ih][iC][PTbin1[count]].Data() );  // puting caption to each pad
			pad2->SetOptionsCaption(20, 43, 0.9,4./5);  // puting caption to each pad
			count++;
		}
	}
	pad2->SetYCaption(8.8/10.);
	//pad1->SetLogY(1);
	pad2->Draw();
	toptex2->Draw();
	//DrawPre(0.20,0.80,0.25,0.84,0.35);
	//pad drown

	count = 0;
	for(int iy=0;iy<ny2;iy++) {
		for(int ix=0;ix<nx2;ix++) {
			TPad* p2 = pad2->GetPad(ix,iy); //include graph to the upper left pad
			p2->cd();
			p2->SetGridx(0);
			p2->SetGridy(1);

			/* Draw 2nd order */
			for( int ic=0;ic<NC;ic++){
				gRAA_L[0][ic][PTbin1[count]]->SetMarkerStyle(OpenMarker[ic]);
				//gRAA_L[0][ic][PTbin1[count]]->SetMarkerColor(gColor[ic]);
				gRAA_L[0][ic][PTbin1[count]]->SetMarkerColor(gColor[ic]);
				gRAA_L[0][ic][PTbin1[count]]->SetMarkerSize(1.1);
				gRAA_L[0][ic][PTbin1[count]]->Draw("p,same");
			} 
			/* Draw 3rd order 
			for( int ic=0;ic<NC;ic++){
				gRAA_L[1][ic][PTbin1[count]]->SetMarkerStyle(gMarker[ic]);
				gRAA_L[1][ic][PTbin1[count]]->SetMarkerColor(gColor[ic]);
				//gRAA_L[1][ic][PTbin1[count]]->SetMarkerColor(gColor[ic]);
				gRAA_L[1][ic][PTbin1[count]]->SetMarkerSize(1.1);
				gRAA_L[1][ic][PTbin1[count]]->Draw("p,same");
			} */
			if( count == 0 ){
				TLegend *leg = new TLegend(0.6,0.5,1.10,0.80,"","brNDC");
				leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.07);
				for( int ic=0;ic<NC;ic++){
					leg->AddEntry(gRAA_L[ih][ic][PTbin1[count]], sCENTs[ic].Data(), "p");
				}
				leg->Draw();
				//alipre->Draw();
			}
			p2->Update();
			count++;
		}
	}
	c_1712->SaveAs("fig4.eps");

	//save graph
	fout->cd();
	for(int ih=0; ih<NH; ih++){
		for(int ic=0; ic<NC; ic++){
			for(int ipt=0; ipt<NPT; ipt++){
				gRAA_Dphi[ih][ic][ipt]->Write();
				gRAA_L[ih][ic][ipt]->Write();
			}
		}
	}
	fout->Close();	

	// drawn done.
}

// function for Get RAA(Pt,Phi) from inclusive RAA
TH1D* GetRAAEPUnfold(TH1D *hidA, TH1D *fit_v2_err , double dphi)
{
	int    NX =  hidA->GetNbinsX();
	TH1D *hid;
	hid = (TH1D*) hidA->Clone();    hid->Reset();
	for (int i = 1;i <= NX;i++){ // for histo bin start from 1 
		double pt    	= hidA->GetBinCenter(i) ;
		double Raa 	 	= hidA->GetBinContent(i) ;
		double eRaa	 	= hidA->GetBinError(i) ;
		double fv2_pt 	= fit_v2_err->GetXaxis()->FindBin( pt );
		double v2  		= fit_v2_err->GetBinContent( fv2_pt );
		double ev2  	= fit_v2_err->GetBinError( fv2_pt );
		double RAAEP 	= Raa * ( 1 + 2*v2* TMath::Cos(2*dphi) );
		double hraa = (Raa+eRaa) * ( 1 + 2*(v2+ev2) * TMath::Cos(2*dphi) );
		double eRAAEP	= eRaa*( 1+ 2* v2 * TMath::Cos(2*dphi) ) + ev2( 2*Raa*TMath::Cos(2*dphi));
		hid->SetBinContent(i, RAAEP );
		hid->SetBinError(i, eRAAEP);
		cout << "error : "  <<i << "\t" << eRaa << "\t" << ev2 << endl;
	}
	return hid;
	// Raa(in) - Raa(out) = 4*v2*Raa
	// Raa(in) = Raa(1+2v2)
	// Raa(out) = Raa(1-2v2)
}

TH1D* GetRAAEPUnfold(TH1D *hidA, TF1 *fit_vn , double dphi, int ih)
{
	int    NX =  hidA->GetNbinsX();
	TH1D *hid;
	hid = (TH1D*) hidA->Clone();    hid->Reset();
	for (int i = 1;i <= NX;i++){ // for histo bin start from 1 
		double pt    	= hidA->GetBinCenter(i) ;
		double Raa 	 	= hidA->GetBinContent(i) ;
		double eRaa	 	= hidA->GetBinError(i) ;
		double vn  		= fit_vn->Eval(pt);
		double evn  	= 0;
		double RAAEP 	= Raa * ( 1 + 2*vn* TMath::Cos(ih*dphi) );
		double eRAAEP	= eRaa*(1+2*vn* TMath::Cos(ih*dphi) ) + evn*( 2*Raa* TMath::Cos(ih*dphi) );
		hid->SetBinContent(i, RAAEP );
		hid->SetBinError(i, eRAAEP);
		//cout << "error : "  <<i << "\t" << eRAAEP << endl;
	}
	return hid;
	// Raa(in) - Raa(out) = 4*v2*Raa
	// Raa(in) = Raa(1+2v2)
	// Raa(out) = Raa(1-2v2)
}



double FlowFitFunctionV2(double *x,double *par){
	static double PI = TMath::Pi();
	return par[0]*(1+2*(par[1]*cos(2*x[0])));
}


double FlowFitFunctionV3(double *x,double *par){
	static double PI = TMath::Pi();
	return par[0]*(1+2*(par[1]*cos(3*x[0])));
}





