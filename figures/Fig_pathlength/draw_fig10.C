void draw_fig10(){

	int OpenMarker[7] = { 24, 25, 26, 27, 28, 30, 32};
	int ClosedMarker[7] = {20, 21, 22, 33, 34, 29, 23 };
	//int gColor[7] = { 1, 2, 4, kCyan+4, 6, 12 8};


	TFile *fin_ALICE = new TFile("ALICE_RAA_Ldep.root","read");
	TFile *fin_PHENIX = new TFile("PHENIX_RAA_Ldep.root","read");
	cout << "load file done" << endl;

	const int NH = 1; // for only 1 harmonics
	const int NC = 6; // 6 centralities
	int Cents[NC+1] = {0, 10, 20, 30, 40, 50, 60};
	
	
	//to compare similar pt bins choose pt in manualy
	int ipt_ALICE = 16;
	double ptALICE_low = 7;
	double ptALICE_high = 8;

	int ipt_PHENIX = 2; // check it 7~8
	double ptPHENIXlow = 7;
	double ptPHENIXhigh = 8;

	// 
	TString Leg_ALICE[NC];
	TString Leg_PHENIX[NC];
	for(int ic=0; ic<NC; ic++){
		Leg_ALICE[ic] = Form("Pb+Pb 2.76TeV %0.1f<P_{T}<%0.1f",  ptALICE_low, ptALICE_high);
		Leg_PHENIX[ic] = Form("Au+Au 200GeV %1.1f<P_{T}<%1.1f",  ptPHENIXlow, ptPHENIXhigh);
	}

	TGraphErrors *gr_ALICE[NH][NC];
	TGraphErrors *gr_PHENIX[NH][NC];

	TMultiGraph *gMultiALICE[NH];
	TMultiGraph *gMultiPHENIX[NH];	

	cout << "Graphe created" << endl;

	for(int ih=0; ih<NH; ih++){
		gMultiALICE[ih] = new TMultiGraph();
		gMultiPHENIX[ih] = new TMultiGraph();
		for(int ic=0; ic<NC; ic++){
			gr_ALICE[ih][ic] = (TGraphErrors*) fin_ALICE->Get(Form("gRAA_Psi%d_L%02d%02d", ih+2, ic, ipt_ALICE)); 
			cout << "get done" << endl;
			gr_PHENIX[ih][ic] = (TGraphErrors*) fin_PHENIX->Get(Form("gRAA_Psi_%02d_L%02d_%02d", ih+2, ic, ipt_PHENIX)); 
			cout << "get done" << endl;
			
			gMultiALICE[ih]->Add( gr_ALICE[ih][ic] );
			cout << ipt_ALICE << ic << endl;
			gMultiPHENIX[ih]->Add ( gr_PHENIX[ih][ic] ) ; 
			cout << ipt_PHENIX << ic << endl;
		}


	}// graph load done.;


	TCanvas *c1 = new TCanvas("c1","c1", 450, 450);
	TH2D *hset = new TH2D("hset", "", 10, 0, 8, 10, 0, 1);
	hset->SetStats(0);
	hset->SetXTitle("Path-Length(fm)");
	hset->SetYTitle("R_{AA}" );
	hset->Draw();

	TF1 *fFitRAAL = new TF1("fFitRAAL", "[0] /  (TMath::Power( (x ) , [1])+ [0]) ",2.5, 7.2 ); 

	TLegend *leg = new TLegend(0.6, 0.5, 0.90, 0.85, "" , "brNDC" ); 
	for(int ih=0; ih<NH; ih++){
		for(int ic=0; ic<NC; ic++){

			gr_ALICE[ih][ic]->SetMarkerStyle(gOMarker[ic]);
			gr_PHENIX[ih][ic]->SetMarkerStyle(gOMarker[ic]);

			gr_ALICE[ih][ic]->SetMarkerColor(kGreen+2);
			gr_PHENIX[ih][ic]->SetMarkerColor(kBlue);

			gr_ALICE[ih][ic]->SetMarkerSize(1.4);
			gr_PHENIX[ih][ic]->SetMarkerSize(1.4);		

			//gr_ALICE[ih][ic]->Draw("p, same");
			//gr_PHENIX[ih][ic]->Draw("p, same"); 

		

		}

		gMultiALICE[ih]->Draw("AP");
		gMultiALICE[ih]->Fit("fFitRAAL","");
		gMultiALICE[ih]->GetXaxis()->SetTitle("Path-Length(fm)");
		gMultiALICE[ih]->GetYaxis()->SetTitle("R_{AA}");

		gMultiPHENIX[ih]->Draw("P");
		gMultiPHENIX[ih]->Fit("fFitRAAL", "");

	}
	
	TLegend *leg = new TLegend(0.48, 0.6, 0.85, 0.88, "" , "brNDC" );
	leg->SetBorderSize(0); leg->SetFillStyle(0);
	TGraphErrors *gr[NC];
	for(int ic=0; ic<NC; ic++){
			gr[ic] = new TGraphErrors();
			gr[ic]->SetMarkerColor(kBlack);
			gr[ic]->SetMarkerStyle(OpenMarker[ic]);
			leg->AddEntry( gr[ic], Form("Centality %02d~%02d", ic*10, ic*10+10 ), "p" );
	}
	leg->AddEntry((TObject*)0, "", "");
	leg->AddEntry( gr_ALICE[0][0], Leg_ALICE[0].Data(), "p" );
	leg->AddEntry( gr_PHENIX[0][0], Leg_PHENIX[0].Data(), "p" );
	leg->Draw()	;

	c1->SaveAs("fig10.eps");

}
