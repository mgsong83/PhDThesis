void draw_fig7(){

	TFile *fin = new TFile("data/Glauber_PathLength23_AuAu_42mb_outside_mean.root","read");


	const int NCent=6;
	const int NphiBin=6;
	float phiBin[NphiBin] = { 7.5, 22.5, 37.5, 52.5, 67.5, 82.5};

	TProfile *hL[NCent];
	TGraphErrors *grL[NCent];
	TGraphErrors *grDraw[NCent];
	for(int ic=0; ic<NCent; ic++){
		hL[ic] = (TProfile*)fin->Get(Form("hL2eccvsDphiRP%02d",ic));
		grL[ic] = new TGraphErrors(hL[ic]);
		grDraw[ic] = new TGraphErrors(NphiBin);
	}
	
	
	for(int ic=0; ic<NCent; ic++){
		for(int iphi=0; iphi<NphiBin; iphi++){
			double _dphi = phiBin[iphi] * TMath::Pi() / 180. ;// change degree<->radian
			double _l = grL[ic]->Eval(_dphi);
			cout<< phiBin[iphi] << "\t" << _l << endl; 

			grDraw[ic]->SetPoint(iphi, phiBin[iphi], _l); 
		}
	}

	TH2D *hset = new TH2D("hset" , "", 90, 0, 90, 100, 1 , 10 );
	hset->SetXTitle("#Delta #phi [degree]");
	hset->SetYTitle("Path-Length [fm]");  
	hset->SetStats(0);
	mc(1,1);
	hset->Draw();
	
	TLegend *leg = new TLegend(0.2, 0.7 , 0.4, 0.9, "" ,"brNDC" );
	leg->SetFillStyle(0); leg->SetTextSize(0.025);

	for(int ic=0; ic<NCent; ic++){
	
		grDraw[ic]->SetMarkerColor(ic+3);
		grDraw[ic]->SetMarkerStyle(ic+20);
		grDraw[ic]->Draw("same p" );	
		leg->AddEntry(grDraw[ic], Form("Cent %02d ~ %02d", ic*10, ic*10+10), "p"); 
	}	
	leg->Draw();

	c1->SaveAs("fig7.eps");

}
