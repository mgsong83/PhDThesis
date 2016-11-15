void draw_fig2(){

	TFile *fin[6];
	TH1D *hL1[6];
	TH1D *hL2[6];
		
	TH2D *hset = new TH2D("hset", "", 100, 0, 2*TMath::Pi(), 100, 0, 12);
	hset->SetStats(0);
	hset->SetXTitle("#Delta#phi");
	hset->SetYTitle("PathLength L");

	mc(1,1);	
	for(int ic=0; ic<6; ic++){
		fin[ic] = new TFile(Form("data/Glauber_PathLength23_PbPb_65mb_outside_mean.root", ic), "read");

		hL1[ic] = (TH1D*)fin[ic]->Get(Form("hL2eccvsDphiRP%02d",ic ) );	
		hL2[ic] = (TH1D*)fin[ic]->Get(Form("hL3eccvsDphiPP%02d",ic ) );	

		hset->SetTitle(Form("Centrality %02d ~ %02d \%", ic*10, ic*10+10 )  );
		hset->Draw();
		hL1[ic]->SetMarkerStyle(24);
		hL1[ic]->SetMarkerSize(0.6);
		hL1[ic]->Draw("same");
		hL2[ic]->SetLineColor(kRed);
		hL2[ic]->SetMarkerColor(kRed);
		hL2[ic]->SetMarkerStyle(20);
		hL2[ic]->SetMarkerSize(0.6);
		hL2[ic]->Draw("same");
		//c1->SaveAs(Form("compareL_%d.pdf", ic)); 
	}

	mc(10,1);
	
	TH2D *href = new TH2D("href","", 10, 0, 2*TMath::Pi(), 10, 1.5, 9);
	href->SetStats(0);
	href->Draw();
	href->SetStats(0);
	href->SetXTitle("#Delta#phi  [rad]");
	href->SetYTitle("PathLength L  [fm]");


	hL1[0]->SetLineColor(kBlue);
	hL1[0]->SetMarkerColor(kBlue);
	hL1[0]->Draw("same");

	hL2[0]->SetLineColor(kBlue);
	hL2[0]->SetMarkerColor(kBlue);
	hL2[0]->Draw("same");

	hL1[4]->SetLineColor(kGreen+1);
	hL1[4]->SetMarkerColor(kGreen+1);
	hL1[4]->Draw("same");

	hL2[4]->SetLineColor(kGreen+1);
	hL2[4]->SetMarkerColor(kGreen+1);
	hL2[4]->Draw("same");

	TLegend *leg = new TLegend(0.58, 0.78 , 0.85, 0.85, "" ,"brNDC" );
	leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
	leg->AddEntry( hL1[0], "Centrality 00~10%", "l" );
	leg->AddEntry( hL1[4], "Centrality 40~50%", "l" );
	leg->Draw();

	c10->SaveAs("fig2.eps");
}
