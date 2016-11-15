void draw_fig9(){



	TFile *fin = new TFile("data/Glauber_PathLength23_PbPb_65mb_outside_mean.root","read");
	TFile *fin2 = new TFile("data/Glauber_PathLength23_AuAu_42mb_outside_mean.root","read");


	const int NCent=6;
	const int NphiBin=6;
	float phiBin[NphiBin] = { 7.5, 22.5, 37.5, 52.5, 67.5, 82.5};

	TProfile *hL[NCent][2];
	TGraphErrors *grL[NCent][2];
	TGraphErrors *grDraw[NCent][2];


	for(int ic=0; ic<NCent; ic++){
		hL[ic][0] = (TProfile*)fin->Get(Form("hL2eccvsDphiRP%02d",ic));
		grL[ic][0] = new TGraphErrors(hL[ic][0]);
		grDraw[ic][0] = new TGraphErrors(NphiBin);

		hL[ic][1] = (TProfile*)fin2->Get(Form("hL2eccvsDphiRP%02d",ic));
		grL[ic][1] = new TGraphErrors(hL[ic][1]);
		grDraw[ic][1] = new TGraphErrors(NphiBin);
	}
cout << "1" << endl;	
	
	for(int ic=0; ic<NCent; ic++){
		for(int iphi=0; iphi<NphiBin; iphi++){
			double _dphi = phiBin[iphi] * TMath::Pi() / 180. ;// change degree<->radian
			double _l = grL[ic][0]->Eval(_dphi);
			grDraw[ic][0]->SetPoint(iphi, phiBin[iphi], _l); 
		}
	}
	for(int ic=0; ic<NCent; ic++){
		for(int iphi=0; iphi<NphiBin; iphi++){
			double _dphi = phiBin[iphi] * TMath::Pi() / 180. ;// change degree<->radian
			double _l = grL[ic][1]->Eval(_dphi);
			grDraw[ic][1]->SetPoint(iphi, phiBin[iphi]+2.2, _l); 
		}
	}

cout << "2" << endl;	

	
	TH2D *hset = new TH2D("hset" , "", 90, 0, 90, 100, 1 , 10 );
	hset->SetXTitle("#Delta #phi [degree]");
	hset->SetYTitle("Path-Length [fm]");  
	hset->SetStats(0);
	mc(1,1);
	hset->Draw();

cout << "3" << endl;	
	TLegend *leg = new TLegend(0.2, 0.7 , 0.4, 0.9, "" ,"brNDC" );
	leg->SetFillStyle(0); leg->SetTextSize(0.025);

	for(int ic=0; ic<NCent; ic++){
	
		grDraw[ic][0]->SetMarkerColor(gColor[ic]);
		grDraw[ic][0]->SetMarkerStyle(gCMarker[2]);
		grDraw[ic][0]->Draw("same p" );	

		grDraw[ic][1]->SetMarkerColor(gColor[ic]);
		grDraw[ic][1]->SetMarkerStyle(gOMarker[2]);
		grDraw[ic][1]->Draw("same p" );	


		leg->AddEntry(grDraw[ic][0], Form("Cent %02d ~ %02d", ic*10, ic*10+10), "p"); 
	}	
	leg->Draw();

	c1->SaveAs("fig9.eps");

}
