void draw_fitresult(){


	TFile *fin[2];
	fin[0] = new TFile("ALICE_RAA_Ldep.root", "read");
	fin[1] = new TFile("PHENIX_RAA_Ldep.root", "read");


	TGraphErrors *grFitresult[2];
	TGraphErrors *grFitChi2NDF[2];

	for(int i=0; i<2; i++){
		grFitresult[i] = (TGraphErrors*) fin[i]->Get("grFitResult");
		grFitChi2NDF[i] = (TGraphErrors*) fin[i]->Get("grFitChi2NDF");
	}// read file and load graph done

	grFitresult[0]->SetMarkerStyle(22);
	grFitresult[0]->SetMarkerColor(kBlue);
	grFitresult[0]->SetMarkerSize(0.5);
	grFitresult[1]->SetMarkerStyle(22);
	grFitresult[1]->SetMarkerColor(kRed);
	grFitresult[1]->SetMarkerSize(0.5);

	grFitChi2NDF[0]->SetMarkerStyle(22);
	grFitChi2NDF[0]->SetMarkerColor(kBlue);
	grFitChi2NDF[0]->SetMarkerSize(0.5);
	grFitChi2NDF[1]->SetMarkerStyle(22);
	grFitChi2NDF[1]->SetMarkerColor(kRed);
	grFitChi2NDF[1]->SetMarkerSize(0.5);



	mc(1,1);
	TH2D *h2 = new TH2D("FitResult","Parameter N vs Pt",100, 0, 20, 10, 0, 3);
	h2->SetStats(0);
	h2->Draw();
	grFitresult[0]->Draw("p, same");
	grFitresult[1]->Draw("p, same");

	mc(2,1);
	TH2D *h3 = new TH2D("FitResult2","Chi2/NDF vs pt",100, 0, 20, 10, -1, 5);
	h3->SetStats(0);
	h3->Draw();
	grFitChi2NDF[0]->Draw("p, same");
	grFitChi2NDF[1]->Draw("p, same");

}
