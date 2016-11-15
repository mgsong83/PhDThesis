// Default Configuration
int gFillStyle = 3002;
int gOMarker[6] = {kFullSquare, kFullCircle, kOpenSquare, kOpenCircle, 32, 27};
int gCMarker[6] = {kFullSquare, kFullCircle, 22, 33, 34, 35};
int gColor[6] = { kBlue, kRed, kBlue, kRed, kRed, kRed};
//int gColor[6] = { kBlue, kBlue, kBlue, kBlue, kBlue, kBlue};
double x_min = 100;
double x_max = 450;
//double y_min = -0.52 * 10e-6;
double y_min = 0 * 10e-6;
//double y_max = 0.1. * 10e-6;
double y_max = 0.47 * 10e-6;
double y_min_ratio = 0.4;
double y_max_ratio = 1.62;

double sc_input = 3.55556e-06;

double leg_y_min = 0.1;
double leg_y_max = 0.3;

TString strSaveFile = "figs_SC_comparison/SC_Comparison_SC42_ToyMC_lowM.eps";
const int ND = 5; // number of Data sets
TString datainfo[ND][3] = {
	"Results/results_Toy_Eta00_lowM_nojet_ToyMCwJet.root", 	"gr_SC_corr_1", 		"SC(4,2) : SP",
	"Results/results_Toy_Eta00_lowM_nojet_ToyMCwJet.root", 	"gr_SC_with_QC_4224", 	"SC(4,2) : QC",
	"Results/results_Toy_Eta00_lowM_jet_ToyMCwJet.root", 		"gr_SC_corr_1", 			"SC(4,2) : SP + PYTHIA Jet",
	"Results/results_Toy_Eta00_lowM_jet_ToyMCwJet.root", 		"gr_SC_with_QC_4224", 	"SC(4,2) : QC + PYTHIA Jet",
	"Results/results_Toy_Eta08_lowM_ToyMCwJet.root", "gr_SC_corr_1", "SC(4,2) : SP + PYTHIA Jet (#Delta#eta>0.8)"
};

double CentBin[8] = { 0,	5,  10, 20, 30, 40, 50, 60};
double CentMult[8] = { 450, 400, 350, 300, 250, 200, 150, 100 };


TFile *fin[ND];
TGraphErrors *gr_SC[ND];
TGraphErrors *gr_SC_ratio[ND]; // ratio to 1st graph 
TF1 *fFitftn;


//=============== Main Code =========
void Comparison_SC42_ToyMC(){
	cout << "start macro" << endl;
	LoadSCResults();
	ConvertC2M();





	fFitftn = new TF1("fFitftn", "pol0");

	TCanvas *c1 = new TCanvas("c1","c1", 600, 700);
	c1->Draw();
	TPad *upPad = new TPad("up","upper pad", 0, 0.3, 1, 1);        
	TPad *lowPad = new TPad("low","lower pad", 0, 0, 1, 0.3);      

	upPad->Draw();
	lowPad->Draw();

	TLine *line_0= new TLine( x_min, 0, x_max, 0);
	TLine *line_1= new TLine( x_min, 1, x_max, 1);

	TGraphErrors *gr_input = new TGraphErrors();
	for(int i=0; i<10; i++){
		gr_input->SetPoint( i, i*100, sc_input);
	}
	gr_input->SetLineColor(kCyan-1);
	gr_input->SetLineWidth(3);


	upPad->cd();
	upPad->SetTopMargin(0.05);
	upPad->SetBottomMargin(0);
	TH2D *href = new TH2D("href", "", 100, x_min, x_max, 100,  y_min, y_max);
	href->GetYaxis()->SetTitle("SC(m,n) ");
	href->GetYaxis()->SetTitleSize(0.04);
	href->GetYaxis()->SetTitleOffset(1.144);
	href->GetYaxis()->SetLabelSize(0.035);
	//	TGaxis::SetExponentOffset(0.01644, -0.0844, "y");
	href->SetStats(0);
	href->Draw();
	line_0->Draw("same l");

	lowPad->cd();
	lowPad->SetTopMargin(0);
	lowPad->SetBottomMargin(0.2);
	TH2D *href_ratio = new TH2D("href_ratio", "" , 100, x_min, x_max, 100, y_min_ratio, y_max_ratio);
	href_ratio->GetYaxis()->SetTitle("Ratio ");
	href_ratio->GetYaxis()->CenterTitle(1);
	href_ratio->GetXaxis()->SetTitle("Multiplicity");
	href_ratio->GetYaxis()->SetTitleSize(0.09);
	href_ratio->GetXaxis()->SetTitleSize(0.08);
	href_ratio->GetYaxis()->SetTitleOffset(0.45);
	href_ratio->GetXaxis()->SetTitleOffset(0.95);
	href_ratio->GetYaxis()->SetLabelSize(0.075);
	href_ratio->GetXaxis()->SetLabelSize(0.08);
	href_ratio->SetStats(0);
	href_ratio->Draw();
	line_1->Draw("same l");
	//
	//Fill Graph
 
	upPad->cd();
	gr_input->Draw("same l");
	TLegend *leg = new  TLegend(0.14, leg_y_min, 0.55, leg_y_max, "SC(3,2) with ToyMC", "brNDC");
	leg->SetFillStyle(0); 
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);
	leg->AddEntry( gr_input, "Input SC value", "l");

	for(int i=0; i<ND; i++){
		Shift_graph_X( gr_SC[i], -1+i/2. );
		gr_SC[i]->Draw("same p");
		leg->AddEntry( gr_SC[i], datainfo[i][2].Data(), "p");
	}
	leg->Draw();

	lowPad->cd();
	for(int i=0; i<ND; i++){
			gr_SC_ratio[i] = CalculateRatio( gr_SC[i] ,  gr_input, 0); // gr2 / gr1
			gr_SC_ratio[i]->SetMarkerColor( gColor[i] );
			gr_SC_ratio[i]->SetMarkerStyle( gOMarker[i] );
			gr_SC_ratio[i]->SetMarkerSize( 1.1 );
			gr_SC_ratio[i]->SetLineColor( gColor[i] );
			gr_SC_ratio[i]->Draw("same p");
		
			gr_SC_ratio[i]->Fit("fFitftn", "n");

	}	


	TLegend *leg_low = new  TLegend(0.15, 0.7, 0.55, 0.9, "", "brNDC");
	leg_low->SetFillStyle(0); 
	leg_low->SetBorderSize(0);
	leg_low->SetTextSize(0.08);
	leg_low->Draw();


	c1->SaveAs( strSaveFile.Data() );

}



// Load SC results //
void LoadSCResults(){
		cout << "read data" << endl;
		for(int i=0; i<ND; i++){
				cout << i << endl;
				fin[i] = new TFile( datainfo[i][0].Data(), "read" );
				fin[i]->Print();
		
				gr_SC[i] = (TGraphErrors*)fin[i]->Get( datainfo[i][1].Data() );
				gr_SC[i]->Print();
		}
		cout << "read done" << endl;
}
// Convert Centrality to Multiplicity 
void ConvertC2M(){
		cout << "Converting cent < = > multi " << endl;
		TGraphErrors *gr_tmp = new TGraphErrors();
		for(int i=0; i<ND; i++){
			int np = gr_SC[i]->GetN();
			for(int ip=0; ip<np; ip++){
				double x = (CentMult[ip] + CentMult[ip+1])/2.;
				double y = gr_SC[i]->GetY()[ip];
				double ex = 0;
				double ey = gr_SC[i]->GetEY()[ip];
			
				gr_tmp->SetPoint( ip, x, y );
				gr_tmp->SetPointError( ip, ex, ey);		
			}

			gr_SC[i] = (TGraphErrors*)gr_tmp->Clone();
			gr_SC[i]->SetMarkerStyle( gOMarker[i] );
			gr_SC[i]->SetMarkerColor( gColor[i] );
		}

}




// CalculateRatio 
TGraphErrors *CalculateRatio( TGraphErrors* gr1, TGraphErrors* gr2, double xshift){
		int NPoint = gr1->GetN();
		TGraphErrors *gr_ratio = new TGraphErrors( NPoint);
		TGraph ger( gr2->GetN(), gr2->GetX(), gr2->GetEY() ); // << Err estimation of gr2
		for(int i=0; i<NPoint; i++){
				double x = gr1->GetX()[i];
				double y1 = gr1->GetY()[i];
				double ey1 = gr1->GetEY()[i];  // << estimation of gr2's y value for x1
				double y2 = gr2->Eval(x);
				double ey2 = ger.Eval(x);
				double ratio = TMath::Abs( y1 / y2 );
				gr_ratio->SetPoint( i, x+xshift, ratio);
				gr_ratio->SetPointError( i, 0, ratio * TMath::Sqrt( ey1*ey1/y1/y1 + ey2*ey2/y2/y2 ));

		}
		return gr_ratio;
}
















