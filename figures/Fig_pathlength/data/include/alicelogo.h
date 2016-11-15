void DrawLogo( Double_t xlow, Double_t ylow, Double_t xup, Double_t yup, char *tex){

   //TImage *img = TImage::Open("~/myrootmac/alice_preliminary.pdf");
   TImage *img = TImage::Open("~/myrootmac/alice_logo_transparent.png");

   if (!img) {
      printf("Could not create an image... exit\n");
      return;
   }
   TPad *l = new TPad("l","l", xlow, ylow, xup,  yup);
   // need tex here
   TLatex *la = new TLatex(0.1, 0, tex);
   la->SetTextSize(1);
   l->Draw();
   l->cd();
   img->Draw();
   la->Draw();
}

void DrawPerfor( Double_t xlow, Double_t ylow, Double_t xup, Double_t yup, Double_t texX, char *datestr=NULL){
  char date[100];
  time_t rawtime = time(NULL);
  strftime(date, 90, "%d/%b/%Y", localtime( &rawtime ));
  TString tex1 = "#color[2]{Performance Plot}";
  TString tex2 = Form("#color[2]{%s}", date) ;
  TImage *img = TImage::Open("~/myrootmac/1004067_01-Icon.jpg");

   if (!img) {
      printf("Could not create an image... exit\n");
      return;
   }
   TPad *l = new TPad("l","l", xlow, ylow, xup,  yup, kWhite);
   l->SetTopMargin(0);
   l->SetLeftMargin(0);
   l->SetBottomMargin(0);
   l->SetFrameFillColor(kWhite);
   l->Draw();
   l->cd();
   TPad *l1 = new TPad("l1", "l1", 0, 0, texX, 1, kWhite);
   l1->Draw();
   l1->cd();
   img->Draw();

   l->cd();
   // need tex here
   TLatex *la1 = new TLatex(texX, 0.5, tex1.Data(););
   la1->SetTextSize(0.3);
   la1->Draw();
   TLatex *la2 = new TLatex(texX, 0.1, tex2.Data(););
   la2->SetTextSize(0.3);
   la2->Draw();
}

void DrawPre( Double_t xlow, Double_t ylow, Double_t xup, Double_t yup, Double_t texX, char *datestr=NULL){
  char date[100];
  time_t rawtime = time(NULL);
  strftime(date, 90, "%d/%b/%Y", localtime( &rawtime ));
  TString tex1 = "#color[2]{Preliminary}";
  TString tex2 = Form("#color[2]{%s}", date) ;
  //TImage *img = TImage::Open("~/myrootmac/alice_logo_transparent.png");
  TASImage *img = new TASImage("~/myrootmac/alice_logo_transparent.png");

   if (!img) {
      printf("Could not create an image... exit\n");
      return;
   }
   TPad *l = new TPad("l","l", xlow, ylow, xup,  yup, kWhite);
   l->SetTopMargin(0);
   l->SetLeftMargin(0);
   l->SetBottomMargin(0);
   l->SetFrameFillColor(kWhite);
   l->Draw();
   l->cd();
   img->Draw();

/*
   l->cd();
   // need tex here
   TLatex *la1 = new TLatex(texX, 0.05, tex1.Data(););
   la1->SetTextSize(0.45);
   la1->Draw();
*/
}
