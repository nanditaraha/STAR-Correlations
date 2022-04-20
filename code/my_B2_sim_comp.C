//Here we compare diff sim of B2 from R2
#include "fluct_common.h"
#include "TString.h"
#include "TLegend.h"
#include "TMath.h"

double ylimits[2];
double xlimits[2];
void RangeFinderH(int, TH1D*);
void RangeFinderG(int, TGraph*);
void RangeFinderGE(int, TGraphErrors*);
void RangeCheck();
TString outfileP = TString("ps/B2_c16_ampt_energy.pdf");

void my_B2_sim_comp(int ds=16){

  TFile *f1 = new TFile(Form("root_files/B2_int_c16_urqmdLN_%d_cuts1.root",ds),"read"); 
  TFile *f2 = new TFile(Form("root_files/B2_int_c16_ampt_%d_cuts1.root",ds),"read"); 
  TFile *f3 = new TFile(Form("root_files/B2_int_c16_hijing_%d_cuts1.root",ds),"read"); 
  
  TGraphErrors *gintB2cp_urqmd[3] = {0};
  TGraphErrors *gintB2cp_ampt[3] = {0};
  TGraphErrors *gintB2cp_hijing[3] = {0};
  TGraphErrors *gintB2cm_urqmd[3] = {0};
  TGraphErrors *gintB2cm_ampt[3] = {0};
  TGraphErrors *gintB2cm_hijing[3] = {0};

  for(int ipart=0; ipart<3;ipart++){
    
    gintB2cp_urqmd[ipart] = (TGraphErrors*)f1->Get(Form("gint_B2cp_%d_urqmd",ipart));
    gintB2cm_urqmd[ipart] = (TGraphErrors*)f1->Get(Form("gint_B2cm_%d_urqmd",ipart));

    gintB2cp_ampt[ipart] = (TGraphErrors*)f2->Get(Form("gint_B2cp_%d_ampt",ipart));
    gintB2cm_ampt[ipart] = (TGraphErrors*)f2->Get(Form("gint_B2cm_%d_ampt",ipart));
    
    gintB2cp_hijing[ipart] = (TGraphErrors*)f3->Get(Form("gint_B2cp_%d_hijing",ipart));
    gintB2cm_hijing[ipart] = (TGraphErrors*)f3->Get(Form("gint_B2cm_%d_hijing",ipart));
  }
    
  gROOT->SetStyle("Modern");
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetPadTopMargin(0.005);
  gStyle->SetPadBottomMargin(0.08);
  gStyle->SetPadLeftMargin(0.12);
  int ican=-1,iline=-1,ivline=-1;
  int itxt=-1;
  //TCanvas *ccan[1000];
  TLatex *txt[1000];
  for (int i=0;i<1000;i++){
    txt[i]      = new TLatex();
    txt[i]      ->SetTextSize(0.06);
    txt[i]      ->SetTextAlign(12);
    txt[i]      ->SetNDC();
  }
  //int ican=0;
  TLine* lhor   = new TLine(0,0,1,1); lhor->SetLineWidth(2); lhor->SetLineColor(16); lhor->SetLineStyle(2);
  TCanvas *ccan    = new TCanvas("","",20,30+20,900,650);
  ccan->Divide(3,2,0.0001,0.0001);
  for (int ipart=0;ipart<3;ipart++){
    ccan->cd(ipart+1);
    RangeFinderG(0,gintB2cp_urqmd[ipart]);
    RangeFinderG(1,gintB2cp_ampt[ipart]);
    RangeFinderG(2,gintB2cp_hijing[ipart]);
    RangeCheck();
    ylimits[0] = -0.1; ylimits[1] = 1.2;
    //if (ipart==2) {ylimits[0] = -0.5, ylimits[1] = 2;}
    gintB2cp_urqmd[ipart]->Draw("AP");
    gintB2cp_urqmd[ipart]->SetMinimum(ylimits[0]);
    gintB2cp_urqmd[ipart]->SetMaximum(ylimits[1]);
    gintB2cp_ampt[ipart]->Draw("P");
    gintB2cp_urqmd[ipart]->Draw("P");
    gintB2cp_hijing[ipart]->Draw("P");
    lhor->DrawLine(0,0.0,16.5,0.0); lhor->DrawLine(0,1.0,16.5,1.0);
    if (ipart==0){
      TLegend *legend1= new TLegend(0.4,0.71,0.7,0.52);
      legend1->SetTextFont(43);
      legend1->SetTextSize(15);
      legend1->SetFillStyle(0);
      legend1->AddEntry(gintB2cp_urqmd[ipart],"Urqmd");
      legend1->AddEntry(gintB2cp_ampt[ipart],"Ampt");
      legend1->AddEntry(gintB2cp_hijing[ipart],"Hijing");
      legend1->Draw();
    }
  }

  for (int ipart=0;ipart<3;ipart++){
    ccan->cd(ipart+4);
    RangeFinderG(0,gintB2cm_urqmd[ipart]);
    RangeFinderG(1,gintB2cm_ampt[ipart]);
    RangeFinderG(2,gintB2cm_hijing[ipart]);
    RangeCheck();
    
    ylimits[0] = -0.1; ylimits[1] = 1.2;
    //if (ipart==1) {ylimits[0] = -1.1; ylimits[1] = 1.4;}
    //if (ipart==2) {ylimits[0] = -0.1, ylimits[1] = 2.5;}
    gintB2cm_urqmd[ipart]->Draw("AP");
    gintB2cm_urqmd[ipart]->SetMinimum(ylimits[0]);
    gintB2cm_urqmd[ipart]->SetMaximum(ylimits[1]);
    gintB2cm_ampt[ipart]->Draw("P");
    gintB2cm_urqmd[ipart]->Draw("P");
    gintB2cm_hijing[ipart]->Draw("P");
    lhor->DrawLine(0,0.0,16.5,0.0); lhor->DrawLine(0,1.0,16.5,1.0);
  }
  ccan->cd(); ccan->Update();
  ccan->Print(outfileP.Data());    
}

void comp_energy(const char* model="ampt",int cent=16){
  const int NFILES = 8;
  const int NCENT = 3; //plotting for cent bins 1,9 and 16
  int ds[NFILES] = {19,20,31,23,25,18,17,16};
  double eneBin[NFILES] = {1,2,3,4,5,6,7,8};
  double centBin[NCENT] = {0,8,15};//centrality bins for 75-80%, 35-40%, 0-5% 
  double eneErr[NFILES] = {0};
  double energy[NFILES] = {7.7,11.5,14.5,19.6,27,39,62.4,200};
  double B2cp_cent[3][3][NFILES]={0};
  double B2cm_cent[3][3][NFILES]={0};
  double B2cp_err[3][3][NFILES]={0};
  double B2cm_err[3][3][NFILES]={0};

  double B2cp_diff[3][3][NFILES]={0};
  double B2cm_diff[3][3][NFILES]={0};
  //double B2cp_diff_err[3][3][NFILES]={0};
  //double B2cm_diff_err[3][3][NFILES]={0};

  TGraphErrors *gintB2cp_model[3][NFILES] = {0};
  TGraphErrors *gintB2cm_model[3][NFILES] = {0};
  TGraphErrors *gsruB2cp_model[3][NFILES] = {0};
  TGraphErrors *gsruB2cm_model[3][NFILES] = {0};
  TGraphErrors *gintB2cp_energy[3][NCENT] = {0};
  TGraphErrors *gintB2cm_energy[3][NCENT] = {0};
  TGraphErrors *gintB2cp_diff[3][NCENT] = {0};
  TGraphErrors *gintB2cm_diff[3][NCENT] = {0};

  int colors[11] = {kBlue,kGreen+2,kRed,1,kOrange+3,kMagenta,kCyan+2,kMagenta+2,kRed-3,15};
  const char* partname[3] = {"#pi#pi","KK","pp"};
  TFile *f1[NFILES] = {0};
  int col = 0;
  
  for (int ifile=0;ifile<NFILES;ifile++){
    f1[ifile] = new TFile(Form("root_files/B2_int_c16_%s_%d_cuts1.root",model,ds[ifile]),"read");
    f1[ifile]->cd();
    for (int ipart=0;ipart<3;ipart++){
      gintB2cp_model[ipart][ifile] = (TGraphErrors*)f1[ifile]->Get(Form("gint_B2cp_%d_%s",ipart,model));
      gintB2cm_model[ipart][ifile] = (TGraphErrors*)f1[ifile]->Get(Form("gint_B2cm_%d_%s",ipart,model));
      gsruB2cp_model[ipart][ifile] = (TGraphErrors*)f1[ifile]->Get(Form("gsru_B2cp_%d_%s",ipart,model));
      gsruB2cm_model[ipart][ifile] = (TGraphErrors*)f1[ifile]->Get(Form("gsru_B2cm_%d_%s",ipart,model));
      
      for (int icent=0;icent<NCENT;icent++){
	double B2int_p = gintB2cp_model[ipart][ifile]->GetPointY(centBin[icent]);
	double B2int_m = gintB2cm_model[ipart][ifile]->GetPointY(centBin[icent]);
	double B2sru_p = gsruB2cp_model[ipart][ifile]->GetPointY(centBin[icent]);
        double B2sru_m = gsruB2cm_model[ipart][ifile]->GetPointY(centBin[icent]);

	double diff_p = abs(B2int_p-B2sru_p);
	double diff_m = abs(B2int_m-B2sru_m);

	B2cp_diff[icent][ipart][ifile] = diff_p;
        B2cm_diff[icent][ipart][ifile] = diff_m;
	
	B2cp_cent[icent][ipart][ifile] = gintB2cp_model[ipart][ifile]->GetPointY(centBin[icent]);
	B2cp_err[icent][ipart][ifile] = gintB2cp_model[ipart][ifile]->GetErrorY(centBin[icent]);

	B2cm_cent[icent][ipart][ifile] = gintB2cm_model[ipart][ifile]->GetPointY(centBin[icent]);   
	B2cm_err[icent][ipart][ifile] = gintB2cm_model[ipart][ifile]->GetErrorY(centBin[icent]);

      }
    }
  }

  for (int ipart=0;ipart<3;ipart++){
    for (int icent=0;icent<3;icent++){
      gintB2cp_energy[ipart][icent] = new TGraphErrors(NFILES,eneBin,B2cp_cent[icent][ipart],eneErr,B2cp_err[icent][ipart]);
      gintB2cp_energy[ipart][icent]->SetName(Form("gint_B2cp_%d_%s",ipart,model));
      gintB2cp_energy[ipart][icent]        ->SetTitle(Form("Integral %s B_{2+};EnergyBin",partname[ipart]));
      gintB2cp_energy[ipart][icent]        ->SetMarkerStyle(20);
      gintB2cp_energy[ipart][icent]        ->SetMarkerColor( colors[icent]);
      gintB2cp_energy[ipart][icent]        ->SetLineColor  ( colors[icent]);
      gintB2cp_energy[ipart][icent]        ->SetLineWidth  ( 2);
      gintB2cm_energy[ipart][icent]        = new TGraphErrors(NFILES,eneBin,B2cm_cent[icent][ipart],eneErr,B2cm_err[icent][ipart]);
      gintB2cm_energy[ipart][icent]->SetName(Form("gint_B2cm_%d_%s",ipart,model));
      gintB2cm_energy[ipart][icent]        ->SetTitle(Form("Integral %s B_{2-};EnergyBin",partname[ipart]));
      gintB2cm_energy[ipart][icent]        ->SetMarkerStyle(20);
      gintB2cm_energy[ipart][icent]        ->SetMarkerColor( colors[icent]);
      gintB2cm_energy[ipart][icent]        ->SetLineColor  ( colors[icent]);
      gintB2cm_energy[ipart][icent]  ->SetLineWidth  ( 2);

      gintB2cp_diff[ipart][icent] = new TGraphErrors(NFILES,eneBin,B2cp_diff[icent][ipart],eneErr,B2cp_err[icent][ipart]);
      gintB2cp_diff[ipart][icent]->SetName(Form("gdiff_B2cp_%d_%s",ipart,model));
      gintB2cp_diff[ipart][icent]        ->SetTitle(Form("Diff of Integral B_{2+} and sum rule %s;EnergyBin",partname[ipart]));
      gintB2cp_diff[ipart][icent]        ->SetMarkerStyle(20);
      gintB2cp_diff[ipart][icent]        ->SetMarkerColor( colors[icent]);
      gintB2cp_diff[ipart][icent]        ->SetLineColor  ( colors[icent]);
      gintB2cp_diff[ipart][icent]        ->SetLineWidth  ( 2);
      gintB2cm_diff[ipart][icent]        = new TGraphErrors(NFILES,eneBin,B2cm_diff[icent][ipart],eneErr,B2cm_err[icent][ipart]);
      gintB2cm_diff[ipart][icent]->SetName(Form("gdiff_B2cm_%d_%s",ipart,model));
      gintB2cm_diff[ipart][icent]        ->SetTitle(Form("Diff of Integral B_{2-} and sum rule %s;EnergyBin",partname[ipart]));
      gintB2cm_diff[ipart][icent]        ->SetMarkerStyle(20);
      gintB2cm_diff[ipart][icent]        ->SetMarkerColor( colors[icent]);
      gintB2cm_diff[ipart][icent]        ->SetLineColor  ( colors[icent]);
      gintB2cm_diff[ipart][icent]  ->SetLineWidth  ( 2);
    }
  }

  //Painting!!
  gROOT->SetStyle("Modern");
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetPadTopMargin(0.005);
  gStyle->SetPadBottomMargin(0.08);
  gStyle->SetPadLeftMargin(0.12);
  int ican=-1,iline=-1,ivline=-1;
  int itxt=-1;

  TLatex *txt[1000];
  for (int i=0;i<1000;i++){
    txt[i]      = new TLatex();
    txt[i]      ->SetTextSize(0.06);
    txt[i]      ->SetTextAlign(12);
    txt[i]      ->SetNDC();
  }
  
  TCanvas *ccan    = new TCanvas("","",20,30+20,900,650);
  ccan->Divide(3,2,0.0001,0.0001);
  for (int ipart=0;ipart<3;ipart++){
    ccan->cd(ipart+1);
    ylimits[0] = -0.1; ylimits[1] = 1.2;
    gintB2cp_energy[ipart][0]->Draw("AP");
    gintB2cp_energy[ipart][0]->SetMinimum(ylimits[0]);
    gintB2cp_energy[ipart][0]->SetMaximum(ylimits[1]);
    gintB2cp_energy[ipart][0]->Draw("P");
    gintB2cp_energy[ipart][1]->Draw("P");
    gintB2cp_energy[ipart][2]->Draw("P");
    if (ipart==0){
      TLegend *legend1= new TLegend(0.4,0.71,0.7,0.52);
      legend1->SetTextFont(43);
      legend1->SetTextSize(15);
      legend1->SetFillStyle(0);
      legend1->AddEntry(gintB2cp_energy[ipart][0],"75-80%");
      legend1->AddEntry(gintB2cp_energy[ipart][1],"35-40%");
      legend1->AddEntry(gintB2cp_energy[ipart][2],"0-5%");
      legend1->Draw();
    }
  }

   for (int ipart=0;ipart<3;ipart++){
    ccan->cd(ipart+4);
    ylimits[0] = -0.1; ylimits[1] = 1.2;
    gintB2cm_energy[ipart][0]->Draw("AP");
    gintB2cm_energy[ipart][0]->SetMinimum(ylimits[0]);
    gintB2cm_energy[ipart][0]->SetMaximum(ylimits[1]);
    gintB2cm_energy[ipart][0]->Draw("P");
    gintB2cm_energy[ipart][1]->Draw("P");
    gintB2cm_energy[ipart][2]->Draw("P");
  }
  ccan->cd(); ccan->Update();
  ccan->Print(outfileP.Data());
  /*
  TCanvas *ccan    = new TCanvas("","",20,30+20,900,650);
  ccan->Divide(3,2,0.0001,0.0001);
  for (int ipart=0;ipart<3;ipart++){
    ccan->cd(ipart+1);
    if(ipart==0)
      ylimits[0] = -0.02; ylimits[1] = 0.13;
    gintB2cp_diff[ipart][0]->Draw("AP");
    gintB2cp_diff[0][0]->SetMinimum(ylimits[0]);
    gintB2cp_diff[0][0]->SetMaximum(ylimits[1]);
    gintB2cp_diff[ipart][0]->Draw("P");
    gintB2cp_diff[ipart][1]->Draw("P");
    gintB2cp_diff[ipart][2]->Draw("P");

    if (ipart==0){
      //TLegend *legend1= new TLegend(0.4,0.71,0.7,0.52);//urqmd
      //TLegend *legend1= new TLegend(0.5,0.91,0.8,0.72);//ampt
      TLegend *legend1= new TLegend(0.45,0.61,0.75,0.42);//hijing
      legend1->SetTextFont(43);
      legend1->SetTextSize(15);
      legend1->SetFillStyle(0);
      legend1->AddEntry(gintB2cp_diff[ipart][0],"75-80%");
      legend1->AddEntry(gintB2cp_diff[ipart][1],"35-40%");
      legend1->AddEntry(gintB2cp_diff[ipart][2],"0-5%");
      legend1->Draw();
    }
  }
  for (int ipart=0;ipart<3;ipart++){
    ccan->cd(ipart+4);
    if(ipart==0)
      ylimits[0] = -0.02; ylimits[1] = 0.13;
    gintB2cm_diff[ipart][0]->Draw("AP");
    gintB2cm_diff[0][0]->SetMinimum(ylimits[0]);
    gintB2cm_diff[0][0]->SetMaximum(ylimits[1]);
    gintB2cm_diff[ipart][0]->Draw("P");
    gintB2cm_diff[ipart][1]->Draw("P");
    gintB2cm_diff[ipart][2]->Draw("P");
  }
  
  
  ccan->cd(); ccan->Update();
  ccan->Print(outfileP.Data());
  */
}

//----------------------------------------------------------------------------
void RangeFinderH(int isteer, TH1D* h){
  bool debug = false;
  if (isteer>=100){
    debug	 = true;
    isteer	-= 100;
    cout<<"debug on... "<<isteer<<endl;
  }
  if (isteer==0){
    ylimits[0]	=  999999;
    ylimits[1]	= -999999;
  }
  if (!h) return;
  if (debug) cout<<"start... "<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
  if (debug) cout<<h->GetName()<<" "<<h->GetEntries()<<endl;
  //if (h->GetEntries()<1){
  //	ylimits[0]	=  0;
  //	ylimits[1]	=  1;
  //	return;
  //}
  int nbx	= h->GetNbinsX();
  for (int ibx=1;ibx<=nbx;ibx++){
    float val	= h->GetBinContent(ibx);
    float vale	= h->GetBinError(ibx);
    //if (val==0) continue;
    if (val-vale<ylimits[0]){ ylimits[0] = val-vale; }	// update lower limit...
    if (val+vale>ylimits[1]){ ylimits[1] = val+vale; }	// update upper limit...
    if (debug) cout<<"check... "<<ibx<<" val="<<val<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
  }
  if (debug) cout<<"FINAL... "<<isteer<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
}
//----------------------------------------------------------------------------
void RangeFinderG(int isteer, TGraph* g){
  bool debug = false;
  if (isteer>=100){
    debug	 = true;
    isteer	-= 100;
    cout<<"debug on... "<<isteer<<endl;
  }
  if (isteer==0){
    ylimits[0]	=  999999;
    ylimits[1]	= -999999;
  }
  //if (debug) cout<<"start... "<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
  //if (debug) cout<<g->GetName()<<" "<<g->GetEntries()<<endl;
  //if (h->GetEntries()<1){
  //	ylimits[0]	=  0;
  //	ylimits[1]	=  1;
  //	return;
  //}
  int nbx	= g->GetN();
  for (int ibx=0;ibx<nbx;ibx++){
    double x,val;
    g->GetPoint(ibx,x,val);
    //if (val==0) continue;
    if (val<ylimits[0]){ ylimits[0] = val; }
    if (val>ylimits[1]){ ylimits[1] = val; }
    if (debug) cout<<"check... "<<ibx<<" val="<<val<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
  }
  if (debug) cout<<"FINAL... "<<isteer<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
}
//----------------------------------------------------------------------------
void RangeFinderGE(int isteer, TGraphErrors* g){
  bool debug = false;
  if (isteer>=100){
    debug	 = true;
    isteer	-= 100;
    cout<<"debug on... "<<isteer<<endl;
  }
  if (isteer==0){
    ylimits[0]	=  999999;
    ylimits[1]	= -999999;
  }
  //if (debug) cout<<"start... "<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
  //if (debug) cout<<g->GetName()<<" "<<g->GetEntries()<<endl;
  //if (h->GetEntries()<1){
  //	ylimits[0]	=  0;
  //	ylimits[1]	=  1;
  //	return;
  //}
  int nbx	= g->GetN();
  for (int ibx=0;ibx<nbx;ibx++){
    double x,val,vale;
    g->GetPoint(ibx,x,val);
    vale	= g->GetErrorY(ibx);
    //if (val==0) continue;
    if (val-vale<ylimits[0]){ ylimits[0] = val-vale; }
    if (val+vale>ylimits[1]){ ylimits[1] = val+vale; }
    if (debug) cout<<"check... "<<ibx<<" val="<<val<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
  }
  if (debug) cout<<"FINAL... "<<isteer<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
}
//----------------------------------------------------------------------------
void RangeCheck(){
  if (ylimits[0]>0){ 
    ylimits[0]	*= 0.9;
  } else if (ylimits[0]==0){
    //		ylimits[0]	= -1;
  } else if (ylimits[0]<0){
    ylimits[0]	*= 1.1;
  }
  if (ylimits[1]>0){ 
    ylimits[1]	*= 1.1;
  } else if (ylimits[1]==0){
    //		ylimits[1]	=    1;
  } else if (ylimits[1]<0){
    ylimits[1]	*= 0.9;
  }
}
//----------------------------------------------------------------------------
