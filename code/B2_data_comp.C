#include "fluct_common.h"
#include "TString.h"
#include "TLegend.h"
#include "TMath.h"

double ylimits[2];
void RangeFinderH(int, TH1D*);
void RangeFinderG(int, TGraph*);
void RangeFinderGE(int, TGraphErrors*);
void RangeCheck();

void B2_data_comp(){
    
  const int NFILES = 2;//3 simulations + 2 data (mixing with 2 pref used) ... will add conv. later
  const int NPid = 3;//pions, kaons and protons...
  //const int NCent = 16;
  
  const char* nfile[NFILES] = {"root_files/B2_mult_m40_data_16_cuts1_W1minvqcut_xp.root",
    "root_files/B2_mult_m40_data_16_cuts1_W1minvqcut_xs.root"
    //"root_files/B2_int_c16_hijing_16_cuts1.root",
    //"root_files/B2_mult_m40_data_24_cuts1_xp.root",
    //"root_files/B2_pub_m40_data_24_cuts1_xp.root"
  };
  

  const char* pid_name[3] = {"pipi","KK","pp"};
  const char* leg_name[NFILES] = {"ds16_xp","ds16_xs"};//"ds16: hijing","ds_24_xp: Pref from mult","ds_24_xp: Pub pref"};
  //const char* leg_name[NFILES] = {"Urqmd","Ampt","Hijing","Data - pref from #LTN#GT","Data - published pref"};
  //int colors[11] = {kBlue,kGreen+2,kRed,1,kOrange+3,kMagenta,kCyan+2,kMagenta+2,kRed-3,15};
  const int color[NFILES] = {kBlue,kGreen+2};//,kRed,kMagenta,kRed-5,1};//kOrange+3,kMagenta,kCyan+2,kMagenta+2,kRed-3,15};
  const int markerStyle[NFILES] = {27,33};//,42,43,34};
  const double markerSize[NFILES] = {1.0,0.6};//,1.0,0.5,1.0};
  const double lineSize[NFILES] = {3,1};//,1.0,0.5,1.0};
  
  TH2D *hB2[NFILES][NCent][NPid]  = {0};
  TH1D *hB2dy[NFILES][NCent][NPid]  = {0};
  TH1D *hB2dphi[NFILES][NCent][NPid]  = {0};
  
  for(int ifile=0;ifile<NFILES;ifile++){
    TFile *f = new TFile(nfile[ifile],"read"); f->cd();
    for (int icent=1; icent<NCent; icent++){
      for (int ipaty=0;ipaty<NPid;ipaty++){  //ipaty here is the part. type. Pions =0, kaons = 1 and protons = 2
	//cout<<"ok1\n";
        hB2[ifile][icent][ipaty] = (TH2D*)f->Get(Form("hB2c_%d_%d",icent,ipaty));
        hB2[ifile][icent][ipaty]->SetName(Form("hB2_%d_%d_%d",ifile,icent,ipaty));
        hB2[ifile][icent][ipaty]->SetStats(0);
	hB2[ifile][icent][ipaty]->GetXaxis()->SetTitleSize(0.06);
	hB2[ifile][icent][ipaty]->GetYaxis()->SetTitleSize(0.05);
	hB2[ifile][icent][ipaty]->GetXaxis()->SetLabelSize(0.05);
        hB2[ifile][icent][ipaty]->GetYaxis()->SetLabelSize(0.05);
	hB2[ifile][icent][ipaty]->GetZaxis()->SetLabelSize(0.05);

        hB2dy[ifile][icent][ipaty] =(TH1D*)f->Get(Form("hB2cdy_%d_%d",icent,ipaty));
        hB2dy[ifile][icent][ipaty]->SetName(Form("hB2dy_%d_%d_%d",ifile,icent,ipaty));
        hB2dy[ifile][icent][ipaty]->SetStats(0);
        hB2dy[ifile][icent][ipaty]->SetLineColor(color[ifile]);
	hB2dy[ifile][icent][ipaty]->SetLineWidth(lineSize[ifile]);
        hB2dy[ifile][icent][ipaty]->SetMarkerColor(color[ifile]);
        hB2dy[ifile][icent][ipaty]->SetMarkerStyle(markerStyle[ifile]);
        hB2dy[ifile][icent][ipaty]->SetMarkerSize(markerSize[ifile]);
	hB2dy[ifile][icent][ipaty]->GetXaxis()->SetTitleSize(0.06);
        hB2dy[ifile][icent][ipaty]->GetXaxis()->SetLabelSize(0.05);
        hB2dy[ifile][icent][ipaty]->GetYaxis()->SetLabelSize(0.05);

	hB2dphi[ifile][icent][ipaty] =(TH1D*)f->Get(Form("hB2cdphi_%d_%d",icent,ipaty));
        hB2dphi[ifile][icent][ipaty]->SetName(Form("hB2dphi_%d_%d_%d",ifile,icent,ipaty));
        hB2dphi[ifile][icent][ipaty]->SetStats(0);
        hB2dphi[ifile][icent][ipaty]->SetLineColor(color[ifile]);
        hB2dphi[ifile][icent][ipaty]->SetMarkerColor(color[ifile]);
	hB2dphi[ifile][icent][ipaty]->SetLineWidth(lineSize[ifile]);
        hB2dphi[ifile][icent][ipaty]->SetMarkerStyle(markerStyle[ifile]);
        hB2dphi[ifile][icent][ipaty]->SetMarkerSize(markerSize[ifile]);
	hB2dphi[ifile][icent][ipaty]->GetXaxis()->SetTitleSize(0.06);
        hB2dphi[ifile][icent][ipaty]->GetXaxis()->SetLabelSize(0.05);
        hB2dphi[ifile][icent][ipaty]->GetYaxis()->SetLabelSize(0.05);

      }
    }
  }

  gStyle->SetPalette(1);
  gROOT->SetStyle("Modern");
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetPadTopMargin(0.005);
  gStyle->SetPadBottomMargin(0.08);
  gStyle->SetPadLeftMargin(0.12);
  int ican=-1,iline=-1,ivline=-1;
  int itxt=-1;
  TCanvas *ccan[10];//[10];
  TLatex *txt[1000];
  for (int i=0;i<1000;i++){
    txt[i]      = new TLatex();
    txt[i]      ->SetTextSize(0.06);
    txt[i]      ->SetTextAlign(12);
    txt[i]      ->SetNDC();
  }
  //int ican=0;
  //TLine* lhor   = new TLine(0,0,1,1); lhor->SetLineWidth(2); lhor->SetLineColor(16); lhor->SetLineStyle(2);
 
  for (int ipart=0;ipart<3;ipart++){
    int can_index=1;
    ccan[ipart]  = new TCanvas("","",20,30+20,750,900);
    ccan[ipart]->Divide(2,2,0.0001,0.0001);
    for (int ifi=0;ifi<NFILES;ifi++){//for data with pub pref and hmult pref. 
      ccan[ipart]->cd(can_index);
      hB2[ifi][16][ipart]->Draw("lego2");
      can_index++;
    }
    TLegend *legend1= new TLegend(0.25,0.8,0.55,0.90);
    legend1->SetTextFont(43);
    legend1->SetTextSize(13);
    legend1->SetBorderSize(0);
    legend1->SetFillStyle(0);//Transparent leg   
    for (int ifi=0;ifi<NFILES;ifi++){//for data with pub pref and hmult pref.
      ccan[ipart]->cd(can_index);
      /* TLegend *legend1= new TLegend(0.25,0.8,0.55,0.90);
      legend1->SetTextFont(43);
      legend1->SetTextSize(13);
      legend1->SetFillStyle(0);//Transparent leg 
      */
      legend1->AddEntry(hB2dy[ifi][16][ipart],leg_name[ifi]);
      if(ifi==0){
	hB2dy[ifi][16][ipart]->Draw();
	/*TLegend *legend1= new TLegend(0.25,0.8,0.55,0.90);
	legend1->SetTextFont(43);
	legend1->SetTextSize(13);
	legend1->SetFillStyle(0);//Transparent leg
	legend1->AddEntry(hB2dy[ifi][16][ipart],leg_name[ifi]);*/
	
      }
      else hB2dy[ifi][16][ipart]->Draw("same");
      legend1->Draw();
    }
    can_index++;
    for (int ifi=0;ifi<NFILES;ifi++){//for data with pub pref and hmult pref.
      ccan[ipart]->cd(can_index);
      //legend1->AddEntry(hB2dphi[ifi][16][ipart],leg_name[ifi]);
      if(ifi==0)
	hB2dphi[ifi][16][ipart]->Draw();
      else hB2dphi[ifi][16][ipart]->Draw("same");
      legend1->Draw();
    }
  }
    
  /*
  //TCanvas *ccan    = new TCanvas("","",20,30+20,900,650);
  ccan->Divide(3,2,0.0001,0.0001);
  TLegend *legend1= new TLegend(0.25,0.8,0.55,0.90);
  legend1->SetTextFont(43);
  legend1->SetTextSize(13);
  legend1->SetFillStyle(0);//Transparent leg 
  for (int ipart=0;ipart<NPid;ipart++){
    for (int ifile=0;ifile<NFILES;ifile++){
      ccan->cd(ipart+1);
      //RangeFinderG(0,hB2dy[ifile][16][ipart]);
      //RangeFinderG(1,hB2dy[ifile][16][ipart]);
      //RangeFinderG(2,hB2dy[ifile][16][ipart]);
      //RangeCheck();
      
      hB2dy[ifile][16][ipart]->GetXaxis()->SetRangeUser(0,0.8);
      ytg_S[ipart] = new TGraph(Form("SP_B2dy_%s.CSV",pid_name[ipart]),"%lg %lg",",");
      ytg_S[ipart]->SetLineColor(color[NFILES]);
      if(ifile==0)
	hB2dy[ifile][16][ipart]->Draw();
      else
	hB2dy[ifile][16][ipart]->Draw("same");
      ytg_S[ipart]->Draw("C");
      if(ipart==0){
      legend1->AddEntry(hB2dy[ifile][16][ipart],leg_name[ifile]);
      }
    }
    if(ipart==0){
      legend1->AddEntry(ytg_S[ipart],leg_name[NFILES]);
      legend1->Draw();
    }
  }
  for (int ipart=0;ipart<3;ipart++){
    for (int ifile=0;ifile<NFILES;ifile++){
      ccan->cd(ipart+4);
      //RangeFinderG(0,hB2dphi[ifile][16][ipart]);
      //RangeFinderG(1,hB2dphi[ifile][16][ipart]);
      //RangeFinderG(2,hB2dphi[ifile][16][ipart]);
      //RangeCheck();
      tg_S[ipart]  = new TGraph(Form("SP_B2dphi_%s.CSV",pid_name[ipart]),"%lg %lg",",");
      tg_S[ipart]->SetLineColor(color[NFILES]);
      hB2dphi[ifile][16][ipart]->GetXaxis()->SetRangeUser(0,180);
      hB2dphi[ifile][16][ipart]->GetYaxis()->SetRangeUser(-0.05,0.5);
      if(ifile==0)
        hB2dphi[ifile][16][ipart]->Draw();
      else
        hB2dphi[ifile][16][ipart]->Draw("same");
      tg_S[ipart]->Draw("C");
    }
  }
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

