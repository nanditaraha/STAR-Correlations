/*This macro compares run 11 and run 10 R2s, mults etc*/
#include "fluct_common.h"
//For my ref: Order of NPairTypes: {pi+pi+, pi-pi-, pi+pi-, pi-pi+, then same for kaons and finally protns}
//NCent = 17 but includes a  min-bias

double GetR2Baseline(TH2D *hmult, bool fDistinguishable);
double drawPlots(int cent);
void printB2_urqmd(int cent=16){

  const int NFILES = 3;
  const int NPID = 3;
  
  TH2D *hB2[NFILES][NCent][NPID]  = {0};
  //TH1D *hB2dy[NFILES][NCent][NPID]  = {0};
  TH1D *hB2dphi[NFILES][NCent][NPID]  = {0};
  /*
  // +ve part of BF
  TH2D *hB2_1[NFILES][NCent][NPID]  = {0};
  TH1D *hB2_1dy[NFILES][NCent][NPID]  = {0};
  TH1D *hB2_1dphi[NFILES][NCent][NPID]  = {0};

  
 // -ve part of BF
  TH2D *hB2_2[NFILES][NCent]  = {0};
  TH1D *hB2_2dy[NFILES][NCent]  = {0};
  TH1D *hB2_2dphi[NFILES][NCent]  = {0};
 */
  //TH2D *hmult[NFILES][NCent][NPairTypes]  = {0};
  //TH1D *hmult1[NFILES][NCent][NPairTypes]  = {0};
  
  const char* partname[3] = {"#pi","K","p"};
  const char* partname1[3] = {"pi","K","p"};
  double dy_Y_min, dy_Y_max;
  double dphi_Y_min, dphi_Y_max;
  TCanvas *c = new TCanvas();//B2 - 2D
  TCanvas *c1 = new TCanvas();//Proj
  //TCanvas *c1 = new TCanvas("c1","c1_1",20,50,900,650);//Projections - B2
  TCanvas *c2 = new TCanvas();//+ve parts of B2 - 2D
  TCanvas *c3 = new TCanvas();//+ve and -ve parts of B2 - projections
  TCanvas *c4 = new TCanvas();//-ve parts of B2 - 2D

  TGraph *tg_S[3] = {0};

  const char* infile[NFILES]={
    "root_files/B2_c4_urqmdLN_16_cuts1.root",
    "root_files/B2_c4_urqmdLN_31_cuts1.root",
    "root_files/B2_c4_urqmdLN_25_cuts1.root"
			      //"BF_m40_urqmd_16_cuts1_xp.root",
			      /*"BF_m40_urqmd_16_cuts1_pions1012.root",
			      "BF_m40_urqmd_16_cuts1_pions1024.root",
			      "BF_m40_urqmd_16_cuts1_pions2012.root",
			      "BF_m40_urqmd_16_cuts1_pions2024.root",
			      "BF_m40_urqmd_16_cuts1_pions2436.root",
			      "BF_m40_urqmd_16_cuts1_pions3648.root"
			      */
  };
  //TFile *f1 = new TFile("/tier2/home/groups/rhi/share/readerroot2021/reader_m40_data_24_cuts1_flipsome_xp.root","read");
  //f->cd();
  
  const int color[NFILES+4] = {kRed,kBlack,kBlue, kOrange,kOrange+2,4,2};
  //const int color[NFILES] = {kBlue,kMagenta,kRed};
  const int markerStyle[NFILES+4] = {27,20,21,22,27,20,21};
  const double markerSize[NFILES+4] = {1,1,1,1,1,1,1};
  //const double* markerSize[NFILES] = {4,20};
  //const char* leg_names[NFILES] = {"Run 10", "Run 11","Run 10 - baseline subtracted","Run 11 - baseline subtracted","Run 10 - URQMD"};
  //const char* leg_names[NFILES] = {"y-bins=06, #varphi-bins=06","y-bins=20, #varphi-bins=24","y-bins=36, #varphi-bins=48"};
  const char* leg_names[NFILES+4] = {"#sqrt{S_{NN}} = 200 GeV","#sqrt{S_{NN}} = 14.5 GeV","#sqrt{S_{NN}} = 27 GeV",
  				   "y-bins=20, #varphi-bins=12","y-bins=20, #varphi-bins=24","y-bins=24, #varphi-bins=36",
  				   "y-bins=36, #varphi-bins=48"
  };
  gStyle->SetPalette(1);

  for (int ipart=0;ipart<3;ipart++){            // pi, then K, then p
    tg_S[ipart]        = new TGraphErrors(Form("Scott_BFdphi_%s.csv",partname1[ipart]),"%lg %lg",",");
    tg_S[ipart]        ->SetTitle(Form("%s B_{2}(#varphi);B_{2}(#varphi)",partname[ipart]));
    tg_S[ipart]        ->SetMarkerStyle(20);
    tg_S[ipart]        ->SetMarkerColor( 4);
    tg_S[ipart]        ->SetLineColor  ( 4);
    tg_S[ipart]        ->SetLineWidth  ( 2);
  }
  
  for(int ifile=0;ifile<NFILES;ifile++){
    TFile *f = new TFile(Form("%s",infile[ifile]),"read");
    f->cd();
    for (int icent=1; icent<NCent; icent++){
      if(icent==16){
	cout<<"OK cent="<<icent<<" \tFile"<<ifile<<endl;
      for (int ipid=0;ipid<NPID;ipid++){
	hB2[ifile][icent][ipid] = (TH2D*)f->Get(Form("hB2c_%d_%d",icent,ipid));
	hB2[ifile][icent][ipid]->SetName(Form("hBF2c_%d_%d_%d",ifile,icent,ipid));
	hB2[ifile][icent][ipid]->SetStats(0);
	//hB2[ifile][icent][ipid]->SetTitle("#pi#pi B_{2}(dy, d#varphi)");
	cout<<"Got B2(dy,dphi)"<<endl;
	//double scale = hB2[ifile][icent][ipid]->GetXaxis()->GetNbins();
	hB2dphi[ifile][icent][ipid] = (TH1D*)f->Get(Form("hB2cdphi_%d_%d",icent,ipid));
	//hB2dphi[ifile][icent][ipid]->Scale(1./scale);
	//hB2dphi[ifile][icent][ipid]->GetXaxis()->SetRangeUser(0,180);
        hB2dphi[ifile][icent][ipid]->SetName(Form("hB2dphi_%d_%d_%d",ifile,icent,ipid));
        hB2dphi[ifile][icent][ipid]->SetStats(0);
        hB2dphi[ifile][icent][ipid]->SetLineColor(color[ifile]);
        hB2dphi[ifile][icent][ipid]->SetMarkerColor(color[ifile]);
        hB2dphi[ifile][icent][ipid]->SetMarkerSize(markerSize[ifile]);
        hB2dphi[ifile][icent][ipid]->SetMarkerStyle(markerStyle[ifile]);
        hB2dphi[ifile][icent][ipid]->GetXaxis()->SetLabelSize(0.04);
        hB2dphi[ifile][icent][ipid]->GetYaxis()->SetLabelSize(0.04);
        hB2dphi[ifile][icent][ipid]->GetXaxis()->SetLabelOffset(0.0);
	cout<<"Got B2(dphi)"<<endl;
	/*
        hB2dphi[ifile][icent]->GetYaxis()->SetTitle("B_{2}(#Delta#varphi)");
        hB2dphi[ifile][icent]->GetXaxis()->SetTitle("#Delta#varphi");
	hB2dy[ifile][icent] =(TH1D*)f->Get(Form("hBFcdy_pipi_%d",icent));
	hB2dy[ifile][icent]->SetName(Form("hB2dy_%d_%d",ifile,icent));
	hB2dy[ifile][icent]->SetStats(0);
	hB2dy[ifile][icent]->SetLineColor(color[ifile]);
	hB2dy[ifile][icent]->SetMarkerColor(color[ifile]);
	hB2dy[ifile][icent]->SetMarkerSize(markerSize[ifile]);
	hB2dy[ifile][icent]->SetMarkerStyle(markerStyle[ifile]);
	hB2dy[ifile][icent]->GetXaxis()->SetLabelSize(0.04);
	hB2dy[ifile][icent]->GetYaxis()->SetLabelSize(0.04);
	hB2dy[ifile][icent]->GetXaxis()->SetLabelOffset(0.0);
	hB2dy[ifile][icent]->GetYaxis()->SetTitle("B_{2}(#Delta y)");
        hB2dy[ifile][icent]->GetXaxis()->SetTitle("#Delta y");

	hB2dphi[ifile][icent] =(TH1D*)f->Get(Form("hBFcdphi_pipi_%d",icent));
	hB2dphi[ifile][icent]->SetName(Form("hB2dphi_%d_%d",ifile,icent,ifile));
	hB2dphi[ifile][icent]->SetStats(0);
	hB2dphi[ifile][icent]->SetLineColor(color[ifile]);
	hB2dphi[ifile][icent]->SetMarkerColor(color[ifile]);
	hB2dphi[ifile][icent]->SetMarkerSize(markerSize[ifile]);
	hB2dphi[ifile][icent]->SetMarkerStyle(markerStyle[ifile]);
	hB2dphi[ifile][icent]->GetXaxis()->SetLabelSize(0.04);
	hB2dphi[ifile][icent]->GetYaxis()->SetLabelSize(0.04);
	hB2dphi[ifile][icent]->GetXaxis()->SetLabelOffset(0.0);
	hB2dphi[ifile][icent]->GetYaxis()->SetTitle("B_{2}(#Delta#varphi)");
        hB2dphi[ifile][icent]->GetXaxis()->SetTitle("#Delta#varphi");

	//B2+
	hB2_1[ifile][icent] = (TH2D*)f->Get(Form("hBFc1_pipi_%d",icent));
        hB2_1[ifile][icent]->SetName(Form("hBF2_1_%d_%d",ifile,icent));
	hB2_1[ifile][icent]->GetXaxis()->SetTitleOffset(1.6);
	hB2_1[ifile][icent]->GetYaxis()->SetTitleOffset(2.4);
        hB2_1[ifile][icent]->SetStats(0);
	hB2_1[ifile][icent]->SetTitle("B_{2}+ #pi#pi");

        hB2_1dy[ifile][icent] =(TH1D*)f->Get(Form("hBFc1dy_pipi_%d",icent));
        hB2_1dy[ifile][icent]->SetName(Form("hB2_1_dy_%d_%d",ifile,icent));
        hB2_1dy[ifile][icent]->SetStats(0);
        hB2_1dy[ifile][icent]->SetLineColor(color[ifile]);
        hB2_1dy[ifile][icent]->SetMarkerColor(color[ifile]);
	hB2_1dy[ifile][icent]->SetMarkerSize(markerSize[ifile]);
        hB2_1dy[ifile][icent]->SetMarkerStyle(markerStyle[ifile]);
        //hB2_1dy[ifile][icent]->GetXaxis()->SetLabelSize(0.04);
        hB2_1dy[ifile][icent]->GetYaxis()->SetLabelSize(0.04);
        //hB2_1dy[ifile][icent]->GetXaxis()->SetLabelOffset(0.0);
        hB2_1dy[ifile][icent]->GetYaxis()->SetTitle("B_{2}(#Delta y)");
        hB2_1dy[ifile][icent]->GetXaxis()->SetTitle("#Delta y");
	hB2_1dy[ifile][icent]->GetZaxis()->SetTitle("");
	//hB2_1dy[ifile][icent]->SetTitle("");

        hB2_1dphi[ifile][icent] =(TH1D*)f->Get(Form("hBFc1dphi_pipi_%d",icent));
        hB2_1dphi[ifile][icent]->SetName(Form("hB2_1_dphi_%d_%d",ifile,icent,ifile));
        hB2_1dphi[ifile][icent]->SetStats(0);
        hB2_1dphi[ifile][icent]->SetLineColor(color[ifile]);
        hB2_1dphi[ifile][icent]->SetMarkerColor(color[ifile]);
        hB2_1dphi[ifile][icent]->SetMarkerStyle(markerStyle[ifile]);
	hB2_1dphi[ifile][icent]->SetMarkerSize(markerSize[ifile]);
        //hB2_1dphi[ifile][icent]->GetXaxis()->SetLabelSize(0.04);
        hB2_1dphi[ifile][icent]->GetYaxis()->SetLabelSize(0.04);
        //hB2_1dphi[ifile][icent]->GetXaxis()->SetLabelOffset(0.0);
        hB2_1dphi[ifile][icent]->GetYaxis()->SetTitle("B_{2}(#Delta#varphi)");
        hB2_1dphi[ifile][icent]->GetXaxis()->SetTitle("#Delta#varphi");
	//hB2_1dphi[ifile][icent]->SetTitle("");
	
	//B2-
	hB2_2[ifile][icent] = (TH2D*)f->Get(Form("hBFc2_pipi_%d",icent));
        hB2_2[ifile][icent]->SetName(Form("hBF2_2_%d_%d",ifile,icent));
	hB2_2[ifile][icent]->GetXaxis()->SetTitleOffset(1.6);
        hB2_2[ifile][icent]->GetYaxis()->SetTitleOffset(2.4);
        hB2_2[ifile][icent]->SetStats(0);
	hB2_2[ifile][icent]->SetTitle("B_{2}- #pi#pi");

        hB2_2dy[ifile][icent] =(TH1D*)f->Get(Form("hBFc2dy_pipi_%d",icent));
        hB2_2dy[ifile][icent]->SetName(Form("hB2_2_dy_%d_%d",ifile,icent));
        hB2_2dy[ifile][icent]->SetStats(0);
        hB2_2dy[ifile][icent]->SetLineColor(color[ifile]);
        hB2_2dy[ifile][icent]->SetMarkerColor(color[ifile]);
	hB2_2dy[ifile][icent]->SetMarkerSize(markerSize[ifile]);
        hB2_2dy[ifile][icent]->SetMarkerStyle(markerStyle[ifile]);
        //hB2_2dy[ifile][icent]->GetXaxis()->SetLabelSize(0.04);
        hB2_2dy[ifile][icent]->GetYaxis()->SetLabelSize(0.04);
        //hB2_2dy[ifile][icent]->GetXaxis()->SetLabelOffset(0.0);
        hB2_2dy[ifile][icent]->GetYaxis()->SetTitle("B_{2}(#Delta y)");
        hB2_2dy[ifile][icent]->GetXaxis()->SetTitle("#Delta y");
	hB2_2dy[ifile][icent]->SetTitle("");

        hB2_2dphi[ifile][icent] =(TH1D*)f->Get(Form("hBFc2dphi_pipi_%d",icent));
        hB2_2dphi[ifile][icent]->SetName(Form("hB2_2_dphi_%d_%d",ifile,icent,ifile));
        hB2_2dphi[ifile][icent]->SetStats(0);
        hB2_2dphi[ifile][icent]->SetLineColor(color[ifile]);
        hB2_2dphi[ifile][icent]->SetMarkerColor(color[ifile]);
	hB2_1dphi[ifile][icent]->SetMarkerSize(markerSize[ifile]);
        hB2_2dphi[ifile][icent]->SetMarkerStyle(markerStyle[ifile]);
        //hB2_2dphi[ifile][icent]->GetXaxis()->SetLabelSize(0.04);
        hB2_2dphi[ifile][icent]->GetYaxis()->SetLabelSize(0.04);
        //hB2_2dphi[ifile][icent]->GetXaxis()->SetLabelOffset(0.0);
        hB2_2dphi[ifile][icent]->GetYaxis()->SetTitle("B_{2}(#Delta#varphi)");
        hB2_2dphi[ifile][icent]->GetXaxis()->SetTitle("#Delta#varphi");
	hB2_2dphi[ifile][icent]->SetTitle("");
	
	
	hmult[ifile][icent] = (TH2D*)f->Get(Form("hmult_%d_%d",icent,));
        hmult[ifile][icent]->SetName(Form("hmult_%d_%d_%d",ifile,icent,));
	hmult[ifile][icent]->SetStats(0);

	hmult1[ifile][icent] = (TH1D*)f->Get(Form("hmult1_%d_%d",icent,));
        hmult1[ifile][icent]->SetName(Form("hmult1_%d_%d_%d",ifile,icent,));
        hmult1[ifile][icent]->SetStats(0);
	hmult1[ifile][icent]->SetLineColor(color[ifile]);
	*/
	//}
      }//end of pid
      } 
    }// end of icent
    
  }// end of ifile
  //pions
    /*  if(pid>=0 && pid<4){
    dy_Y_min = 0.002; 
    dy_Y_max = 0.0055; 
    dphi_Y_min = 0.0018;
    dphi_Y_max = 0.0073;
  }
  //LS Kaons
  else if(pid>3 && pid<6){
    dy_Y_min = -0.006;
    dy_Y_max = 0.006;
    dphi_Y_min = -0.001;
    dphi_Y_max = 0.005;
  }
  //ULS Kaons
  else if(pid>5 && pid<8){
    dy_Y_min = -0.001;
    dy_Y_max = 0.006;
    dphi_Y_min = 0.002;
    dphi_Y_max = 0.011;
  }
  //LS protons
  else if(pid>7 && pid<10){
    dy_Y_min = -0.025;
    dy_Y_max = 0.002;
    dphi_Y_min = -0.025;
    dphi_Y_max = -0.003;
  }  
  //ULS protons
  else if(pid+2>9 && pid+2<12){
    dy_Y_min = 0.0;
    dy_Y_max = 0.012;
    dphi_Y_min = 0.0;
    dphi_Y_max = 0.012;
  }

  //Drawing B2(dy,dphi)
  */
  c->Divide(3,1);  
  //int k=1;
  
  for(int i=0;i<3;i++){
    c->cd(i+1);
    TLegend *leg= new TLegend(0.20,0.35,0.6,0.25);
    leg->SetTextFont(43);                                                                                                                                      
    leg->SetTextSize(12);                                                                                                                                      
    //leg->SetBorderSize(0);                                                                                                                                     
    
    for(int j=0;j<NFILES;j++){
      if(j==0)
	hB2dphi[j][cent][i]->Draw();
      else hB2dphi[j][cent][i]->Draw("same");
      cout<<"file# = "<<j<<endl;
      if(i==0){
	leg->AddEntry(hB2dphi[j][cent][i],leg_names[j]);
	leg->Draw();
      }
    }
  }
      //hB2[j][cent]->Draw("surf1");
      //k++;
  
  //End drawing
  
  
  //Drawing B2(dy) and B2(dphi)
  /*
  c1->Divide(2,1,0.0001,0.0001);
  c1->cd(1);//dy compt on 1st section of canvas    
  TMultiGraph *mg1 = new TMultiGraph();
  mg1->Add(ytg_J);
  ytg_J->Draw("alp");
  //mg1->Draw("LP");
  TLegend *leg= new TLegend(0.20,0.35,0.6,0.25);
  leg->SetTextFont(43);
  leg->SetTextSize(12);
  leg->SetBorderSize(0);
  for(int j=0;j<NFILES;j++){
    if(j==0)
      hB2dy[j][cent]->Draw();
    else hB2dy[j][cent]->Draw("same");
    leg->SetTextFont(43);
    leg->SetTextSize(15);
    leg->SetBorderSize(0);
    leg->AddEntry(hB2dphi[j][cent],leg_names[j]);
    
  }
  leg->Draw("same");
  */
  c1->Divide(3,1,0.0001,0.0001);
  TLegend *leg= new TLegend(0.20,0.35,0.6,0.25);
    leg->SetTextFont(43);                                                                                                                         
    leg->SetTextSize(15);
  for(int pid=0;pid<3;pid++){
    TMultiGraph *mg = new TMultiGraph();
    c1->cd(pid+1);//dphi compt on 1st section of canvas    
    mg->Add(tg_S[pid]);
    tg_S[pid]->Draw("alp");
    mg->Draw("LP");
    for(int j=0;j<1;j++){
      //if(j==0)
      //hB2dphi[j][cent][pid]->Draw();
      //else
      hB2dphi[j][cent][pid]->Draw("same");
      if(pid==0){ leg->AddEntry(hB2dphi[j][cent][pid],"UrQMD");
	leg->AddEntry(tg_S[pid],"Scott's Plot");
	leg->Draw();
      }
			       
      //hB2dphi[j][cent][pid]->Draw();
    }
  }
  /*
  for(int pid=0;pid<3;pid++){
    // TMultiGraph *mg = new TMultiGraph();
    c1->cd(pid+4);//dphi compt on 1st section of canvas
    for(int j=0;j<1;j++){
      if(j==0)                                                                                                                                               
	hB2dphi[j][cent][pid]->Draw();                                                                                                                    
      else                                                                                                                                                   
	hB2dphi[j][cent][pid]->Draw("same");

    }
    }
  */
  //Drawing each part!! B2+ and B2-
  /*
  c2->Divide(NFILES,1);
  int l=1;
  //Draw B2+
  for(int j=0;j<NFILES;j++){
    c2->cd(l);
    hB2_1[j][cent]->Draw("surf1");
    l++;
  }
  //Draw B2-

  c4->Divide(NFILES,1,0.001,0.001);
  int m=1;
  for(int j=0;j<NFILES;j++){
    c4->cd(m);
    hB2_2[j][cent]->Draw("surf1");
    m++;
  }

  c3->Divide(2,2,0.0001,0.0001);
  c3->cd(1);
  TLegend *leg1= new TLegend(0.16,0.37,0.55,0.12);
  leg1->SetTextFont(43);
  leg1->SetTextSize(14);
  leg1->SetBorderSize(0);
  for(int j=0;j<NFILES;j++){
    if(j==0)
      hB2_1dy[j][cent]->Draw();
    else hB2_1dy[j][cent]->Draw("same");
    //leg1->AddEntry(hB2dphi[j][cent],leg_names[j]);
    
  }
  //leg1->Draw("same");
  c3->cd(2);
  for(int j=0;j<NFILES;j++){
    if(j==0)
      hB2_1dphi[j][cent]->Draw();
    else hB2_1dphi[j][cent]->Draw("same");
  }
  c3->cd(3);
  for(int j=0;j<NFILES;j++){
    if(j==0)
      hB2_2dy[j][cent]->Draw();
    else hB2_2dy[j][cent]->Draw("same");
  }
  c3->cd(4);
  for(int j=0;j<NFILES;j++){
    if(j==0)
      hB2_2dphi[j][cent]->Draw();
    else hB2_2dphi[j][cent]->Draw("same");
  }

  //End drawing
  */ 
  
}

double GetR2Baseline(TH2D *hmult, bool fDistinguishable){
  double result   = 0.;
  if (!hmult){ cout<<"CalcRc::GetR2Baseline -- no hmult... exit"<<endl; exit(0); }
  double nent     = (double)hmult->GetEntries();

  TH2D *hmultwork = (TH2D*)hmult->Clone("hmultwork");
  hmultwork               ->Scale(1./nent);       // now its a probability density                                                         
  double sumnum           = 0;
  double sumden1          = 0;
  double sumden2          = 0;
  TAxis* multaxisx        = (TAxis*)hmultwork->GetXaxis();
  TAxis* multaxisy        = (TAxis*)hmultwork->GetYaxis();
  TH1D* hmultworkx        = (TH1D*)hmultwork->ProjectionX();
  double meanmult         = hmultworkx->GetMean();
  if (!fDistinguishable){
    for (int ibx=1;ibx<=multaxisx->GetNbins();ibx++){
      double n1i       = multaxisx->GetBinCenter(ibx);
      double pi        = hmultworkx->GetBinContent(ibx);
      sumnum          += pi*n1i*(n1i-1.);
      sumden1         += pi*n1i;
    }
    if ( sumnum>0.0 && sumden1>0.0 )
      result = (sumnum/sumden1/sumden1) - 1.0;       // R2integral 
  } else if (fDistinguishable){
    for (int ibx=1;ibx<=multaxisx->GetNbins();ibx++){
      for (int iby=1;iby<=multaxisy->GetNbins();iby++){
	double n1i       = multaxisx->GetBinCenter(ibx);
	double n2i       = multaxisy->GetBinCenter(iby);
	double pi        = hmultwork->GetBinContent(ibx,iby);
	sumnum          += pi*n1i*n2i;
	sumden1         += pi*n1i;//Mean mult f1=<n_a> - 1st part. 
	sumden2         += pi*n2i;//Mean mult f1=<n_a> - 2nd part.
      }
    }
    if (sumnum>0.0 && sumden1>0.0 && sumden2>0 )
      result = (sumnum/sumden1/sumden2) - 1.0;

  }       // end fDistinguishable... 
  delete hmultwork; hmultwork=0;
  return result;
}
