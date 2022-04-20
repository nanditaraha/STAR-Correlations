#include "fluct_common.h"
#include "TString.h"
#include "TLegend.h"
#include "TMath.h"

double ylimits[2];
void RangeFinderH(int, TH1D*);
void RangeFinderG(int, TGraph*);
void RangeFinderGE(int, TGraphErrors*);
void RangeCheck();

void B2_sys_errors(){
    
  const int NFILES = 16;//3 simulations + 2 data (mixing with 2 pref used) ... will add conv. later
  const int NPid = 3;//pions, kaons and protons...
    
  const char* pid_name[3] = {"pipi","KK","pp"};
  const char* leg_name[NFILES] = {"Cuts1","Cuts2","Cuts3","Cuts4","Cuts5","Cuts6","Cuts7","Cuts8","Cuts9","Cuts10","Cuts11","Cuts12","Cuts13","Cuts14","Cuts15","Cuts16"};
  //const char* leg_name[NFILES] = {"Urqmd","Ampt","Hijing","Data - pref from #LTN#GT","Data - published pref"};
  //int colors[11] = {kBlue,kGreen+2,kRed,1,kOrange+3,kMagenta,kCyan+2,kMagenta+2,kRed-3,15};
  //running out of colors :)
  const int color[NFILES] = {kBlue,kGreen,kRed,kMagenta,kRed-5,1,kOrange+3,kMagenta,kCyan+2,kMagenta+2,kRed-3,15,kGreen+2,kBlack,kBlue+2,kOrange};
  const int markerStyle[NFILES] = {43,33,42,27,34,4,44,45,34,25,21,40,41,29,20,43};
  const double markerSize[NPid] = {0.7,1.,1.};//,1.0,0.5,1.0};
  const double lineSize[NFILES] = {3,1};//,1.0,0.5,1.0};
  
  TH2D *hB2[NFILES][NCent][NPid]  = {0};
  TH1D *hB2dy[NFILES][NCent][NPid]  = {0};
  TH1D *hB2dphi[NFILES][NCent][NPid]  = {0};
  TH1D *hB2dy_avg[NCent][NPid]  = {0};
  TH1D *hB2dphi_avg[NCent][NPid]  = {0};
  //Histogram for pion B2(dy) vs. cuts, 0-5% centrality only for a test (bin 3)
  TH1D *hB2dy_bin = new TH1D("Bin_hist","#sigma_{i} of B_{2}(dy) vs. cuts 0-5% for a bin",NFILES,0.5,NFILES+0.5);
  
  TH1D *hB2dy_sys[NFILES][NCent][NPid]  = {0};
  TH1D *hB2dphi_sys[NFILES][NCent][NPid]  = {0};
  
  for(int ifile=0;ifile<NFILES;ifile++){
    TFile *f = new TFile(Form("root_files/B2_W1_c100_data_16_cuts%d_W1.root",ifile+1),"read"); f->cd();
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
	hB2dy[ifile][icent][ipaty]->SetMarkerColor(color[ifile]);
	//hB2dy[ifile][icent][ipaty]->SetLineWidth(lineSize[ifile]);
	hB2dy[ifile][icent][ipaty]->SetMarkerStyle(markerStyle[ifile]);
        //hB2dy[ifile][icent][ipaty]->SetMarkerSize(markerSize[ifile]);
	hB2dy[ifile][icent][ipaty]->GetXaxis()->SetTitleSize(0.06);
        hB2dy[ifile][icent][ipaty]->GetXaxis()->SetLabelSize(0.05);
        hB2dy[ifile][icent][ipaty]->GetYaxis()->SetLabelSize(0.05);
	
	hB2dphi[ifile][icent][ipaty] =(TH1D*)f->Get(Form("hB2cdphi_%d_%d",icent,ipaty));
        hB2dphi[ifile][icent][ipaty]->SetName(Form("hB2dphi_%d_%d_%d",ifile,icent,ipaty));
        hB2dphi[ifile][icent][ipaty]->SetStats(0);
        hB2dphi[ifile][icent][ipaty]->SetLineColor(color[ifile]);
        hB2dphi[ifile][icent][ipaty]->SetMarkerColor(color[ifile]);
	//hB2dphi[ifile][icent][ipaty]->SetLineWidth(lineSize[ifile]);
        hB2dphi[ifile][icent][ipaty]->SetMarkerStyle(markerStyle[ifile]);
        //hB2dphi[ifile][icent][ipaty]->SetMarkerSize(markerSize[ifile]);
	hB2dphi[ifile][icent][ipaty]->GetXaxis()->SetTitleSize(0.06);
        hB2dphi[ifile][icent][ipaty]->GetXaxis()->SetLabelSize(0.05);
        hB2dphi[ifile][icent][ipaty]->GetYaxis()->SetLabelSize(0.05);
    
	//Booking hists for sys-errors
	int nbins_dy = hB2dy[ifile][icent][ipaty]->GetNbinsX();
	double min_dy = hB2dy[ifile][icent][ipaty]->GetXaxis()->GetXmin();
	double max_dy = hB2dy[ifile][icent][ipaty]->GetXaxis()->GetXmax();
        hB2dy_sys[ifile][icent][ipaty] = new TH1D(Form("hB2dy_sys%d_%d",icent,ipaty),"B_{2}(dy)",nbins_dy,min_dy,max_dy);
        hB2dy_sys[ifile][icent][ipaty]        ->SetMarkerStyle(20);
        hB2dy_sys[ifile][icent][ipaty]        ->SetMarkerColor( kRed);
        hB2dy_sys[ifile][icent][ipaty]        ->SetLineColor  ( kRed);
        //hB2dy_sys[ifile][icent][ipaty]        ->SetLineWidth  ( 2);
	hB2dy_sys[ifile][icent][ipaty]->SetMarkerSize(0.7);
	hB2dy_sys[ifile][icent][ipaty]->GetXaxis()->SetTitleSize(0.06);
        hB2dy_sys[ifile][icent][ipaty]->GetXaxis()->SetLabelSize(0.05);
        hB2dy_sys[ifile][icent][ipaty]->GetYaxis()->SetLabelSize(0.05);
	hB2dy_sys[ifile][icent][ipaty]->GetXaxis()->SetTitle("#Deltay");

        int nbins_dphi = hB2dphi[ifile][icent][ipaty]->GetNbinsX();
	double min_dphi = hB2dphi[ifile][icent][ipaty]->GetXaxis()->GetXmin();
        double max_dphi = hB2dphi[ifile][icent][ipaty]->GetXaxis()->GetXmax();	
        hB2dphi_sys[ifile][icent][ipaty] = new TH1D(Form("hB2dphi_sys%d_%d",icent,ipaty),"B_{2}(d#phi)",nbins_dphi,min_dphi,max_dphi);
        hB2dphi_sys[ifile][icent][ipaty]        ->SetMarkerStyle(20);
        hB2dphi_sys[ifile][icent][ipaty]        ->SetMarkerColor( kRed);
        hB2dphi_sys[ifile][icent][ipaty]        ->SetLineColor  ( kRed);
        //hB2dphi_sys[ifile][icent][ipaty]        ->SetLineWidth  ( 2);
	hB2dphi_sys[ifile][icent][ipaty]->SetMarkerSize(0.7);
	hB2dphi_sys[ifile][icent][ipaty]->GetXaxis()->SetTitleSize(0.06);
        hB2dphi_sys[ifile][icent][ipaty]->GetXaxis()->SetLabelSize(0.05);
        hB2dphi_sys[ifile][icent][ipaty]->GetYaxis()->SetLabelSize(0.05);
	hB2dphi_sys[ifile][icent][ipaty]->GetXaxis()->SetTitle("#Delta#phi");
      }
    }
  }

  //Booking avg hists:
  for (int icent=1; icent<NCent; icent++){
      for (int ipaty=0;ipaty<NPid;ipaty++){
	int nbins_dy = hB2dy[0][icent][ipaty]->GetNbinsX();
        double min_dy = hB2dy[0][icent][ipaty]->GetXaxis()->GetXmin();
        double max_dy = hB2dy[0][icent][ipaty]->GetXaxis()->GetXmax();
        hB2dy_avg[icent][ipaty] = new TH1D(Form("hB2dy_avg%d_%d",icent,ipaty),"B_{2}(dy) - Average of all cuts",nbins_dy,min_dy,max_dy);
        hB2dy_avg[icent][ipaty]        ->SetMarkerStyle(20);
        hB2dy_avg[icent][ipaty]        ->SetMarkerColor( kBlue);
        hB2dy_avg[icent][ipaty]        ->SetLineColor  ( kBlue);
        hB2dy_avg[icent][ipaty]        ->SetLineWidth  ( 2);
	hB2dy_avg[icent][ipaty]->GetXaxis()->SetTitle("dy");
        hB2dy_avg[icent][ipaty]->SetMarkerSize(0.7);
  
	int nbins_dphi = hB2dphi[0][icent][ipaty]->GetNbinsX();
	double min_dphi = hB2dphi[0][icent][ipaty]->GetXaxis()->GetXmin();
	double max_dphi = hB2dphi[0][icent][ipaty]->GetXaxis()->GetXmax();
	hB2dphi_avg[icent][ipaty] = new TH1D(Form("hB2dphi_avg%d_%d",icent,ipaty),"B_{2}(d#phi) - Average of all cuts",nbins_dphi,min_dphi,max_dphi);
	hB2dphi_avg[icent][ipaty]        ->SetMarkerStyle(20);
	hB2dphi_avg[icent][ipaty]        ->SetMarkerColor( kBlue);
	hB2dphi_avg[icent][ipaty]        ->SetLineColor  ( kBlue);
	hB2dphi_avg[icent][ipaty]        ->SetLineWidth  ( 2);
	hB2dphi_avg[icent][ipaty]->GetXaxis()->SetTitle("d#phi");
	hB2dphi_avg[icent][ipaty]->SetMarkerSize(0.7);
      }
  }
  const int NbinX = hB2dy[1][1][1]->GetXaxis()->GetNbins();
  const int NbinY = hB2dphi[1][1][1]->GetXaxis()->GetNbins();
  //Stat errors and values//
  double dy_stat_x[NFILES][NCent][NPid][NbinX];//bin centers
  //the x-axis errors are zero for all values (dy, dphi) both stat and sys. So the one below will be used for all four.
  //double err_dy_x[NFILES][NCent][NPid][NbinY];//Made for all files and is zero
  double sys_err_dy[NCent][NPid][NbinX];  //Actual errors in bin contents - sigma of 14 files and same for all. 
  //double dy_sys_x[NFILES][NCent][NPid][NbinX];//bin centers
  double dy_sys_y[NFILES][NCent][NPid][NbinX];//bin contents
  
  double sys_err_dphi[NCent][NPid][NbinY];
  //double dphi_sys_x[NFILES][NCent][NPid][NbinY];
  double dphi_sys_y[NFILES][NCent][NPid][NbinY];
  
  //Calculating sys error for B2(dy) for each bin - which is the sigma of B2(dy or dphi) for 14 cut sets
  for (int icent=1; icent<NCent; icent++){
    for (int ipaty=0;ipaty<NPid;ipaty++){
      for (int ibinx=0;ibinx<NbinX;ibinx++){
	double dy_var=0, dy_avg=0;
	for(int ifile=0;ifile<NFILES;ifile++)
	  dy_avg=hB2dy[ifile][icent][ipaty]->GetBinError(ibinx+1)/NFILES;
	for(int ifile=0;ifile<NFILES;ifile++){
	  //Finding average quadratures of errors...
	  //dy_var+=hB2dy[ifile][icent][ipaty]->GetBinError(ibinx+1)*hB2dy[ifile][icent][ipaty]->GetBinError(ibinx+1)/NFILES;
	  //finding sigma...
	  dy_var+=(hB2dy[ifile][icent][ipaty]->GetBinError(ibinx+1)-dy_avg)*
	  (hB2dy[ifile][icent][ipaty]->GetBinError(ibinx+1)-dy_avg)/NFILES;
	  /*
	  if(ibinx==2 && icent==16 && ipaty==0) {//hB2dy_bin->Fill(hB2dy[ifile][icent][ipaty]->GetBinContent(ibinx+1));
	    hB2dy_bin->SetBinContent(ifile+1,hB2dy[ifile][icent][ipaty]->GetBinContent(ibinx+1));
	    hB2dy_bin->SetBinError(ifile+1,hB2dy[ifile][icent][ipaty]->GetBinError(ibinx+1));
	    cout<<hB2dy_bin->GetBinContent(ifile+1)<<endl;
	    }*/
	}
	//if(ibinx==1 && icent==16 && ipaty==0) cout<<hB2dy[ifile][icent][ipaty]->GetBinError(ibinx+1)<<endl;        	
	sys_err_dy[icent][ipaty][ibinx] = sqrt(dy_var);
	if(ibinx==2 && icent==16 && ipaty==0) cout<<"Sigma="<<sys_err_dy[icent][ipaty][ibinx]<<endl;
	for(int ifile=0;ifile<NFILES;ifile++){
	  if(ibinx==2 && icent==16 && ipaty==0) {//hB2dy_bin->Fill(hB2dy[ifile][icent][ipaty]->GetBinContent(ibinx+1));
            //hB2dy_bin->SetBinContent(ifile+1,hB2dy[ifile][icent][ipaty]->GetBinContent(ibinx+1));
            hB2dy_bin->SetBinContent(ifile+1,hB2dy[ifile][icent][ipaty]->GetBinError(ibinx+1));
	    hB2dy_bin->SetBinError(ifile+1,sys_err_dy[icent][ipaty][ibinx]);
            cout<<hB2dy_bin->GetBinContent(ifile+1)<<endl;
	  }
	}
      }//end of B2(dy) errors

      //Similarly for B2(dphi)
      for (int ibiny=0;ibiny<NbinY;ibiny++){
	double dphi_avg=0; double dphi_var=0;
	for(int ifile=0;ifile<NFILES;ifile++)  dphi_avg+= hB2dphi[ifile][icent][ipaty]->GetBinContent(ibiny+1)/NFILES;
	for(int ifile=0;ifile<NFILES;ifile++)
	  dphi_var+=hB2dphi[ifile][icent][ipaty]->GetBinError(ibiny+1)*hB2dphi[ifile][icent][ipaty]->GetBinError(ibiny+1)/NFILES;
	//dphi_var+=(hB2dphi[ifile][icent][ipaty]->GetBinContent(ibiny+1)-dphi_avg)*
	//(hB2dphi[ifile][icent][ipaty]->GetBinContent(ibiny+1)-dphi_avg)/NFILES;
	sys_err_dphi[icent][ipaty][ibiny] = sqrt(dphi_var);
      }
    }
  }

  //Adding the errors and points to sys-graph  
  for(int ifile=0;ifile<NFILES;ifile++){
    for (int icent=1; icent<NCent; icent++){
      for (int ipaty=0;ipaty<NPid;ipaty++){
	double e_dy1[NbinX],e_dy[NbinX];
	for (int ibinx=0;ibinx<NbinX;ibinx++){	  
	  //err_dy_x[ifile][icent][ipaty][ibinx] = 0;
	  //dy_sys_x[ifile][icent][ipaty][ibinx] = hB2dy[ifile][icent][ipaty]->GetBinCenter(ibinx+1);
	  //dy_sys_y[ifile][icent][ipaty][ibinx] = hB2dy[ifile][icent][ipaty]->GetBinContent(ibinx+1);	    
	  hB2dy_sys[ifile][icent][ipaty]->SetBinContent(ibinx+1,hB2dy[ifile][icent][ipaty]->GetBinContent(ibinx+1));
	  hB2dy_sys[ifile][icent][ipaty]->SetBinError(ibinx+1,sys_err_dy[icent][ipaty][ibinx]);
	  //B2dy_sys[ifile][icent][ipaty]->SetPointError(ibinx,err_dy_x[ifile][icent][ipaty][ibinx],sys_err_dy[icent][ipaty][ibinx]);
	}
	//Same for dphi
	for (int ibiny=0;ibiny<NbinY;ibiny++){
	  //err_dy_x[ifile][icent][ipaty][ibiny]=0;
	  //dphi_sys_x[ifile][icent][ipaty][ibiny] = hB2dphi[ifile][icent][ipaty]->GetBinCenter(ibiny+1);
	  //dphi_sys_y[ifile][icent][ipaty][ibiny] = hB2dphi[ifile][icent][ipaty]->GetBinContent(ibiny+1);
	 hB2dphi_sys[ifile][icent][ipaty]->SetBinContent(ibiny+1,hB2dphi[ifile][icent][ipaty]->GetBinContent(ibiny+1));
	 hB2dphi_sys[ifile][icent][ipaty]->SetBinError(ibiny+1,sys_err_dphi[icent][ipaty][ibiny]);
	}
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
  
  for (int ipart=0;ipart<3;ipart++){
    int can_index=1;
    ccan[ipart]  = new TCanvas("","",20,30+20,750,900);
    ccan[ipart]->Divide(2,2,0.0001,0.0001);
    ccan[ipart]->cd(can_index);
    hB2[0][16][ipart]->Draw("lego2");
    can_index++;ccan[ipart]->cd(can_index);
    hB2[0][16][ipart]->Draw("colz");
    can_index++;ccan[ipart]->cd(can_index);

    //hB2dy[0][16][ipart]->SetBarWidth(1);
    //hB2dy[0][16][ipart]->SetBarOffset(0.25);
    //hB2dy[0][16][ipart]->SetLineColor(kRed);

    //To draw gray boxes around sys errors, use default fill style and use a fill color of kGray!
    //hB2dy_sys[0][16][ipart]->SetFillColor(kGray); 
    
    //This fill style just draws boxes around the errors
    hB2dy[0][16][ipart]->SetFillStyle(0);
    hB2dy_sys[0][16][ipart]->SetFillStyle(0);
    hB2dphi_sys[0][16][ipart]->SetFillStyle(0);
    hB2dphi[0][16][ipart]->SetFillStyle(0);
    //hB2dy_sys[0][16][ipart]->Draw();
    
    hB2dy[0][16][ipart]->Draw("e1,x0");
    hB2dy_sys[0][16][ipart]->Draw("same,e2");
    //hB2dy[0][16][ipart]->Draw("same");
    TLegend *legend1= new TLegend(0.3,0.9,0.6,0.75);
    legend1->SetTextFont(43);
    legend1->SetTextSize(15);
    legend1->SetFillStyle(0);
    legend1->AddEntry(hB2dy_sys[0][16][ipart],"Systematic Error");
    legend1->AddEntry(hB2dy[0][16][ipart],"Statistical Error");
    legend1->Draw(); 
    can_index++;ccan[ipart]->cd(can_index);
    //gB2dphi_stat[0][16][ipart]->Draw("AP");
    
    hB2dphi[0][16][ipart]->Draw("e1,x0");
    hB2dphi_sys[0][16][ipart]->Draw("same,e2");

  }
  
  //Drawing all cuts for comparison - just for looking - not used for calculations!!
  /*
  for (int ipart=0;ipart<3;ipart++){
    int can_index=1;
    ccan[ipart]  = new TCanvas("","",0,0,820,550);
    ccan[ipart]->SetWindowSize(820 + (820 - ccan[ipart]->GetWw()), 550 + (550 - 100 - ccan[ipart]->GetWh()));
    ccan[ipart]->Divide(2,1);//,0.0001,0.0001);

    TLegend *legend= new TLegend(0.6,0.9,0.6,0.75);
    legend->SetTextFont(43);
    legend->SetTextSize(15);
    legend->SetFillStyle(0);
    
    TLegend *legend1= new TLegend(0.3,0.9,0.6,0.75);
    legend1->SetTextFont(43);
    legend1->SetTextSize(15);
    legend1->SetFillStyle(0);

    double dy_min[3]={0.05,0.03,0.06};
    double dy_max[3]={0.3,0.126,0.28};
    double dphi_min[3]={-0.1,-0.07,0.0};
    double dphi_max[3]={0.4,0.25,0.3};
    
    
    //Drawing all 16 files overlaid
    for(int ifile=0;ifile<NFILES;ifile++){
      ccan[ipart]->cd(can_index);
      hB2dy[ifile][16][ipart]->GetYaxis()->SetRangeUser(dy_min[ipart],dy_max[ipart]);
      if(ifile==0)
	hB2dy[ifile][16][ipart]->Draw();
      else hB2dy[ifile][16][ipart]->Draw("same");
      if(ipart==0){
	if(ifile<NFILES/2)
	  legend->AddEntry(hB2dy[ifile][16][ipart],leg_name[ifile]);
	else legend1->AddEntry(hB2dy[ifile][16][ipart],leg_name[ifile]);	
      }
    }
    if(ipart==0){
      legend->Draw();
      legend1->Draw();
    }
    
    can_index++; ccan[ipart]->cd(can_index);
    for(int ifile=0;ifile<NFILES;ifile++){
      hB2dphi[ifile][16][ipart]->GetYaxis()->SetRangeUser(dphi_min[ipart],dphi_max[ipart]);
      if(ifile==0)
	hB2dphi[ifile][16][ipart]->Draw();
      else hB2dphi[ifile][16][ipart]->Draw("same");
      
    }
  }//end of overlaid files
    
  //Drawing Avg. hist ... - just looking not used anywhere in calculations!!
  //int can_index=1;
  auto cc  = new TCanvas("","",20,30+20,1150,820);
  cc->Divide(3,2,0.0001,0.0001);
  for (int ipart=0;ipart<3;ipart++){
    for(int ifile=0;ifile<NFILES;ifile++){
      hB2dy_avg[16][ipart]->Add(hB2dy[ifile][16][ipart]);
      hB2dphi_avg[16][ipart]->Add(hB2dphi[ifile][16][ipart]);
    }
    hB2dy_avg[16][ipart]->Scale(1/NFILES);
    hB2dphi_avg[16][ipart]->Scale(1/NFILES);
  }
  for (int ipart=0;ipart<3;ipart++){
    cc->cd(ipart+1);
    hB2dy_avg[16][ipart]->Draw();
  }
 for (int ipart=0;ipart<3;ipart++){
    cc->cd(ipart+4);
    hB2dphi_avg[16][ipart]->Draw();
  }
  */
 //An illustration showing a bin dist for all 16 files!
  hB2dy_bin->Draw("e1");
  //TFile *fout   = new TFile("test.root","recreate");
  //hB2dy_bin->Write();
  cout<<"Bin rms="<<hB2dy_bin->GetRMS(2)<<endl;

  cout<<"B2dphi ratio: bin 1"<<endl;
  for(int i=0;i<3;i++)
    cout<<hB2dphi_sys[0][16][i]->GetBinError(1)/hB2dphi[0][16][i]->GetBinError(1)<<endl;
  cout<<"B2dphi ratio: bin 5"<<endl;
  for(int i=0;i<3;i++)
    cout<<hB2dphi_sys[0][16][i]->GetBinError(5)/hB2dphi[0][16][i]->GetBinError(5)<<endl;
  cout<<"B2dphi ratio: bin 10"<<endl;
  for(int i=0;i<3;i++)
    cout<<hB2dphi_sys[0][16][i]->GetBinError(10)/hB2dphi[0][16][i]->GetBinError(10)<<endl;
  cout<<"B2dphi ratio: bin 15"<<endl;
  for(int i=0;i<3;i++)
    cout<<hB2dphi_sys[0][16][i]->GetBinError(15)/hB2dphi[0][16][i]->GetBinError(15)<<endl;
  cout<<"B2dphi ratio: bin 19"<<endl;
  for(int i=0;i<3;i++)
    cout<<hB2dphi_sys[0][16][i]->GetBinError(19)/hB2dphi[0][16][i]->GetBinError(19)<<endl;
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

