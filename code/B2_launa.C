#include "fluct_common.h"
#include "TString.h"
#include "TLegend.h"
#include "TMath.h"

double ylimits[2];
void RangeFinderH(int, TH1D*);
void RangeFinderG(int, TGraph*);
void RangeFinderGE(int, TGraphErrors*);
void RangeCheck();
bool fDistinguishable;
double multinfo[10];
void GetNPairs(TH2D* h2); 
TGraph *gxs_pim,*gxs_pip,*gxs_km,*gxs_kp,*gxs_pm,*gxs_pp;
TGraph *gxs_pime,*gxs_pipe,*gxs_kme,*gxs_kpe,*gxs_pme,*gxs_ppe;
double xsinfo[2];
void GetdNdy(double c,int k);
void BuildB2(double p1, double p1e,
	     double p2, double p2e,
	     TH2D* r2a,TH2D* r2b, TH2D* r2c, TH2D* r2d,
	     TH2D* b2a,TH2D* b2b, TH2D* b2c);
void GenericProject(TH2D* h2,TString steerav,TString steerax,TH1D* h1);


//----------------------------------------------------------------------------------------
//	To plot Balance Functions	(for data)
//		
//	This version:
//		-using more stable equation for B2 uncertainty propagation
//		-uses either real data or simulation data input (set SIMINPUT at top)
//			if siminput, will get prefactors from hy_particles_kept2
//		-calculates also the "sum rules" for the B2 integrals as cross-check 
//----------------------------------------------------------------------------------------


void B2(){


	//DATA
//	TString* onefile	= new TString("/Users/launaasllanaj/Desktop/analysis2022/readerroot2021_22/reader_m40_data_24_cuts1_xp.root");	//kDS==24 (run-11 auau 200 GeV) 
 	TString* onefile	= new TString("~/RHIC/root_files/reader_m40_data_24_cuts1_xp.root");	 

	
	//SIMULATION!!
//	TString* onefile	= new TString("/Users/launaasllanaj/Desktop/analysis2022/readerroot2021_22/reader_c16_urqmdLN_16_cuts1.root");
  	
  	
  	int kCS			= 1;		// cuts set 
  	bool SIMINPUT	= false;	// set to TRUE if simulation input

	int colors[11]		= {1,kBlue,kBlue+1,kGreen+2,kCyan+1,kCyan+2,kMagenta+2,kRed-3,15,kRed+2,kRed};
  	//
  	//
  	TString ofstring	= TString(onefile->Data());
  	int ind		= ofstring.Index("reader_") + 6; 
  	ofstring.Remove(0,ind);
  	ind		= ofstring.Index(".root"); 
  	ofstring.Remove(ind,ind+5);
	TString dir			= TString("../ps/");
  	TString base		= TString("B2") + ofstring;
  	TString outfile		= dir + base + TString(".ps");
  	TString outfileO	= dir + base + TString(".ps(");
  	TString outfileC	= dir + base + TString(".ps]");
  	TString outfileP	= dir + base + TString(".pdf");
  	dir					= TString("../root/");
  	TString outfileroot	= dir + base + TString(".root");
  	cout<<"psfile name   = "<<outfile.Data()<<endl;
	cout<<"rootfile name = "<<outfileroot.Data()<<endl;
  	//
  	//
	TH2D* hmult[NCent][NPairTypes]		= {0};
  	TH1D* hmult1[NCent][NPairTypes]		= {0};
  	TH1D* hmult2[NCent][NPairTypes]		= {0};
  	TH2D* hrho2[NCent][NPairTypes]		= {0};
  	TH2D* hrho1rho1[NCent][NPairTypes]	= {0};
  	TH2D* hR2[NCent][NPairTypes]		= {0};
  	//
  	double aintrho1rho1[NCent][NPairTypes]	= {0};
  	double aintrho1rho1e[NCent][NPairTypes]	= {0};
  	double aintrho2[NCent][NPairTypes]		= {0};
  	double aintrho2e[NCent][NPairTypes]		= {0};
  	//
  	TH2D* hic_rho1rho1_[NPairTypes]		= {0};
  	TH2D* hic_rho2_[NPairTypes]			= {0};
  	TH2D* hic_rho1rho1_rat[NPairTypes]	= {0};
  	TH2D* hic_rho2_rat[NPairTypes]		= {0};
  	const int anbm	= 500;
  	double ax1m		=     0.001;
  	double ax2m		= 10000.001;
  	double logax1m	= TMath::Log10(ax1m);
  	double logax2m	= TMath::Log10(ax2m);
  	double abwm		= (logax2m-logax1m)/((double)anbm);
  	double xbm[anbm+1];
  	xbm[0]			= ax1m;
  	for (int i=1;i<=anbm;i++){
    	xbm[i]		= TMath::Power(10,logax1m+i*abwm);
  	}
  	for (int ipaty=0;ipaty<NPairTypes;ipaty++){
    	hic_rho1rho1_[ipaty]	= new TH2D(Form("hic_rho1rho1_%d",ipaty),
					   Form("hic_rho1rho1_%d",ipaty),anbm,xbm,anbm,xbm);
    	hic_rho2_[ipaty] 		= new TH2D(Form("hic_rho2_%d",ipaty),
					   Form("hic_rho2_%d",ipaty),anbm,xbm,anbm,xbm);
    	hic_rho1rho1_rat[ipaty]	= new TH2D(Form("hic_rho1rho1_rat%d",ipaty),
					   Form("hic_rho1rho1_rat%d",ipaty),anbm,xbm,400,0,2);
    	hic_rho2_rat[ipaty] 	= new TH2D(Form("hic_rho2_rat%d",ipaty),
					   Form("hic_rho2_rat%d",ipaty),anbm,xbm,400,0,2);
  	}
  	//
  	TH2D *hR2c[NCent][NPairTypes]		= {0};	//Corrected R2
  	TH1D *hR2cdy[NCent][NPairTypes]		= {0};	//Corrected R2dy
  	TH1D *hR2cdphi[NCent][NPairTypes]	= {0};	//Corrected R2dphi
  	TH2D *hrho2c[NCent][NPairTypes]		= {0};	//Corrected rho2
  	TH1D *hrho2cdy[NCent][NPairTypes]	= {0};	//Corrected rho2dy
  	TH1D *hrho2cdphi[NCent][NPairTypes]	= {0};	//Corrected rho2dphi
  	//
  	TH2D *hB2cp[NCent][3]		= {0};		// B2+
  	TH2D *hB2cm[NCent][3]		= {0};		// B2-
  	TH2D *hB2c[NCent][3]		= {0};		// B2
  	TH1D *hB2cpdy[NCent][3]		= {0};		// B2+dy
  	TH1D *hB2cmdy[NCent][3]		= {0};		// B2-dy
  	TH1D *hB2cdy[NCent][3]		= {0};		// B2dy
  	TH1D *hB2cpdphi[NCent][3]	= {0};		// B2+dphi
  	TH1D *hB2cmdphi[NCent][3]	= {0};		// B2-dphi
  	TH1D *hB2cdphi[NCent][3]	= {0};		// B2dphi
  	//
  	double aintB2cp[NCent][3]	= {0};		// B2+ integral
  	double aintB2cm[NCent][3]	= {0};		// B2- integral
  	double aintB2c[NCent][3]	= {0};		// B2 integral
 	double aintB2cpe[NCent][3]	= {0};		// B2+ integral uncertainty
  	double aintB2cme[NCent][3]	= {0};		// B2- integral uncertainty
  	double aintB2ce[NCent][3]	= {0};		// B2 integral uncertainty
    //
    double armsB2cdyp[NCent][3]	= {0};		// B2+ dy rms
  	double armsB2cdym[NCent][3]	= {0};		// B2- dy rms
  	double armsB2cdy[NCent][3]	= {0};		// B2  dy rms
 	double armsB2cdype[NCent][3]= {0};		// B2+ dy rms uncertainty
  	double armsB2cdyme[NCent][3]= {0};		// B2- dy rms uncertainty
  	double armsB2cdye[NCent][3] = {0};		// B2  dy rms uncertainty
	double armsB2cdphip[NCent][3]= {0};		// B2+ dphi rms
  	double armsB2cdphim[NCent][3]= {0};		// B2- dphi rms
  	double armsB2cdphi[NCent][3] = {0};		// B2  dphi  rms
 	double armsB2cdphipe[NCent][3]= {0};	// B2+ dphi rms uncertainty
  	double armsB2cdphime[NCent][3]= {0};	// B2- dphi rms uncertainty
  	double armsB2cdphie[NCent][3] = {0};	// B2  dphi rms uncertainty
	//
  	double nevents[NCent][NPairTypes] 	= {0};
  	double npairs[NCent][NPairTypes] 	= {0};
  	double mult[NCent][NPARTICLENAMES]	= {0};		//Prefactors from hy_particles_kept (SIM data only!)
  	bool   multdone[NCent][NPARTICLENAMES]	= {0};
	//
	TH1D*  hy_particles_kept2[NPARTICLENAMES][NCent] = {0};		//LOOK AT INDEXING HERE (originates in reader_loop)!!!!!!!!
  	//
//   TH1D *hB2cp_dy[NCent][3]	= {0};   
//   TH1D *hB2cm_dy[NCent][3]	= {0};   
//   TH1D *hB2c_dy[NCent][3]		= {0};   
// 	TH1D *hB2cp_dphi[NCent][3]	= {0};   
// 	TH1D *hB2cm_dphi[NCent][3]	= {0};   
// 	TH1D *hB2c_dphi[NCent][3]	= {0};   

  	//
  	//
  	//======== Sum Rule Graphs...
  	double C2int[NCent][NPairTypes]	= {0};	//Integral(C_2) obtained from hmult 
  	TGraphErrors* gintB2cp_cent[3]	= {0};	//Integral obtained from B2+
  	TGraphErrors* gintB2cm_cent[3]	= {0};	//Integral obtained from B2-
  	TGraphErrors* gintB2c_cent[3]	= {0};	//Integral obtained from B2
  	TGraphErrors* gsruB2cp_cent[3]	= {0};	//Integral from sum rule B2+
  	TGraphErrors* gsruB2cm_cent[3]	= {0};	//Integral from sum rule B2-
  	TGraphErrors* gsruB2c_cent[3]	= {0};	//Integral from sum rule B2
	//
	TGraphErrors* grmsB2cp_dy_cent[3]	= {0};	
	TGraphErrors* grmsB2cm_dy_cent[3]	= {0};	
	TGraphErrors* grmsB2c_dy_cent[3]	= {0};	
	TGraphErrors* grmsB2cp_dphi_cent[3]	= {0};	
	TGraphErrors* grmsB2cm_dphi_cent[3]	= {0};	
	TGraphErrors* grmsB2c_dphi_cent[3]	= {0};	
	//
	//
	//
  	const char* partname[3] = {"#pi#pi","KK","pp"};
  	//
  	for (int ipart=0;ipart<3;ipart++){		// pi, then K, then p
    	gintB2cp_cent[ipart]	= new TGraphErrors();
    	gintB2cp_cent[ipart]	->SetTitle(Form("Integral %s B_{2+};CentralityBin",partname[ipart]));
    	gintB2cp_cent[ipart]	->SetMarkerStyle(20);
    	gintB2cp_cent[ipart]	->SetMarkerColor( 4);
    	gintB2cp_cent[ipart]	->SetLineColor  ( 4);
    	gintB2cp_cent[ipart]	->SetLineWidth  ( 2);
    	//
    	gintB2cm_cent[ipart]	= new TGraphErrors();
    	gintB2cm_cent[ipart]	->SetTitle(Form("Integral %s B_{2-};CentralityBin",partname[ipart]));
    	gintB2cm_cent[ipart]	->SetMarkerStyle(20);
    	gintB2cm_cent[ipart]	->SetMarkerColor( 4);
    	gintB2cm_cent[ipart]	->SetLineColor  ( 4);
    	gintB2cm_cent[ipart]	->SetLineWidth  ( 2);
    	//
    	gintB2c_cent[ipart]	= new TGraphErrors();
    	gintB2c_cent[ipart]	->SetTitle(Form("Integral %s B_{2};CentralityBin",partname[ipart]));
    	gintB2c_cent[ipart]	->SetMarkerStyle(20);
    	gintB2c_cent[ipart]	->SetMarkerColor( 4);
    	gintB2c_cent[ipart]	->SetLineColor  ( 4);
    	gintB2c_cent[ipart]	->SetLineWidth  ( 2);
    	//
    	//
    	gsruB2cp_cent[ipart]	= new TGraphErrors();
    	gsruB2cp_cent[ipart]	->SetTitle(Form("Sum Rule %s B_{2+};CentralityBin",partname[ipart]));
    	gsruB2cp_cent[ipart]	->SetMarkerStyle(21);
    	gsruB2cp_cent[ipart]	->SetMarkerColor(kGreen+2);
    	gsruB2cp_cent[ipart]	->SetLineColor  (kGreen+2);
    	gsruB2cp_cent[ipart]	->SetLineWidth  ( 2);
    	//
    	gsruB2cm_cent[ipart]	= new TGraphErrors();
    	gsruB2cm_cent[ipart]	->SetTitle(Form("Sum Rule %s B_{2-};CentralityBin",partname[ipart]));
    	gsruB2cm_cent[ipart]	->SetMarkerStyle(21);
    	gsruB2cm_cent[ipart]	->SetMarkerColor(kGreen+2);
    	gsruB2cm_cent[ipart]	->SetLineColor  (kGreen+2);
    	gsruB2cm_cent[ipart]	->SetLineWidth  ( 2);
  	    //
  	    gsruB2c_cent[ipart]	= new TGraphErrors();
    	gsruB2c_cent[ipart]	->SetTitle(Form("Sum Rule %s B_{2-};CentralityBin",partname[ipart]));
    	gsruB2c_cent[ipart]	->SetMarkerStyle(21);
    	gsruB2c_cent[ipart]	->SetMarkerColor(kGreen+2);
    	gsruB2c_cent[ipart]	->SetLineColor  (kGreen+2);
    	gsruB2c_cent[ipart]	->SetLineWidth  ( 2);
		//
		//
    	grmsB2cp_dy_cent[ipart]	= new TGraphErrors();
    	grmsB2cp_dy_cent[ipart]	->SetTitle(Form("RMS %s #LTB_{2+}#GT vs. #Deltay;CentralityBin",partname[ipart]));
    	grmsB2cp_dy_cent[ipart]	->SetMarkerStyle(20);
    	grmsB2cp_dy_cent[ipart]	->SetMarkerColor( 2);
    	grmsB2cp_dy_cent[ipart]	->SetLineColor  ( 2);
    	grmsB2cp_dy_cent[ipart]	->SetLineWidth  ( 2);
		//
    	grmsB2cm_dy_cent[ipart]	= new TGraphErrors();
    	grmsB2cm_dy_cent[ipart]	->SetTitle(Form("RMS %s #LTB_{2-}#GT vs. #Deltay;CentralityBin",partname[ipart]));
    	grmsB2cm_dy_cent[ipart]	->SetMarkerStyle(20);
    	grmsB2cm_dy_cent[ipart]	->SetMarkerColor( 2);
    	grmsB2cm_dy_cent[ipart]	->SetLineColor  ( 2);
    	grmsB2cm_dy_cent[ipart]	->SetLineWidth  ( 2);
		//
    	grmsB2c_dy_cent[ipart]	= new TGraphErrors();
    	grmsB2c_dy_cent[ipart]	->SetTitle(Form("RMS %s #LTB_{2}#GT vs. #Deltay;CentralityBin",partname[ipart]));
    	grmsB2c_dy_cent[ipart]	->SetMarkerStyle(20);
    	grmsB2c_dy_cent[ipart]	->SetMarkerColor( 2);
    	grmsB2c_dy_cent[ipart]	->SetLineColor  ( 2);
    	grmsB2c_dy_cent[ipart]	->SetLineWidth  ( 2);
		//
    	grmsB2cp_dphi_cent[ipart]	= new TGraphErrors();
    	grmsB2cp_dphi_cent[ipart]	->SetTitle(Form("RMS %s #LTB_{2+}#GT vs. #Delta#varphi;CentralityBin",partname[ipart]));
    	grmsB2cp_dphi_cent[ipart]	->SetMarkerStyle(20);
    	grmsB2cp_dphi_cent[ipart]	->SetMarkerColor( 2);
    	grmsB2cp_dphi_cent[ipart]	->SetLineColor  ( 2);
    	grmsB2cp_dphi_cent[ipart]	->SetLineWidth  ( 2);
		//
    	grmsB2cm_dphi_cent[ipart]	= new TGraphErrors();
    	grmsB2cm_dphi_cent[ipart]	->SetTitle(Form("RMS %s #LTB_{2-}#GT vs.#Delta#varphi;CentralityBin",partname[ipart]));
    	grmsB2cm_dphi_cent[ipart]	->SetMarkerStyle(20);
    	grmsB2cm_dphi_cent[ipart]	->SetMarkerColor( 2);
    	grmsB2cm_dphi_cent[ipart]	->SetLineColor  ( 2);
    	grmsB2cm_dphi_cent[ipart]	->SetLineWidth  ( 2);
		//
    	grmsB2c_dphi_cent[ipart]	= new TGraphErrors();
    	grmsB2c_dphi_cent[ipart]	->SetTitle(Form("RMS %s #LTB_{2}#GT vs.#Delta#varphi ;CentralityBin",partname[ipart]));
    	grmsB2c_dphi_cent[ipart]	->SetMarkerStyle(20);
    	grmsB2c_dphi_cent[ipart]	->SetMarkerColor( 2);
    	grmsB2c_dphi_cent[ipart]	->SetLineColor  ( 2);
    	grmsB2c_dphi_cent[ipart]	->SetLineWidth  ( 2);
		//

  	}
  	//
  	GetdNdy(-1,0);	//initialize xsec routine (not used if sim input, only used if data input)
  	//
  	//
  	//
  	TString* infile	= onefile;
  	cout<<"Opening "<<infile->Data()<<endl;
  	TFile *f		= new TFile(infile->Data(),"read"); f->cd();
	//
 	//for (int ipaty=0;ipaty<NPairTypes;ipaty++){
 	for (int ipaty=0;ipaty<12;ipaty++){		//NPairTypes=12
    	int kPID1			= PairTypes_Info[ipaty][0];
    	int kPID2			= PairTypes_Info[ipaty][1];
    	fDistinguishable	= true;
    	if (kPID1==kPID2) fDistinguishable	= false;
    	//
    	int	ir2;
    	for (int icent=1;icent<NCent;icent++){
      		TString* fstri	= new TString(Form("_%d_%d",icent,ipaty));
      		TString* fstro	= new TString(Form("_%d_%d",icent,ipaty));	// same! (now that "ie" loop removed...
      		//
      		hmult[icent][ipaty]		= (TH2D*)f->Get(Form("hmult%s",fstri->Data()));
      		//cout<<icent<<" "<<ipaty<<" "<<ie<<" "<<hmult2D[icent][ipaty]<<endl;
      		hmult[icent][ipaty]		->SetDirectory(0);
      		hmult[icent][ipaty]		->SetName(Form("hmult%s",fstro->Data()));
      		//
      		hmult1[icent][ipaty]	= (TH1D*)hmult[icent][ipaty]->ProjectionX();
      		hmult1[icent][ipaty]	->SetDirectory(0);
      		hmult1[icent][ipaty]	->SetName(Form("hmult1%s",fstro->Data()));
      		hmult1[icent][ipaty]	->SetLineWidth(2);
      		hmult1[icent][ipaty]	->SetLineColor(1);
			//
      		hmult2[icent][ipaty]	= (TH1D*)hmult[icent][ipaty]->ProjectionY();
      		hmult2[icent][ipaty]	->SetDirectory(0);
      		hmult2[icent][ipaty]	->SetName(Form("hmult2%s",fstro->Data()));
      		hmult2[icent][ipaty]	->SetLineWidth(2);
      		hmult2[icent][ipaty]	->SetLineColor(2);
      		//
      		GetNPairs(hmult[icent][ipaty]);
     		 //cout<<"icent="<<icent<<" ipaty="<<ipaty<<" ie="<<ie<<"\t NPairs="<<NPairs<<endl;
      		double thisNent			= multinfo[0];
      		double thisf11			= multinfo[1];
      		double thisf12			= multinfo[2];
      		double thisf2			= multinfo[3];
      		nevents[icent][ipaty]	= thisNent;		//nevents in ensemble
      		npairs[icent][ipaty]	= thisf2;		//npairs actual = f2
      		C2int[icent][ipaty]		= multinfo[4];	//C2 integral from hmult
      		//
      		//======== Get info from hy_particles_kept 
      		//======== NOTICE hyparticles kept is indexed "the other way" !!!!! (binning from in reader_loop)
			if(SIMINPUT){		//for hy_kept histograms 
				if (!multdone[icent][kPID1] && thisNent){	// use particle1 in this ipaty
	  				hy_particles_kept2[kPID1][icent] = (TH1D*)f->Get( Form("hy_particles_kept2_%d_%d",kPID1,icent) );
	  				mult[icent][kPID1] = hy_particles_kept2[kPID1][icent]->Integral()/thisNent;
					if (mult[icent][kPID1]) multdone[icent][kPID1]	= true;
	  					//<<icent<<" "<<ipaty<<" kPID1="<<kPID1<<" \t "<<hy_particles_kept2[icent][kPID1]->Integral()/thisNent<<endl;
	  					//} else {
	  					//cout<<icent<<" "<<ipaty<<" kPID1="<<kPID1<<" \t already done..."<<endl;
				}
				if (!multdone[icent][kPID2] && thisNent){	//use particle2 in this ipaty
	  				hy_particles_kept2[kPID2][icent] = (TH1D*)f->Get( Form("hy_particles_kept2_%d_%d",kPID2,icent) );
	 				mult[icent][kPID2] = hy_particles_kept2[kPID2][icent]->Integral()/thisNent;
	  				if (mult[icent][kPID2]) multdone[icent][kPID2]	= true;
	  					//cout<<icent<<" "<<ipaty<<" kPID2="<<kPID2<<" \t "<<hy_particles_kept2[icent][kPID2]->Integral()/thisNent<<endl;
					} else {
	  					//cout<<icent<<" "<<ipaty<<" kPID2="<<kPID1<<" \t already done..."<<endl;
				}
    		}	//======= end SIMPUT
    		//
      		//
      		ir2			= 1;		// (dy,dphi) only here! 
      		//
      		hrho1rho1[icent][ipaty]	= (TH2D*)f->Get(Form("hrho1rho1%d%s",ir2,fstri->Data()));
      		hrho1rho1[icent][ipaty]	->SetDirectory(0);
      		hrho1rho1[icent][ipaty]	->SetName(Form("hrho1rho1%s",fstro->Data()));
      		hrho2[icent][ipaty]		= (TH2D*)f->Get(Form("hrho2%d%s",ir2,fstri->Data()));
      		hrho2[icent][ipaty]		->SetDirectory(0);
      		hrho2[icent][ipaty]		->SetName(Form("hrho2%s",fstro->Data()));
      		hR2[icent][ipaty]		= (TH2D*)f->Get(Form("hR2%d%s",ir2,fstri->Data()));
      		hR2[icent][ipaty]		->SetDirectory(0);
      		hR2[icent][ipaty]		->SetName(Form("hR2%s",fstro->Data()));
      		//
      		//--------------------------------------------------------------
     		//---- Scale by 1/binwidths...
     		//double DYNB	= hrho1rho1[icent][ipaty]->GetNbinsX();
     		//double DYL	= hrho1rho1[icent][ipaty]->GetXaxis()->GetXmin();
     		//double DYU	= hrho1rho1[icent][ipaty]->GetXaxis()->GetXmax();
     		//double DYBW	= (DYU-DYL)/DYNB;
     		//double DPHINB	= hrho1rho1[icent][ipaty]->GetNbinsY();
     		//double DPHIL	= hrho1rho1[icent][ipaty]->GetYaxis()->GetXmin();
     		//double DPHIU	= hrho1rho1[icent][ipaty]->GetYaxis()->GetXmax();
     		//double DPHIBW	= (DPHIU-DPHIL)/DPHINB;		// degrees
     		//hrho1rho1[icent][ipaty]	->Scale(1./DYBW/DPHIBW);
     		//hrho2[icent][ipaty]		->Scale(1./DYBW/DPHIBW);
     		//
     		// if (icent==16&&ie==0){
     		// cout<<ipaty
     		// 	<<"\t dy: \t"<<hrho1rho1[icent][ipaty]->GetNbinsX()
     		// 	<<" "<<hrho1rho1[icent][ipaty]->GetXaxis()->GetXmin()
     		// 	<<" "<<hrho1rho1[icent][ipaty]->GetXaxis()->GetXmax()
     		// 	<<"\t dphi: \t"<<hrho1rho1[icent][ipaty]->GetNbinsY()
     		// 	<<" "<<hrho1rho1[icent][ipaty]->GetYaxis()->GetXmin()
     		// 	<<" "<<hrho1rho1[icent][ipaty]->GetYaxis()->GetXmax()
     		// 	<<endl;
     		// }
     		//
     		//------------------------------------------------------------------
     		//
      		if (!SIMINPUT){				// real data, so crossing-corrected hists DO exist...
				hrho2c[icent][ipaty]	= (TH2D*)f->Get(Form("hrho2c%s",fstri->Data()));
				hrho2c[icent][ipaty]	->SetDirectory(0);
				hrho2c[icent][ipaty]	->SetName(Form("hrho2c%s",fstro->Data()));
				hR2c[icent][ipaty]		= (TH2D*)f->Get(Form("hR2c%s",fstri->Data()));
				hR2c[icent][ipaty]		->SetDirectory(0);
				hR2c[icent][ipaty]		->SetName(Form("hR2c%s",fstro->Data()));
				hR2cdy[icent][ipaty]	= (TH1D*)f->Get(Form("hR2cdy%s",fstri->Data()));
				hR2cdy[icent][ipaty]	->SetDirectory(0);
				hR2cdy[icent][ipaty]	->SetName(Form("hR2cdy%s",fstro->Data()));
				hR2cdphi[icent][ipaty]	= (TH1D*)f->Get(Form("hR2cdphi%s",fstri->Data()));
				hR2cdphi[icent][ipaty]	->SetDirectory(0);
				hR2cdphi[icent][ipaty]	->SetName(Form("hR2cdphi%s",fstro->Data()));
      		} else if (SIMINPUT){		// sim data, so crossing-corrected hists DO NOT exist...
				hrho2c[icent][ipaty]	= (TH2D*)f->Get(Form("hrho21%s",fstri->Data()));
				hrho2c[icent][ipaty]	->SetDirectory(0);
				hrho2c[icent][ipaty]	->SetName(Form("hrhoc%s",fstro->Data()));
				hR2c[icent][ipaty]		= (TH2D*)f->Get(Form("hR21%s",fstri->Data()));
				hR2c[icent][ipaty]		->SetDirectory(0);
				hR2c[icent][ipaty]		->SetName(Form("hR2c%s",fstro->Data()));
				hR2cdy[icent][ipaty]	= (TH1D*)f->Get(Form("hR2dy%s",fstri->Data()));
				hR2cdy[icent][ipaty]	->SetDirectory(0);
				hR2cdy[icent][ipaty]	->SetName(Form("hR2cdy%s",fstro->Data()));
				hR2cdphi[icent][ipaty]	= (TH1D*)f->Get(Form("hR2dphi%s",fstri->Data()));
				hR2cdphi[icent][ipaty]	->SetDirectory(0);
				hR2cdphi[icent][ipaty]	->SetName(Form("hR2cdphi%s",fstro->Data()));
      		}
      		//
      		//========= Form rho2 and rho1rho1 integrals...
      		//			NOTE: the units here are COUNTS.
      		//				  thus, you just integrate like a physicist. just add up all bin contents.
      		//			  OR, you can integrate like a mathematician IF you scale each bin content with 1/BWx/BWy first!
      		//
      		double ainte	= 0;
      		//aintrho1rho1[icent][ipaty]= hrho1rho1[icent][ipaty]->IntegralAndError(1,0,1,0,ainte,"width");				
      		aintrho1rho1[icent][ipaty]	= hrho1rho1[icent][ipaty]->IntegralAndError(1,0,1,0,ainte,"");				
      		aintrho1rho1e[icent][ipaty]	= ainte;
      		//aintrho2[icent][ipaty]	= hrho2[icent][ipaty]->IntegralAndError(1,0,1,0,ainte,"width");				
      		aintrho2[icent][ipaty]		= hrho2[icent][ipaty]->IntegralAndError(1,0,1,0,ainte,"");				
      		aintrho2e[icent][ipaty]		= ainte;
      		//
      		hic_rho1rho1_[ipaty]	->Fill(thisf11*thisf12,aintrho1rho1[icent][ipaty]);
      		hic_rho2_[ipaty]		->Fill(thisf2         ,aintrho2[icent][ipaty]    );
      		hic_rho1rho1_rat[ipaty]	->Fill(thisf11*thisf12,aintrho1rho1[icent][ipaty]/thisf11/thisf12);
      		hic_rho2_rat[ipaty]		->Fill(thisf2         ,aintrho2[icent][ipaty]/thisf2);
//       		if (icent==16){
// 				cout<<ipaty<<" \t "
// 	    			<<thisf11*thisf12<<" "
// 	   				<<aintrho1rho1[icent][ipaty]<<" +- "<<aintrho1rho1e[icent][ipaty]<<" \t "
// 	    			<<thisf2<<" "
// 	    			<<aintrho2[icent][ipaty]<<" +- "<<aintrho2e[icent][ipaty]<<" \t "
// 	    			<<aintrho1rho1[icent][ipaty]/thisf11/thisf12<<" \t "				
// 	    			<<aintrho2[icent][ipaty]/thisf2
// 	    			<<endl;
//      		}
    
    
    	} //end cent
  	} //end paty
	//
	
	
	
	
	//==== Set up empty hists for B2s these will have the same binning as the appropriate R2s
	//
	for (int ipart=0;ipart<3;ipart++){		// pi, then K, then p
    	for (int icent=1;icent<NCent;icent++){
      		TString* fstri	= new TString(Form("_%d_%d",icent,ipart));
      		TString* fstro	= new TString(Form("_%d_%d",icent,ipart));		// same! (no ie loop anymore)
      		int ipaty		= 4*ipart;		// clone the LSpos to make the B2
      		//
      		hB2cp[icent][ipart]		= (TH2D*)hR2c[icent][ipaty]->Clone();	// B2+
      		hB2cm[icent][ipart]		= (TH2D*)hR2c[icent][ipaty]->Clone();	// B2-
      		hB2c[icent][ipart]		= (TH2D*)hR2c[icent][ipaty]->Clone();	// B2
      		hB2cp[icent][ipart]		->SetName(Form("hB2cp%s",fstro->Data()));
      		hB2cm[icent][ipart]		->SetName(Form("hB2cm%s",fstro->Data()));
      		hB2c[icent][ipart]		->SetName(Form("hB2c%s",fstro->Data()));
      		hB2cp[icent][ipart]		->SetDirectory(0);
      		hB2cm[icent][ipart]		->SetDirectory(0);
      		hB2c[icent][ipart]		->SetDirectory(0);
      		hB2cp[icent][ipart]		->Reset();
      		hB2cm[icent][ipart]		->Reset();
      		hB2c[icent][ipart]		->Reset();
      		//
      		hB2cpdy[icent][ipart]	= (TH1D*)hR2cdy[icent][ipaty]->Clone();	// B2+ dy
      		hB2cmdy[icent][ipart]	= (TH1D*)hR2cdy[icent][ipaty]->Clone();	// B2- dy
      		hB2cdy[icent][ipart]	= (TH1D*)hR2cdy[icent][ipaty]->Clone();	// B2 dy
      		hB2cpdy[icent][ipart]	->SetName(Form("hB2cpdy%s",fstro->Data()));
      		hB2cmdy[icent][ipart]	->SetName(Form("hB2cmdy%s",fstro->Data()));
      		hB2cdy[icent][ipart]	->SetName(Form("hB2cdy%s" ,fstro->Data()));
      		hB2cpdy[icent][ipart]	->SetDirectory(0);
      		hB2cmdy[icent][ipart]	->SetDirectory(0);
      		hB2cdy[icent][ipart]	->SetDirectory(0);
   			//
      		hB2cpdphi[icent][ipart]	= (TH1D*)hR2cdphi[icent][ipaty]->Clone();	// B2+ dphi
      		hB2cmdphi[icent][ipart]	= (TH1D*)hR2cdphi[icent][ipaty]->Clone();	// B2- dphi
      		hB2cdphi[icent][ipart]	= (TH1D*)hR2cdphi[icent][ipaty]->Clone();	// B2 dphi
      		hB2cpdphi[icent][ipart]	->SetName(Form("hB2cpdphi%s",fstro->Data()));
      		hB2cmdphi[icent][ipart]	->SetName(Form("hB2cmdphi%s",fstro->Data()));
      		hB2cdphi[icent][ipart]	->SetName(Form("hB2cdphi%s" ,fstro->Data()));
      		hB2cpdphi[icent][ipart]	->SetDirectory(0);
      		hB2cmdphi[icent][ipart]	->SetDirectory(0);
      		hB2cdphi[icent][ipart]	->SetDirectory(0);
      		//
      		hB2cpdy[icent][ipart]	->SetLineColor(kBlue);
      		hB2cmdy[icent][ipart]	->SetLineColor(kBlue);
      		hB2cdy[icent][ipart]	->SetLineColor(kBlue);
      		hB2cpdphi[icent][ipart]	->SetLineColor(kBlue);
      		hB2cmdphi[icent][ipart]	->SetLineColor(kBlue);
      		hB2cdphi[icent][ipart]	->SetLineColor(kBlue);
      		hB2cpdy[icent][ipart]	->SetMarkerColor(kBlue);
      		hB2cmdy[icent][ipart]	->SetMarkerColor(kBlue);
      		hB2cdy[icent][ipart]	->SetMarkerColor(kBlue);
      		hB2cpdphi[icent][ipart]	->SetMarkerColor(kBlue);
      		hB2cmdphi[icent][ipart]	->SetMarkerColor(kBlue);
      		hB2cdphi[icent][ipart]	->SetMarkerColor(kBlue);
      		//
      		//==== Set new Titles for all three B2s and all six Projections
      		TString thistit			= (TString)hB2c[icent][ipart]->GetTitle();
      		int thisind				= thistit.Index("R",1,TString::kExact);
      		thistit.Replace(thisind,1,"B");
      		//cout<<icent<<" "<<ipart<<" "<<ie<<" \t"<<thistit.Data()<<endl;
      		hB2c[icent][ipart]	->SetTitle(thistit.Data());
      		thisind					= thistit.Index("2",1,TString::kExact);
      		thistit.Insert(thisind+1,"+",1);
      		hB2cp[icent][ipart]	->SetTitle(thistit.Data());
      		//cout<<icent<<" "<<ipart<<" "<<ie<<" \t"<<thistit.Data()<<endl;
      		thisind					= thistit.Index("2+",1,TString::kExact);
      		thistit.Replace(thisind,2,"2-");
      		hB2cm[icent][ipart]	->SetTitle(thistit.Data());
      		//cout<<icent<<" "<<ipart<<" "<<ie<<" \t"<<thistit.Data()<<endl;
      		//
      		thistit					= (TString)hB2cdy[icent][ipart]->GetTitle();
      		thisind					= thistit.Index("R",1,TString::kExact);
      		thistit.Replace(thisind,1,"B");
      		hB2cdy[icent][ipart]->SetTitle(thistit.Data());
      		thisind					= thistit.Index("2",1,TString::kExact);
      		thistit.Insert(thisind+1,"+",1);
      		hB2cpdy[icent][ipart]->SetTitle(thistit.Data());
      		thisind					= thistit.Index("2+",1,TString::kExact);
      		thistit.Replace(thisind,2,"2-");
      		hB2cmdy[icent][ipart]->SetTitle(thistit.Data());
      		//
      		thistit					= (TString)hB2cdphi[icent][ipart]->GetTitle();
      		thisind					= thistit.Index("R",1,TString::kExact);
      		thistit.Replace(thisind,1,"B");
      		hB2cdphi[icent][ipart]->SetTitle(thistit.Data());
      		thisind					= thistit.Index("2",1,TString::kExact);
      		thistit.Insert(thisind+1,"+",1);
      		hB2cpdphi[icent][ipart]->SetTitle(thistit.Data());
      		thisind					= thistit.Index("2+",1,TString::kExact);
      		thistit.Replace(thisind,2,"2-");
      		hB2cmdphi[icent][ipart]->SetTitle(thistit.Data());
      		//
      		//
      		double cpct	= (16-icent)*5. + 2.5;	// center of cent bin in percent
      		int kp,km;							
      		if (ipart==0){ kp=0; km=7; }		// pi+, pi-
      		if (ipart==1){ kp=1; km=8; }		// K+, K-
      		if (ipart==2){ kp=2; km=9; }		// p, pbar
      		//
      		//__________________________________________________________________
      		//
      		//===== Prefactors come from hy_particles_kept for simulation input!
      		//
      		double prefp,prefpe,prefm,prefme;
      		if (!SIMINPUT){	
      						//***** real data input...
				GetdNdy(cpct,kp);			// mult pos
				prefp 	= xsinfo[0]; 		
				prefpe	= xsinfo[1];		
				GetdNdy(cpct,km);			// mult neg
				prefm	= xsinfo[0]; 		
				prefme	= xsinfo[1];		
		
				
				
// 				GetdNdy(cpct,kp);			// mult pos
// 				prefp 	= xsinfo[0]*0.6; 		
// 				prefpe	= xsinfo[1]*0.6;		
// 				GetdNdy(cpct,km);			// mult neg
// 				prefm	= xsinfo[0]*0.6; 		
// 				prefme	= xsinfo[1]*0.6;		
		

// 				prefp 	= hmult1[icent][ipaty]->GetMean(); 		//
// 				prefm 	= hmult1[icent][ipaty+1]->GetMean(); 		//
//  				prefpe 	= hmult1[icent][ipaty]->GetMeanError(); 		//
//  				prefme 	= hmult1[icent][ipaty+1]->GetMeanError(); 		//
 								
      		} else {					//**** simulation input...
				prefp	=      mult[icent][kp] ;  	// from hy_particles_kept2
				prefpe	= sqrt(mult[icent][kp]);  	// from hy_particles_kept2
				prefm	=      mult[icent][km] ;  	// from hy_particles_kept2
				prefme	= sqrt(mult[icent][km]);  	// from hy_particles_kept2
      		}
      		if (prefp<=0.||prefm<=0.){
				cout<<"PREF NEG \t"<<icent<<" "<<kp<<" "<<km<<" "<<prefp<<" "<<prefm<<endl;
			continue; 
      		}
      		//
      		//
      		//__________________________________________________________________
			//
      		//======= Build B2...
      		BuildB2(prefp,prefpe,				//  INPUT <N+> and uncertainty
			prefm,prefme,				//  INPUT <N-> and uncertainty
			hR2c[icent][0+ipart*4],		//  INPUT R2 LSpos
			hR2c[icent][1+ipart*4],		//  INPUT R2 LSneg
			hR2c[icent][2+ipart*4],		//  INPUT R2 ULSpm
			hR2c[icent][3+ipart*4],		//  INPUT R2 ULSmp
			hB2cp[icent][ipart],		// OUTPUT B2+
			hB2cm[icent][ipart],		// OUTPUT B2-
			hB2c[icent][ipart] 	);		// OUTPUT B2
		//
      		double aNBdy	= (double) hB2c[icent][ipart]->GetNbinsX();
      		double aNBdphi	= (double) hB2c[icent][ipart]->GetNbinsY();
      		double aNB		= aNBdy * aNBdphi;
		//
		//hB2cp[icent][ipart]	->Scale(1./aNB);	// Final B2+
		//hB2cm[icent][ipart]	->Scale(1./aNB);	// Final B2-
		//hB2c[icent][ipart]	->Scale(1./aNB);	// Final B2
		
      
      		//__________________________________________________________________
			//
      		//===== B2 integrals...
      		double ainte	= 0;
      		aintB2cp[icent][ipart]	 = hB2cp[icent][ipart]->IntegralAndError(1,0,1,0,ainte,"")/aNB;				
      		aintB2cpe[icent][ipart]	 = ainte/aNB;
      		aintB2cm[icent][ipart]	 = hB2cm[icent][ipart]->IntegralAndError(1,0,1,0,ainte,"")/aNB;				
      		aintB2cme[icent][ipart]	 = ainte/aNB;
      		aintB2c[icent][ipart]	 = hB2c[icent][ipart] ->IntegralAndError(1,0,1,0,ainte,"")/aNB;				
      		aintB2ce[icent][ipart]	 = ainte/aNB;
      		//cout<<icent<<" "<<ipart<<" "<<ie<<"\t NB="<<aNB<<" \t integral: "
      		//	  <<aintB2c[icent][ipart]<<" +/- "<<aintB2cpe[icent][ipart]<<endl;
      	
      
      		//__________________________________________________________________
			//
        	//===== B2 RMS...
//         	hB2cp_dy[icent][ipart] = (TH1D*)hB2cp[icent][ipart]->ProjectionX();	
// 			hB2cp_dy[icent][ipart]->SetDirectory(0);
// 			hB2cm_dy[icent][ipart] = (TH1D*)hB2cm[icent][ipart]->ProjectionX();	
// 			hB2cm_dy[icent][ipart]->SetDirectory(0);
//         	hB2c_dy[icent][ipart] = (TH1D*)hB2c[icent][ipart]->ProjectionX();	
// 			hB2c_dy[icent][ipart]->SetDirectory(0);
// 
// 	        hB2cp_dphi[icent][ipart] = (TH1D*)hB2cp[icent][ipart]->ProjectionY();	
// 			hB2cp_dphi[icent][ipart]->SetDirectory(0);
// 			hB2cm_dphi[icent][ipart] = (TH1D*)hB2cm[icent][ipart]->ProjectionY();	
// 			hB2cm_dphi[icent][ipart]->SetDirectory(0);
//         	hB2c_dphi[icent][ipart] = (TH1D*)hB2c[icent][ipart]->ProjectionY();	
// 			hB2c_dphi[icent][ipart]->SetDirectory(0);



			//
    		
      		
      		//__________________________________________________________________
      		//
      		//====== SUM RULES for B2 integrals...
      		double intC2pp		= C2int[icent][0+ipart*4];
      		double intC2mm		= C2int[icent][1+ipart*4];
      		double intC2pm		= C2int[icent][2+ipart*4];
      		double intC2mp		= C2int[icent][3+ipart*4];
      		double sumrule_B2p	= (1./prefp)*(intC2pm - intC2pp);
      		double sumrule_B2m	= (1./prefm)*(intC2mp - intC2mm);
      		double sumrule_B2	= (sumrule_B2p + sumrule_B2m)/2.;
      		//cout<<"icent="<<icent<<" ipart="<<ipart<<" \t "
      		//	  <<aintB2cp[icent][ipart]<<" "<<sumrule_B2p<<" \t "
      		//	  <<aintB2cm[icent][ipart]<<" "<<sumrule_B2m<<" \t "
      		//	  <<endl;
      	
      		int kntmp			= gintB2cp_cent[ipart]->GetN();
      		gsruB2cp_cent[ipart]	->SetPoint     (kntmp,(double)icent+0.1,sumrule_B2p);
      		gintB2cp_cent[ipart]	->SetPoint     (kntmp,(double)icent-0.1,aintB2cp[icent][ipart]);
      		gintB2cp_cent[ipart]	->SetPointError(kntmp,(double)   0.    ,aintB2cpe[icent][ipart]);
      		gsruB2cm_cent[ipart]	->SetPoint     (kntmp,(double)icent+0.1,sumrule_B2m);
      		gintB2cm_cent[ipart]	->SetPoint     (kntmp,(double)icent-0.1,aintB2cm[icent][ipart]);
      		gintB2cm_cent[ipart]	->SetPointError(kntmp,(double)   0.    ,aintB2cme[icent][ipart]);
       		gsruB2c_cent[ipart]	->SetPoint     (kntmp,(double)icent+0.1,sumrule_B2);
      		gintB2c_cent[ipart]	->SetPoint     (kntmp,(double)icent-0.1,aintB2c[icent][ipart]);
      		gintB2c_cent[ipart]	->SetPointError(kntmp,(double)   0.    ,aintB2ce[icent][ipart]);
      		//cout<<icent<<" \t"	<<aintB2cp[icent][ipart]<<" "<<aintB2cpe[icent][ipart]<<" "<<sumrule_B2p<<" "
      		//	  <<aintB2cm[icent][ipart]<<" "<<aintB2cme[icent][ipart]<<" "<<sumrule_B2m<<" "
      		//	  <<endl;
		      	
      	  	//__________________________________________________________________
			//
      		//==== Fill B2 projections...
      		GenericProject(hB2cp[icent][ipart],"Avg","X",hB2cpdy[icent][ipart]);
      		GenericProject(hB2cm[icent][ipart],"Avg","X",hB2cmdy[icent][ipart]);
      		GenericProject(hB2c[icent][ipart] ,"Avg","X",hB2cdy[icent][ipart] );
      		GenericProject(hB2cp[icent][ipart],"Avg","Y",hB2cpdphi[icent][ipart]);
      		GenericProject(hB2cm[icent][ipart],"Avg","Y",hB2cmdphi[icent][ipart]);
      		GenericProject(hB2c[icent][ipart] ,"Avg","Y",hB2cdphi[icent][ipart] );

		armsB2cdyp[icent][ipart]        = hB2cpdy[icent][ipart]->GetRMS();
                armsB2cdype[icent][ipart]       = hB2cpdy[icent][ipart]->GetRMSError();
                armsB2cdym[icent][ipart]        = hB2cmdy[icent][ipart]->GetRMS();
                armsB2cdyme[icent][ipart]       = hB2cmdy[icent][ipart]->GetRMSError();
                armsB2cdy[icent][ipart]         = hB2cdy[icent][ipart]->GetRMS();
                armsB2cdye[icent][ipart]        = hB2cdy[icent][ipart]->GetRMSError();
		
		armsB2cdphip[icent][ipart]       = hB2cpdphi[icent][ipart]->GetRMS();
                armsB2cdphipe[icent][ipart]      = hB2cpdphi[icent][ipart]->GetRMSError();
                armsB2cdphim[icent][ipart]      = hB2cmdphi[icent][ipart]->GetRMS();
                armsB2cdphime[icent][ipart]     = hB2cmdphi[icent][ipart]->GetRMSError();
                armsB2cdphi[icent][ipart]        = hB2cdphi[icent][ipart]->GetRMS();
                armsB2cdphie[icent][ipart]       = hB2cdphi[icent][ipart]->GetRMSError();
		
		int krmsmp_dy                   = grmsB2cp_dy_cent[ipart]->GetN();
                grmsB2cp_dy_cent[ipart] ->SetPoint     (krmsmp_dy,(double)icent-0.1,armsB2cdyp[icent][ipart]);
                grmsB2cp_dy_cent[ipart] ->SetPointError(krmsmp_dy,(double)   0.    ,armsB2cdype[icent][ipart]);
                grmsB2cm_dy_cent[ipart] ->SetPoint     (krmsmp_dy,(double)icent-0.1,armsB2cdym[icent][ipart]);
                grmsB2cm_dy_cent[ipart] ->SetPointError(krmsmp_dy,(double)   0.    ,armsB2cdyme[icent][ipart]);
                grmsB2c_dy_cent[ipart]  ->SetPoint     (krmsmp_dy,(double)icent-0.1,armsB2cdy[icent][ipart]);
                grmsB2c_dy_cent[ipart]  ->SetPointError(krmsmp_dy,(double)   0.    ,armsB2cdye[icent][ipart]);

		int krmsmp_dphi                 = grmsB2cp_dphi_cent[ipart]->GetN();
                grmsB2cp_dphi_cent[ipart]       ->SetPoint     (krmsmp_dphi,(double)icent-0.1,armsB2cdphip[icent][ipart]);
                grmsB2cp_dphi_cent[ipart]       ->SetPointError(krmsmp_dphi,(double)   0.    ,armsB2cdphipe[icent][ipart]);
                grmsB2cm_dphi_cent[ipart]       ->SetPoint     (krmsmp_dphi,(double)icent-0.1,armsB2cdphim[icent][ipart]);
                grmsB2cm_dphi_cent[ipart]       ->SetPointError(krmsmp_dphi,(double)   0.    ,armsB2cdphime[icent][ipart]);
                grmsB2c_dphi_cent[ipart]        ->SetPoint     (krmsmp_dphi,(double)icent-0.1,armsB2cdphi[icent][ipart]);
                grmsB2c_dphi_cent[ipart]        ->SetPointError(krmsmp_dphi,(double)   0.    ,armsB2cdphie[icent][ipart]);

		
		//==== Fill B2 projections...
			//GenericProject(hB2cp[icent][ipart],"Add","X",hB2cpdy[icent][ipart]);
			//GenericProject(hB2cm[icent][ipart],"Add","X",hB2cmdy[icent][ipart]);
			//GenericProject(hB2c[icent][ipart] ,"Add","X",hB2cdy[icent][ipart] );
			//GenericProject(hB2cp[icent][ipart],"Add","Y",hB2cpdphi[icent][ipart]);
			//GenericProject(hB2cm[icent][ipart],"Add","Y",hB2cmdphi[icent][ipart]);
			//GenericProject(hB2c[icent][ipart] ,"Add","Y",hB2cdphi[icent][ipart] );
      
      
      		//===== scale by 1/bin-widths...
      		//double DYNB	= hB2c[icent][ipart]->GetNbinsX();
      		//double DYL	= hB2c[icent][ipart]->GetXaxis()->GetXmin();
      		//double DYU	= hB2c[icent][ipart]->GetXaxis()->GetXmax();
      		//double DYBW	= (DYU-DYL)/DYNB;
      		//double DPHINB	= hB2c[icent][ipart]->GetNbinsY();
      		//double DPHIL	= hB2c[icent][ipart]->GetYaxis()->GetXmin();
      		//double DPHIU	= hB2c[icent][ipart]->GetYaxis()->GetXmax();
      		//double DPHIBW	= (DPHIU-DPHIL)/DPHINB;		// degrees
      		////cout<<DYBW<<" "<<DPHIBW<<endl;
      		//hB2cpdy[icent][ipart]->Scale(1./DPHIBW);
      		//hB2cmdy[icent][ipart]->Scale(1./DPHIBW);
      		//hB2cdy[icent][ipart] ->Scale(1./DPHIBW);
      		//hB2cpdphi[icent][ipart]->Scale(1./DYBW);
      		//hB2cmdphi[icent][ipart]->Scale(1./DYBW);
      		//hB2cdphi[icent][ipart] ->Scale(1./DYBW);
    		//
    
    
    
    	}	// end icent
  	}	// end ipart
	//
  	f->Close();

  	//========================== done calculating balance functions! 

  	
  	//_____________________________________________________
  	// prepare to paint...
  	//
  	gROOT->SetStyle("Modern");
  	gStyle->SetOptStat(0);
  	gStyle->SetPadRightMargin(0.14);
  	gStyle->SetPadTopMargin(0.005);
  	gStyle->SetPadBottomMargin(0.08);
  	gStyle->SetPadLeftMargin(0.12);
  	int ican=-1,iline=-1,ivline=-1; 
  	int itxt=-1;
  	TCanvas *ccan[1000];
  	TLatex *txt[1000]; 
  	for (int i=0;i<1000;i++){
    	txt[i]	= new TLatex();
    	txt[i]	->SetTextSize(0.06);
    	txt[i]	->SetTextAlign(12);
    	txt[i]	->SetNDC();
  	}
  	//_____________________________________________________
	//
	//
  	
  	
  	
  	
  	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// rho1rho1 integrals by ipaty...

/*	
	++ican; ccan[ican]	= new TCanvas(Form("ccan%d",ican),Form("ccan%d",ican),20*ican,30+20*ican,900,650);
	ccan[ican]->cd(); ccan[ican]->Divide(5,5,0.0001,0.0001);
	for (int ipaty=0;ipaty<NPairTypes;ipaty++){
  		ccan[ican]->cd(1+ipaty);
  		gPad->SetLogx(1);
  		gPad->SetLogy(1);
  		hic_rho1rho1_[ipaty]->SetMarkerStyle(20);
  		hic_rho1rho1_[ipaty]->SetMarkerSize(0.5);
  		hic_rho1rho1_[ipaty]->SetMarkerColor(4);
  		hic_rho1rho1_[ipaty]->Draw("P");
  		//cout<<hic_rho1rho1_[ipaty]->Integral()<<endl;
  	}
  	ccan[ican]->cd(); ccan[ican]->Update();
  	ccan[ican]->Print(outfileO.Data());	
*/

  	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  	// rho2 integrals by ipaty...

/*
  	++ican; ccan[ican]	= new TCanvas(Form("ccan%d",ican),Form("ccan%d",ican),20*ican,30+20*ican,900,650);
  	ccan[ican]->cd(); ccan[ican]->Divide(5,5,0.0001,0.0001);
  	for (int ipaty=0;ipaty<NPairTypes;ipaty++){
  		ccan[ican]->cd(1+ipaty);
  		gPad->SetLogx(1);
  		gPad->SetLogy(1);
  		hic_rho2_[ipaty]->SetMarkerStyle(20);
  		hic_rho2_[ipaty]->SetMarkerSize(0.5);
  		hic_rho2_[ipaty]->SetMarkerColor(4);
  		hic_rho2_[ipaty]->Draw("P");
  	}
  	ccan[ican]->cd(); ccan[ican]->Update();
  	ccan[ican]->Print(outfile.Data());	
 */ 
 
  	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  	// rho1rho1 and rho2 integrals normalized to factorial moments...
 
 /* 	
  	++ican; ccan[ican]	= new TCanvas(Form("ccan%d",ican),Form("ccan%d",ican),20*ican,30+20*ican,900,650);
  	ccan[ican]->cd(); ccan[ican]->Divide(5,5,0.0001,0.0001);
  	for (int ipaty=0;ipaty<NPairTypes;ipaty++){
  		ccan[ican]->cd(1+ipaty);
  		gPad->SetLogx(1);
  		hic_rho1rho1_rat[ipaty]->SetMarkerStyle(24);
  		hic_rho1rho1_rat[ipaty]->SetMarkerSize(0.5);
  		hic_rho1rho1_rat[ipaty]->SetMarkerColor(4);
  		hic_rho1rho1_rat[ipaty]->Draw("P");
  		hic_rho2_rat[ipaty]->SetMarkerStyle(25);
  		hic_rho2_rat[ipaty]->SetMarkerSize(0.5);
  		hic_rho2_rat[ipaty]->SetMarkerColor(6);
  		hic_rho2_rat[ipaty]->Draw("P same");
  		//cout<<hic_rho1rho1_[ipaty]->Integral()<<endl;
  	}
  	ccan[ican]->cd(); ccan[ican]->Update();
  	ccan[ican]->Print(outfileO.Data());	
*/

  	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  	// R2s and B2s and dy, dphi projections...
  
	cout<<"Painting central collisions only..."<<endl;
	for (int icent=16;icent<=16;icent++){			//!!!!!!!!!!!! PAINT CENTRAL ONLY FOR NOW......
	
    	//-------- Page 1 in series is R2 info...
		for (int ipart=0;ipart<3;ipart++){
    		TString particlename;
			if (ipart==0){ particlename = TString("PIONS");   } else
			if (ipart==1){ particlename = TString("KAONS");   } else
	  		if (ipart==2){ particlename = TString("PROTONS"); }
      		//
      		++ican; ccan[ican]	= new TCanvas(Form("ccan%d",ican),Form("ccan%d",ican),20*ican,30+20*ican,1200,900);
      		ccan[ican]->cd(); ccan[ican]->Divide(4,4,0.0001,0.0001);
			
			//---- Row 1
      		ccan[ican]->cd(1); 
      		gPad->SetLogy(1);
      		hmult1[icent][ipart*4+2]->Draw();
      		hmult2[icent][ipart*4+2]->Draw("same");
      		ccan[ican]->cd(2); 
      		++itxt; txt[itxt]->SetTextSize(0.09); 
      		txt[itxt]->DrawLatex(0.01,0.9,Form("%s",base.Data()));						
      		++itxt; txt[itxt]->SetTextSize(0.09); 
      		txt[itxt]->DrawLatex(0.01,0.8,Form("Particle = %d, %s",ipart,particlename.Data()));
      		++itxt; txt[itxt]->SetTextSize(0.09); 
      		txt[itxt]->DrawLatex(0.01,0.7,Form("Centrality = %s",centnames17[icent]));
		
			//---- Row 2
      		for (int kpa=0;kpa<4;kpa++){
				ccan[ican]->cd(kpa+5); 
				int kpaty	= ipart*4 + kpa;
				hR2c[icent][kpaty]->Draw("lego2");
      		}
      		
			//---- Row 3
      		for (int kpa=0;kpa<4;kpa++){
				ccan[ican]->cd(kpa+9); 
				int kpaty	= ipart*4 + kpa;
				hR2cdy[icent][kpaty]->Draw("E");
 				//hR2cdy[icent][kpaty]->Draw("hist same");
      		}
			
			//---- Row 4
      		for (int kpa=0;kpa<4;kpa++){
				ccan[ican]->cd(kpa+13); 
				int kpaty	= ipart*4 + kpa;
				hR2cdphi[icent][kpaty]->Draw("E");
				//hR2cdphi[icent][kpaty]->Draw("hist same");
      		}
      
    	ccan[ican]->cd(); ccan[ican]->Update();
  		ccan[ican]->Print(outfileO.Data());	
      
      
		//-------- Page 2 in series is B2 DISTRIBUTION info...
    	++ican; ccan[ican]	= new TCanvas(Form("ccan%d",ican),Form("ccan%d",ican),20*ican,30+20*ican,1200,900);
    	ccan[ican]->cd(); ccan[ican]->Divide(3,4,0.0001,0.0001);
    	//
    
    	double n_dy = hB2c[icent][ipart]->GetNbinsX();
    	double n_dphi = hB2c[icent][ipart]->GetNbinsY();

    	ccan[ican]->cd(1); 
    	hB2cp[icent][ipart]->Draw("lego2");
    	++itxt; txt[itxt]->SetTextAlign(32);
    	txt[itxt]->DrawLatex(0.99,0.90,Form("Int: %.2f#pm%.2f",aintB2cp[icent][ipart],aintB2cpe[icent][ipart]));
     	//txt[itxt]->DrawLatex(0.99,0.90,Form("Int: %.2f#pm%.2f",aintB2cp[icent][ipart]/aNB[icent][ipart],aintB2cpe[icent][ipart]/aNB[icent][ipart]));
    	//
    	ccan[ican]->cd(2); 
    	hB2cm[icent][ipart]->Draw("lego2");
    	++itxt; txt[itxt]->SetTextAlign(32);
    	txt[itxt]->DrawLatex(0.99,0.90,Form("Int: %.2f#pm%.2f",aintB2cm[icent][ipart],aintB2cme[icent][ipart]));
    	//
    	ccan[ican]->cd(3); 
    	hB2c[icent][ipart]->Draw("lego2");
    	++itxt; txt[itxt]->SetTextAlign(32);
    	txt[itxt]->DrawLatex(0.99,0.90,Form("Int: %.2f#pm%.2f",aintB2c[icent][ipart],aintB2ce[icent][ipart]));
    	//
    	ccan[ican]->cd(4); 
    	hB2cp[icent][ipart]->Draw("colz");
    	//
    	ccan[ican]->cd(5); 
    	hB2cm[icent][ipart]->Draw("colz");
    	//
    	ccan[ican]->cd(6); 
    	hB2c[icent][ipart]->Draw("colz");
    	//
    	ccan[ican]->cd(7); 
    	RangeFinderH(0,hB2cpdy[icent][ipart]);
    	RangeCheck();
    	hB2cpdy[icent][ipart]->SetMinimum(ylimits[0]);
    	hB2cpdy[icent][ipart]->SetMaximum(ylimits[1]);
    	hB2cpdy[icent][ipart]->Draw("E");
    	//hB2cpdy[icent][ipart]->Draw("hist same");
    	++itxt; txt[itxt]->DrawLatex(0.7,0.88,Form("INT: %.2f",hB2cpdy[icent][ipart]->Integral()/n_dy));
    	//++itxt; txt[itxt]->DrawLatex(0.7,0.82,Form("Int: %.2f",hB2cpdy[icent][ipart]->Integral("width")));
    	//
    	ccan[ican]->cd(8); 
      	RangeFinderH(0,hB2cmdy[icent][ipart]);
      	RangeCheck();
      	hB2cmdy[icent][ipart]->SetMinimum(ylimits[0]);
      	hB2cmdy[icent][ipart]->SetMaximum(ylimits[1]);
      	hB2cmdy[icent][ipart]->Draw("E");
      	//hB2cmdy[icent][ipart]->Draw("hist same");
      	++itxt; txt[itxt]->DrawLatex(0.7,0.88,Form("INT: %.2f",hB2cmdy[icent][ipart]->Integral()/n_dy));
      	//++itxt; txt[itxt]->DrawLatex(0.7,0.82,Form("Int: %.2f",hB2cmdy[icent][ipart]->Integral("width")));
      	//
      	ccan[ican]->cd(9); 
      	RangeFinderH(0,hB2cdy[icent][ipart]);
      	RangeCheck();
      	hB2cdy[icent][ipart]->SetMinimum(ylimits[0]);
      	hB2cdy[icent][ipart]->SetMaximum(ylimits[1]);
      	hB2cdy[icent][ipart]->Draw("E");
      	//hB2cdy[icent][ipart]->Draw("hist same");
      	++itxt; txt[itxt]->DrawLatex(0.7,0.88,Form("INT: %.2f",hB2cdy[icent][ipart]->Integral()/n_dy));
      	//++itxt; txt[itxt]->DrawLatex(0.7,0.82,Form("Int: %.2f",hB2cdy[icent][ipart]->Integral("width")));
      	//
      	ccan[ican]->cd(10); 
      	RangeFinderH(0,hB2cpdphi[icent][ipart]);
      	RangeCheck();
      	hB2cpdphi[icent][ipart]->SetMinimum(ylimits[0]);
      	hB2cpdphi[icent][ipart]->SetMaximum(ylimits[1]);
      	hB2cpdphi[icent][ipart]->Draw("E");
      	//hB2cpdphi[icent][ipart]->Draw("hist same");
      	++itxt; txt[itxt]->DrawLatex(0.7,0.88,Form("INT: %.2f",hB2cpdphi[icent][ipart]->Integral()/n_dphi));
      	//++itxt; txt[itxt]->DrawLatex(0.7,0.82,Form("Int: %.2f",hB2cpdphi[icent][ipart]->Integral("width")));
      	//
      	ccan[ican]->cd(11); 
      	RangeFinderH(0,hB2cmdphi[icent][ipart]);
      	RangeCheck();
      	hB2cmdphi[icent][ipart]->SetMinimum(ylimits[0]);
      	hB2cmdphi[icent][ipart]->SetMaximum(ylimits[1]);
      	hB2cmdphi[icent][ipart]->Draw("E");
      	//hB2cmdphi[icent][ipart]->Draw("hist same");
      	++itxt; txt[itxt]->DrawLatex(0.7,0.88,Form("INT: %.2f",hB2cmdphi[icent][ipart]->Integral()/n_dphi));
      	//++itxt; txt[itxt]->DrawLatex(0.7,0.82,Form("Int: %.2f",hB2cmdphi[icent][ipart]->Integral("width")));
      	//
      	ccan[ican]->cd(12); 
      	RangeFinderH(0,hB2cdphi[icent][ipart]);
      	RangeCheck();
      	hB2cdphi[icent][ipart]->SetMinimum(ylimits[0]);
      	hB2cdphi[icent][ipart]->SetMaximum(ylimits[1]);
      	hB2cdphi[icent][ipart]->Draw("E");
   		//hB2cdphi[icent][ipart]->Draw("hist same");
      	++itxt; txt[itxt]->DrawLatex(0.7,0.88,Form("INT: %.2f",hB2cdphi[icent][ipart]->Integral()/n_dphi));
      	//++itxt; txt[itxt]->DrawLatex(0.7,0.82,Form("Int: %.2f",hB2cdphi[icent][ipart]->Integral("width")));
      
      	ccan[ican]->cd(); ccan[ican]->Update();
     	ccan[ican]->Print(outfile.Data());	
      
      
		//-------- Page 3 in series is B2 DISTRIBUTION info...
      
      	TLine* lhor	= new TLine(0,0,1,1); lhor->SetLineWidth(2); lhor->SetLineColor(16); lhor->SetLineStyle(2);

       
        //-------- B2+, B2-  and B2 integrals, versus centrality, for each particle pair(one page per energy)
  		//
  		++ican; ccan[ican]	= new TCanvas(Form("ccan%d",ican),Form("ccan%d",ican),20*ican,30+20*ican,1200,600);
  		ccan[ican]->cd(); ccan[ican]->Divide(3,2,0.0001,0.0001);
  		//for (int ipart=0;ipart<3;ipart++){
    		ccan[ican]->cd(1);
    		RangeFinderG(0,gintB2cp_cent[ipart]);
    		RangeCheck();
    		gintB2cp_cent[ipart]->Draw("AP");
    		if (ipart==0){ ylimits[0] = -0.1; }
    		if (ipart==1){ ylimits[0] = -0.5; }
    		if (ipart==2){ ylimits[0] = -0.5; }
    	
       		if (ipart==0){ ylimits[1] = 1.7; }
    		if (ipart==1){ ylimits[1] = 1.9; }
    		if (ipart==2){ ylimits[1] = 1.9; }
 	
    		gintB2cp_cent[ipart]->SetMinimum(ylimits[0]);
    		gintB2cp_cent[ipart]->SetMaximum(ylimits[1]);
    		lhor->DrawLine(0,0.0,16.5,0.0); lhor->DrawLine(0,1.0,16.5,1.0); 
     		//if (ipart==0){
       			//++itxt; txt[itxt]->SetTextSize(0.04); txt[itxt]->SetTextAlign(22); txt[itxt]->SetTextFont(42); 
       			//txt[itxt]->DrawLatex(0.5,0.91,Form("%s",base.Data()));						
     		//}
    		//
    		ccan[ican]->cd(2);
    		RangeFinderG(0,gintB2cm_cent[ipart]);
    		RangeCheck();
    		gintB2cm_cent[ipart]->Draw("AP");
    		if (ipart==0){ ylimits[0] = -0.1; }
    		if (ipart==1){ ylimits[0] = -0.5; }
    		if (ipart==2){ ylimits[0] = -1.; }
    	
    		if (ipart==0){ ylimits[1] = 1.7; }
    		if (ipart==1){ ylimits[1] = 1.9; }
    		if (ipart==2){ ylimits[1] = 1.9; }

    		gintB2cm_cent[ipart]->SetMinimum(ylimits[0]);
    		gintB2cm_cent[ipart]->SetMaximum(ylimits[1]);
    		lhor->DrawLine(0,0.0,16.5,0.0); lhor->DrawLine(0,1.0,16.5,1.0); 
    		//
     		ccan[ican]->cd(3);
    		RangeFinderG(0,gintB2c_cent[ipart]);
    		RangeCheck();
    		gintB2c_cent[ipart]->Draw("AP");
    		if (ipart==0){ ylimits[0] = -0.1; }
    		if (ipart==1){ ylimits[0] = -0.5; }
    		if (ipart==2){ ylimits[0] = -1.; }
    	
    		if (ipart==0){ ylimits[1] = 1.7; }
    		if (ipart==1){ ylimits[1] = 1.9; }
    		if (ipart==2){ ylimits[1] = 1.9; }

    		gintB2c_cent[ipart]->SetMinimum(ylimits[0]);
    		gintB2c_cent[ipart]->SetMaximum(ylimits[1]);
    		lhor->DrawLine(0,0.0,16.5,0.0); lhor->DrawLine(0,1.0,16.5,1.0); 
			//
  		
  		ccan[ican]->cd(); ccan[ican]->Update();
  		ccan[ican]->Print(outfile.Data());	

    	//-------- Page 4 in series is B2 DISTRIBUTION info...
      
        //---- B2+, B2-  and B2 RMS values, versus centrality, for each particle pair(one page per energy)
  		//
  		++ican; ccan[ican]	= new TCanvas(Form("ccan%d",ican),Form("ccan%d",ican),20*ican,30+20*ican,1200,700);
  		ccan[ican]->cd(); ccan[ican]->Divide(3,2,0.0001,0.0001);

    		ccan[ican]->cd(1);
    		RangeFinderG(0,grmsB2cp_dy_cent[ipart]);
    		RangeCheck();
    		grmsB2cp_dy_cent[ipart]->Draw("AP");
    		if (ipart==0){ ylimits[0] = -0.2; }
    		if (ipart==1){ ylimits[0] = 0.; }
    		if (ipart==2){ ylimits[0] = 0.; }
    	
       		if (ipart==0){ ylimits[1] = 1.2; }
    		if (ipart==1){ ylimits[1] = 1.; }
    		if (ipart==2){ ylimits[1] = 1.5; }
 	
    		grmsB2cp_dy_cent[ipart]->SetMinimum(ylimits[0]);
    		grmsB2cp_dy_cent[ipart]->SetMaximum(ylimits[1]);
			//
    		ccan[ican]->cd(2);
    		RangeFinderG(0,grmsB2cm_dy_cent[ipart]);
    		RangeCheck();
    		grmsB2cm_dy_cent[ipart]->Draw("AP");
    		if (ipart==0){ ylimits[0] = -0.2; }
    		if (ipart==1){ ylimits[0] = 0.; }
    		if (ipart==2){ ylimits[0] = 0.; }
    	
    		if (ipart==0){ ylimits[1] = 1.2; }
    		if (ipart==1){ ylimits[1] = 1.; }
    		if (ipart==2){ ylimits[1] = 1.5; }
			//
    		grmsB2cm_dy_cent[ipart]->SetMinimum(ylimits[0]);
    		grmsB2cm_dy_cent[ipart]->SetMaximum(ylimits[1]);

     		ccan[ican]->cd(3);
    		RangeFinderG(0,grmsB2c_dy_cent[ipart]);
    		RangeCheck();
    		grmsB2c_dy_cent[ipart]->Draw("AP");
    		if (ipart==0){ ylimits[0] = -0.2; }
    		if (ipart==1){ ylimits[0] = 0.; }
    		if (ipart==2){ ylimits[0] = 0.; }
    	
    		if (ipart==0){ ylimits[1] = 1.2; }
    		if (ipart==1){ ylimits[1] = 1.; }
    		if (ipart==2){ ylimits[1] = 1.5; }

    		grmsB2c_dy_cent[ipart]->SetMinimum(ylimits[0]);
    		grmsB2c_dy_cent[ipart]->SetMaximum(ylimits[1]);

    		ccan[ican]->cd(4);
    		RangeFinderG(0,grmsB2cp_dphi_cent[ipart]);
    		RangeCheck();
    		grmsB2cp_dphi_cent[ipart]->Draw("AP");
//     	if (ipart==0){ ylimits[0] = -0.2; }
//     	if (ipart==1){ ylimits[0] = 0.; }
//     	if (ipart==2){ ylimits[0] = 0.; }
//     	
//        	if (ipart==0){ ylimits[1] = 1.2; }
//     	if (ipart==1){ ylimits[1] = 1.; }
//     	if (ipart==2){ ylimits[1] = 1.5; }
//  	
//     	grmsB2cp_dphi_cent[ipart]->SetMinimum(ylimits[0]);
//     	grmsB2cp_dphi_cent[ipart]->SetMaximum(ylimits[1]);

    		ccan[ican]->cd(5);
    		RangeFinderG(0,grmsB2cm_dphi_cent[ipart]);
    		RangeCheck();
    		grmsB2cm_dphi_cent[ipart]->Draw("AP");
//     		if (ipart==0){ ylimits[0] = -0.2; }
//     		if (ipart==1){ ylimits[0] = 0.; }
//     		if (ipart==2){ ylimits[0] = 0.; }
//     	
//     		if (ipart==0){ ylimits[1] = 1.2; }
//     		if (ipart==1){ ylimits[1] = 1.; }
//     		if (ipart==2){ ylimits[1] = 1.5; }
// 
//     		grmsB2cm_dphi_cent[ipart]->SetMinimum(ylimits[0]);
//     		grmsB2cm_dphi_cent[ipart]->SetMaximum(ylimits[1]);

     		ccan[ican]->cd(6);
    		RangeFinderG(0,grmsB2c_dphi_cent[ipart]);
    		RangeCheck();
    		grmsB2c_dphi_cent[ipart]->Draw("AP");
//     		if (ipart==0){ ylimits[0] = -0.2; }
//     		if (ipart==1){ ylimits[0] = 0.; }
//     		if (ipart==2){ ylimits[0] = 0.; }
//     	
//     		if (ipart==0){ ylimits[1] = 1.2; }
//     		if (ipart==1){ ylimits[1] = 1.; }
//     		if (ipart==2){ ylimits[1] = 1.5; }
// 
//     		grmsB2c_dphi_cent[ipart]->SetMinimum(ylimits[0]);
//     		grmsB2c_dphi_cent[ipart]->SetMaximum(ylimits[1]);


  		ccan[ican]->cd(); ccan[ican]->Update();
  		ccan[ican]->Print(outfile.Data());	

      

    	}	// end ipart pi, k, p
 
  	}	// end icent


  	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*	
TLine* lhor	= new TLine(0,0,1,1); lhor->SetLineWidth(2); lhor->SetLineColor(16); lhor->SetLineStyle(2);

  	//===== compare B2+, B2- integrals to sum rule values, versus centrality, for each particle (one page per energy)
  	++ican; ccan[ican]	= new TCanvas(Form("ccan%d",ican),Form("ccan%d",ican),20*ican,30+20*ican,900,600);
  	ccan[ican]->cd(); ccan[ican]->Divide(3,2,0.0001,0.0001);
  	for (int ipart=0;ipart<3;ipart++){
    	ccan[ican]->cd(1+ipart);
    	RangeFinderG(0,gintB2cp_cent[ipart]);
    	RangeFinderG(1,gsruB2cp_cent[ipart]);
    	RangeCheck();
    	gintB2cp_cent[ipart]->Draw("AP");
    	if (ipart==0){ ylimits[0] = -0.1; }
    	if (ipart==1){ ylimits[0] = -0.5; }
    	if (ipart==2){ ylimits[0] = -1.; }
    	gintB2cp_cent[ipart]->SetMinimum(ylimits[0]);
    	gintB2cp_cent[ipart]->SetMaximum(ylimits[1]);
    	gsruB2cp_cent[ipart]->Draw("PL");
    	gintB2cp_cent[ipart]->Draw("P");
    	lhor->DrawLine(0,0.0,16.5,0.0); lhor->DrawLine(0,1.0,16.5,1.0); 
    	if (ipart==0){
      		++itxt; txt[itxt]->SetTextSize(0.04); txt[itxt]->SetTextAlign(22); txt[itxt]->SetTextFont(42); 
      		txt[itxt]->DrawLatex(0.5,0.91,Form("%s",base.Data()));						
    	}
  	}
  	for (int ipart=0;ipart<3;ipart++){
    	ccan[ican]->cd(4+ipart);
    	RangeFinderG(0,gintB2cm_cent[ipart]);
    	RangeFinderG(1,gsruB2cm_cent[ipart]);
    	RangeCheck();
    	gintB2cm_cent[ipart]->Draw("AP");
    	if (ipart==0){ ylimits[0] = -0.1; }
    	if (ipart==1){ ylimits[0] = -0.5; }
    	if (ipart==2){ ylimits[0] = -1.; }
    	gintB2cm_cent[ipart]->SetMinimum(ylimits[0]);
    	gintB2cm_cent[ipart]->SetMaximum(ylimits[1]);
    	gsruB2cm_cent[ipart]->Draw("PL");
    	gintB2cm_cent[ipart]->Draw("P");
    	lhor->DrawLine(0,0.0,16.5,0.0); lhor->DrawLine(0,1.0,16.5,1.0); 
  	}
 
 
  	ccan[ican]->cd(); ccan[ican]->Update();
  	ccan[ican]->Print(outfile.Data());	
*/
  	
  	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



	//========== close-out...
  	ccan[ican]->Print(outfileC.Data());	
  	//
  	//
  	//==== build the root file...
  	cout<<"Writing "<<outfileroot.Data()<<endl;
  	TFile *fout	= new TFile(outfileroot.Data(),"recreate");
  	fout->cd();
  	//
  	for (int ip=0;ip<NPARTICLENAMES;ip++){			// particles (fluct_common)
    	for (int icent=1;icent<NCent;icent++){
    		if (hy_particles_kept2[ip][icent]){ 
				hy_particles_kept2[ip][icent]->Write();		// indexing originates in reader_loop
      		}
    	}	// end icent
  	}	// end ip
  	//
  	for (int ipaty=0;ipaty<12;ipaty++){			// pi,K,p only (leave out all nuclei pair info)
    	for (int icent=1;icent<NCent;icent++){
      		hmult[icent][ipaty]		->Write();
      		hmult1[icent][ipaty]	->Write();
      		hmult2[icent][ipaty]	->Write();
    	}	// end icent
  	}	// end ipaty
  	//
  	for (int ipart=0;ipart<3;ipart++){
    	for (int icent=1;icent<NCent;icent++){
      		hB2cp[icent][ipart]		->Write();
      		hB2cm[icent][ipart]		->Write();
      		hB2c[icent][ipart]		->Write();
      		hB2cpdy[icent][ipart]	->Write();
      		hB2cmdy[icent][ipart]	->Write();
      		hB2cdy[icent][ipart]	->Write();
      		hB2cpdphi[icent][ipart]	->Write();
      		hB2cmdphi[icent][ipart]	->Write();
      		hB2cdphi[icent][ipart]	->Write();
    	}	// end icent
  	}	// end ipart
  
  	fout->Close();	
	//
	//
  	//======= build and scp the pdf file...
	if (ican>=0){
    	TString sexec;
    	sexec	 = TString(Form("ps2pdf %s %s",outfile.Data(),outfileP.Data()));
    	cout<<sexec.Data()<<endl;
    	gSystem->Exec(sexec.Data());
    	sexec	 = TString(Form("/bin/rm %s",outfile.Data()));
    	cout<<sexec.Data()<<endl;
    	gSystem->Exec(sexec.Data());
    	//sexec	 = TString(Form("/usr/bin/scp %s wjllope@rhic22.physics.wayne.edu:/Library/WebServer/WebPages/files/",outfileP.Data()));
    	//sexec	 = TString(Form("cp %s .",outfileP.Data()));
    	cout<<sexec.Data()<<endl;
    	gSystem->Exec(sexec.Data());
  	}








}


//**************** end of macro







//========================================================================================
void BuildB2(double prefp, double prefpe,
	     double prefm, double prefme,
	     TH2D* R2pp, TH2D* R2mm, TH2D* R2pm, TH2D* R2mp,
	     TH2D* B2p,  TH2D* B2m,  TH2D* B2  ){
  
  // 	double DYNB		= R2pp->GetNbinsX();
  // 	double DYL		= R2pp->GetXaxis()->GetXmin();
  // 	double DYU		= R2pp->GetXaxis()->GetXmax();
  // 	double DYBW		= (DYU-DYL)/DYNB;
  // 	double DPHINB	= R2pp->GetNbinsY();
  // 	double DPHIL	= R2pp->GetYaxis()->GetXmin();
  // 	double DPHIU	= R2pp->GetYaxis()->GetXmax();
  // 	double DPHIBW	= (DPHIU-DPHIL)/DPHINB;		// degrees
  //
  //cout<<"PREF pos: "<<prefp<<" +- "<<prefpe<<" \t neg: "<<prefm<<" +- "<<prefme<<endl;


	for (int ibx=1;ibx<=R2pp->GetNbinsX();ibx++){
    	for (int iby=1;iby<=R2pp->GetNbinsY();iby++){
      		//
      		double valpp	 = R2pp->GetBinContent(ibx,iby);
      		double valppe	 = R2pp->GetBinError  (ibx,iby);
      		double valmm	 = R2mm->GetBinContent(ibx,iby);
      		double valmme	 = R2mm->GetBinError  (ibx,iby);
      		double valpm	 = R2pm->GetBinContent(ibx,iby);
      		double valpme	 = R2pm->GetBinError  (ibx,iby);
      		double valmp	 = R2mp->GetBinContent(ibx,iby);
      		double valmpe	 = R2mp->GetBinError  (ibx,iby);
      		//
      		double b2p		= prefm*valpm - prefp*valpp;
      		double b2pe2	= pow(prefme*valpm,2) + pow(prefm*valpme,2)
    		+ pow(prefpe*valpp,2) + pow(prefp*valppe,2);
      		double b2pe		= sqrt(b2pe2);
      		//
      		double b2m		= prefp*valmp - prefm*valmm;
      		double b2me2	= pow(prefpe*valmp,2) + pow(prefp*valmpe,2)
			+ pow(prefme*valmm,2) + pow(prefm*valmme,2);
     	 	double b2me		= sqrt(b2me2);
      		//
      		double b2e2		 = pow(b2pe,2) + pow(b2me,2);
      		double b2		 = (b2p + b2m)/2.;
      		double b2e		 = sqrt(b2e2)/2.;
      		//
      		//cout<<ibx<<" "<<iby<<" \t "<<valppe<<" "<<valmme<<" "<<valpme<<" "<<valmpe<<" \t "<<b2pe<<endl;
      		//cout<<b2p<<" "<<b2pe<<" "<<b2m<<" "<<b2me<<" "<<b2<<" "<<b2e<<" "<<endl;
      		//cout<<ibx<<" "<<iby<<" \t "<<b2pe<<" "<<b2me<<" \t "<<b2e<<" "<<endl;
      		//
      		B2p	->SetBinContent(ibx,iby,b2p );
     		B2m	->SetBinContent(ibx,iby,b2m );
      		B2	->SetBinContent(ibx,iby,b2  );
      		B2p	->SetBinError  (ibx,iby,b2pe);
      		B2m	->SetBinError  (ibx,iby,b2me);
      		B2	->SetBinError  (ibx,iby,b2e );
      		//
    	}
  	}
  	
	 		//R2pp->Scale(1./DYBW/DPHIBW);
 	 		//R2mm->Scale(1./DYBW/DPHIBW);
  	 		//R2pm->Scale(1./DYBW/DPHIBW);
 	 		//R2mp->Scale(1./DYBW/DPHIBW);
	 		//B2p->Scale(1./DYBW/DPHIBW);
  	 		//B2m->Scale(1./DYBW/DPHIBW);
  	 		//B2->Scale(1./DYBW/DPHIBW);

}

//========================================================================================
void GetdNdy(double cpct,int kpid){
	double xs  = 0.;
  	double xse = 0.;
  	//
  	//======= 200 GeV
  	const int xs_n			= 9;
  	double xs_cpct[xs_n]	= {75,65,55,45,35,25,15,7.5,2.5};
  	double xs_pim[xs_n]		= {10.9,21.1,36.3,58.9,89.6,136,196,261,327};
 	double xs_pime[xs_n]	= { 0.8, 1.6, 2.8, 4.5, 6.8, 10, 15, 20, 25};
  	double xs_pip[xs_n]		= {10.8,21.1,36.2,58.7,89.2,135,194,257,322};
  	double xs_pipe[xs_n]	= { 0.8, 1.6, 2.7, 4.5, 6.8, 10, 15, 20, 25};
  	double xs_km[xs_n]		= {1.38,2.89,5.19,8.37,13.2,19.7,28.7,39.8,49.5};
  	double xs_kme[xs_n]		= {0.13,0.26,0.47,0.78, 1.3, 2.0, 3.1, 4.6, 6.2};
  	double xs_kp[xs_n]		= {1.41,2.98,5.40,8.69,13.6,20.5,30.0,40.8,51.3};
  	double xs_kpe[xs_n]		= {0.13,0.27,0.49,0.81, 1.3, 2.0, 3.2, 4.7, 6.5};
  	double xs_pm[xs_n]		= {0.915,1.84,3.16,4.93,7.46,11.2,15.7,21.4,26.7};
  	double xs_pme[xs_n]		= {0.081,0.16,0.29,0.46,0.72, 1.1, 1.7, 2.5, 3.4};
  	double xs_pp[xs_n]		= {1.09,2.20,3.88,6.17,9.30,14.4,20.1,28.2,34.7};
  	double xs_ppe[xs_n]		= {0.10,0.20,0.35,0.57,0.89, 1.4, 2.2, 3.3, 4.4};
  	//
  	if (cpct<0.){	//---- call with cpercent<0 to initialize!
    	gxs_pim		= new TGraph(xs_n,xs_cpct,xs_pim );
    	gxs_pip		= new TGraph(xs_n,xs_cpct,xs_pip );
    	gxs_km		= new TGraph(xs_n,xs_cpct,xs_km  );
    	gxs_kp		= new TGraph(xs_n,xs_cpct,xs_kp  );
    	gxs_pm		= new TGraph(xs_n,xs_cpct,xs_pm  );
    	gxs_pp		= new TGraph(xs_n,xs_cpct,xs_pp  );
    	gxs_pime	= new TGraph(xs_n,xs_cpct,xs_pime);
    	gxs_pipe	= new TGraph(xs_n,xs_cpct,xs_pipe);
    	gxs_kme		= new TGraph(xs_n,xs_cpct,xs_kme );
    	gxs_kpe		= new TGraph(xs_n,xs_cpct,xs_kpe );
    	gxs_pme		= new TGraph(xs_n,xs_cpct,xs_pme );
    	gxs_ppe		= new TGraph(xs_n,xs_cpct,xs_ppe );
  	}
  	//
  	//======= if no xsec for this particle, return zero...
	if (kpid<0||(kpid%7)>3){ return xs;     } else
	if (kpid==0){ xs = gxs_pip->Eval(cpct); xse = gxs_pipe->Eval(cpct); } else	// pi+
    if (kpid==1){ xs = gxs_kp ->Eval(cpct); xse = gxs_pime->Eval(cpct); } else	// K+
	if (kpid==2){ xs = gxs_pp ->Eval(cpct); xse = gxs_kpe ->Eval(cpct); } else	// p
	if (kpid==7){ xs = gxs_pim->Eval(cpct); xse = gxs_kme ->Eval(cpct); } else	// pi-
	if (kpid==8){ xs = gxs_km ->Eval(cpct); xse = gxs_ppe ->Eval(cpct); } else	// K-
	if (kpid==9){ xs = gxs_pm ->Eval(cpct); xse = gxs_pme ->Eval(cpct); } 		// pbar
  	xsinfo[0] = xs;
  	xsinfo[1] = xse;
  	return;
  	//
	
};


//========================================================================================
//
//	Service routine to project any 2D R2 onto a 1D R2 projection
//	Called from Loop right after SetFLip (crossing correction)...
//	This routine assumes that the 1D is properly binned for the assumed projection! 
//
//	examples:
//			GenericProject(hR2dydphi,"Avg","Y",h1)		--> h1 is R2(dphi) from (dy,dphi), crossing correctable
//			GenericProject(hR2dydphi,"Avg","X",h1)		--> h1 is R2(dy) from (dy,dphi), crossing correctable
//			GenericProject(hR2y1y2,"Avg","D",h1)		--> h1 is R2(dy) from (y1,y2), NOT crossing correctable
//			GenericProject(hR2y1y2,"Avg","A",h1)		--> h1 is R2(<y>) from (y1,y2), NOT crossing correctable
//
//	now includes protection on h1 (should be empty at call!) !!
//
//
void GenericProject(TH2D* h2,TString steerav,TString steerax,TH1D* h1){
	//
  	if (!h1){ cout<<"reader::GenericProject -- no h1 pointer, cannot project!"<<endl; exit(0); }
  	if (!h2){ cout<<"reader::GenericProject -- no h2 pointer, cannot project!"<<endl; exit(0); }
  	bool doadd=false, doavg=false;
  	doadd	= steerav.Contains("Add",TString::kIgnoreCase);	
  	doavg	= steerav.Contains("Avg",TString::kIgnoreCase);	
  	if (!doadd && !doavg){ cout<<"reader::GenericProject steering issue."<<endl; exit(0); }
  	if ( doadd &&  doavg){ cout<<"reader::GenericProject steering issue."<<endl; exit(0); }
  	int kProjAxis	= -1;
  	if (steerax.Contains("X",TString::kIgnoreCase)) kProjAxis = 0;	// project onto x axis
  	if (steerax.Contains("Y",TString::kIgnoreCase)) kProjAxis = 1;	// project onto y axis
  	if (steerax.Contains("A",TString::kIgnoreCase)) kProjAxis = 2;	// project onto axis-avg axis
  	if (steerax.Contains("D",TString::kIgnoreCase)) kProjAxis = 3;	// project onto axis-diff axis
  	int kR2			=  0;
  	TString histtit	= (TString)h2->GetTitle();
  	if (histtit.Contains("y_{1},y_{2}",TString::kIgnoreCase)) kR2 = 0;	// h2 is (y1,y2);
  	if (histtit.Contains("dy,d#phi"   ,TString::kIgnoreCase)) kR2 = 1;	// h2 is (dy,dphi);
  	if (histtit.Contains("dy,dq"      ,TString::kIgnoreCase)) kR2 = 2;	// h2 is (dy,dq);
  	//kR2	= 1;	// hardwired in this B2 macro
  	//
  	h1			->Reset();
  	h1->SetMaximum(1.0);
  	h1->SetMinimum(-1.0);
  	TH1D *h1e	= (TH1D*)h1->Clone("h1e");
  	TH1D *h1N	= (TH1D*)h1->Clone("h1N");
  	//
  	if (kProjAxis==0){		// project onto X-axis of h2
    	for (int ibx=1;ibx<=h2->GetXaxis()->GetNbins();ibx++){
      		for (int iby=1;iby<=h2->GetYaxis()->GetNbins();iby++){
				double x,y,a,d,val,vale;
					x	= h2->GetXaxis()->GetBinCenter(ibx);
					val		= h2->GetBinContent(ibx,iby);
					vale	= h2->GetBinError(ibx,iby);
					h1		->Fill(x, val);
					h1e		->Fill(x, vale*vale);
					h1N		->Fill(x, 1.0);
      			}	// iby
    	}// ibx
    
  	} else if (kProjAxis==1){		// project onto Y-axis of h2
	for (int ibx=1;ibx<=h2->GetXaxis()->GetNbins();ibx++){
    	for (int iby=1;iby<=h2->GetYaxis()->GetNbins();iby++){
			double x,y,a,d,val,vale;
				y		= h2->GetYaxis()->GetBinCenter(iby);
				val		= h2->GetBinContent(ibx,iby);
				vale	= h2->GetBinError(ibx,iby);
				h1		->Fill(y, val);
				h1e		->Fill(y, vale*vale);
				h1N		->Fill(y, 1.0);
      		}	// iby
    	}// ibx
  	} else if (kProjAxis==2){		// project onto avg-axis of h2
    for (int ibx=1;ibx<=h2->GetXaxis()->GetNbins();ibx++){
    	for (int iby=1;iby<=h2->GetYaxis()->GetNbins();iby++){
			double x,y,a,d,val,vale;
				a		= (h2->GetXaxis()->GetBinCenter(ibx) + h2->GetYaxis()->GetBinCenter(iby))/2.;
				val		= h2->GetBinContent(ibx,iby);
				vale	= h2->GetBinError(ibx,iby);
				h1		->Fill(a, val);
				h1e		->Fill(a, vale*vale);
				h1N		->Fill(a, 1.0);
      		}// iby
    	}// ibx
  	} else if (kProjAxis==3){		// project onto diff-axis of h2
    for (int ibx=1;ibx<=h2->GetXaxis()->GetNbins();ibx++){
    	for (int iby=1;iby<=h2->GetYaxis()->GetNbins();iby++){
			double x,y,a,d,val,vale;
				d		= h2->GetXaxis()->GetBinCenter(ibx) - h2->GetYaxis()->GetBinCenter(iby);
				val		= h2->GetBinContent(ibx,iby);
				vale	= h2->GetBinError(ibx,iby);
				h1		->Fill(d, val);
				h1e		->Fill(d, vale*vale);
				h1N		->Fill(d, 1.0);
      		}	// iby
    	}	// ibx
  	}
  	//
  	//---- now finally fill the projection plot h1
  	for (int ibi=1;ibi<=h1->GetNbinsX();ibi++){
    	double v,e,n;
    		v	=       h1	->GetBinContent(ibi);
    		e	= sqrt( h1e	->GetBinContent(ibi));
    		n	=       h1N	->GetBinContent(ibi);
    		if (n>0){
      			if (doavg){
					h1->SetBinContent(ibi,v/n);					
					h1->SetBinError(ibi,e/n);
      			} else if (doadd){
					h1->SetBinContent(ibi,v);					
					h1->SetBinError(ibi,e);
      			}
    		}
  	}
  	//
  	delete h1e;
  	delete h1N;
  	//
}


//========================================================================================
//
//======= get factorial moments...
//
void GetNPairs(TH2D *hmult){
  for (int i=0;i<10;i++){ multinfo[i] = 0.; }
  //
  double result 	= -1.;
  if (!hmult){ cout<<"no hmult... exit"<<endl; exit(0); }
  //
  double nent 	= (double)hmult->GetEntries();
  if (nent==0){ 
    //cout<<"no entries..."<<endl; 
    return result; 
  }
  TH2D *hmultwork		= (TH2D*)hmult->Clone("hmultwork");
  hmultwork			->Scale(1./nent);
  TAxis* multaxisx	= (TAxis*)hmultwork->GetXaxis();
  TAxis* multaxisy	= (TAxis*)hmultwork->GetYaxis();
  TH1D* hmultworkx	= (TH1D*)hmultwork->ProjectionX();
  TH1D* hmultworky	= (TH1D*)hmultwork->ProjectionY();
  multinfo[0]			= nent;
  //
  if (!fDistinguishable){
    for (int ibx=1;ibx<=multaxisx->GetNbins();ibx++){
      double n1i	 = multaxisx->GetBinCenter(ibx);
      double pi	 = hmultworkx->GetBinContent(ibx);		
      multinfo[1]	+= pi*n1i;			// f1 
      multinfo[3]	+= pi*n1i*(n1i-1.);	// f2 
    }
    for (int iby=1;iby<=multaxisy->GetNbins();iby++){
      double n2i	 = multaxisy->GetBinCenter(iby);
      double pi	 = hmultworky->GetBinContent(iby);		
      multinfo[2]	+= pi*n2i;			// f1, same as multinfo[1]
    }
  } else if (fDistinguishable){
    for (int ibx=1;ibx<=multaxisx->GetNbins();ibx++){
      for (int iby=1;iby<=multaxisy->GetNbins();iby++){
	double n1i	 = multaxisx->GetBinCenter(ibx);
	double n2i	 = multaxisy->GetBinCenter(iby);
	double pi	 = hmultwork->GetBinContent(ibx,iby);		
	multinfo[1]	+= pi*n1i;		// f1 particle 1
	multinfo[2]	+= pi*n2i;		// f1 particle 2
	multinfo[3]	+= pi*n1i*n2i;	// f2 
      }
    }
  }	// end fDistinguishable...
	//
  multinfo[4]	= multinfo[3] - multinfo[1]*multinfo[2];	// correlator integral
  //
  delete hmultwork;  hmultwork =0;
  delete hmultworkx; hmultworkx=0;
  delete hmultworky; hmultworky=0;
}

//========================================================================================
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
//========================================================================================
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
//========================================================================================
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
//========================================================================================
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
//========================================================================================


