//
//---- fluct_common.h
//---- main definitions used by many codes and macros...
//
	static const int 		NSpecies	= 7;
	static const double 	Species_mass[NSpecies] = {0.13956995,0.493677,0.93827231,
												 	  1.87561339,2.80925,2.80923,3.727417};
	static const char 	*SpeciesNames[NSpecies]   		  = {"#pi", "K", "p",  "d",  "t","^{3}He","#alpha"};
	//
//	static const double Species_ptmin[NSpecies]           = {  0.2, 0.2, 0.4,  0.8,  1.2,     1.2,      1.6};	
//	static const double Species_ptmax[NSpecies]           = {  2.0, 2.0, 2.0,  4.0,  4.0,     4.0,      4.0};	
//	static const double Species_pmax_dedx[NSpecies]       = {  0.6, 0.6, 0.9,  1.6,  1.8,     4.0,      4.0};	
//	static const double Species_pmax_tof[NSpecies]        = {  1.7, 1.7, 3.0,  6.0,  9.0,     9.0,     12.0};	
	static const double Species_ptmin[NSpecies]           = {  0.2, 0.2, 0.4,  0.2,  0.2,     0.2,      0.2};	
	static const double Species_ptmax[NSpecies]           = {  8.0, 8.0, 8.0,  8.0,  8.0,     8.0,      8.0};	
	static const double Species_pmax_dedx[NSpecies]       = {  8.0, 8.0, 8.0,  8.0,  8.0,     8.0,      8.0};	
	static const double Species_pmax_tof[NSpecies]        = {  8.0, 8.0, 8.0,  8.0,  8.0,     8.0,      8.0};	
	static const double Species_pmax_dedx_FXT[NSpecies]   = {  1.3, 0.6, 8.0,  8.0,  8.0,     8.0,      8.0};	
	static const double Species_pmax_tof_FXT[NSpecies]    = {  1.3, 0.6, 8.0,  8.0,  8.0,     8.0,      8.0};	
//
	static const double Species_Charge[NSpecies]          = {    1,   1,   1,    1,    1,       2,        2};	
	static const double Species_Ddedxcut[NSpecies]        = {    0,   0,   0,    1,    1,       3,        2};	
	static const double Species_Z[NSpecies]               = {    0,   0,   1,    1,    1,       2,        2};	
	static const double Species_A[NSpecies]               = {    0,   0,   1,    2,    3,       3,        4};	
	static const double Species_bichdiffcutu[NSpecies]    = {    0, 0.0, 0.0,  1.0,  1.1,     1.0,      1.0};	
	static const double Species_bichdiffcutd[NSpecies]    = {    0, 0.0, 0.0,  1.0,  1.1,     1.0,      1.0};	
	//
	static const double Species_ynbfxt[NSpecies]	= {  28,  24,  22,  16,  18,  18,  16};	
	static const double Species_ylfxt[NSpecies]		= {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};	
	static const double Species_yufxt[NSpecies]		= { 0.4, 0.2, 0.1,-0.2,-0.1,-0.1,-0.2};	
	//
	//	ParticleID = SpeciesID for chg>0, and =SpeciesID+NSpecies for chg<0
	//
	static const int NPARTICLENAMES					= 14;
	static const int kParticleIDPionPlus			=  0;
	static const int kParticleIDKaonPlus			=  1;
	static const int kParticleIDProtonPlus			=  2;
	static const int kParticleIDDeuteronPlus		=  3;
	static const int kParticleIDTritonPlus			=  4;
	static const int kParticleID3HePlus				=  5;
	static const int kParticleID4HePlus				=  6;
	static const int kParticleIDPionMinus			=  7;	
	static const int kParticleIDKaonMinus			=  8;
	static const int kParticleIDProtonMinus			=  9;
	static const int kParticleIDDeuteronMinus		= 10;
	static const int kParticleIDTritonMinus			= 11;
	static const int kParticleID3HeMinus			= 12;
	static const int kParticleID4HeMinus			= 13;
	//
	static const char *ParticleIDNames[NPARTICLENAMES] = {"#pi+", "K+", "p+",  "d+",  "t+","^{3}He+","#alpha+",
											   		 	  "#pi-", "K-", "p-",  "d-",  "t-","^{3}He-","#alpha-" };
	//
	static const int ParticlePDGids[NPARTICLENAMES]	= { 211, 321, 2212, 1000010020, 1000010030, 1000020030, 1000020040,
	                                                   -211,-321,-2212,-1000010020,-1000010030,-1000020030,-1000020040 };
	//
	//	PairTypes:	specifies 2 particles to calculate R2 for... 
	//			PairTypes[0][...] is ParticleID for particle 1 in this pair
	//			PairTypes[1][...] is ParticleID for particle 2 in this pair
	//			PairTypes[2][...] is kCentralityStyle for this pair: =2 refmult2, =3 refmult3
	//
// 	static const int NPairTypes		=  2;
// 	static const int PairTypes_Info[NPairTypes][3]= {	{ kParticleIDPionPlus,     kParticleIDPionPlus    ,2 },
// 														{ kParticleIDProtonPlus,   kParticleIDProtonPlus  ,3 } };
	//
// 	static const int NPairTypes		= 3;
// 	static const int PairTypes_Info[NPairTypes][3]= {{ kParticleIDPionPlus,     kParticleIDPionPlus    ,2 },
// 													 { kParticleIDPionPlus,     kParticleIDPionMinus   ,2 },
//  													 { kParticleIDProtonPlus,   kParticleIDProtonPlus  ,3 } };
	//	
// 	static const int NPairTypes		= 4;
// 	static const int PairTypes_Info[NPairTypes][3]= {	{ kParticleIDPionPlus,     kParticleIDPionPlus    ,2 },		// pi
// 														{ kParticleIDPionMinus,    kParticleIDPionMinus   ,2 },		// pi
// 														{ kParticleIDPionPlus,     kParticleIDPionMinus   ,2 },		// pi
// 														{ kParticleIDPionMinus,    kParticleIDPionPlus    ,2 } };	// pi
	//
	//---- ReaderRoot_2021
	static const int NPairTypes		= 12;
	static const int PairTypes_Info[NPairTypes][3]= {	{ kParticleIDPionPlus,     kParticleIDPionPlus    ,2 },		// pi
														{ kParticleIDPionMinus,    kParticleIDPionMinus   ,2 },		// pi
														{ kParticleIDPionPlus,     kParticleIDPionMinus   ,2 },		// pi
														{ kParticleIDPionMinus,    kParticleIDPionPlus    ,2 },		// pi
														{ kParticleIDKaonPlus,     kParticleIDKaonPlus    ,2 },		// K
														{ kParticleIDKaonMinus,    kParticleIDKaonMinus   ,2 },		// K
														{ kParticleIDKaonPlus,     kParticleIDKaonMinus   ,2 },		// K
														{ kParticleIDKaonMinus,    kParticleIDKaonPlus    ,2 },		// K
														{ kParticleIDProtonPlus,   kParticleIDProtonPlus  ,2 },		// p
														{ kParticleIDProtonMinus,  kParticleIDProtonMinus ,2 },		// p
														{ kParticleIDProtonPlus,   kParticleIDProtonMinus ,2 },		// p
														{ kParticleIDProtonMinus,  kParticleIDProtonPlus  ,2 } };	// p
	
	//---- ReaderRoot_2020
// 	static const int NPairTypes		= 21;
// 	static const int PairTypes_Info[NPairTypes][3]= {	{ kParticleIDPionPlus,     kParticleIDPionPlus    ,2 },		// pi
// 														{ kParticleIDPionMinus,    kParticleIDPionMinus   ,2 },		// pi
// 														{ kParticleIDPionPlus,     kParticleIDPionMinus   ,2 },		// pi
// 														{ kParticleIDPionMinus,    kParticleIDPionPlus    ,2 },		// pi
// 														{ kParticleIDKaonPlus,     kParticleIDKaonPlus    ,2 },		// K
// 														{ kParticleIDKaonMinus,    kParticleIDKaonMinus   ,2 },		// K
// 														{ kParticleIDKaonPlus,     kParticleIDKaonMinus   ,2 },		// K
// 														{ kParticleIDKaonMinus,    kParticleIDKaonPlus    ,2 },		// K
// 														{ kParticleIDProtonPlus,   kParticleIDProtonPlus  ,3 },		// p
// 														{ kParticleIDProtonMinus,  kParticleIDProtonMinus ,3 },		// p
// 														{ kParticleIDProtonPlus,   kParticleIDProtonMinus ,3 },		// p
// 														{ kParticleIDProtonMinus,  kParticleIDProtonPlus  ,3 },		// p
// 														{ kParticleIDProtonPlus,   kParticleIDDeuteronPlus,3 },		// pd
// 														{ kParticleIDProtonPlus,   kParticleIDTritonPlus  ,3 },		// pt
// 														{ kParticleIDProtonPlus,   kParticleID3HePlus     ,3 },		// ph
// 														{ kParticleIDDeuteronPlus, kParticleIDDeuteronPlus,3 },		// dd
// 														{ kParticleIDDeuteronPlus, kParticleIDTritonPlus  ,3 },		// dt
// 														{ kParticleIDDeuteronPlus, kParticleID3HePlus     ,3 },		// dh
// 														{ kParticleIDTritonPlus,   kParticleIDTritonPlus  ,3 },		// tt
// 														{ kParticleIDTritonPlus,   kParticleID3HePlus     ,3 },		// th
// 														{ kParticleID3HePlus,      kParticleID3HePlus     ,3 } };	// hh
	//
	//
	static const int NR2TYPES	= 3;	// used in reader to get hists. Also defined in class CalcRm - ensure consistency!!!
	//
	static const bool DUMPDATA			= true;
	//
	static const bool DODETADEPENDENCE	= false;
	//
static const int NE					=  10;   //No. of datasets
	static const float roots[NE] 		= { 7.7 ,  11.5,   14.5 ,   19.6,    27,    39,   62.4,   200,   3.0,  200  };
	static const char* energies[NE]		= {"7.7", "11.5", "14.5", "19.6",  "27",  "39", "62.4", "200", "3.0", "200" };
	static const char* estrings[NE]		= {"007", "011",   "015",  "019", "027", "039",  "062", "20010", "003", "20011" };
	static const int   edsid[NE]		= {  19 ,   20 ,     31 ,    23 ,   25 ,   18 ,    17 ,    16,    48,   24  };
	//
	//static const char* npartfiles[NE]	= {"007","011","019","019","027","039","062","200run11","3.0fxt","009"};
	//
	static const int NMODEL = 2;
	static const char* modnames[NMODEL] = {"urqmd","pythia"};	
	//
	static const int NCent			=  17;
	static const int NCentUse		=  10;
	static const int NPart			=   8;
	static const int NBIN			=   2;
	static const int NCUM			=   9;
	static const int NSAM			= 100;		// just needed for input filenames
	static const int NSAMfluct		= 100;		// just needed for input filenames	
	//
	static const char *particleids[NPart]	= {"#pi^{+}", "#pi^{-}", 
	                            	 			"K^{+}"  , "K^{-}"  , 
	                            	 			"p"      , "#bar{p}",
	                            	 			"e^{+}"  , "e^{-}"};
	static const char* centnames17[NCent]	= {"MB","75-80%","70-75%","65-70%","60-65%",
								  				"55-60%","50-55%","45-50%","40-45%","35-40%",
								  				"30-35%","25-30%","20-25%","15-20%","10-15%",
												"5-10%","0-5%"};
	static const char* centnames[NCentUse]	= {"MB","70-80%","60-70%","50-60%","40-50%",
								  				"30-40%","20-30%","10-20%","5-10%","0-5%"};
	static const char* centstrings[NCentUse]= {"MB","7080","6070","5060","4050",
								  				"3040","2030","1020","0510","0005"};
	static const char* centnamessim[12]		= {"MB","100-90%","90-80%","70-80%","60-70%","50-60%","40-50%",
								  				"30-40%","20-30%","10-20%","5-10%","0-5%"};
	//
	static const char* centnamesSD0[NCentUse]	= {"MB","70-80%","60-70%","50-60%","40-50%",
								  					"30-40%","20-30%","10-20%","5-10%","0-5%"};
	static const char* centnamesSD1[NCentUse]	= {"MB","30-35%","25-30%","20-25%","15-20%",
								  					"10-15%","7.5-10%","5-7.5%","2.5-5%","0-2.5%"};
	static const char* centnamesSD2[NCentUse]	= {"MB","10-12.5%","8.75-10%","7.5-8.75%","6.25-7.5%",
								  					"5-6.25%","3.75-5%","2.5-3.75%","1.25-2.5%","0-1.25%"};
	//
	//
	//
// 	static const char *particlenames[NPart]	= {"#pi^{+}+#pi^{-}", "#pi^{+}-#pi^{-}",
// 	                            	   			"K^{+}+K^{-}"    , "K^{+}-K^{-}"    ,
// 	                            	   			"p+#bar{p}"      , "p-#bar{p}"      };	                               
// 	static const char* groupnames[4]		= {"tot","net","pos","neg"};
// 	static const char* binnames[NBIN]		= {"SampSing","Data"};
// 	static const char* cumnames[NCUM]		= {"C_{1}","C_{2}","C_{3}","C_{4}",
// 												"R_{12}","R_{31}",
// 												"R_{32}","R_{42}","R_{43}"};
	//
// 	static const int NSET				=  12;
// 	static const char* setstrings[NSET]	= { "netp_TPC",				// 0
// 											"netq_TPC",				// 1
// 											"netk_TPC",				// 2
// 											"netpi_TPC",			// 3
// 											"netp_TPCmTOF",			// 4
// 											"netq_TPCmTOF",			// 5
// 											"netk_TPCmTOF",			// 6
// 											"netpi_TPCmTOF",		// 7
// 											"netp_TPCpTOF",			// 8
// 											"netq_TPCpTOF",			// 9
// 											"netk_TPCpTOF",			//10
// 											"netpi_TPCpTOF",		//11
// 											};
// 	static const char* setnames[NSET]	= { "p, |y|<0.5, RM3",					// 0
// 											"q, |#eta|<0.5, RM2",				// 1
// 											"k, |#eta|<0.5, RM2",				// 2
// 											"#pi, |#eta|<0.5, RM2",				// 3
// 											"p, |y|<0.5, RM3, mTOF",			// 4
// 											"q, |#eta|<0.5, RM2, mTOF",			// 5
// 											"k, |#eta|<0.5, RM2, mTOF",			// 6
// 											"#pi, |#eta|<0.5, RM2, mTOF",		// 7
// 											"p, |y|<0.5, RM3, pTOF",			// 8
// 											"q, |#eta|<0.5, RM2, pTOF",			// 9
// 											"k, |#eta|<0.5, RM2, pTOF",			//10
// 											"#pi, |#eta|<0.5, RM2, pTOF"		//11
// 											};
	//
//	static const int	deta_iecmax[NSET]	= {7,4,4,4,7,4,4,4,7,4,4,4};
	//
	//---- efftype: 1=netp, 2=netq 
//	static const int	efftype[NSET]		= {1,2,3,4,1,2,3,4,1,2,3,4};
	//
	//---- centvarid: 1=RM, 2=RM2, 3=RM3, 4=EMCE 
//	static const int	centvarids[NSET]			= {3,2,2,2,3,2,2,2,3,2,2,2};
//	static const int 	NCENTVARTYPES				= 4;
//	static const char*	centvarnames[NCENTVARTYPES]	= {"refmult","refmult2","refmult3","emcsume"};
	//
// 	static const int NUDYN_NPROFILES=  16; 
// 	static const char* nudyn_profnames[NUDYN_NPROFILES] = {
// 		"#LTN_{A}#GT",
// 		"#LTN_{B}#GT",
// 		"#LTN_{A}(N_{A}-1)#GT",
// 		"#LTN_{B}(N_{B}-1)#GT",
// 		"#LTN_{A}N_{B}#GT",
// 		"nab2",
// 		"na2b",
// 		"nab3",
// 		"na3b",
// 		"na2b2",
// 		"na2",
// 		"na3",
// 		"na4",
// 		"nb2",
// 		"nb3",
// 		"nb4"
// 	};
	//
	//
	//
// 	static const int NSETSIM				=  8;
// 	static const char* setstringssim[NSETSIM]	= { "netp",			// 0
// 											"netq",					// 1
// 											"netk",					// 2
// 											"netpi",				// 3
// 											"netp_tpceff",			// 4
// 											"netq_tpceff",			// 5
// 											"netk_tpceff",			// 6
// 											"netpi_tpceff"			// 7
// 											};
// 	static const char* setnamessim[NSETSIM]	= { "p, |y|<0.5, RM3",				// 0
// 											"q, |#eta|<0.5, RM2",				// 1
// 											"k, |#eta|<0.5, RM2",				// 2
// 											"#pi, |#eta|<0.5, RM2",				// 3
// 											"p, |y|<0.5, RM3, TPCeff",			// 4
// 											"q, |#eta|<0.5, RM2, TPCeff",		// 5
// 											"k, |#eta|<0.5, RM2, TPCeff",		// 6
// 											"#pi, |#eta|<0.5, RM2, TPCeff"		// 7
// 											};
// 	static const int centvaridssim[NSETSIM]	= {3,2,2,2,3,2,2,2};
//
//---- the end
//


