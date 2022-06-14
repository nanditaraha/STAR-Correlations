//
//	PairTypes:	specifies 2 particles to calculate R2 for... 
//			PairTypes[0][...] is ParticleID for particle 1 in this pair
//			PairTypes[1][...] is ParticleID for particle 2 in this pair
//			PairTypes[2][...] is kCentralityStyle for this pair: =2 refmult2, =3 refmult3
//
//---- PairTypes Mixed
//
static const int NPairTypes		= 24;
static const int PairTypes_Info[NPairTypes][3]= {	
						// pi-K
						{ kParticleIDPionPlus,     kParticleIDKaonPlus    ,2 },	// 
						{ kParticleIDPionMinus,    kParticleIDKaonMinus   ,2 },	// 
						{ kParticleIDPionPlus,     kParticleIDKaonMinus   ,2 },	// 
						{ kParticleIDPionMinus,    kParticleIDKaonPlus    ,2 },	// 
						// pi-p
						{ kParticleIDPionPlus,     kParticleIDProtonPlus  ,2 },	// 
						{ kParticleIDPionMinus,    kParticleIDProtonMinus ,2 },	// 
						{ kParticleIDPionPlus,     kParticleIDProtonMinus ,2 },	// 
						{ kParticleIDPionMinus,    kParticleIDProtonPlus  ,2 },	//
						// k-pi
						{ kParticleIDKaonPlus,     kParticleIDPionPlus    ,2 },	// 
						{ kParticleIDKaonMinus,    kParticleIDPionMinus   ,2 },	// 
						{ kParticleIDKaonPlus,     kParticleIDPionMinus   ,2 },	// 
						{ kParticleIDKaonMinus,    kParticleIDPionPlus    ,2 },	// 
						// k-p
						{ kParticleIDKaonPlus,     kParticleIDProtonPlus  ,2 },	// 
						{ kParticleIDKaonMinus,    kParticleIDProtonMinus ,2 },	// 
						{ kParticleIDKaonPlus,     kParticleIDProtonMinus ,2 },	// 
						{ kParticleIDKaonMinus,    kParticleIDProtonPlus  ,2 },	//
						// p-pi
						{ kParticleIDProtonPlus,   kParticleIDPionPlus    ,2 },	// 
						{ kParticleIDProtonMinus,  kParticleIDPionMinus   ,2 },	// 
						{ kParticleIDProtonPlus,   kParticleIDPionMinus   ,2 },	// 
						{ kParticleIDProtonMinus,  kParticleIDPionPlus    ,2 },	//
						// p-K
						{ kParticleIDProtonPlus,   kParticleIDKaonPlus    ,2 },	// 
						{ kParticleIDProtonMinus,  kParticleIDKaonMinus   ,2 },	// 
						{ kParticleIDProtonPlus,   kParticleIDKaonMinus   ,2 },	// 
						{ kParticleIDProtonMinus,  kParticleIDKaonPlus    ,2 }	// 
					};
//
//---- End PairTypes Mixed
