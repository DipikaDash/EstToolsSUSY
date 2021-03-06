#ifndef ESTTOOLS_LMPARAMETERS_HH_
#define ESTTOOLS_LMPARAMETERS_HH_

#include "../utils/EstHelper.hh"

namespace EstTools{

const TString inputdir = "../../macros/run/plots_19_01_05_0LSR";
const TString outputdir = ".";

const TString datadir = ".";
const TString datadir_2017 = "../plots_19_01_05_0LSR17";
const TString lumistr = "35.9";
const TString lumistr_2017 = "41.37";
//const TString lumistr = "41.37";

TString getLumi(){return lumistr(TRegexp("[0-9]+.[0-9]"));}

// lumi and base weight
const TString wgtvar = lumistr+"*weight*truePUWeight*btagWeight*topptWeight*sdMVAWeight*resTopWeight";
//const TString wgtvar = lumistr+"*weight*PUScale2017Reco*btagWeight";
const TString wgtvar_2017 = lumistr_2017+"*weight*truePUWeight*btagWeight";
//const TString wgtvar = lumistr+"*weight*topptWeight*truePUWeight*btagWeight";

// photon trigger eff.
const TString phowgt = wgtvar;
//const TString phowgt = wgtvar + "*qcdRespTailWeight";

// No Lepton SF
//const TString lepvetowgt = wgtvar;
//const TString lepselwgt  = wgtvar;
//const TString vetoes = " && nvetolep==0 && nvetotau==0";

// Tag-and-Probe Lepton SF
const TString lepvetowgt = wgtvar + "*(leptnpweightLM*lepvetoweightLM*(njets<5 || nbjets<1) + leptnpweightHM*lepvetoweightHM*(njets>=5 && nbjets>=1))";
const TString lepselwgt  = wgtvar + "*(leptnpweightLM*(njets<5 || nbjets<1) + leptnpweightHM*(njets>=5 && nbjets>=1))";
const TString lepvetowgt_2017 = wgtvar_2017 + "*(leptnpweightLM*lepvetoweightLM*(njets<5 || nbjets<1) + leptnpweightHM*lepvetoweightHM*(njets>=5 && nbjets>=1))";
const TString lepselwgt_2017  = wgtvar_2017 + "*(leptnpweightLM*(njets<5 || nbjets<1) + leptnpweightHM*(njets>=5 && nbjets>=1))";
const TString vetoes = " && nvetolep==0 && (nvetotau==0 || (ismc && npromptgentau>0))";

// 1LCR Lepton SF
//const TString lepvetowgt = wgtvar + "*lepvetoweight";
//const TString lepselwgt  = wgtvar + "*lepselweight";
//const TString vetoes = " && ((nvetolep==0 && nvetotau==0) || (ismc && (ngoodgenele>0 || ngoodgenmu>0 || npromptgentau>0)))";

// 1Lep LLB method
bool ADD_LEP_TO_MET = false;
bool ICHEPCR = false;
bool SPLITTF = true; // split TF to CR-SR and SR-extrapolation
bool isDphiCut = true; //Adding bool to add dphiCut > 0.7
const TString revert_vetoes = " && nvetolep>0 && mtlepmet<100";

// MET+LEP LL method
//bool ADD_LEP_TO_MET = true;
const TString lepcrsel = " && nvetolep==1 && mtlepmet<100 && origmet>100";
//const TString lepcrsel = " && nprimlep==1 && mtlepmet<100 && origmet>100";

// lepton trigger eff.
//const TString trigLepCR = " && (passtrige || passtrigmu)";
//const TString onelepcrwgt  = lepselwgt + "*trigEleWeight*trigMuWeight";

const TString trigLepCR = " && passtriglepOR";
const TString onelepcrwgt  = lepselwgt;

// qcd weights
const TString qcdwgt = wgtvar + "*qcdRespTailWeight";
//const TString qcdwgt = wgtvar;
const TString qcdvetowgt = lepvetowgt + "*qcdRespTailWeight";
//const TString qcdvetowgt = lepvetowgt;

// signal weights
const TString sigwgt = lepvetowgt + "*btagFastSimWeight*isrWeightTight*(0.85*(nivf>=1) + 1.0*(nivf==0))";
//const TString sigwgt = lepvetowgt;

// triggers
const TString trigSR = " && (passmetmht || ismc)";
const TString trigPhoCR = " && passtrigphoOR && origmet<200";
const TString phoBadEventRemoval = " && (!(lumi==189375 && event==430170481) && !(lumi==163479 && event==319690728) && !(lumi==24214 && event==55002562) && !(lumi==12510 && event==28415512) && !(lumi==16662 && event==32583938) && !(lumi==115657 && event==226172626) && !(lumi==149227 && event==431689582) && !(lumi==203626 && event==398201606))";
const TString trigDiLepCR = " && passtrigdilepOR && dileppt>200";
const TString datasel = " && passjson && (passmetfilters || process==10) && j1chEnFrac>0.1 && j1chEnFrac<0.99 && (origmet/calomet<5)";
//const TString datasel = " && (passmetfilters || process==10) && j1chEnFrac>0.1 && j1chEnFrac<0.99 && (origmet/calomet<5)";
//const TString datasel = " && passjetid && passjson && (passmetfilters || process==10) && j1chEnFrac>0.1 && j1chEnFrac<0.99";
const TString qcdSpikeRemovals = " && (!(lumi==40062 && event==91000735))";
const TString dphi_invert = " && (dphij1met<0.1 || dphij2met<0.1 || dphij3met<0.1)";
const TString dphi_cut = " && ( ((mtcsv12met<175 && nsdtop==0 && nsdw==0 && nrestop==0) && (dphij1met>0.5 && dphij2met>0.15 && dphij3met>0.15)) || (!(mtcsv12met<175 && nsdtop==0 && nsdw==0 && nrestop==0) && (dphij1met>0.5 && dphij2met>0.5 && dphij3met>0.5 && dphij4met>0.5)) )"; // ( ((passLM) && dPhiLM) || ((!passLM) && dPhiHM) )

// ------------------------------------------------------------------------
// search regions and control regions

const TString baseline = "met>250 && njets>=2";

std::map<TString, TString> cutMap = []{
    // Underscore "_" not allowed in the names!!!
    std::map<TString, TString> cmap = {
        {"lmNoDPhi",  "ak8isrpt>200 && dphiisrmet>2 && nsdtop==0 && nsdw==0 && nrestop==0 && metovsqrtht>10"},
        {"dPhiLM",    "dphij1met>0.5 && dphij2met>0.15 && dphij3met>0.15"},
        {"hmNoDPhi",  "njets>=5 && nbjets>=1"},
        {"dPhiHM",    "dphij1met>0.5 && dphij2met>0.5 && dphij3met>0.5 && dphij4met>0.5"},
        {"invertDPhi","(dphij1met<0.1 || dphij2met<0.1 || dphij3met<0.1)"},

        {"nb0",       "nbjets==0"},
        {"nb1",       "nbjets==1"},
        {"nb2",       "nbjets>=2"},
        {"nbeq2",     "nbjets==2"},
        {"nb3",       "nbjets>=3"},
        {"nivf0",     "nivf==0"},
        {"nivf1",     "nivf>=1"},
        {"lowptisr",  "ak8isrpt>300 && ak8isrpt<500"},
        {"highptisr", "ak8isrpt>500"},
        {"nj2to5",    "njets>=2 && njets<=5"},
        {"nj6",       "njets>=6"},
        {"nj7",       "njets>=7"},
        {"lowmtb",    "mtcsv12met<175"},
        {"highmtb",   "mtcsv12met>175"},
        {"lowptb",    "csvj1pt<40"},
        {"medptb",    "csvj1pt>40 && csvj1pt<70"},
        {"highptb",   "csvj1pt>70"},
        {"lowptb12",  "(csvj1pt+csvj2pt)<80"},
        {"medptb12",  "(csvj1pt+csvj2pt)>80 && (csvj1pt+csvj2pt)<140"},
        {"highptb12", "(csvj1pt+csvj2pt)>140"},
        {"nt0",       "nsdtop==0"},
        {"nt1",       "nsdtop==1"},
        {"nt2",       "nsdtop==2"},
        {"ntgeq1",    "nsdtop>=1"},
        {"nw0",       "nsdw==0"},
        {"nw1",       "nsdw==1"},
        {"nw2",       "nsdw==2"},
        {"nwgeq1",    "nsdw>=1"},
        {"nrt0",      "nrestop==0"},
        {"nrt1",      "nrestop==1"},
        {"nrt2",      "nrestop==2"},
        {"nrtgeq1",   "nrestop>=1"},
	{"nrtntnwgeq3", "(nsdtop+nrestop+nsdw) >= 3"},
	{"htlt1000",    "ht<1000"},	
	{"ht1000to1500",    "ht>1000 && ht<1500"},	
	{"htgt1500",    "ht>1500"},	
	{"htlt1300",    "ht<1300"},	
	{"htgt1300",    "ht>1300"},	
	{"dphitop",     "dphisrt12>0.7"},
    };

    cmap["lm"] = createCutString("lmNoDPhi_dPhiLM", cmap);
    cmap["hm"] = createCutString("hmNoDPhi_dPhiHM", cmap);
    return cmap;
}();


std::map<TString, TString> labelMap{
  {"lowptisr", R"($300<\ptisr<500$\,GeV)"},
  {"ntgeq1", R"($\nt\geq1$)"},
  {"nt2", R"($\nt=2$)"},
  {"nivf0", R"($\nsv=0$)"},
  {"nivf1", R"($\nsv\geq1$)"},
  {"nw2", R"($\nw=2$)"},
  {"nj2to5", R"($2\leq\nj\leq5$)"},
  {"nb2", R"($\nb\geq2$)"},
  {"nbeq2", R"($\nb=2$)"},
  {"nb3", R"($\nb\geq3$)"},
  {"nb1", R"($\nb=1$)"},
  {"nb0", R"($\nb=0$)"},
  {"nrt2", R"($\nrt=2$)"},
  {"highptisr", R"($\ptisr>500$\,GeV)"},
  {"nj7", R"($\nj\geq7$)"},
  {"highptb", R"($\ptb>70$\,GeV)"},
  {"hm", R"(high \dm)"},
  {"nw0", R"($\nw=0$)"},
  {"nwgeq1", R"($\nw\geq1$)"},
  {"nw1", R"($\nw=1$)"},
  {"nrt0", R"($\nrt=0$)"},
  {"nrt1", R"($\nrt=1$)"},
  {"lowptb", R"($\ptb<40$\,GeV)"},
  {"medptb", R"($40<\ptb<70$\,GeV)"},
  {"nt0", R"($\nt=0$)"},
  {"lm", R"(low \dm)"},
  {"lowptb12", R"($\ptbonetwo<80$\,GeV)"},
  {"highptb12", R"($\ptbonetwo>140$\,GeV)"},
  {"lowmtb", R"($\mtb<175$~\GeV)"},
  {"highmtb", R"($\mtb>175$~\GeV)"},
  {"nt1", R"($\nt=1$)"},
  {"medptb12", R"($80<\ptbonetwo<140$\,GeV)"},
  {"nrtgeq1", R"($\nrt\geq1$)"},
  {"nj6", R"($\nj\geq6$)"},
  {"nrtntnwgeq3", R"($\nt+\nrt+\nw >= 3$)"},
  {"htlt1000",    R"($\ht<1000$)"},	
  {"ht1000to1500",R"($1000<\ht<1500$)"},	
  {"htgt1500",    R"($\ht>1500$)"},	
  {"htlt1300",    R"($\ht<1300$)"},	
  {"htgt1300",    R"($\ht>1300$)"},	
  {"dphitop",     R"($\dphi>0.7$)"},
};

std::map<TString, TString> plotLabelMap{
  {"lowptisr", R"(300 < p_{T}(ISR) < 500)"},
  {"ntgeq1", R"(N_{t} #geq 1)"},
  {"nt2", R"(N_{t} = 2)"},
  {"nivf0", R"(N_{SV} = 0)"},
  {"nivf1", R"(N_{SV} #geq 1)"},
  {"nw2", R"(N_{W} = 2)"},
  {"nj2to5", R"(2 #leq N_{j} #leq 5)"},
  {"nb2", R"(N_{b} #geq 2)"},
  {"nbeq2", R"(N_{b} = 2)"},
  {"nb3", R"(N_{b} #geq 3)"},
  {"nb1", R"(N_{b} = 1)"},
  {"nb0", R"(N_{b} = 0)"},
  {"nrt2", R"(N_{res} = 2)"},
  {"highptisr", R"(p_{T}(ISR) > 500)"},
  {"nj7", R"(N_{j} #geq 7)"},
  {"highptb", R"(p_{T}(b) > 70)"},
  {"hm", R"(high #Deltam)"},
  {"nw0", R"(N_{W} = 0)"},
  {"nwgeq1", R"(N_{W} #geq 1)"},
  {"nw1", R"(N_{W} = 1)"},
  {"nrt0", R"(N_{res} = 0)"},
  {"nrt1", R"(N_{res} = 1)"},
  {"lowptb", R"(p_{T}(b) < 40)"},
  {"medptb", R"(40 < p_{T}(b) < 70)"},
  {"nt0", R"(N_{t} = 0)"},
  {"lm", R"(low #Deltam)"},
  {"lowptb12", R"(p_{T}(b_{12}) < 80)"},
  {"highptb12", R"(p_{T}(b_{12}) > 140)"},
  {"lowmtb", R"(M_{T}(b_{1,2},#vec{p}_{T}^{miss}) < 175)"},
  {"highmtb", R"(M_{T}(b_{1,2},#vec{p}_{T}^{miss}) > 175)"},
  {"nt1", R"(N_{t} = 1)"},
  {"medptb12", R"(80 < p_{T}(b_{12}) < 140)"},
  {"nrtgeq1", R"(N_{res} #geq 1)"},
  {"nj6", R"(N_{j} #geq 6)"},
  {"nrtntnwgeq3", R"(N_{t}+N_{res}+N_{W} #geq 3)"},
  {"htlt1000",    R"(H_{T}<1000)"},	
  {"ht1000to1500",R"(1000<H_{T}<1500)"},	
  {"htgt1500",    R"(H_{T}>1500)"},	
  {"htlt1300",    R"(H_{T}<1300)"},	
  {"htgt1300",    R"(H_{T}>1300)"},	
  {"dphitop",     R"(#Delta#phi(t_{1,2}, #slash{E}_{T})>0.7)"},
};


//std::vector<TString> srbins{
//  //---------- low deltaM ----------
//  // 0b, 0 or >=1 ivf
//  "lm_nb0_nivf0_highptisr_nj2to5",
//  "lm_nb0_nivf0_highptisr_nj6",
//  "lm_nb0_nivf1_highptisr_nj2to5",
//  "lm_nb0_nivf1_highptisr_nj6",
//
//  // 1b, 0 or >=1 ivf
//  "lm_nb1_nivf0_lowmtb_lowptisr_lowptb",
//  "lm_nb1_nivf0_lowmtb_lowptisr_medptb",
//  "lm_nb1_nivf0_lowmtb_highptisr_lowptb",
//  "lm_nb1_nivf0_lowmtb_highptisr_medptb",
//  // ---
//  "lm_nb1_nivf1_lowmtb_lowptb",
//
//  // 2b
//  "lm_nb2_lowmtb_lowptisr_lowptb12",
//  "lm_nb2_lowmtb_lowptisr_medptb12",
//  "lm_nb2_lowmtb_lowptisr_highptb12_nj7",
//  "lm_nb2_lowmtb_highptisr_lowptb12",
//  "lm_nb2_lowmtb_highptisr_medptb12",
//  "lm_nb2_lowmtb_highptisr_highptb12_nj7",
//  //---------- low deltaM ----------
//
//
//  //---------- high deltaM ----------
//  // low mtb
//  "hm_nb1_lowmtb_nj7_nrtgeq1",
//  "hm_nb2_lowmtb_nj7_nrtgeq1",
//
//  // high mtb
//  // 0 taggged
//  "hm_nb1_highmtb_nj7_nt0_nrt0_nw0",
//  "hm_nb2_highmtb_nj7_nt0_nrt0_nw0",
//
//	// nb1
//  // 1 tagged
//  "hm_nb1_highmtb_ntgeq1_nrt0_nw0",
//  "hm_nb1_highmtb_nt0_nrtgeq1_nw0",
//  // 1+1
//  "hm_nb1_highmtb_ntgeq1_nrt0_nwgeq1",
//  "hm_nb1_highmtb_nt0_nrtgeq1_nwgeq1",
//
//	// nb2
//  // 1 tagged
//  "hm_nb2_highmtb_nt1_nrt0_nw0",
//  "hm_nb2_highmtb_nt0_nrt1_nw0",
//  "hm_nb2_highmtb_nt0_nrt0_nw1",
//
//  // 1+1
//  "hm_nb2_highmtb_nt1_nrt0_nw1",
//  "hm_nb2_highmtb_nt0_nrt1_nw1",
//  "hm_nb2_highmtb_nt1_nrt1_nw0",
//
//  // 2
//  "hm_nb2_highmtb_nt2_nrt0_nw0",
//  "hm_nb2_highmtb_nt0_nrt2_nw0",
//  "hm_nb2_highmtb_nt0_nrt0_nw2",
//  //---------- high deltaM ----------
//
//};

std::vector<TString> srbins{
  //---------- low deltaM ----------
  // 0b, 0 or >=1 ivf
  "lm_nb0_nivf0_highptisr_nj2to5",
  "lm_nb0_nivf0_highptisr_nj6",
  "lm_nb0_nivf1_highptisr_nj2to5",
  "lm_nb0_nivf1_highptisr_nj6",

  // 1b, 0 or >=1 ivf
  "lm_nb1_nivf0_lowmtb_lowptisr_lowptb",
  "lm_nb1_nivf0_lowmtb_lowptisr_medptb",
  "lm_nb1_nivf0_lowmtb_highptisr_lowptb",
  "lm_nb1_nivf0_lowmtb_highptisr_medptb",
  // ---
  "lm_nb1_nivf1_lowmtb_lowptb",

  // 2b
  "lm_nb2_lowmtb_lowptisr_lowptb12",
  "lm_nb2_lowmtb_lowptisr_medptb12",
  "lm_nb2_lowmtb_lowptisr_highptb12_nj7",
  "lm_nb2_lowmtb_highptisr_lowptb12",
  "lm_nb2_lowmtb_highptisr_medptb12",
  "lm_nb2_lowmtb_highptisr_highptb12_nj7",
  //---------- low deltaM ----------


  //---------- high deltaM ----------
  // low mtb
  "hm_nb1_lowmtb_nj7_nrtgeq1",
  "hm_nb2_lowmtb_nj7_nrtgeq1",

  // high mtb
  // 0 taggged
  "hm_nb1_highmtb_nj7_nt0_nrt0_nw0",
  "hm_nb2_highmtb_nj7_nt0_nrt0_nw0",

	// nb1
  // 1 tagged
  "hm_nb1_highmtb_ntgeq1_nrt0_nw0_htlt1000",
  "hm_nb1_highmtb_ntgeq1_nrt0_nw0_ht1000to1500",
  "hm_nb1_highmtb_ntgeq1_nrt0_nw0_htgt1500",
  "hm_nb1_highmtb_nt0_nrt0_nwgeq1_htlt1300",
  "hm_nb1_highmtb_nt0_nrt0_nwgeq1_htgt1300",
  "hm_nb1_highmtb_nt0_nrtgeq1_nw0_htlt1000",
  "hm_nb1_highmtb_nt0_nrtgeq1_nw0_ht1000to1500",
  "hm_nb1_highmtb_nt0_nrtgeq1_nw0_htgt1500",

  // 1+1
  "hm_nb1_highmtb_ntgeq1_nrt0_nwgeq1",
  "hm_nb1_highmtb_ntgeq1_nrtgeq1_nw0",
  "hm_nb1_highmtb_nt0_nrtgeq1_nwgeq1",

	// nb2
  // 1 tagged
  "hm_nbeq2_highmtb_nt1_nrt0_nw0_htlt1000",
  "hm_nbeq2_highmtb_nt1_nrt0_nw0_ht1000to1500",
  "hm_nbeq2_highmtb_nt1_nrt0_nw0_htgt1500",
  "hm_nbeq2_highmtb_nt0_nrt0_nw1_htlt1300",
  "hm_nbeq2_highmtb_nt0_nrt0_nw1_htgt1300",
  "hm_nbeq2_highmtb_nt0_nrt1_nw0_htlt1000",
  "hm_nbeq2_highmtb_nt0_nrt1_nw0_ht1000to1500",
  "hm_nbeq2_highmtb_nt0_nrt1_nw0_htgt1500",

  // 1+1eq
  "hm_nbeq2_highmtb_nt1_nrt0_nw1",
  "hm_nbeq2_highmtb_nt0_nrt1_nw1",
  "hm_nbeq2_highmtb_nt1_nrt1_nw0_htlt1300",
  "hm_nbeq2_highmtb_nt1_nrt1_nw0_htgt1300",

  // 2
  "hm_nbeq2_highmtb_nt2_nrt0_nw0_htlt1300",
  "hm_nbeq2_highmtb_nt2_nrt0_nw0_htgt1300",
  "hm_nbeq2_highmtb_nt0_nrt2_nw0_htlt1300",
  "hm_nbeq2_highmtb_nt0_nrt2_nw0_htgt1300",
  "hm_nbeq2_highmtb_nt0_nrt0_nw2",

  // 3
  "hm_nbeq2_highmtb_nrtntnwgeq3_htlt1300",
  "hm_nbeq2_highmtb_nrtntnwgeq3_htgt1300",

	// nb3
  //1 tagged
  "hm_nb3_highmtb_nt1_nrt0_nw0_htlt1000",
  "hm_nb3_highmtb_nt1_nrt0_nw0_ht1000to1500",
  "hm_nb3_highmtb_nt1_nrt0_nw0_htgt1500",
  "hm_nb3_highmtb_nt0_nrt0_nw1",
  "hm_nb3_highmtb_nt0_nrt1_nw0_htlt1000",
  "hm_nb3_highmtb_nt0_nrt1_nw0_ht1000to1500",
  "hm_nb3_highmtb_nt0_nrt1_nw0_htgt1500",

  //1+1
  "hm_nb3_highmtb_nt1_nrt0_nw1",
  "hm_nb3_highmtb_nt0_nrt1_nw1",
  "hm_nb3_highmtb_nt1_nrt1_nw0_htlt1300",
  "hm_nb3_highmtb_nt1_nrt1_nw0_htgt1300",

  //2
  "hm_nb3_highmtb_nt2_nrt0_nw0_htlt1300",
  "hm_nb3_highmtb_nt2_nrt0_nw0_htgt1300",
  "hm_nb3_highmtb_nt0_nrt2_nw0_htlt1300",
  "hm_nb3_highmtb_nt0_nrt2_nw0_htgt1300",
  "hm_nb3_highmtb_nt0_nrt0_nw2",

  // 3
  "hm_nb3_highmtb_nrtntnwgeq3",
  
  //---------- high deltaM ----------

};

std::vector<TString> srbinsDphiCut{
  //---------- low deltaM ----------
  // 0b, 0 or >=1 ivf
  "lm_nb0_nivf0_highptisr_nj2to5",
  "lm_nb0_nivf0_highptisr_nj6",
  "lm_nb0_nivf1_highptisr_nj2to5",
  "lm_nb0_nivf1_highptisr_nj6",

  // 1b, 0 or >=1 ivf
  "lm_nb1_nivf0_lowmtb_lowptisr_lowptb",
  "lm_nb1_nivf0_lowmtb_lowptisr_medptb",
  "lm_nb1_nivf0_lowmtb_highptisr_lowptb",
  "lm_nb1_nivf0_lowmtb_highptisr_medptb",
  // ---
  "lm_nb1_nivf1_lowmtb_lowptb",

  // 2b
  "lm_nb2_lowmtb_lowptisr_lowptb12",
  "lm_nb2_lowmtb_lowptisr_medptb12",
  "lm_nb2_lowmtb_lowptisr_highptb12_nj7",
  "lm_nb2_lowmtb_highptisr_lowptb12",
  "lm_nb2_lowmtb_highptisr_medptb12",
  "lm_nb2_lowmtb_highptisr_highptb12_nj7",
  //---------- low deltaM ----------


  //---------- high deltaM ----------
  // low mtb
  "hm_nb1_lowmtb_nj7_nrtgeq1_dphitop",
  "hm_nb2_lowmtb_nj7_nrtgeq1_dphitop",

  // high mtb
  // 0 taggged
  "hm_nb1_highmtb_nj7_nt0_nrt0_nw0",
  "hm_nb2_highmtb_nj7_nt0_nrt0_nw0",

	// nb1
  // 1 tagged
  "hm_nb1_highmtb_ntgeq1_nrt0_nw0_htlt1000_dphitop",
  "hm_nb1_highmtb_ntgeq1_nrt0_nw0_ht1000to1500_dphitop",
  "hm_nb1_highmtb_ntgeq1_nrt0_nw0_htgt1500_dphitop",
  "hm_nb1_highmtb_nt0_nrt0_nwgeq1_htlt1300_dphitop",
  "hm_nb1_highmtb_nt0_nrt0_nwgeq1_htgt1300_dphitop",
  "hm_nb1_highmtb_nt0_nrtgeq1_nw0_htlt1000_dphitop",
  "hm_nb1_highmtb_nt0_nrtgeq1_nw0_ht1000to1500_dphitop",
  "hm_nb1_highmtb_nt0_nrtgeq1_nw0_htgt1500_dphitop",

  // 1+1
  "hm_nb1_highmtb_ntgeq1_nrt0_nwgeq1_dphitop",
  "hm_nb1_highmtb_ntgeq1_nrtgeq1_nw0_dphitop",
  "hm_nb1_highmtb_nt0_nrtgeq1_nwgeq1_dphitop",

	// nb2
  // 1 tagged
  "hm_nbeq2_highmtb_nt1_nrt0_nw0_htlt1000_dphitop",
  "hm_nbeq2_highmtb_nt1_nrt0_nw0_ht1000to1500_dphitop",
  "hm_nbeq2_highmtb_nt1_nrt0_nw0_htgt1500_dphitop",
  "hm_nbeq2_highmtb_nt0_nrt0_nw1_htlt1300_dphitop",
  "hm_nbeq2_highmtb_nt0_nrt0_nw1_htgt1300_dphitop",
  "hm_nbeq2_highmtb_nt0_nrt1_nw0_htlt1000_dphitop",
  "hm_nbeq2_highmtb_nt0_nrt1_nw0_ht1000to1500_dphitop",
  "hm_nbeq2_highmtb_nt0_nrt1_nw0_htgt1500_dphitop",

  // 1+1eq
  "hm_nbeq2_highmtb_nt1_nrt0_nw1_dphitop",
  "hm_nbeq2_highmtb_nt0_nrt1_nw1_dphitop",
  "hm_nbeq2_highmtb_nt1_nrt1_nw0_htlt1300_dphitop",
  "hm_nbeq2_highmtb_nt1_nrt1_nw0_htgt1300_dphitop",

  // 2
  "hm_nbeq2_highmtb_nt2_nrt0_nw0_htlt1300_dphitop",
  "hm_nbeq2_highmtb_nt2_nrt0_nw0_htgt1300_dphitop",
  "hm_nbeq2_highmtb_nt0_nrt2_nw0_htlt1300_dphitop",
  "hm_nbeq2_highmtb_nt0_nrt2_nw0_htgt1300_dphitop",
  "hm_nbeq2_highmtb_nt0_nrt0_nw2_dphitop",

  // 3
  "hm_nbeq2_highmtb_nrtntnwgeq3_htlt1300_dphitop",
  "hm_nbeq2_highmtb_nrtntnwgeq3_htgt1300_dphitop",

	// nb3
  //1 tagged
  "hm_nb3_highmtb_nt1_nrt0_nw0_htlt1000_dphitop",
  "hm_nb3_highmtb_nt1_nrt0_nw0_ht1000to1500_dphitop",
  "hm_nb3_highmtb_nt1_nrt0_nw0_htgt1500_dphitop",
  "hm_nb3_highmtb_nt0_nrt0_nw1_dphitop",
  "hm_nb3_highmtb_nt0_nrt1_nw0_htlt1000_dphitop",
  "hm_nb3_highmtb_nt0_nrt1_nw0_ht1000to1500_dphitop",
  "hm_nb3_highmtb_nt0_nrt1_nw0_htgt1500_dphitop",

  //1+1
  "hm_nb3_highmtb_nt1_nrt0_nw1_dphitop",
  "hm_nb3_highmtb_nt0_nrt1_nw1_dphitop",
  "hm_nb3_highmtb_nt1_nrt1_nw0_htlt1300_dphitop",
  "hm_nb3_highmtb_nt1_nrt1_nw0_htgt1300_dphitop",

  //2
  "hm_nb3_highmtb_nt2_nrt0_nw0_htlt1300_dphitop",
  "hm_nb3_highmtb_nt2_nrt0_nw0_htgt1300_dphitop",
  "hm_nb3_highmtb_nt0_nrt2_nw0_htlt1300_dphitop",
  "hm_nb3_highmtb_nt0_nrt2_nw0_htgt1300_dphitop",
  "hm_nb3_highmtb_nt0_nrt0_nw2_dphitop",

  // 3
  "hm_nb3_highmtb_nrtntnwgeq3_dphitop",
  
  //---------- high deltaM ----------

};

std::map<TString, TString> srcuts = []{
    std::map<TString, TString> cuts;
    
    for (auto name : srbins){
      if((!name.Contains("lm") && !name.Contains("nj7_nt0_nrt0_nw0")) && isDphiCut) name+="_dphitop";
      std::cout << "name: " << name << std::endl;
      cuts[name] = createCutString(name, cutMap);
    }
    return cuts;
}();

std::map<TString, TString> srlabels = []{
    std::map<TString, TString> cmap;
    for (auto s: srbins){
      if((!s.Contains("lm") && !s.Contains("nj7_nt0_nrt0_nw0")) && isDphiCut) s+="_dphitop";
      std::cout << "s: " << s << std::endl;
      cmap[s] = s;
    }
    return cmap;
}();

std::map<TString, std::vector<int>> srMETbins{
  //---------- low deltaM ----------
  // 0b, 0 or >=1 ivf
  {"lm_nb0_nivf0_highptisr_nj2to5",              {     450, 550, 650, 750, 1000}},
  {"lm_nb0_nivf0_highptisr_nj6",                 {     450, 550, 650, 750, 1000}},
  {"lm_nb0_nivf1_highptisr_nj2to5",              {     450, 550, 650, 750, 1000}},
  {"lm_nb0_nivf1_highptisr_nj6",                 {     450, 550, 650, 750, 1000}},

  // 1b, 0 or >=1 ivf
  {"lm_nb1_nivf0_lowmtb_lowptisr_lowptb",        {300, 400, 500, 600, 1000}},
  {"lm_nb1_nivf0_lowmtb_lowptisr_medptb",        {300, 400, 500, 600, 1000}},
  {"lm_nb1_nivf0_lowmtb_highptisr_lowptb",            {450, 550, 650, 750, 1000}},
  {"lm_nb1_nivf0_lowmtb_highptisr_medptb",            {450, 550, 650, 750, 1000}},

  {"lm_nb1_nivf1_lowmtb_lowptb",                 {300, 400, 500, 1000}},

  // 2b
  {"lm_nb2_lowmtb_lowptisr_lowptb12",            {300, 400, 500, 1000}},
  {"lm_nb2_lowmtb_lowptisr_medptb12",            {300, 400, 500, 1000}},
  {"lm_nb2_lowmtb_lowptisr_highptb12_nj7",       {300, 400, 500, 1000}},
  {"lm_nb2_lowmtb_highptisr_lowptb12",                {450, 550, 650, 1000}},
  {"lm_nb2_lowmtb_highptisr_medptb12",                {450, 550, 650, 1000}},
  {"lm_nb2_lowmtb_highptisr_highptb12_nj7",           {450, 550, 650, 1000}},
  //---------- low deltaM ----------


  //---------- high deltaM ----------
  // low mtb
  {"hm_nb1_lowmtb_nj7_nrtgeq1",                    {250, 300, 400, 500, 1000}},
  {"hm_nb2_lowmtb_nj7_nrtgeq1",                    {250, 300, 400, 500, 1000}},

  // high mtb
  // 0 taggged
  {"hm_nb1_highmtb_nj7_nt0_nrt0_nw0",              {250, 350, 450, 550, 1000}},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw0",              {250, 350, 450, 550, 1000}},

  // nb1
  {"hm_nb1_highmtb_ntgeq1_nrt0_nw0_htlt1000",      {250, 550, 650, 1000}},
  {"hm_nb1_highmtb_ntgeq1_nrt0_nw0_ht1000to1500",  {250, 550, 650, 1000}},
  {"hm_nb1_highmtb_ntgeq1_nrt0_nw0_htgt1500",      {250, 550, 650, 1000}},
  {"hm_nb1_highmtb_nt0_nrt0_nwgeq1_htlt1300",      {250, 350, 450, 1000}},
  {"hm_nb1_highmtb_nt0_nrt0_nwgeq1_htgt1300",      {250, 350, 450, 1000}},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nw0_htlt1000",      {250, 350, 450, 550, 650, 1000}},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nw0_ht1000to1500",  {250, 350, 450, 550, 650, 1000}},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nw0_htgt1500",      {250, 350, 450, 550, 650, 1000}},

  // 1+1
  {"hm_nb1_highmtb_ntgeq1_nrt0_nwgeq1",            {250, 550, 1000}},
  {"hm_nb1_highmtb_ntgeq1_nrtgeq1_nw0",            {250, 550, 1000}},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nwgeq1",            {250, 550, 1000}},

  // nb2
  {"hm_nbeq2_highmtb_nt1_nrt0_nw0_htlt1000",       {250, 550, 650, 1000}},
  {"hm_nbeq2_highmtb_nt1_nrt0_nw0_ht1000to1500",   {250, 550, 650, 1000}},
  {"hm_nbeq2_highmtb_nt1_nrt0_nw0_htgt1500",       {250, 550, 650, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt0_nw1_htlt1300",       {250, 350, 450, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt0_nw1_htgt1300",       {250, 350, 450, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt1_nw0_htlt1000",       {250, 350, 450, 550, 650, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt1_nw0_ht1000to1500",   {250, 350, 450, 550, 650, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt1_nw0_htgt1500",       {250, 350, 450, 550, 650, 1000}},

  // 1+1
  {"hm_nbeq2_highmtb_nt1_nrt0_nw1",                {250, 550, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt1_nw1",                {250, 550, 1000}},
  {"hm_nbeq2_highmtb_nt1_nrt1_nw0_htlt1300",       {250, 350, 450, 1000}},
  {"hm_nbeq2_highmtb_nt1_nrt1_nw0_htgt1300",       {250, 350, 450, 1000}},

  // 2
  {"hm_nbeq2_highmtb_nt2_nrt0_nw0_htlt1300",       {250, 450, 600, 1000}},
  {"hm_nbeq2_highmtb_nt2_nrt0_nw0_htgt1300",       {250, 450, 600, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt2_nw0_htlt1300",       {250, 450, 600, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt2_nw0_htgt1300",       {250, 450, 600, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt0_nw2",                {250, 1000}},

  // 3
  {"hm_nbeq2_highmtb_nrtntnwgeq3_htlt1300",        {250, 400, 1000}},
  {"hm_nbeq2_highmtb_nrtntnwgeq3_htgt1300",        {250, 400, 1000}},

	// nb3
  //1 tagged
  {"hm_nb3_highmtb_nt1_nrt0_nw0_htlt1000",         {250, 350, 550, 1000}},
  {"hm_nb3_highmtb_nt1_nrt0_nw0_ht1000to1500",     {250, 350, 550, 1000}},
  {"hm_nb3_highmtb_nt1_nrt0_nw0_htgt1500",         {250, 350, 550, 1000}},
  {"hm_nb3_highmtb_nt0_nrt0_nw1",                  {250, 350, 550, 1000}},
  {"hm_nb3_highmtb_nt0_nrt1_nw0_htlt1000",         {250, 350, 550, 1000}},
  {"hm_nb3_highmtb_nt0_nrt1_nw0_ht1000to1500",     {250, 350, 550, 1000}},
  {"hm_nb3_highmtb_nt0_nrt1_nw0_htgt1500",         {250, 350, 550, 1000}},

  //1+1
  {"hm_nb3_highmtb_nt1_nrt0_nw1",                  {250, 500, 1000}},
  {"hm_nb3_highmtb_nt0_nrt1_nw1",                  {250, 500, 1000}},
  {"hm_nb3_highmtb_nt1_nrt1_nw0_htlt1300",         {250, 350, 500, 1000}},
  {"hm_nb3_highmtb_nt1_nrt1_nw0_htgt1300",         {250, 350, 500, 1000}},

  //2
  {"hm_nb3_highmtb_nt2_nrt0_nw0_htlt1300",         {250, 500, 1000}},
  {"hm_nb3_highmtb_nt2_nrt0_nw0_htgt1300",         {250, 500, 1000}},
  {"hm_nb3_highmtb_nt0_nrt2_nw0_htlt1300",         {250, 500, 1000}},
  {"hm_nb3_highmtb_nt0_nrt2_nw0_htgt1300",         {250, 500, 1000}},
  {"hm_nb3_highmtb_nt0_nrt0_nw2",                  {250, 1000}},

  // 3
  {"hm_nb3_highmtb_nrtntnwgeq3",                   {250, 400, 1000}},
  
  //---------- high deltaM ----------

};

std::map<TString, std::vector<int>> srMETDphiTopbins{
  //---------- low deltaM ----------
  // 0b, 0 or >=1 ivf
  {"lm_nb0_nivf0_highptisr_nj2to5",              {     450, 550, 650, 750, 1000}},
  {"lm_nb0_nivf0_highptisr_nj6",                 {     450, 550, 650, 750, 1000}},
  {"lm_nb0_nivf1_highptisr_nj2to5",              {     450, 550, 650, 750, 1000}},
  {"lm_nb0_nivf1_highptisr_nj6",                 {     450, 550, 650, 750, 1000}},

  // 1b, 0 or >=1 ivf
  {"lm_nb1_nivf0_lowmtb_lowptisr_lowptb",        {300, 400, 500, 600, 1000}},
  {"lm_nb1_nivf0_lowmtb_lowptisr_medptb",        {300, 400, 500, 600, 1000}},
  {"lm_nb1_nivf0_lowmtb_highptisr_lowptb",            {450, 550, 650, 750, 1000}},
  {"lm_nb1_nivf0_lowmtb_highptisr_medptb",            {450, 550, 650, 750, 1000}},

  {"lm_nb1_nivf1_lowmtb_lowptb",                 {300, 400, 500, 1000}},

  // 2b
  {"lm_nb2_lowmtb_lowptisr_lowptb12",            {300, 400, 500, 1000}},
  {"lm_nb2_lowmtb_lowptisr_medptb12",            {300, 400, 500, 1000}},
  {"lm_nb2_lowmtb_lowptisr_highptb12_nj7",       {300, 400, 500, 1000}},
  {"lm_nb2_lowmtb_highptisr_lowptb12",                {450, 550, 650, 1000}},
  {"lm_nb2_lowmtb_highptisr_medptb12",                {450, 550, 650, 1000}},
  {"lm_nb2_lowmtb_highptisr_highptb12_nj7",           {450, 550, 650, 1000}},
  //---------- low deltaM ----------


  //---------- high deltaM ----------
  // low mtb
  {"hm_nb1_lowmtb_nj7_nrtgeq1_dphitop",                    {250, 300, 400, 500, 1000}},
  {"hm_nb2_lowmtb_nj7_nrtgeq1_dphitop",                    {250, 300, 400, 500, 1000}},

  // high mtb
  // 0 taggged
  {"hm_nb1_highmtb_nj7_nt0_nrt0_nw0",              {250, 350, 450, 550, 1000}},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw0",              {250, 350, 450, 550, 1000}},

  // nb1
  {"hm_nb1_highmtb_ntgeq1_nrt0_nw0_htlt1000_dphitop",      {250, 550, 650, 1000}},
  {"hm_nb1_highmtb_ntgeq1_nrt0_nw0_ht1000to1500_dphitop",  {250, 550, 650, 1000}},
  {"hm_nb1_highmtb_ntgeq1_nrt0_nw0_htgt1500_dphitop",      {250, 550, 650, 1000}},
  {"hm_nb1_highmtb_nt0_nrt0_nwgeq1_htlt1300_dphitop",      {250, 350, 450, 1000}},
  {"hm_nb1_highmtb_nt0_nrt0_nwgeq1_htgt1300_dphitop",      {250, 350, 450, 1000}},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nw0_htlt1000_dphitop",      {250, 350, 450, 550, 650, 1000}},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nw0_ht1000to1500_dphitop",  {250, 350, 450, 550, 650, 1000}},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nw0_htgt1500_dphitop",      {250, 350, 450, 550, 650, 1000}},

  // 1+1
  {"hm_nb1_highmtb_ntgeq1_nrt0_nwgeq1_dphitop",            {250, 550, 1000}},
  {"hm_nb1_highmtb_ntgeq1_nrtgeq1_nw0_dphitop",            {250, 550, 1000}},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nwgeq1_dphitop",            {250, 550, 1000}},

  // nb2
  {"hm_nbeq2_highmtb_nt1_nrt0_nw0_htlt1000_dphitop",       {250, 550, 650, 1000}},
  {"hm_nbeq2_highmtb_nt1_nrt0_nw0_ht1000to1500_dphitop",   {250, 550, 650, 1000}},
  {"hm_nbeq2_highmtb_nt1_nrt0_nw0_htgt1500_dphitop",       {250, 550, 650, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt0_nw1_htlt1300_dphitop",       {250, 350, 450, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt0_nw1_htgt1300_dphitop",       {250, 350, 450, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt1_nw0_htlt1000_dphitop",       {250, 350, 450, 550, 650, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt1_nw0_ht1000to1500_dphitop",   {250, 350, 450, 550, 650, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt1_nw0_htgt1500_dphitop",       {250, 350, 450, 550, 650, 1000}},

  // 1+1
  {"hm_nbeq2_highmtb_nt1_nrt0_nw1_dphitop",                {250, 550, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt1_nw1_dphitop",                {250, 550, 1000}},
  {"hm_nbeq2_highmtb_nt1_nrt1_nw0_htlt1300_dphitop",       {250, 350, 450, 1000}},
  {"hm_nbeq2_highmtb_nt1_nrt1_nw0_htgt1300_dphitop",       {250, 350, 450, 1000}},

  // 2
  {"hm_nbeq2_highmtb_nt2_nrt0_nw0_htlt1300_dphitop",       {250, 450, 600, 1000}},
  {"hm_nbeq2_highmtb_nt2_nrt0_nw0_htgt1300_dphitop",       {250, 450, 600, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt2_nw0_htlt1300_dphitop",       {250, 450, 600, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt2_nw0_htgt1300_dphitop",       {250, 450, 600, 1000}},
  {"hm_nbeq2_highmtb_nt0_nrt0_nw2_dphitop",                {250, 1000}},

  // 3
  {"hm_nbeq2_highmtb_nrtntnwgeq3_htlt1300_dphitop",        {250, 400, 1000}},
  {"hm_nbeq2_highmtb_nrtntnwgeq3_htgt1300_dphitop",        {250, 400, 1000}},

	// nb3
  //1 tagged
  {"hm_nb3_highmtb_nt1_nrt0_nw0_htlt1000_dphitop",         {250, 350, 550, 1000}},
  {"hm_nb3_highmtb_nt1_nrt0_nw0_ht1000to1500_dphitop",     {250, 350, 550, 1000}},
  {"hm_nb3_highmtb_nt1_nrt0_nw0_htgt1500_dphitop",         {250, 350, 550, 1000}},
  {"hm_nb3_highmtb_nt0_nrt0_nw1_dphitop",                  {250, 350, 550, 1000}},
  {"hm_nb3_highmtb_nt0_nrt1_nw0_htlt1000_dphitop",         {250, 350, 550, 1000}},
  {"hm_nb3_highmtb_nt0_nrt1_nw0_ht1000to1500_dphitop",     {250, 350, 550, 1000}},
  {"hm_nb3_highmtb_nt0_nrt1_nw0_htgt1500_dphitop",         {250, 350, 550, 1000}},

  //1+1
  {"hm_nb3_highmtb_nt1_nrt0_nw1_dphitop",                  {250, 500, 1000}},
  {"hm_nb3_highmtb_nt0_nrt1_nw1_dphitop",                  {250, 500, 1000}},
  {"hm_nb3_highmtb_nt1_nrt1_nw0_htlt1300_dphitop",         {250, 350, 500, 1000}},
  {"hm_nb3_highmtb_nt1_nrt1_nw0_htgt1300_dphitop",         {250, 350, 500, 1000}},

  //2
  {"hm_nb3_highmtb_nt2_nrt0_nw0_htlt1300_dphitop",         {250, 500, 1000}},
  {"hm_nb3_highmtb_nt2_nrt0_nw0_htgt1300_dphitop",         {250, 500, 1000}},
  {"hm_nb3_highmtb_nt0_nrt2_nw0_htlt1300_dphitop",         {250, 500, 1000}},
  {"hm_nb3_highmtb_nt0_nrt2_nw0_htgt1300_dphitop",         {250, 500, 1000}},
  {"hm_nb3_highmtb_nt0_nrt0_nw2_dphitop",                  {250, 1000}},

  // 3
  {"hm_nb3_highmtb_nrtntnwgeq3_dphitop",                   {250, 400, 1000}},
  
  //---------- high deltaM ----------

};

// normalization cuts for Rz
std::map<TString, TString> normMap{
  {"lm_nb0_nivf0", createCutString("lmNoDPhi_nb0_nivf0", cutMap)},
  {"lm_nb0_nivf1", createCutString("lmNoDPhi_nb0_nivf1", cutMap)},
  {"lm_nb1_nivf0", createCutString("lmNoDPhi_nb1_nivf0", cutMap)},
  {"lm_nb1_nivf1", createCutString("lmNoDPhi_nb1_nivf1", cutMap)},
  {"lm_nb2",       createCutString("lmNoDPhi_nb2", cutMap)},

  {"hm_nb1",       createCutString("hmNoDPhi_nb1", cutMap)},
  {"hm_nb2",       createCutString("hmNoDPhi_nb2", cutMap)},
};

// normalize photon to Data after baseline+this cut to calc Sgamma
std::map<TString, TString> phoNormMap = []{
  if (ICHEPCR) return normMap;
  else return std::map<TString, TString>{
    {"lm_nb0", createCutString("lmNoDPhi_nb0", cutMap)},
    {"lm_nb1", createCutString("lmNoDPhi_nb1", cutMap)},
    {"lm_nb2", createCutString("lmNoDPhi_nb2", cutMap)},

    {"hm_nb1", createCutString("hmNoDPhi_nb1", cutMap)},
    {"hm_nb2", createCutString("hmNoDPhi_nb2", cutMap)},
  };
}();

//std::map<TString, TString> phoNormMap = normMap;

std::map<TString, TString> phocrMapping{
  //---------- low deltaM ----------
  // 0b, 0 or >=1 ivf
  {"lm_nb0_nivf0_highptisr_nj2to5",              "lm_nb0_highptisr_nj2to5"},
  {"lm_nb0_nivf0_highptisr_nj6",                 "lm_nb0_highptisr_nj6"},
  {"lm_nb0_nivf1_highptisr_nj2to5",              "lm_nb0_highptisr_nj2to5"},
  {"lm_nb0_nivf1_highptisr_nj6",                 "lm_nb0_highptisr_nj6"},

  // 1b, 0 or >=1 ivf
  {"lm_nb1_nivf0_lowmtb_lowptisr_lowptb",        "lm_nb1_lowmtb_lowptisr_lowptb"},
  {"lm_nb1_nivf0_lowmtb_lowptisr_medptb",        "lm_nb1_lowmtb_lowptisr_medptb"},
  {"lm_nb1_nivf0_lowmtb_highptisr_lowptb",       "lm_nb1_lowmtb_highptisr_lowptb"},
  {"lm_nb1_nivf0_lowmtb_highptisr_medptb",       "lm_nb1_lowmtb_highptisr_medptb"},
  {"lm_nb1_nivf1_lowmtb_lowptb",                 "lm_nb1_lowmtb_lowptb"},

  // 2b
  {"lm_nb2_lowmtb_lowptisr_lowptb12",            "lm_nb2_lowmtb_lowptisr_lowptb12"},
  {"lm_nb2_lowmtb_lowptisr_medptb12",            "lm_nb2_lowmtb_lowptisr_medptb12"},
  {"lm_nb2_lowmtb_lowptisr_highptb12_nj7",       "lm_nb2_lowmtb_lowptisr_highptb12_nj7"},
  {"lm_nb2_lowmtb_highptisr_lowptb12",           "lm_nb2_lowmtb_highptisr_lowptb12"},
  {"lm_nb2_lowmtb_highptisr_medptb12",           "lm_nb2_lowmtb_highptisr_medptb12"},
  {"lm_nb2_lowmtb_highptisr_highptb12_nj7",      "lm_nb2_lowmtb_highptisr_highptb12_nj7"},
  //---------- low deltaM ----------


  //---------- high deltaM ----------
  // low mtb
  {"hm_nb1_lowmtb_nj7_nrtgeq1",        "hm_nb1_lowmtb_nj7"},
  {"hm_nb2_lowmtb_nj7_nrtgeq1",        "hm_nb2_lowmtb_nj7"},

  // high mtb
  // 0 taggged
  {"hm_nb1_highmtb_nj7_nt0_nrt0_nw0",  "hm_nb1_highmtb_nj7"},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw0",  "hm_nb2_highmtb_nj7"},

  // nb1
  // 1
  {"hm_nb1_highmtb_ntgeq1_nrt0_nw0",   "hm_nb1_highmtb"},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nw0",   "hm_nb1_highmtb"},
  // 1+1
  {"hm_nb1_highmtb_ntgeq1_nrt0_nwgeq1",  "hm_nb1_highmtb"},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nwgeq1",  "hm_nb1_highmtb"},

  // nb2
  //1
  {"hm_nb2_highmtb_nt1_nrt0_nw0",      "hm_nb2_highmtb"},
  {"hm_nb2_highmtb_nt0_nrt1_nw0",      "hm_nb2_highmtb"},
  {"hm_nb2_highmtb_nt0_nrt0_nw1",      "hm_nb2_highmtb"},

  // 1+1
  {"hm_nb2_highmtb_nt1_nrt0_nw1",        "hm_nb2_highmtb"},
  {"hm_nb2_highmtb_nt0_nrt1_nw1",        "hm_nb2_highmtb"},
  {"hm_nb2_highmtb_nt1_nrt1_nw0",        "hm_nb2_highmtb"},

  // 2
  {"hm_nb2_highmtb_nt2_nrt0_nw0",        "hm_nb2_highmtb"},
  {"hm_nb2_highmtb_nt0_nrt2_nw0",        "hm_nb2_highmtb"},
  {"hm_nb2_highmtb_nt0_nrt0_nw2",        "hm_nb2_highmtb"},
  //---------- high deltaM ----------

};


std::map<TString, TString> phocrCuts = []{
    std::map<TString, TString> cuts;
    for (auto sr2cr : phocrMapping)
      cuts[sr2cr.first] = createCutString(sr2cr.second, cutMap);
    return cuts;
}();

std::map<TString, TString> phocrlabels = phocrMapping;

std::map<TString, std::vector<int>> phocrMETbins = srMETbins;


std::map<TString, TString> lepcrMapping {
  //---------- low deltaM ----------
  // 0b, 0 or >=1 ivf
  {"lm_nb0_nivf0_highptisr_nj2to5",              "lm_nb0_nivf0_highptisr_nj2to5"},
  {"lm_nb0_nivf0_highptisr_nj6",                 "lm_nb0_nivf0_highptisr_nj6"},
  {"lm_nb0_nivf1_highptisr_nj2to5",              "lm_nb0_nivf1_highptisr_nj2to5"},
  {"lm_nb0_nivf1_highptisr_nj6",                 "lm_nb0_nivf1_highptisr_nj6"},

  // 1b, 0 or >=1 ivf
  {"lm_nb1_nivf0_lowmtb_lowptisr_lowptb",        "lm_nb1_nivf0_lowmtb_lowptisr_lowptb"},
  {"lm_nb1_nivf0_lowmtb_lowptisr_medptb",        "lm_nb1_nivf0_lowmtb_lowptisr_medptb"},
  {"lm_nb1_nivf0_lowmtb_highptisr_lowptb",       "lm_nb1_nivf0_lowmtb_highptisr_lowptb"},
  {"lm_nb1_nivf0_lowmtb_highptisr_medptb",       "lm_nb1_nivf0_lowmtb_highptisr_medptb"},
  {"lm_nb1_nivf1_lowmtb_lowptb",                 "lm_nb1_nivf1_lowmtb_lowptb"},

  // 2b
  {"lm_nb2_lowmtb_lowptisr_lowptb12",            "lm_nb2_lowmtb_lowptisr_lowptb12"},
  {"lm_nb2_lowmtb_lowptisr_medptb12",            "lm_nb2_lowmtb_lowptisr_medptb12"},
  {"lm_nb2_lowmtb_lowptisr_highptb12_nj7",       "lm_nb2_lowmtb_lowptisr_highptb12_nj7"},
  {"lm_nb2_lowmtb_highptisr_lowptb12",           "lm_nb2_lowmtb_highptisr_lowptb12"},
  {"lm_nb2_lowmtb_highptisr_medptb12",           "lm_nb2_lowmtb_highptisr_medptb12"},
  {"lm_nb2_lowmtb_highptisr_highptb12_nj7",      "lm_nb2_lowmtb_highptisr_highptb12_nj7"},
  //---------- low deltaM ----------


  //---------- high deltaM ----------
  // low mtb
  {"hm_nb1_lowmtb_nj7_nrtgeq1",                  	"hm_nb1_lowmtb_nj7"},
  {"hm_nb2_lowmtb_nj7_nrtgeq1",			 	"hm_nb2_lowmtb_nj7"},

  // high mtb
  // 0 taggged
  {"hm_nb1_highmtb_nj7_nt0_nrt0_nw0",		 	"hm_nb1_highmtb_nj7"},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw0",		 	"hm_nb2_highmtb_nj7"},

	// nb1
  // 1 tagged
  {"hm_nb1_highmtb_ntgeq1_nrt0_nw0_htlt1000",		"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_ntgeq1_nrt0_nw0_ht1000to1500",	"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_ntgeq1_nrt0_nw0_htgt1500",		"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_nt0_nrt0_nwgeq1_htlt1300",		"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_nt0_nrt0_nwgeq1_htgt1300",		"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nw0_htlt1000",		"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nw0_ht1000to1500",	"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nw0_htgt1500",		"hm_nb1_highmtb"},

  // 1+1
  {"hm_nb1_highmtb_ntgeq1_nrt0_nwgeq1",		 	"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_ntgeq1_nrtgeq1_nw0",		 	"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nwgeq1",		 	"hm_nb1_highmtb"},

	// nb2
  // 1 tagged
  {"hm_nbeq2_highmtb_nt1_nrt0_nw0_htlt1000",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt1_nrt0_nw0_ht1000to1500",	"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt1_nrt0_nw0_htgt1500",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt0_nw1_htlt1300",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt0_nw1_htgt1300",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt1_nw0_htlt1000",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt1_nw0_ht1000to1500",	"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt1_nw0_htgt1500",		"hm_nbeq2_highmtb"},

  // 1+1eq
  {"hm_nbeq2_highmtb_nt1_nrt0_nw1",		 	"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt1_nw1",		 	"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt1_nrt1_nw0_htlt1300",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt1_nrt1_nw0_htgt1300",		"hm_nbeq2_highmtb"},

  // 2
  {"hm_nbeq2_highmtb_nt2_nrt0_nw0_htlt1300",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt2_nrt0_nw0_htgt1300",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt2_nw0_htlt1300",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt2_nw0_htgt1300",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt0_nw2",		 	"hm_nbeq2_highmtb"},

  // 3
  {"hm_nbeq2_highmtb_nrtntnwgeq3_htlt1300",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nrtntnwgeq3_htgt1300",		"hm_nbeq2_highmtb"},

	// nb3
  //1 tagged
  {"hm_nb3_highmtb_nt1_nrt0_nw0_htlt1000",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt1_nrt0_nw0_ht1000to1500",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt1_nrt0_nw0_htgt1500",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt0_nrt0_nw1",		 	"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt0_nrt1_nw0_htlt1000",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt0_nrt1_nw0_ht1000to1500",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt0_nrt1_nw0_htgt1500",		"hm_nb3_highmtb"},

  //1+1
  {"hm_nb3_highmtb_nt1_nrt0_nw1",		 	"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt0_nrt1_nw1",		 	"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt1_nrt1_nw0_htlt1300",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt1_nrt1_nw0_htgt1300",		"hm_nb3_highmtb"},

  //2
  {"hm_nb3_highmtb_nt2_nrt0_nw0_htlt1300",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt2_nrt0_nw0_htgt1300",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt0_nrt2_nw0_htlt1300",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt0_nrt2_nw0_htgt1300",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt0_nrt0_nw2",		 	"hm_nb3_highmtb"},

  // 3
  {"hm_nb3_highmtb_nrtntnwgeq3",		 	"hm_nb3_highmtb"},
  
  //---------- high deltaM ----------
};

std::map<TString, TString> lepcrDphiTopMapping {
  //---------- low deltaM ----------
  // 0b, 0 or >=1 ivf
  {"lm_nb0_nivf0_highptisr_nj2to5",              "lm_nb0_nivf0_highptisr_nj2to5"},
  {"lm_nb0_nivf0_highptisr_nj6",                 "lm_nb0_nivf0_highptisr_nj6"},
  {"lm_nb0_nivf1_highptisr_nj2to5",              "lm_nb0_nivf1_highptisr_nj2to5"},
  {"lm_nb0_nivf1_highptisr_nj6",                 "lm_nb0_nivf1_highptisr_nj6"},

  // 1b, 0 or >=1 ivf
  {"lm_nb1_nivf0_lowmtb_lowptisr_lowptb",        "lm_nb1_nivf0_lowmtb_lowptisr_lowptb"},
  {"lm_nb1_nivf0_lowmtb_lowptisr_medptb",        "lm_nb1_nivf0_lowmtb_lowptisr_medptb"},
  {"lm_nb1_nivf0_lowmtb_highptisr_lowptb",       "lm_nb1_nivf0_lowmtb_highptisr_lowptb"},
  {"lm_nb1_nivf0_lowmtb_highptisr_medptb",       "lm_nb1_nivf0_lowmtb_highptisr_medptb"},
  {"lm_nb1_nivf1_lowmtb_lowptb",                 "lm_nb1_nivf1_lowmtb_lowptb"},

  // 2b
  {"lm_nb2_lowmtb_lowptisr_lowptb12",            "lm_nb2_lowmtb_lowptisr_lowptb12"},
  {"lm_nb2_lowmtb_lowptisr_medptb12",            "lm_nb2_lowmtb_lowptisr_medptb12"},
  {"lm_nb2_lowmtb_lowptisr_highptb12_nj7",       "lm_nb2_lowmtb_lowptisr_highptb12_nj7"},
  {"lm_nb2_lowmtb_highptisr_lowptb12",           "lm_nb2_lowmtb_highptisr_lowptb12"},
  {"lm_nb2_lowmtb_highptisr_medptb12",           "lm_nb2_lowmtb_highptisr_medptb12"},
  {"lm_nb2_lowmtb_highptisr_highptb12_nj7",      "lm_nb2_lowmtb_highptisr_highptb12_nj7"},
  //---------- low deltaM ----------


  //---------- high deltaM ----------
  // low mtb
  {"hm_nb1_lowmtb_nj7_nrtgeq1_dphitop",                  	"hm_nb1_lowmtb_nj7"},
  {"hm_nb2_lowmtb_nj7_nrtgeq1_dphitop",			 	"hm_nb2_lowmtb_nj7"},

  // high mtb
  // 0 taggged
  {"hm_nb1_highmtb_nj7_nt0_nrt0_nw0",		 	"hm_nb1_highmtb_nj7"},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw0",		 	"hm_nb2_highmtb_nj7"},

	// nb1
  // 1 tagged
  {"hm_nb1_highmtb_ntgeq1_nrt0_nw0_htlt1000_dphitop",		"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_ntgeq1_nrt0_nw0_ht1000to1500_dphitop",	"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_ntgeq1_nrt0_nw0_htgt1500_dphitop",		"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_nt0_nrt0_nwgeq1_htlt1300_dphitop",		"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_nt0_nrt0_nwgeq1_htgt1300_dphitop",		"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nw0_htlt1000_dphitop",		"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nw0_ht1000to1500_dphitop",	"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nw0_htgt1500_dphitop",		"hm_nb1_highmtb"},

  // 1+1
  {"hm_nb1_highmtb_ntgeq1_nrt0_nwgeq1_dphitop",		 	"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_ntgeq1_nrtgeq1_nw0_dphitop",		 	"hm_nb1_highmtb"},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nwgeq1_dphitop",		 	"hm_nb1_highmtb"},

	// nb2
  // 1 tagged
  {"hm_nbeq2_highmtb_nt1_nrt0_nw0_htlt1000_dphitop",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt1_nrt0_nw0_ht1000to1500_dphitop",	"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt1_nrt0_nw0_htgt1500_dphitop",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt0_nw1_htlt1300_dphitop",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt0_nw1_htgt1300_dphitop",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt1_nw0_htlt1000_dphitop",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt1_nw0_ht1000to1500_dphitop",	"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt1_nw0_htgt1500_dphitop",		"hm_nbeq2_highmtb"},

  // 1+1eq
  {"hm_nbeq2_highmtb_nt1_nrt0_nw1_dphitop",		 	"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt1_nw1_dphitop",		 	"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt1_nrt1_nw0_htlt1300_dphitop",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt1_nrt1_nw0_htgt1300_dphitop",		"hm_nbeq2_highmtb"},

  // 2
  {"hm_nbeq2_highmtb_nt2_nrt0_nw0_htlt1300_dphitop",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt2_nrt0_nw0_htgt1300_dphitop",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt2_nw0_htlt1300_dphitop",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt2_nw0_htgt1300_dphitop",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nt0_nrt0_nw2_dphitop",		 	"hm_nbeq2_highmtb"},

  // 3
  {"hm_nbeq2_highmtb_nrtntnwgeq3_htlt1300_dphitop",		"hm_nbeq2_highmtb"},
  {"hm_nbeq2_highmtb_nrtntnwgeq3_htgt1300_dphitop",		"hm_nbeq2_highmtb"},

	// nb3
  //1 tagged
  {"hm_nb3_highmtb_nt1_nrt0_nw0_htlt1000_dphitop",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt1_nrt0_nw0_ht1000to1500_dphitop",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt1_nrt0_nw0_htgt1500_dphitop",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt0_nrt0_nw1_dphitop",		 	"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt0_nrt1_nw0_htlt1000_dphitop",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt0_nrt1_nw0_ht1000to1500_dphitop",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt0_nrt1_nw0_htgt1500_dphitop",		"hm_nb3_highmtb"},

  //1+1
  {"hm_nb3_highmtb_nt1_nrt0_nw1_dphitop",		 	"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt0_nrt1_nw1_dphitop",		 	"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt1_nrt1_nw0_htlt1300_dphitop",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt1_nrt1_nw0_htgt1300_dphitop",		"hm_nb3_highmtb"},

  //2
  {"hm_nb3_highmtb_nt2_nrt0_nw0_htlt1300_dphitop",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt2_nrt0_nw0_htgt1300_dphitop",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt0_nrt2_nw0_htlt1300_dphitop",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt0_nrt2_nw0_htgt1300_dphitop",		"hm_nb3_highmtb"},
  {"hm_nb3_highmtb_nt0_nrt0_nw2_dphitop",		 	"hm_nb3_highmtb"},

  // 3
  {"hm_nb3_highmtb_nrtntnwgeq3_dphitop",		 	"hm_nb3_highmtb"},
  
  //---------- high deltaM ----------
};


std::map<TString, TString> lepcrCuts = []{
    std::map<TString, TString> cuts;
    for (auto sr2cr : lepcrMapping){
      TString sr2crCut = sr2cr.first;
      if((!sr2cr.first.Contains("lm") && !sr2cr.first.Contains("nj7_nt0_nrt0_nw0")) && isDphiCut) sr2crCut+="_dphitop";
      cuts[sr2crCut] = createCutString(sr2cr.second, cutMap);
    }
    return cuts;
}();

std::map<TString, TString> lepcrlabels = lepcrMapping;
std::map<TString, TString> lepcrDphiToplabels = lepcrDphiTopMapping;
std::map<TString, std::vector<int>> lepcrMETbins = srMETbins;
std::map<TString, std::vector<int>> lepcrMETDphiTopbins = srMETDphiTopbins;

// qcd-cr: inverted dPhi cut applied on CR samples now
std::map<TString, TString> qcdcrMapping =[]{
  auto crmap = lepcrMapping;
  for (auto &s : crmap){
    s.second.ReplaceAll("lm_", "lmNoDPhi_");
    s.second.ReplaceAll("hm_", "hmNoDPhi_");
  }
  return crmap;
}();
std::map<TString, TString> qcdcrCuts = []{
    std::map<TString, TString> cuts;
    for (auto sr2cr : qcdcrMapping)
      cuts[sr2cr.first] = createCutString(sr2cr.second, cutMap);
    return cuts;
}();
std::map<TString, TString> qcd1to1crCuts = []{
    std::map<TString, TString> cuts;
    for (auto name : srbins){
      TString crname = name;
      crname.ReplaceAll("lm_", "lmNoDPhi_");
      crname.ReplaceAll("hm_", "hmNoDPhi_");
      cuts[name] = createCutString(crname, cutMap);
    }
    return cuts;
}();
std::map<TString, TString> qcdcrlabels = lepcrlabels;
std::map<TString, std::vector<int>> qcdcrMETbins {
  //---------- low deltaM ----------
  // 0b, 0 or >=1 ivf
  {"lm_nb0_nivf0_highptisr_nj2to5",              {     450, 550, 650, 750, 1000}},
  {"lm_nb0_nivf0_highptisr_nj6",                 {     450, 550, 650, 750, 1000}},
  {"lm_nb0_nivf1_highptisr_nj2to5",              {     450, 550, 650,      1000}}, // merge last 2 bins
  {"lm_nb0_nivf1_highptisr_nj6",                 {     450,                1000}}, // merge last 3 bins

  // 1b, 0 or >=1 ivf
  {"lm_nb1_nivf0_lowmtb_lowptisr_lowptb",        {300, 400,           1000}}, // merge last 3 bins
  {"lm_nb1_nivf0_lowmtb_lowptisr_medptb",        {300, 400,           1000}}, // merge last 3 bins
  {"lm_nb1_nivf0_lowmtb_highptisr_lowptb",            {450, 550, 650,      1000}}, // merge last 2 bins
  {"lm_nb1_nivf0_lowmtb_highptisr_medptb",            {450, 550, 650,      1000}}, // merge last 2 bins

  {"lm_nb1_nivf1_lowmtb_lowptb",                 {300, 400,      1000}}, // merge last 2 bins

  // 2b
  {"lm_nb2_lowmtb_lowptisr_lowptb12",            {300,           1000}}, // merge all 3 bins
  {"lm_nb2_lowmtb_lowptisr_medptb12",            {300, 400,      1000}}, // merge last 2 bins
  {"lm_nb2_lowmtb_lowptisr_highptb12_nj7",       {300, 400,      1000}}, // merge last 2 bins
  {"lm_nb2_lowmtb_highptisr_lowptb12",                {450,           1000}}, // merge all 3 bins
  {"lm_nb2_lowmtb_highptisr_medptb12",                {450,           1000}}, // merge all 3 bins
  {"lm_nb2_lowmtb_highptisr_highptb12_nj7",           {450, 550,      1000}}, // merge last 2 bins
  //---------- low deltaM ----------


  //---------- high deltaM ----------
  // low mtb
  {"hm_nb1_lowmtb_nj7_nrtgeq1",        {250, 300, 400, 500, 1000}},
  {"hm_nb2_lowmtb_nj7_nrtgeq1",        {250, 300, 400, 500, 1000}},

  // high mtb
  // 0 taggged
  {"hm_nb1_highmtb_nj7_nt0_nrt0_nw0",  {250, 350, 450, 550, 1000}},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw0",  {250, 350, 450, 550, 1000}},

  // nb1
  // 1
  {"hm_nb1_highmtb_ntgeq1_nrt0_nw0",                  {550, 650, 1000}},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nw0",   {250, 350, 450, 550, 650, 1000}},
  // 1+1
  {"hm_nb1_highmtb_ntgeq1_nrt0_nwgeq1",               {550, 1000}},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nwgeq1",{250, 350, 450, 550, 1000}},

  // nb2
  // 1
  {"hm_nb2_highmtb_nt1_nrt0_nw0",                     {550, 650, 1000}},
  {"hm_nb2_highmtb_nt0_nrt0_nw1",      {250, 350, 450, 550, 650, 1000}},
  {"hm_nb2_highmtb_nt0_nrt1_nw0",      {250, 350, 450, 550, 650, 1000}},

  // 1+1
  {"hm_nb2_highmtb_nt1_nrt0_nw1",                     {550, 1000}},
  {"hm_nb2_highmtb_nt0_nrt1_nw1",      {250, 350, 450, 550, 1000}},
  {"hm_nb2_highmtb_nt1_nrt1_nw0",      {250, 350, 450, 1000}},

  // 2
  {"hm_nb2_highmtb_nt2_nrt0_nw0",      {250, 1000}},
  {"hm_nb2_highmtb_nt0_nrt2_nw0",      {250, 1000}},
  {"hm_nb2_highmtb_nt0_nrt0_nw2",      {250, 1000}},
  //---------- high deltaM ----------
};


map<TString, Category> srCatMap(){
  map<TString, Category> cmap;
  for (auto &name : srbins){
    auto nameMet = name;
    if((!name.Contains("lm") && !name.Contains("nj7_nt0_nrt0_nw0")) && isDphiCut) name+="_dphitop";
    cmap[name] = Category(name, srcuts.at(name), srlabels.at(name), BinInfo("met", "#slash{E}_{T}", srMETbins.at(nameMet), "GeV"));
  }
  return cmap;
}

map<TString, Category> phoCatMap(){
  map<TString, Category> cmap;
  const auto &cuts = ICHEPCR ? srcuts: phocrCuts;
  const auto &labels = ICHEPCR ? srlabels: phocrlabels;
  for (auto &name : srbins){
    cmap[name] = Category(name, cuts.at(name), labels.at(name), BinInfo("met", "#slash{E}_{T}^{#gamma}", phocrMETbins.at(name), "GeV"));
  }
  return cmap;
}

map<TString, Category> lepCatMap(){
  TString varlabel = ADD_LEP_TO_MET ? "p_{T}^{W}" : "#slash{E}_{T}";
  map<TString, Category> cmap;
  //map<TString, TString> labels;
  const auto &cuts = ICHEPCR ? srcuts: lepcrCuts;
  const auto &labels = ICHEPCR ? srlabels: lepcrlabels;
  //if(isDphiCut) labels = ICHEPCR ? srlabels: lepcrDphiToplabels;
  for (auto &name : srbins){
    //if((!name.Contains("lm") && !name.Contains("nj7_nt0_nrt0_nw0")) && isDphiCut) name+="_dphitop";
    std::cout << "name: " << name << std::endl;
    std::cout << "cutName: " << cuts.at(name) << std::endl;
    std::cout << "labelName: " << labels.at(name) << std::endl;
    //std::cout << "METName: " << lepcrMETDphiTopbins.at(name).first << std::endl;
    if(isDphiCut) cmap[name] = Category(name, cuts.at(name), labels.at(name), BinInfo("met", varlabel, lepcrMETbins.at(name), "GeV"));
    else          cmap[name] = Category(name, cuts.at(name), labels.at(name), BinInfo("met", varlabel, lepcrMETDphiTopbins.at(name), "GeV"));
  }
  return cmap;
}

map<TString, Category> qcdCatMap(){
  map<TString, Category> cmap;
  const auto &cuts = ICHEPCR ? qcd1to1crCuts: qcdcrCuts;
  const auto &labels = ICHEPCR ? srlabels: qcdcrlabels;
  for (auto &name : srbins){
    cmap[name] = Category(name, cuts.at(name), labels.at(name), BinInfo("met", "#slash{E}_{T}", qcdcrMETbins.at(name), "GeV"));
  }
  return cmap;
}

map<TString, Category> zllCatMap{
  {"on-z",  Category("on-z",  "dilepmass > 80 && dilepmass < 100",                      "on Z",   BinInfo("met", "#slash{E}_{T}^{ll}", vector<double>{200, 1000}, "GeV"))},
  {"off-z", Category("off-z", "dilepmass > 50 && (dilepmass < 80 || dilepmass > 100)",  "off Z",  BinInfo("met", "#slash{E}_{T}^{ll}", vector<double>{200, 1000}, "GeV"))}
};


// ------------------------------------------------------------------------
// samples

BaseConfig phoConfig(){
  BaseConfig     config;

  config.inputdir = inputdir;
  config.outputdir = outputdir+"/znunu";
  config.header = "#sqrt{s} = 13 TeV, "+lumistr+" fb^{-1}";

  config.addSample("singlepho",   "Data",           datadir+"/photoncr/singlepho",  "1.0",  datasel + trigPhoCR);

  config.addSample("photon",      "Photon",         "photoncr/photon",     phowgt, datasel + trigPhoCR + phoBadEventRemoval);
//  config.addSample("photon",      "Photon",         "photoncr/gjets",      phowgt, datasel + trigPhoCR);
  config.addSample("znunu-sr",    "Z#rightarrow#nu#nu",   "sr/znunu",      lepvetowgt, datasel + trigSR + vetoes);

  config.sel = baseline;
  config.categories = srbins;
  config.catMaps = srCatMap();
  config.crCatMaps = phoCatMap();

  return config;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BaseConfig zllConfig(){
  BaseConfig     config;

  config.inputdir = inputdir;
  config.outputdir = outputdir+"/zllcr";
  config.header = "#sqrt{s} = 13 TeV, "+lumistr+" fb^{-1}";

  config.addSample("dyll",      "Z#rightarrowll+jets",    "zllcr/z-soup",                    lepselwgt, datasel + trigDiLepCR);
  config.addSample("ttbar",     "t#bar{t}",               "zllcr/t-soup",                    lepselwgt, datasel + trigDiLepCR);
  config.addSample("doublelep", "Data",                   datadir+"/zllcr/doublelep",       "1.0",     datasel + trigDiLepCR);

  config.sel = baseline;
  config.catMaps = zllCatMap;
  for (auto &c : zllCatMap) config.categories.push_back(c.first);

  return config;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BaseConfig lepConfig(){
  BaseConfig     config;

  config.inputdir = inputdir;
  config.outputdir = outputdir+"/LLB";
  config.header = "#sqrt{s} = 13 TeV, "+lumistr+" fb^{-1}";

  // samples for cr categories
  if (ADD_LEP_TO_MET){
    config.addSample("singlelep",   "Data",          datadir+"/lepcr/singlelep",       "1.0",     datasel + trigLepCR + lepcrsel);
    config.addSample("ttbar",       "t#bar{t}",      "lepcr/ttbar",           onelepcrwgt, datasel + trigLepCR + lepcrsel);
    config.addSample("wjets",       "W+jets",        "lepcr/wjets",           onelepcrwgt, datasel + trigLepCR + lepcrsel);
    config.addSample("tW",          "tW",            "lepcr/tW",              onelepcrwgt, datasel + trigLepCR + lepcrsel);
    config.addSample("ttW",         "ttW",           "lepcr/ttW",             onelepcrwgt, datasel + trigLepCR + lepcrsel);
//    config.addSample("qcd",         "QCD",           "lepcr/qcd",             onelepcrwgt, datasel + trigLepCR + lepcrsel);
  }else{
    config.addSample("singlelep",   "Data",          datadir+"/sr/met",    "1.0",          datasel + trigSR + revert_vetoes);
    config.addSample("singlelep-2017",   "Data",     datadir_2017+"/sr/met",    "1.0",     datasel + trigSR + revert_vetoes);
    config.addSample("ttbar",       "t#bar{t}",      "sr/ttbar",           lepselwgt,      datasel + trigSR + revert_vetoes);
    config.addSample("wjets",       "W+jets",        "sr/wjets",           lepselwgt,      datasel + trigSR + revert_vetoes);
    config.addSample("ttbar-2017",  "t#bar{t}",      datadir_2017+"/sr/ttbar",  lepselwgt_2017, datasel + trigSR + revert_vetoes);
    config.addSample("wjets-2017",  "W+jets",        datadir_2017+"/sr/wjets",  lepselwgt_2017, datasel + trigSR + revert_vetoes);
    config.addSample("tW",          "tW",            "sr/tW",              lepselwgt,      datasel + trigSR + revert_vetoes);
    config.addSample("ttW",         "ttW",           "sr/ttW",             lepselwgt,      datasel + trigSR + revert_vetoes);
//    config.addSample("qcd",         "QCD",           "sr/qcd",             lepselwgt, datasel + trigSR + revert_vetoes);
  }

  // samples for sr categories
  config.addSample("ttbar-sr",       "t#bar{t}",      "sr/ttbar",                lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("wjets-sr",       "W+jets",        "sr/wjets",                lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("ttbar-sr-2017",  "t#bar{t}",      datadir_2017+"/sr/ttbar",  lepvetowgt_2017, datasel + trigSR + vetoes);
  config.addSample("wjets-sr-2017",  "W+jets",        datadir_2017+"/sr/wjets",  lepvetowgt_2017, datasel + trigSR + vetoes);
  config.addSample("tW-sr",          "tW",            "sr/tW",                   lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("ttW-sr",         "ttW",           "sr/ttW",                  lepvetowgt, datasel + trigSR + vetoes);
//  config.addSample("qcd-sr",         "QCD",           "qcd",                     lepvetowgt, datasel + trigSR + vetoes);
//  config.addSample("rare-sr",        "Rare",          "sr/rare",                 lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("ttZ-sr",         "ttZ",           "sr/ttZ",                  lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("diboson-sr",     "Diboson",       "sr/diboson",              lepvetowgt, datasel + trigSR + vetoes);

  // samples for splitting the TF (optional, see l.splitTF)
  if (SPLITTF){
    config.addSample("ttbar-sr-int",       "t#bar{t}",      "sr/ttbar",           lepvetowgt, datasel + trigSR + vetoes);
    config.addSample("wjets-sr-int",       "W+jets",        "sr/wjets",           lepvetowgt, datasel + trigSR + vetoes);
    config.addSample("ttbar-sr-int-2017",  "t#bar{t}",      datadir_2017+"/sr/ttbar",  lepvetowgt_2017, datasel + trigSR + vetoes);
    config.addSample("wjets-sr-int-2017",  "W+jets",        datadir_2017+"/sr/wjets",  lepvetowgt_2017, datasel + trigSR + vetoes);
    config.addSample("tW-sr-int",          "tW",            "sr/tW",              lepvetowgt, datasel + trigSR + vetoes);
    config.addSample("ttW-sr-int",         "ttW",           "sr/ttW",             lepvetowgt, datasel + trigSR + vetoes);
    config.addSample("ttZ-sr-int",         "ttZ",           "sr/ttZ",             lepvetowgt, datasel + trigSR + vetoes);
    config.addSample("diboson-sr-int",     "Diboson",       "sr/diboson",         lepvetowgt, datasel + trigSR + vetoes);
  }

  config.sel = baseline;
  if(isDphiCut) config.categories = srbinsDphiCut;
  else          config.categories = srbins;
  config.catMaps = srCatMap();
  config.crCatMaps = lepCatMap();

  return config;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BaseConfig srConfig(){
  BaseConfig     config;

  config.inputdir = inputdir;
  config.outputdir = outputdir+"/testSR";
  config.header = "#sqrt{s} = 13 TeV, "+lumistr+" fb^{-1}";

  config.addSample("ttbar",       "t#bar{t}",      "sr/ttbar",        lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("wjets",       "W+jets",        "sr/wjets",        lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("znunu",       "Z#rightarrow#nu#nu", "sr/znunu",   lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("qcd",         "QCD",           "sr/qcd-sr",       lepvetowgt, datasel + trigSR + vetoes + qcdSpikeRemovals);
  config.addSample("tW",          "tW",            "sr/tW",           lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("ttW",         "ttW",           "sr/ttW",          lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("ttZ",         "ttZ",           "sr/ttZ",          lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("diboson",     "Diboson",       "sr/diboson",      lepvetowgt, datasel + trigSR + vetoes);

//  config.addSample("T2fbd_500_420", "T2fbd(500,420)", "signals/T2fbd_500_420",  sigwgt, datasel + trigSR + vetoes);
  config.addSample("T2fbd_500_450", "T2fbd(500,450)", "signals/T2fbd_500_450",  sigwgt, datasel + trigSR + vetoes);
//  config.addSample("T2fbd_500_480", "T2fbd(500,480)", "signals/T2fbd_500_480",  sigwgt, datasel + trigSR + vetoes);
//  config.addSample("T2cc_500_490",  "T2cc(500,490)",  "signals/T2cc_500_490",   sigwgt, datasel + trigSR + vetoes);

//  config.addSample("T2tt_450_250",  "T2tt(450,250)",  "signals/T2tt_450_250",   sigwgt, datasel + trigSR + vetoes);
//  config.addSample("T2tt_400_300",  "T2tt(400,300)",  "signals/T2tt_400_300",   sigwgt, datasel + trigSR + vetoes); //FIXME
    config.addSample("T2tt_400_300",  "T2tt(425,300)",  "signals/T2tt_425_300",   sigwgt, datasel + trigSR + vetoes); //FIXME
  config.addSample("T2tt_700_400",  "T2tt(700,400)",  "signals/T2tt_700_400",   sigwgt, datasel + trigSR + vetoes);
  config.addSample("T2tt_1000_1",   "T2tt(1000,1)",   "signals/T2tt_1000_1",    sigwgt, datasel + trigSR + vetoes);
//  config.addSample("T2tt_1100_1",   "T2tt(1100,1)",   "signals/T2tt_1100_1",    sigwgt, datasel + trigSR + vetoes);
//  config.addSample("T2bW_850_1",    "T2bW(850,1)",    "signals/T2bW_850_1",     sigwgt, datasel + trigSR + vetoes);
//  config.addSample("T2bW_550_350",  "T2bW(550,350)",  "signals/T2bW_550_350",   sigwgt, datasel + trigSR + vetoes);


  COLOR_MAP["T2fbd_500_420"] = kRed;
  COLOR_MAP["T2fbd_500_450"] = kBlue;
  COLOR_MAP["T2fbd_500_480"] = kBlack;
  COLOR_MAP["T2cc_500_490"]  = kMagenta;

  COLOR_MAP["T2tt_450_250"]  = kOrange;
  COLOR_MAP["T2tt_700_400"]  = kCyan+2;
  COLOR_MAP["T2tt_1000_1"]   = kViolet+2;
  COLOR_MAP["T2bW_550_350"]  = kSpring-9;
  COLOR_MAP["T2bW_850_1"]    = kPink+2;

  config.sel = baseline;
  config.categories = srbins;
  config.catMaps = srCatMap();

  return config;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BaseConfig qcdConfig(){
  BaseConfig     config;

  config.inputdir = inputdir;
  config.outputdir = outputdir+"/QCD";
  config.header = "#sqrt{s} = 13 TeV, "+lumistr+" fb^{-1}";

  // qcdcr
  config.addSample("data-cr",     "Data",             datadir+"/sr/met",  "1.0",      datasel + trigSR + vetoes + dphi_invert);
  config.addSample("qcd-cr",      "QCD",              "sr/qcd-cr",       qcdwgt,      datasel + trigSR + dphi_invert + qcdSpikeRemovals);

  config.addSample("qcd-withveto-cr",  "QCD",         "sr/qcd-cr",       qcdvetowgt,  datasel + trigSR + vetoes + dphi_invert + qcdSpikeRemovals);
  config.addSample("qcd-withveto-sr",  "QCD",         "sr/qcd-sr",       qcdvetowgt,  datasel + trigSR + vetoes + qcdSpikeRemovals);

  // qcdcr: other bkg subtraction
  config.addSample("ttbar-cr",       "t#bar{t}",      "sr/ttbar",        lepvetowgt,  datasel + trigSR + vetoes + dphi_invert);
  config.addSample("wjets-cr",       "W+jets",        "sr/wjets",        lepvetowgt,  datasel + trigSR + vetoes + dphi_invert);
  config.addSample("tW-cr",          "tW",            "sr/tW",           lepvetowgt,  datasel + trigSR + vetoes + dphi_invert);
  config.addSample("ttW-cr",         "ttW",           "sr/ttW",          lepvetowgt,  datasel + trigSR + vetoes + dphi_invert);
  config.addSample("znunu-cr",       "Znunu",         "sr/znunu",        lepvetowgt,  datasel + trigSR + vetoes + dphi_invert);

  // onelepcr: norm correction for other bkg subtraction
  config.addSample("data-norm",      "Data",          datadir+"/sr/met", "1.0",       datasel + trigSR + revert_vetoes + dphi_cut);
  config.addSample("ttbar-norm",     "t#bar{t}",      "sr/ttbar",        lepselwgt,   datasel + trigSR + revert_vetoes + dphi_cut);
  config.addSample("wjets-norm",     "W+jets",        "sr/wjets",        lepselwgt,   datasel + trigSR + revert_vetoes + dphi_cut);
  config.addSample("tW-norm",        "tW",            "sr/tW",           lepselwgt,   datasel + trigSR + revert_vetoes + dphi_cut);
  config.addSample("ttW-norm",       "ttW",           "sr/ttW",          lepselwgt,   datasel + trigSR + revert_vetoes + dphi_cut);
  config.addSample("qcd-norm",       "QCD",           "sr/qcd-sr",       lepselwgt,   datasel + trigSR + revert_vetoes + dphi_cut + qcdSpikeRemovals);

  // qcdsr
  config.addSample("qcd-sr",         "QCD",           "sr/qcd-sr",       qcdwgt,      datasel + trigSR + qcdSpikeRemovals);

  // samples for splitting the TF, optional (see z.splitTF)
  if (SPLITTF){
    config.addSample("qcd-sr-int",     "QCD",           "sr/qcd-sr",       qcdwgt,      datasel + trigSR + dphi_cut + qcdSpikeRemovals);
  }

  config.sel = baseline;
  config.categories = srbins;
  config.catMaps = srCatMap();    // sr DPhi cut applied
  config.crCatMaps = qcdCatMap(); // no DPhi cut in category def

  return config;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BaseConfig sigConfig(){
  BaseConfig     config;

  config.inputdir = inputdir;
  config.outputdir = outputdir+"/sig";
  config.header = "#sqrt{s} = 13 TeV, "+lumistr+" fb^{-1}";

  config.addSample("data-sr",        "Data",             datadir+"/sr/met",                    "1.0",  datasel + trigSR + vetoes);

//  config.addSample("T2fbd_375_355",  "T2-4bd(375,355)",  "sig/T2fbd_375_355",  sigwgt, datasel + trigSR + vetoes);
//  config.addSample("T2fbd_375_325",  "T2-4bd(375,325)",  "sig/T2fbd_375_325",  sigwgt, datasel + trigSR + vetoes);
//  config.addSample("T2fbd_375_295",  "T2-4bd(375,295)",  "sig/T2fbd_375_295",  sigwgt, datasel + trigSR + vetoes);
//
  config.sel = baseline;
  config.categories = srbins;
  config.catMaps = srCatMap();

  return config;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

map<TString, BinInfo> varDict {
  {"norm",      BinInfo("met", "#slash{E}_{T}", vector<int>{0, 1000}, "GeV")},
  {"met",       BinInfo("met", "#slash{E}_{T}", vector<int>{250, 350, 450, 550, 650, 750, 1000}, "GeV")},
  {"metgx",       BinInfo("met", "#slash{E}_{T}^{(#gamma)}", vector<int>{250, 350, 450, 550, 650, 850}, "GeV")},
  {"metzg",       BinInfo("met", "#slash{E}_{T}^{#gamma/ll}", vector<int>{250, 350, 450, 550, 650, 850}, "GeV")},
  {"origmet",   BinInfo("origmet", "Original #slash{E}_{T}", 20, 0, 500, "GeV")},
  {"njets",     BinInfo("njets", "N_{j}", 8, -0.5, 7.5)},
  {"njl",       BinInfo("njl", "N_{j}^{ISR}", 4, 0.5, 4.5)},
  {"nlbjets",   BinInfo("nlbjets", "N_{B}^{loose}", 5, -0.5, 4.5)},
  {"nbjets",    BinInfo("nbjets",  "N_{B}^{medium}", 5, -0.5, 4.5)},
  {"dphij1met", BinInfo("dphij1met", "#Delta#phi(j_{1},#slash{E}_{T})", 30, 0, 3)},
  {"dphij2met", BinInfo("dphij2met", "#Delta#phi(j_{2},#slash{E}_{T})", 30, 0, 3)},
  {"metovsqrtht",BinInfo("metovsqrtht", "#slash{E}_{T}/#sqrt{H_{T}}", 10, 0, 20)},
  {"dphiisrmet",BinInfo("dphiisrmet", "#Delta#phi(j_{1}^{ISR},#slash{E}_{T})", vector<double>{0, 2, 3})},
  {"dphiisrmet_fine",BinInfo("dphiisrmet", "#Delta#phi(j_{1}^{ISR},#slash{E}_{T})", 12, 0, 3)},
  {"mtcsv12met",BinInfo("mtcsv12met", "min(m_{T}(b_{1},#slash{E}_{T}),m_{T}(b_{2},#slash{E}_{T}))", 6, 0, 300)},
  {"leptonpt",  BinInfo("leptonpt", "p_{T}^{lep} [GeV]", 12, 0, 600)},
  {"leptonptovermet",  BinInfo("leptonpt/met", "p_{T}^{lep}/#slash{E}_{T}", 20, 0, 1.)},
  {"ak8isrpt",  BinInfo("ak8isrpt", "p_{T}(ISR) [GeV]",  6, 200, 800)},
  {"csvj1pt",   BinInfo("csvj1pt", "p_{T}(b_{1}) [GeV]", 8, 20, 100)},
  {"ptb12",     BinInfo("csvj1pt+csvj2pt", "p_{T}(b_{1})+p_{T}(b_{2}) [GeV]", 8, 40, 200)},
  {"dphilepisr",  BinInfo("dphilepisr", "#Delta#phi(lep, j_{1}^{ISR})", 30, 0, 3)},
  {"drlepisr",  BinInfo("drlepisr", "#DeltaR(lep, j_{1}^{ISR})", 25, 0, 5)},
};

}
#endif /* ESTTOOLS_LMPARAMETERS_HH_ */
