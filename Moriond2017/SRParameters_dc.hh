#ifndef ESTTOOLS_LMPARAMETERS_HH_
#define ESTTOOLS_LMPARAMETERS_HH_

#include "../utils/EstHelper.hh"

namespace EstTools{

const TString inputdir = "/data/hqu/trees/20170206_lepSF";
const TString outputdir = "/tmp/hqu/plots/20170206_lepSF/unblind_36p8ifb_datacard";

const TString datadir = ".";
const TString lumistr = "36.8";

TString getLumi(){return lumistr(TRegexp("[0-9]+.[0-9]"));}

// lumi and base weight
const TString wgtvar = lumistr+"*weight*truePUWeight*btagWeight*topptWeight*sdMVAWeight*resTopWeight";
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
const TString vetoes = " && nvetolep==0 && (nvetotau==0 || (ismc && npromptgentau>0))";

// 1LCR Lepton SF
//const TString lepvetowgt = wgtvar + "*lepvetoweight";
//const TString lepselwgt  = wgtvar + "*lepselweight";
//const TString vetoes = " && ((nvetolep==0 && nvetotau==0) || (ismc && (ngoodgenele>0 || ngoodgenmu>0 || npromptgentau>0)))";

// sr weight w/o top/W SF
const TString lepvetowgt_no_wtopsf = lumistr+"*weight*truePUWeight*btagWeight*topptWeight*(leptnpweightLM*lepvetoweightLM*(njets<5 || nbjets<1) + leptnpweightHM*lepvetoweightHM*(njets>=5 && nbjets>=1))";

// 1Lep LLB method
bool ADD_LEP_TO_MET = false;
bool ICHEPCR = false;
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
const TString sigwgt = lepvetowgt + "*btagFastSimWeight*isrWeightTight";
const TString siglepselwgt = lepselwgt + "*btagFastSimWeight*isrWeightTight";

// triggers
const TString trigSR = " && (passmetmht || ismc)";
const TString trigPhoCR = " && passtrigphoOR && origmet<200";
const TString phoBadEventRemoval = " && (!(lumi==189375 && event==430170481) && !(lumi==163479 && event==319690728) && !(lumi==24214 && event==55002562) && !(lumi==12510 && event==28415512) && !(lumi==16662 && event==32583938) && !(lumi==115657 && event==226172626) && !(lumi==149227 && event==431689582) && !(lumi==203626 && event==398201606))";
const TString trigDiLepCR = " && passtrigdilepOR && dileppt>200";
const TString datasel = " && passjson && (passmetfilters || process==10) && j1chEnFrac>0.1 && j1chEnFrac<0.99 && (origmet/calomet<5)";
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
        {"nivf0",     "nivf==0"},
        {"nivf1",     "nivf>=1"},
        {"lowptisr",  "ak8isrpt>300 && ak8isrpt<500"},
        {"highptisr", "ak8isrpt>500"},
        {"nj2to5",    "njets>=2 && njets<=5"},
        {"nj6",       "njets>=6"},
        {"nj5to6",    "njets>=5 && njets<=6"},
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
        {"nt2",       "nsdtop>=2"},
        {"ntgeq1",    "nsdtop>=1"},
        {"nw0",       "nsdw==0"},
        {"nw1",       "nsdw==1"},
        {"nw2",       "nsdw>=2"},
        {"nwgeq1",    "nsdw>=1"},
        {"nrt0",      "nrestop==0"},
        {"nrt1",      "nrestop==1"},
        {"nrt2",      "nrestop>=2"},
        {"nrtgeq1",   "nrestop>=1"},
    };

    cmap["lm"] = createCutString("lmNoDPhi_dPhiLM", cmap);
    cmap["hm"] = createCutString("hmNoDPhi_dPhiHM", cmap);
    return cmap;
}();

vector<TString> signals {
  "T2tt_150_1",
  "T2tt_150_25",
  "T2tt_150_50",
  "T2tt_150_63",
  "T2tt_167_1",
  "T2tt_175_1",
  "T2tt_175_25",
  "T2tt_175_50",
  "T2tt_175_75",
  "T2tt_175_88",
  "T2tt_183_1",
  "T2tt_192_25",
  "T2tt_200_1",
  "T2tt_200_25",
  "T2tt_200_50",
  "T2tt_200_75",
  "T2tt_200_100",
  "T2tt_200_113",
  "T2tt_208_25",
  "T2tt_217_50",
  "T2tt_225_1",
  "T2tt_225_25",
  "T2tt_225_50",
  "T2tt_225_75",
  "T2tt_225_100",
  "T2tt_225_125",
  "T2tt_225_138",
  "T2tt_233_50",
  "T2tt_242_75",
  "T2tt_250_25",
  "T2tt_250_50",
  "T2tt_250_75",
  "T2tt_250_100",
  "T2tt_250_125",
  "T2tt_250_150",
  "T2tt_250_163",
  "T2tt_258_75",
  "T2tt_267_100",
  "T2tt_275_1",
  "T2tt_275_25",
  "T2tt_275_50",
  "T2tt_275_75",
  "T2tt_275_100",
  "T2tt_275_125",
  "T2tt_275_150",
  "T2tt_275_175",
  "T2tt_275_188",
  "T2tt_283_100",
  "T2tt_292_125",
  "T2tt_300_1",
  "T2tt_300_25",
  "T2tt_300_50",
  "T2tt_300_75",
  "T2tt_300_100",
  "T2tt_300_125",
  "T2tt_300_150",
  "T2tt_300_175",
  "T2tt_300_200",
  "T2tt_300_213",
  "T2tt_308_125",
  "T2tt_317_150",
  "T2tt_325_75",
  "T2tt_325_100",
  "T2tt_325_125",
  "T2tt_325_150",
  "T2tt_325_200",
  "T2tt_333_150",
  "T2tt_342_175",
  "T2tt_350_125",
  "T2tt_350_150",
  "T2tt_350_175",
  "T2tt_350_200",
  "T2tt_350_225",
  "T2tt_350_250",
  "T2tt_350_263",
  "T2tt_358_175",
  "T2tt_367_200",
  "T2tt_375_75",
  "T2tt_375_125",
  "T2tt_375_175",
  "T2tt_375_225",
  "T2tt_392_225",
  "T2tt_400_175",
  "T2tt_400_275",
  "T2tt_408_225",
  "T2tt_425_125",
  "T2tt_425_150",
  "T2tt_425_175",
  "T2tt_425_200",
  "T2tt_425_225",
  "T2tt_425_250",
  "T2tt_425_275",
  "T2tt_425_300",
  "T2tt_433_250",
  "T2tt_442_275",
  "T2tt_450_1",
  "T2tt_450_50",
  "T2tt_450_100",
  "T2tt_450_150",
  "T2tt_450_175",
  "T2tt_450_200",
  "T2tt_450_225",
  "T2tt_450_250",
  "T2tt_450_275",
  "T2tt_450_300",
  "T2tt_450_325",
  "T2tt_450_363",
  "T2tt_458_275",
  "T2tt_475_175",
  "T2tt_475_200",
  "T2tt_475_225",
  "T2tt_475_250",
  "T2tt_475_275",
  "T2tt_475_300",
  "T2tt_475_325",
  "T2tt_475_350",
  "T2tt_475_375",
  "T2tt_492_325",
  "T2tt_500_1",
  "T2tt_500_50",
  "T2tt_500_100",
  "T2tt_500_150",
  "T2tt_500_200",
  "T2tt_500_225",
  "T2tt_500_250",
  "T2tt_500_275",
  "T2tt_500_300",
  "T2tt_500_325",
  "T2tt_500_350",
  "T2tt_500_375",
  "T2tt_500_400",
  "T2tt_500_413",
  "T2tt_517_350",
  "T2tt_525_225",
  "T2tt_525_250",
  "T2tt_525_275",
  "T2tt_525_300",
  "T2tt_525_325",
  "T2tt_525_350",
  "T2tt_525_375",
  "T2tt_525_400",
  "T2tt_525_425",
  "T2tt_525_438",
  "T2tt_533_350",
  "T2tt_542_375",
  "T2tt_550_1",
  "T2tt_550_50",
  "T2tt_550_100",
  "T2tt_550_150",
  "T2tt_550_200",
  "T2tt_550_250",
  "T2tt_550_275",
  "T2tt_550_300",
  "T2tt_550_325",
  "T2tt_550_350",
  "T2tt_550_375",
  "T2tt_550_400",
  "T2tt_550_425",
  "T2tt_550_450",
  "T2tt_558_375",
  "T2tt_567_400",
  "T2tt_575_275",
  "T2tt_575_300",
  "T2tt_575_325",
  "T2tt_575_350",
  "T2tt_575_375",
  "T2tt_575_400",
  "T2tt_575_425",
  "T2tt_575_450",
  "T2tt_575_475",
  "T2tt_575_488",
  "T2tt_583_400",
  "T2tt_592_425",
  "T2tt_600_1",
  "T2tt_600_50",
  "T2tt_600_100",
  "T2tt_600_150",
  "T2tt_600_200",
  "T2tt_600_250",
  "T2tt_600_300",
  "T2tt_600_325",
  "T2tt_600_350",
  "T2tt_600_375",
  "T2tt_600_400",
  "T2tt_600_425",
  "T2tt_600_450",
  "T2tt_600_475",
  "T2tt_600_500",
  "T2tt_600_513",
  "T2tt_608_425",
  "T2tt_617_450",
  "T2tt_625_325",
  "T2tt_625_350",
  "T2tt_625_375",
  "T2tt_625_400",
  "T2tt_625_425",
  "T2tt_625_450",
  "T2tt_625_475",
  "T2tt_625_500",
  "T2tt_625_525",
  "T2tt_625_538",
  "T2tt_633_450",
  "T2tt_642_475",
  "T2tt_650_1",
  "T2tt_650_50",
  "T2tt_650_100",
  "T2tt_650_150",
  "T2tt_650_200",
  "T2tt_650_250",
  "T2tt_650_300",
  "T2tt_650_350",
  "T2tt_650_375",
  "T2tt_650_400",
  "T2tt_650_425",
  "T2tt_650_450",
  "T2tt_650_475",
  "T2tt_650_500",
  "T2tt_650_525",
  "T2tt_650_550",
  "T2tt_650_563",
  "T2tt_658_475",
  "T2tt_667_500",
  "T2tt_675_375",
  "T2tt_675_400",
  "T2tt_675_425",
  "T2tt_675_450",
  "T2tt_675_475",
  "T2tt_675_500",
  "T2tt_675_525",
  "T2tt_675_550",
  "T2tt_675_575",
  "T2tt_675_588",
  "T2tt_683_500",
  "T2tt_692_525",
  "T2tt_700_1",
  "T2tt_700_50",
  "T2tt_700_100",
  "T2tt_700_150",
  "T2tt_700_200",
  "T2tt_700_250",
  "T2tt_700_300",
  "T2tt_700_350",
  "T2tt_700_400",
  "T2tt_700_425",
  "T2tt_700_450",
  "T2tt_700_475",
  "T2tt_700_500",
  "T2tt_700_525",
  "T2tt_700_550",
  "T2tt_700_575",
  "T2tt_700_600",
  "T2tt_700_613",
  "T2tt_708_525",
  "T2tt_717_550",
  "T2tt_725_425",
  "T2tt_725_450",
  "T2tt_725_475",
  "T2tt_725_500",
  "T2tt_725_525",
  "T2tt_725_550",
  "T2tt_725_575",
  "T2tt_725_600",
  "T2tt_725_625",
  "T2tt_725_638",
  "T2tt_733_550",
  "T2tt_742_575",
  "T2tt_750_1",
  "T2tt_750_50",
  "T2tt_750_100",
  "T2tt_750_150",
  "T2tt_750_200",
  "T2tt_750_250",
  "T2tt_750_300",
  "T2tt_750_350",
  "T2tt_750_400",
  "T2tt_750_450",
  "T2tt_750_475",
  "T2tt_750_500",
  "T2tt_750_525",
  "T2tt_750_550",
  "T2tt_750_575",
  "T2tt_750_600",
  "T2tt_750_625",
  "T2tt_750_650",
  "T2tt_758_575",
  "T2tt_767_600",
  "T2tt_775_475",
  "T2tt_775_500",
  "T2tt_775_525",
  "T2tt_775_550",
  "T2tt_775_575",
  "T2tt_775_600",
  "T2tt_775_625",
  "T2tt_775_650",
  "T2tt_783_600",
  "T2tt_792_625",
  "T2tt_800_1",
  "T2tt_800_50",
  "T2tt_800_100",
  "T2tt_800_150",
  "T2tt_800_200",
  "T2tt_800_250",
  "T2tt_800_300",
  "T2tt_800_350",
  "T2tt_800_400",
  "T2tt_800_450",
  "T2tt_800_500",
  "T2tt_800_525",
  "T2tt_800_550",
  "T2tt_800_575",
  "T2tt_800_600",
  "T2tt_800_625",
  "T2tt_800_650",
  "T2tt_808_625",
  "T2tt_817_650",
  "T2tt_825_525",
  "T2tt_825_550",
  "T2tt_825_575",
  "T2tt_825_600",
  "T2tt_825_625",
  "T2tt_825_650",
  "T2tt_833_650",
  "T2tt_850_1",
  "T2tt_850_50",
  "T2tt_850_100",
  "T2tt_850_150",
  "T2tt_850_200",
  "T2tt_850_250",
  "T2tt_850_300",
  "T2tt_850_350",
  "T2tt_850_400",
  "T2tt_850_450",
  "T2tt_850_500",
  "T2tt_850_550",
  "T2tt_850_575",
  "T2tt_850_600",
  "T2tt_850_625",
  "T2tt_850_650",
  "T2tt_875_575",
  "T2tt_875_600",
  "T2tt_875_625",
  "T2tt_875_650",
  "T2tt_900_1",
  "T2tt_900_50",
  "T2tt_900_100",
  "T2tt_900_150",
  "T2tt_900_200",
  "T2tt_900_250",
  "T2tt_900_300",
  "T2tt_900_350",
  "T2tt_900_400",
  "T2tt_900_450",
  "T2tt_900_500",
  "T2tt_900_550",
  "T2tt_900_600",
  "T2tt_900_625",
  "T2tt_900_650",
  "T2tt_925_625",
  "T2tt_925_650",
  "T2tt_950_1",
  "T2tt_950_50",
  "T2tt_950_100",
  "T2tt_950_150",
  "T2tt_950_200",
  "T2tt_950_250",
  "T2tt_950_300",
  "T2tt_950_350",
  "T2tt_950_400",
  "T2tt_950_450",
  "T2tt_950_500",
  "T2tt_950_550",
  "T2tt_950_600",
  "T2tt_950_650",
  "T2tt_1000_1",
  "T2tt_1000_50",
  "T2tt_1000_100",
  "T2tt_1000_150",
  "T2tt_1000_200",
  "T2tt_1000_250",
  "T2tt_1000_300",
  "T2tt_1000_350",
  "T2tt_1000_400",
  "T2tt_1000_450",
  "T2tt_1000_500",
  "T2tt_1000_550",
  "T2tt_1000_600",
  "T2tt_1000_650",
  "T2tt_1050_1",
  "T2tt_1050_50",
  "T2tt_1050_100",
  "T2tt_1050_150",
  "T2tt_1050_200",
  "T2tt_1050_250",
  "T2tt_1050_300",
  "T2tt_1050_350",
  "T2tt_1050_400",
  "T2tt_1050_450",
  "T2tt_1050_500",
  "T2tt_1050_550",
  "T2tt_1050_600",
  "T2tt_1050_650",
  "T2tt_1100_1",
  "T2tt_1100_50",
  "T2tt_1100_100",
  "T2tt_1100_150",
  "T2tt_1100_200",
  "T2tt_1100_250",
  "T2tt_1100_300",
  "T2tt_1100_350",
  "T2tt_1100_400",
  "T2tt_1100_450",
  "T2tt_1100_500",
  "T2tt_1100_550",
  "T2tt_1100_600",
  "T2tt_1100_650",
  "T2tt_1150_1",
  "T2tt_1150_50",
  "T2tt_1150_100",
  "T2tt_1150_150",
  "T2tt_1150_200",
  "T2tt_1150_250",
  "T2tt_1150_300",
  "T2tt_1150_350",
  "T2tt_1150_400",
  "T2tt_1150_450",
  "T2tt_1150_500",
  "T2tt_1150_550",
  "T2tt_1150_600",
  "T2tt_1150_650",
  "T2tt_1200_1",
  "T2tt_1200_50",
  "T2tt_1200_100",
  "T2tt_1200_150",
  "T2tt_1200_200",
  "T2tt_1200_250",
  "T2tt_1200_300",
  "T2tt_1200_350",
  "T2tt_1200_400",
  "T2tt_1200_450",
  "T2tt_1200_500",
  "T2tt_1200_550",
  "T2tt_1200_600",
  "T2tt_1200_650",
  "T2fbd_250_170",
  "T2fbd_250_180",
  "T2fbd_250_190",
  "T2fbd_250_200",
  "T2fbd_250_210",
  "T2fbd_250_220",
  "T2fbd_250_230",
  "T2fbd_250_240",
  "T2fbd_275_195",
  "T2fbd_275_205",
  "T2fbd_275_215",
  "T2fbd_275_225",
  "T2fbd_275_235",
  "T2fbd_275_245",
  "T2fbd_275_255",
  "T2fbd_275_265",
  "T2fbd_300_220",
  "T2fbd_300_230",
  "T2fbd_300_240",
  "T2fbd_300_250",
  "T2fbd_300_260",
  "T2fbd_300_270",
  "T2fbd_300_280",
  "T2fbd_300_290",
  "T2fbd_325_245",
  "T2fbd_325_255",
  "T2fbd_325_265",
  "T2fbd_325_275",
  "T2fbd_325_285",
  "T2fbd_325_295",
  "T2fbd_325_305",
  "T2fbd_325_315",
  "T2fbd_350_270",
  "T2fbd_350_280",
  "T2fbd_350_290",
  "T2fbd_350_300",
  "T2fbd_350_310",
  "T2fbd_350_320",
  "T2fbd_350_330",
  "T2fbd_350_340",
  "T2fbd_375_295",
  "T2fbd_375_305",
  "T2fbd_375_315",
  "T2fbd_375_325",
  "T2fbd_375_335",
  "T2fbd_375_345",
  "T2fbd_375_355",
  "T2fbd_375_365",
  "T2fbd_400_320",
  "T2fbd_400_330",
  "T2fbd_400_340",
  "T2fbd_400_350",
  "T2fbd_400_360",
  "T2fbd_400_370",
  "T2fbd_400_380",
  "T2fbd_400_390",
  "T2fbd_425_345",
  "T2fbd_425_355",
  "T2fbd_425_365",
  "T2fbd_425_375",
  "T2fbd_425_385",
  "T2fbd_425_395",
  "T2fbd_425_405",
  "T2fbd_425_415",
  "T2fbd_450_370",
  "T2fbd_450_380",
  "T2fbd_450_390",
  "T2fbd_450_400",
  "T2fbd_450_410",
  "T2fbd_450_420",
  "T2fbd_450_430",
  "T2fbd_450_440",
  "T2fbd_475_395",
  "T2fbd_475_405",
  "T2fbd_475_415",
  "T2fbd_475_425",
  "T2fbd_475_435",
  "T2fbd_475_445",
  "T2fbd_475_455",
  "T2fbd_475_465",
  "T2fbd_500_420",
  "T2fbd_500_430",
  "T2fbd_500_440",
  "T2fbd_500_450",
  "T2fbd_500_460",
  "T2fbd_500_470",
  "T2fbd_500_480",
  "T2fbd_500_490",
  "T2fbd_525_445",
  "T2fbd_525_455",
  "T2fbd_525_465",
  "T2fbd_525_475",
  "T2fbd_525_485",
  "T2fbd_525_495",
  "T2fbd_525_505",
  "T2fbd_525_515",
  "T2fbd_550_470",
  "T2fbd_550_480",
  "T2fbd_550_490",
  "T2fbd_550_500",
  "T2fbd_550_510",
  "T2fbd_550_520",
  "T2fbd_550_530",
  "T2fbd_550_540",
  "T2fbd_575_495",
  "T2fbd_575_505",
  "T2fbd_575_515",
  "T2fbd_575_525",
  "T2fbd_575_535",
  "T2fbd_575_545",
  "T2fbd_575_555",
  "T2fbd_575_565",
  "T2fbd_600_520",
  "T2fbd_600_530",
  "T2fbd_600_540",
  "T2fbd_600_550",
  "T2fbd_600_560",
  "T2fbd_600_570",
  "T2fbd_600_580",
  "T2fbd_600_590",
  "T2fbd_625_545",
  "T2fbd_625_555",
  "T2fbd_625_565",
  "T2fbd_625_575",
  "T2fbd_625_585",
  "T2fbd_625_595",
  "T2fbd_625_605",
  "T2fbd_625_615",
  "T2fbd_650_570",
  "T2fbd_650_580",
  "T2fbd_650_590",
  "T2fbd_650_600",
  "T2fbd_650_610",
  "T2fbd_650_620",
  "T2fbd_650_630",
  "T2fbd_650_640",
  "T2fbd_675_595",
  "T2fbd_675_605",
  "T2fbd_675_615",
  "T2fbd_675_625",
  "T2fbd_675_635",
  "T2fbd_675_645",
  "T2fbd_675_655",
  "T2fbd_675_665",
  "T2fbd_700_620",
  "T2fbd_700_630",
  "T2fbd_700_640",
  "T2fbd_700_650",
  "T2fbd_700_660",
  "T2fbd_700_670",
  "T2fbd_700_680",
  "T2fbd_700_690",
  "T2fbd_725_645",
  "T2fbd_725_655",
  "T2fbd_725_665",
  "T2fbd_725_675",
  "T2fbd_725_685",
  "T2fbd_725_695",
  "T2fbd_725_705",
  "T2fbd_725_715",
  "T2fbd_750_670",
  "T2fbd_750_680",
  "T2fbd_750_690",
  "T2fbd_750_700",
  "T2fbd_750_710",
  "T2fbd_750_720",
  "T2fbd_750_730",
  "T2fbd_750_740",
  "T2fbd_775_695",
  "T2fbd_775_705",
  "T2fbd_775_715",
  "T2fbd_775_725",
  "T2fbd_775_735",
  "T2fbd_775_745",
  "T2fbd_775_755",
  "T2fbd_775_765",
  "T2fbd_800_720",
  "T2fbd_800_730",
  "T2fbd_800_740",
  "T2fbd_800_750",
  "T2fbd_800_760",
  "T2fbd_800_770",
  "T2fbd_800_780",
  "T2fbd_800_790"
};

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

  // 5-6j
  // nb1
  // 1 tagged
  "hm_nb1_highmtb_nj5to6_ntgeq1_nrt0_nw0",
  "hm_nb1_highmtb_nj5to6_nt0_nrtgeq1_nw0",
  // 1+1
  "hm_nb1_highmtb_nj5to6_ntgeq1_nrt0_nwgeq1",
  "hm_nb1_highmtb_nj5to6_nt0_nrtgeq1_nwgeq1",

  // nb2
  // 1 tagged
  "hm_nb2_highmtb_nj5to6_nt1_nrt0_nw0",
  "hm_nb2_highmtb_nj5to6_nt0_nrt1_nw0",
  "hm_nb2_highmtb_nj5to6_nt0_nrt0_nw1",

  // 1+1
  "hm_nb2_highmtb_nj5to6_nt1_nrt0_nw1",
  "hm_nb2_highmtb_nj5to6_nt0_nrt1_nw1",
  "hm_nb2_highmtb_nj5to6_nt1_nrt1_nw0",

  // 2
  "hm_nb2_highmtb_nj5to6_nt2_nrt0_nw0",
  "hm_nb2_highmtb_nj5to6_nt0_nrt2_nw0",
  "hm_nb2_highmtb_nj5to6_nt0_nrt0_nw2",

  // >=7j
	// nb1
  // 1 tagged
  "hm_nb1_highmtb_nj7_ntgeq1_nrt0_nw0",
  "hm_nb1_highmtb_nj7_nt0_nrtgeq1_nw0",
  // 1+1
  "hm_nb1_highmtb_nj7_ntgeq1_nrt0_nwgeq1",
  "hm_nb1_highmtb_nj7_nt0_nrtgeq1_nwgeq1",

	// nb2
  // 1 tagged
  "hm_nb2_highmtb_nj7_nt1_nrt0_nw0",
  "hm_nb2_highmtb_nj7_nt0_nrt1_nw0",
  "hm_nb2_highmtb_nj7_nt0_nrt0_nw1",

  // 1+1
  "hm_nb2_highmtb_nj7_nt1_nrt0_nw1",
  "hm_nb2_highmtb_nj7_nt0_nrt1_nw1",
  "hm_nb2_highmtb_nj7_nt1_nrt1_nw0",

  // 2
  "hm_nb2_highmtb_nj7_nt2_nrt0_nw0",
  "hm_nb2_highmtb_nj7_nt0_nrt2_nw0",
  "hm_nb2_highmtb_nj7_nt0_nrt0_nw2",
  //---------- high deltaM ----------

};

std::map<TString, TString> srcuts = []{
    std::map<TString, TString> cuts;
    for (auto name : srbins)
      cuts[name] = createCutString(name, cutMap);
    return cuts;
}();

std::map<TString, TString> srlabels = []{
    std::map<TString, TString> cmap;
    for (auto s: srbins) cmap[s] = s;
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
  {"hm_nb1_lowmtb_nj7_nrtgeq1",        {250, 300, 400, 500, 1000}},
  {"hm_nb2_lowmtb_nj7_nrtgeq1",        {250, 300, 400, 500, 1000}},

  // high mtb
  // 0 taggged
  {"hm_nb1_highmtb_nj7_nt0_nrt0_nw0",         {250, 350, 450, 550, 650, 1000}}, // {250, 350, 450, 550, 1000}},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw0",         {250, 350, 450, 550, 650, 1000}}, // {250, 350, 450, 550, 1000}},

  // 5-6j
  // nb1
  {"hm_nb1_highmtb_nj5to6_ntgeq1_nrt0_nw0",                  {550, 650, 1000}},
  {"hm_nb1_highmtb_nj5to6_nt0_nrtgeq1_nw0",   {250, 350, 450, 550, 650, 1000}},

  // 1+1
  {"hm_nb1_highmtb_nj5to6_ntgeq1_nrt0_nwgeq1",               {550, 650, 1000}}, // {550, 1000}},
  {"hm_nb1_highmtb_nj5to6_nt0_nrtgeq1_nwgeq1",{250, 350, 450, 550, 650, 1000}}, // {250, 350, 450, 550, 1000}},

  // nb2
  {"hm_nb2_highmtb_nj5to6_nt1_nrt0_nw0",                     {550, 650, 1000}},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt1_nw0",      {250, 350, 450, 550, 650, 1000}},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt0_nw1",      {250, 350, 450, 550, 650, 1000}},

  // 1+1
  {"hm_nb2_highmtb_nj5to6_nt1_nrt0_nw1",                     {550, 650, 1000}}, // {550, 1000}},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt1_nw1",      {250, 350, 450, 550, 650, 1000}}, // {250, 350, 450, 550, 1000}},
  {"hm_nb2_highmtb_nj5to6_nt1_nrt1_nw0",      {250, 350, 450, 550, 650, 1000}}, // {250, 350, 450, 1000}},

  // 2
  {"hm_nb2_highmtb_nj5to6_nt2_nrt0_nw0",      {250, 350, 450, 550, 650, 1000}}, // {250, 1000}},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt2_nw0",      {250, 350, 450, 550, 650, 1000}}, // {250, 1000}},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt0_nw2",      {250, 350, 450, 550, 650, 1000}}, // {250, 1000}},

  // >=7j
  // nb1
  {"hm_nb1_highmtb_nj7_ntgeq1_nrt0_nw0",                  {550, 650, 1000}},
  {"hm_nb1_highmtb_nj7_nt0_nrtgeq1_nw0",   {250, 350, 450, 550, 650, 1000}},

  // 1+1
  {"hm_nb1_highmtb_nj7_ntgeq1_nrt0_nwgeq1",               {550, 650, 1000}}, // {550, 1000}},
  {"hm_nb1_highmtb_nj7_nt0_nrtgeq1_nwgeq1",{250, 350, 450, 550, 650, 1000}}, // {250, 350, 450, 550, 1000}},

  // nb2
  {"hm_nb2_highmtb_nj7_nt1_nrt0_nw0",                     {550, 650, 1000}},
  {"hm_nb2_highmtb_nj7_nt0_nrt1_nw0",      {250, 350, 450, 550, 650, 1000}},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw1",      {250, 350, 450, 550, 650, 1000}},

  // 1+1
  {"hm_nb2_highmtb_nj7_nt1_nrt0_nw1",                     {550, 650, 1000}}, // {550, 1000}},
  {"hm_nb2_highmtb_nj7_nt0_nrt1_nw1",      {250, 350, 450, 550, 650, 1000}}, // {250, 350, 450, 550, 1000}},
  {"hm_nb2_highmtb_nj7_nt1_nrt1_nw0",      {250, 350, 450, 550, 650, 1000}}, // {250, 350, 450, 1000}},

  // 2
  {"hm_nb2_highmtb_nj7_nt2_nrt0_nw0",      {250, 350, 450, 550, 650, 1000}}, // {250, 1000}},
  {"hm_nb2_highmtb_nj7_nt0_nrt2_nw0",      {250, 350, 450, 550, 650, 1000}}, // {250, 1000}},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw2",      {250, 350, 450, 550, 650, 1000}}, // {250, 1000}},
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

  // 5-6j
  // nb1
  // 1
  {"hm_nb1_highmtb_nj5to6_ntgeq1_nrt0_nw0",   "hm_nb1_highmtb_nj5to6"},
  {"hm_nb1_highmtb_nj5to6_nt0_nrtgeq1_nw0",   "hm_nb1_highmtb_nj5to6"},
  // 1+1
  {"hm_nb1_highmtb_nj5to6_ntgeq1_nrt0_nwgeq1",  "hm_nb1_highmtb_nj5to6"},
  {"hm_nb1_highmtb_nj5to6_nt0_nrtgeq1_nwgeq1",  "hm_nb1_highmtb_nj5to6"},

  // nb2
  //1
  {"hm_nb2_highmtb_nj5to6_nt1_nrt0_nw0",      "hm_nb2_highmtb_nj5to6"},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt1_nw0",      "hm_nb2_highmtb_nj5to6"},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt0_nw1",      "hm_nb2_highmtb_nj5to6"},

  // 1+1
  {"hm_nb2_highmtb_nj5to6_nt1_nrt0_nw1",        "hm_nb2_highmtb_nj5to6"},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt1_nw1",        "hm_nb2_highmtb_nj5to6"},
  {"hm_nb2_highmtb_nj5to6_nt1_nrt1_nw0",        "hm_nb2_highmtb_nj5to6"},

  // 2
  {"hm_nb2_highmtb_nj5to6_nt2_nrt0_nw0",        "hm_nb2_highmtb_nj5to6"},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt2_nw0",        "hm_nb2_highmtb_nj5to6"},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt0_nw2",        "hm_nb2_highmtb_nj5to6"},


  // >=7j
  // nb1
  // 1
  {"hm_nb1_highmtb_nj7_ntgeq1_nrt0_nw0",   "hm_nb1_highmtb_nj7"},
  {"hm_nb1_highmtb_nj7_nt0_nrtgeq1_nw0",   "hm_nb1_highmtb_nj7"},
  // 1+1
  {"hm_nb1_highmtb_nj7_ntgeq1_nrt0_nwgeq1",  "hm_nb1_highmtb_nj7"},
  {"hm_nb1_highmtb_nj7_nt0_nrtgeq1_nwgeq1",  "hm_nb1_highmtb_nj7"},

  // nb2
  //1
  {"hm_nb2_highmtb_nj7_nt1_nrt0_nw0",      "hm_nb2_highmtb_nj7"},
  {"hm_nb2_highmtb_nj7_nt0_nrt1_nw0",      "hm_nb2_highmtb_nj7"},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw1",      "hm_nb2_highmtb_nj7"},

  // 1+1
  {"hm_nb2_highmtb_nj7_nt1_nrt0_nw1",        "hm_nb2_highmtb_nj7"},
  {"hm_nb2_highmtb_nj7_nt0_nrt1_nw1",        "hm_nb2_highmtb_nj7"},
  {"hm_nb2_highmtb_nj7_nt1_nrt1_nw0",        "hm_nb2_highmtb_nj7"},

  // 2
  {"hm_nb2_highmtb_nj7_nt2_nrt0_nw0",        "hm_nb2_highmtb_nj7"},
  {"hm_nb2_highmtb_nj7_nt0_nrt2_nw0",        "hm_nb2_highmtb_nj7"},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw2",        "hm_nb2_highmtb_nj7"},
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
  {"hm_nb1_lowmtb_nj7_nrtgeq1",          "hm_nb1_lowmtb_nj7"},
  {"hm_nb2_lowmtb_nj7_nrtgeq1",          "hm_nb2_lowmtb_nj7"},

  // high mtb
  // 0 taggged
  {"hm_nb1_highmtb_nj7_nt0_nrt0_nw0",    "hm_nb1_highmtb_nj7"},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw0",    "hm_nb2_highmtb_nj7"},

  // 5-6j
  // nb1
  // 1
  {"hm_nb1_highmtb_nj5to6_ntgeq1_nrt0_nw0",     "hm_nb1_highmtb_nj5to6"},
  {"hm_nb1_highmtb_nj5to6_nt0_nrtgeq1_nw0",     "hm_nb1_highmtb_nj5to6"},
  // 1+1
  {"hm_nb1_highmtb_nj5to6_ntgeq1_nrt0_nwgeq1",  "hm_nb1_highmtb_nj5to6"},
  {"hm_nb1_highmtb_nj5to6_nt0_nrtgeq1_nwgeq1",  "hm_nb1_highmtb_nj5to6"},

  // nb2
  // 1
  {"hm_nb2_highmtb_nj5to6_nt1_nrt0_nw0",        "hm_nb2_highmtb_nj5to6"},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt1_nw0",        "hm_nb2_highmtb_nj5to6"},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt0_nw1",        "hm_nb2_highmtb_nj5to6"},
  // 1+1
  {"hm_nb2_highmtb_nj5to6_nt1_nrt0_nw1",        "hm_nb2_highmtb_nj5to6"},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt1_nw1",        "hm_nb2_highmtb_nj5to6"},
  {"hm_nb2_highmtb_nj5to6_nt1_nrt1_nw0",        "hm_nb2_highmtb_nj5to6"},

  // 2
  {"hm_nb2_highmtb_nj5to6_nt2_nrt0_nw0",        "hm_nb2_highmtb_nj5to6"},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt2_nw0",        "hm_nb2_highmtb_nj5to6"},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt0_nw2",        "hm_nb2_highmtb_nj5to6"},

  // >=7j
  // nb1
  // 1
  {"hm_nb1_highmtb_nj7_ntgeq1_nrt0_nw0",     "hm_nb1_highmtb_nj7"},
  {"hm_nb1_highmtb_nj7_nt0_nrtgeq1_nw0",     "hm_nb1_highmtb_nj7"},
  // 1+1
  {"hm_nb1_highmtb_nj7_ntgeq1_nrt0_nwgeq1",  "hm_nb1_highmtb_nj7"},
  {"hm_nb1_highmtb_nj7_nt0_nrtgeq1_nwgeq1",  "hm_nb1_highmtb_nj7"},

  // nb2
  // 1
  {"hm_nb2_highmtb_nj7_nt1_nrt0_nw0",        "hm_nb2_highmtb_nj7"},
  {"hm_nb2_highmtb_nj7_nt0_nrt1_nw0",        "hm_nb2_highmtb_nj7"},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw1",        "hm_nb2_highmtb_nj7"},
  // 1+1
  {"hm_nb2_highmtb_nj7_nt1_nrt0_nw1",        "hm_nb2_highmtb_nj7"},
  {"hm_nb2_highmtb_nj7_nt0_nrt1_nw1",        "hm_nb2_highmtb_nj7"},
  {"hm_nb2_highmtb_nj7_nt1_nrt1_nw0",        "hm_nb2_highmtb_nj7"},

  // 2
  {"hm_nb2_highmtb_nj7_nt2_nrt0_nw0",        "hm_nb2_highmtb_nj7"},
  {"hm_nb2_highmtb_nj7_nt0_nrt2_nw0",        "hm_nb2_highmtb_nj7"},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw2",        "hm_nb2_highmtb_nj7"},
  //---------- high deltaM ----------

};


std::map<TString, TString> lepcrCuts = []{
    std::map<TString, TString> cuts;
    for (auto sr2cr : lepcrMapping)
      cuts[sr2cr.first] = createCutString(sr2cr.second, cutMap);
    return cuts;
}();

std::map<TString, TString> lepcrlabels = lepcrMapping;
std::map<TString, std::vector<int>> lepcrMETbins = srMETbins;

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
  {"hm_nb1_highmtb_nj7_nt0_nrt0_nw0",         {250, 350, 450, 550, 650, 1000}}, // {250, 350, 450, 550, 1000}},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw0",         {250, 350, 450, 550, 650, 1000}}, // {250, 350, 450, 550, 1000}},

  // 5-6j
  // nb1
  {"hm_nb1_highmtb_nj5to6_ntgeq1_nrt0_nw0",                  {550, 650, 1000}},
  {"hm_nb1_highmtb_nj5to6_nt0_nrtgeq1_nw0",   {250, 350, 450, 550, 650, 1000}},

  // 1+1
  {"hm_nb1_highmtb_nj5to6_ntgeq1_nrt0_nwgeq1",               {550, 650, 1000}}, // {550, 1000}},
  {"hm_nb1_highmtb_nj5to6_nt0_nrtgeq1_nwgeq1",{250, 350, 450, 550, 650, 1000}}, // {250, 350, 450, 550, 1000}},

  // nb2
  {"hm_nb2_highmtb_nj5to6_nt1_nrt0_nw0",                     {550, 650, 1000}},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt1_nw0",      {250, 350, 450, 550, 650, 1000}},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt0_nw1",      {250, 350, 450, 550, 650, 1000}},

  // 1+1
  {"hm_nb2_highmtb_nj5to6_nt1_nrt0_nw1",                     {550, 650, 1000}}, // {550, 1000}},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt1_nw1",      {250, 350, 450, 550, 650, 1000}}, // {250, 350, 450, 550, 1000}},
  {"hm_nb2_highmtb_nj5to6_nt1_nrt1_nw0",      {250, 350, 450, 550, 650, 1000}}, // {250, 350, 450, 1000}},

  // 2
  {"hm_nb2_highmtb_nj5to6_nt2_nrt0_nw0",      {250, 350, 450, 550, 650, 1000}}, // {250, 1000}},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt2_nw0",      {250, 350, 450, 550, 650, 1000}}, // {250, 1000}},
  {"hm_nb2_highmtb_nj5to6_nt0_nrt0_nw2",      {250, 350, 450, 550, 650, 1000}}, // {250, 1000}},

  // >=7j
  // nb1
  {"hm_nb1_highmtb_nj7_ntgeq1_nrt0_nw0",                  {550, 650, 1000}},
  {"hm_nb1_highmtb_nj7_nt0_nrtgeq1_nw0",   {250, 350, 450, 550, 650, 1000}},

  // 1+1
  {"hm_nb1_highmtb_nj7_ntgeq1_nrt0_nwgeq1",               {550, 650, 1000}}, // {550, 1000}},
  {"hm_nb1_highmtb_nj7_nt0_nrtgeq1_nwgeq1",{250, 350, 450, 550, 650, 1000}}, // {250, 350, 450, 550, 1000}},

  // nb2
  {"hm_nb2_highmtb_nj7_nt1_nrt0_nw0",                     {550, 650, 1000}},
  {"hm_nb2_highmtb_nj7_nt0_nrt1_nw0",      {250, 350, 450, 550, 650, 1000}},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw1",      {250, 350, 450, 550, 650, 1000}},

  // 1+1
  {"hm_nb2_highmtb_nj7_nt1_nrt0_nw1",                     {550, 650, 1000}}, // {550, 1000}},
  {"hm_nb2_highmtb_nj7_nt0_nrt1_nw1",      {250, 350, 450, 550, 650, 1000}}, // {250, 350, 450, 550, 1000}},
  {"hm_nb2_highmtb_nj7_nt1_nrt1_nw0",      {250, 350, 450, 550, 650, 1000}}, // {250, 350, 450, 1000}},

  // 2
  {"hm_nb2_highmtb_nj7_nt2_nrt0_nw0",      {250, 350, 450, 550, 650, 1000}}, // {250, 1000}},
  {"hm_nb2_highmtb_nj7_nt0_nrt2_nw0",      {250, 350, 450, 550, 650, 1000}}, // {250, 1000}},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw2",      {250, 350, 450, 550, 650, 1000}}, // {250, 1000}},
  //---------- high deltaM ----------
};


map<TString, Category> srCatMap(){
  map<TString, Category> cmap;
  for (auto &name : srbins){
    cmap[name] = Category(name, srcuts.at(name), srlabels.at(name), BinInfo("met", "#slash{E}_{T}", srMETbins.at(name), "GeV"));
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
  const auto &cuts = ICHEPCR ? srcuts: lepcrCuts;
  const auto &labels = ICHEPCR ? srlabels: lepcrlabels;
  for (auto &name : srbins){
    cmap[name] = Category(name, cuts.at(name), labels.at(name), BinInfo("met", varlabel, lepcrMETbins.at(name), "GeV"));
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
  config.crMapping = phocrMapping;

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

  if (ADD_LEP_TO_MET){
    config.addSample("singlelep",   "Data",          datadir+"/lepcr/singlelep",       "1.0",     datasel + trigLepCR + lepcrsel);
    config.addSample("ttbar",       "t#bar{t}",      "lepcr/ttbar",           onelepcrwgt, datasel + trigLepCR + lepcrsel);
    config.addSample("wjets",       "W+jets",        "lepcr/wjets",           onelepcrwgt, datasel + trigLepCR + lepcrsel);
    config.addSample("tW",          "tW",            "lepcr/tW",              onelepcrwgt, datasel + trigLepCR + lepcrsel);
    config.addSample("ttW",         "ttW",           "lepcr/ttW",             onelepcrwgt, datasel + trigLepCR + lepcrsel);
//    config.addSample("qcd",         "QCD",           "lepcr/qcd",             onelepcrwgt, datasel + trigLepCR + lepcrsel);
  }else{
    config.addSample("singlelep",   "Data",          datadir+"/sr/met",    "1.0",     datasel + trigSR + revert_vetoes);
    config.addSample("ttbar",       "t#bar{t}",      "sr/ttbar",           lepselwgt, datasel + trigSR + revert_vetoes);
    config.addSample("wjets",       "W+jets",        "sr/wjets",           lepselwgt, datasel + trigSR + revert_vetoes);
    config.addSample("tW",          "tW",            "sr/tW",              lepselwgt, datasel + trigSR + revert_vetoes);
    config.addSample("ttW",         "ttW",           "sr/ttW",             lepselwgt, datasel + trigSR + revert_vetoes);
//    config.addSample("qcd",         "QCD",           "sr/qcd",             lepselwgt, datasel + trigSR + revert_vetoes);

    // add signals in lepcr
    for (const auto &sig : signals){
      config.addSample("lepcr_"+sig,   sig,          "signals/"+sig,    siglepselwgt, datasel + trigSR + revert_vetoes);
    }
  }

  config.addSample("ttbar-sr",       "t#bar{t}",      "sr/ttbar",           lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("wjets-sr",       "W+jets",        "sr/wjets",           lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("tW-sr",          "tW",            "sr/tW",              lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("ttW-sr",         "ttW",           "sr/ttW",             lepvetowgt, datasel + trigSR + vetoes);
//  config.addSample("qcd-sr",         "QCD",           "qcd",                lepvetowgt, datasel + trigSR + vetoes);
//  config.addSample("rare-sr",        "Rare",          "sr/rare",            lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("ttZ-sr",         "ttZ",           "sr/ttZ",             lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("diboson-sr",     "Diboson",       "sr/diboson",         lepvetowgt, datasel + trigSR + vetoes);

  config.sel = baseline;
  config.categories = srbins;
  config.catMaps = srCatMap();
  config.crCatMaps = lepCatMap();
  config.crMapping = lepcrMapping;

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

  config.sel = baseline;
  config.categories = srbins;
  config.catMaps = srCatMap();    // sr DPhi cut applied
  config.crCatMaps = qcdCatMap(); // no DPhi cut in category def
  config.crMapping = qcdcrMapping;

  return config;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// --------------
std::vector<TString> mergedSRbins{
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
  "hm_nb1_highmtb_ntgeq1_nrt0_nw0",
  "hm_nb1_highmtb_nt0_nrtgeq1_nw0",
  // 1+1
  "hm_nb1_highmtb_ntgeq1_nrt0_nwgeq1",
  "hm_nb1_highmtb_nt0_nrtgeq1_nwgeq1",

  // nb2
  // 1 tagged
  "hm_nb2_highmtb_nt1_nrt0_nw0",
  "hm_nb2_highmtb_nt0_nrt1_nw0",
  "hm_nb2_highmtb_nt0_nrt0_nw1",

  // 1+1
  "hm_nb2_highmtb_nt1_nrt0_nw1",
  "hm_nb2_highmtb_nt0_nrt1_nw1",
  "hm_nb2_highmtb_nt1_nrt1_nw0",

  // 2
  "hm_nb2_highmtb_nt2_nrt0_nw0",
  "hm_nb2_highmtb_nt0_nrt2_nw0",
  "hm_nb2_highmtb_nt0_nrt0_nw2",
  //---------- high deltaM ----------
};

std::map<TString, TString> mergedSRcuts = []{
    std::map<TString, TString> cuts;
    for (auto name : mergedSRbins)
      cuts[name] = createCutString(name, cutMap);
    return cuts;
}();

std::map<TString, TString> mergedSRlabels = []{
    std::map<TString, TString> cmap;
    for (auto s: mergedSRbins) cmap[s] = s;
    return cmap;
}();

std::map<TString, std::vector<int>> mergedSRMETbins{
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
  {"hm_nb1_lowmtb_nj7_nrtgeq1",        {250, 300, 400, 500, 1000}},
  {"hm_nb2_lowmtb_nj7_nrtgeq1",        {250, 300, 400, 500, 1000}},

  // high mtb
  // 0 taggged
  {"hm_nb1_highmtb_nj7_nt0_nrt0_nw0",  {250, 350, 450, 550, 1000}},
  {"hm_nb2_highmtb_nj7_nt0_nrt0_nw0",  {250, 350, 450, 550, 1000}},

  // nb1
  {"hm_nb1_highmtb_ntgeq1_nrt0_nw0",                  {550, 650, 1000}},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nw0",   {250, 350, 450, 550, 650, 1000}},

  // 1+1
  {"hm_nb1_highmtb_ntgeq1_nrt0_nwgeq1",               {550, 1000}},
  {"hm_nb1_highmtb_nt0_nrtgeq1_nwgeq1",{250, 350, 450, 550, 1000}},

  // nb2
  {"hm_nb2_highmtb_nt1_nrt0_nw0",                     {550, 650, 1000}},
  {"hm_nb2_highmtb_nt0_nrt1_nw0",      {250, 350, 450, 550, 650, 1000}},
  {"hm_nb2_highmtb_nt0_nrt0_nw1",      {250, 350, 450, 550, 650, 1000}},

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

map<TString, Category> mergedSRCatMap(){
  map<TString, Category> cmap;
  for (auto &name : mergedSRbins){
    cmap[name] = Category(name, mergedSRcuts.at(name), mergedSRlabels.at(name), BinInfo("met", "#slash{E}_{T}", mergedSRMETbins.at(name), "GeV"));
  }
  return cmap;
}

BaseConfig srConfig(){
  BaseConfig     config;

  config.inputdir = inputdir;
  config.outputdir = outputdir+"/sig";
  config.header = "#sqrt{s} = 13 TeV, "+lumistr+" fb^{-1}";

  config.addSample("data",        "Data",          datadir+"/sr/met",    "1.0",      datasel + trigSR + vetoes);
  config.addSample("ttZ",         "ttZ",           "sr/ttZ",             lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("diboson",     "Diboson",       "sr/diboson",         lepvetowgt, datasel + trigSR + vetoes);

  config.sel = baseline;
  config.categories = mergedSRbins;
  config.catMaps = mergedSRCatMap();

  return config;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BaseConfig signalConfig(TString signal){
  BaseConfig     config;

  config.inputdir = inputdir;
  config.outputdir = outputdir+"/sig";

  // signals in sr
  config.addSample(signal,          signal,           "signals/"+signal,    sigwgt,      datasel + trigSR + vetoes);

  config.sel = baseline;
  config.categories = mergedSRbins;
  config.catMaps = mergedSRCatMap();

  return config;
}

BaseConfig lepcrSignalConfig(TString signal){
  BaseConfig     config;

  config.inputdir = inputdir;
  config.outputdir = outputdir+"/LLB";

  // add signals in lepcr
  config.addSample("lepcr_"+signal,   signal,          "signals/"+signal,   siglepselwgt, datasel + trigSR + revert_vetoes);

  config.sel = baseline;
  config.categories = srbins;
  config.catMaps = srCatMap();
  config.crCatMaps = lepCatMap();
  config.crMapping = lepcrMapping;

  return config;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
map<std::string, std::string> makeBinMap(TString control_region){
  map<std::string, vector<TString>> results; // srbinname_met -> [(sr_sub1, cr_sub1), ...]

  const auto &crMapping = control_region=="phocr" ? phocrMapping : (control_region=="lepcr" ? lepcrMapping : qcdcrMapping);

  const auto &merged_srCatMap = mergedSRCatMap();
  const auto &split_srCatMap = srCatMap();

  for (const auto &merged_cat_name : mergedSRbins){
    vector<TString> categories_to_process; // get the categories to consider
    if (merged_cat_name.BeginsWith("lm") || merged_cat_name.Contains("lowmtb")) continue;
    if (merged_cat_name.Contains("nt0_nrt0_nw0")) categories_to_process.push_back(merged_cat_name);
    else {
      categories_to_process.push_back(TString(merged_cat_name).ReplaceAll("highmtb", "highmtb_nj7"));
      categories_to_process.push_back(TString(merged_cat_name).ReplaceAll("highmtb", "highmtb_nj5to6"));
    }

    // create the sr-to-cr map
    const auto& merged_bin = merged_srCatMap.at(merged_cat_name);
    for (const auto &split_cat_name : categories_to_process){
      // loop over all categories to consider, e.g., 5-6j, >=7j
      const auto &split_bin = split_srCatMap.at(split_cat_name);
      for (unsigned ibin=0; ibin<merged_bin.bin.nbins; ++ibin){
        std::string mergedsr_binname = ("bin_"+merged_cat_name+"_"+merged_bin.bin.binnames.at(ibin)).Data();
        if (merged_bin.bin.plotbins.at(ibin+1) == split_bin.bin.plotbins.at(ibin+1)) {
          if (merged_cat_name.Contains("nt0_nrt0_nw0")) continue;
          else{
            // no splitting in MET: merged in nj
            auto splitsrbinname = "bin_"+split_cat_name+"_"+split_bin.bin.binnames.at(ibin);
            auto crbinname = "bin_"+control_region+"_"+TString(crMapping.at(split_cat_name)).ReplaceAll("NoDPhi_","_")+"_"+split_bin.bin.binnames.at(ibin);
            results[mergedsr_binname]; // touch it: initialize it if not, otherwise should append (5-6j, and >=7j)
            results[mergedsr_binname].push_back("<"+splitsrbinname+">*("+crbinname+")");
          }
        }else{
          // also merge in MET
          for (unsigned icr = ibin; icr < split_bin.bin.nbins; ++icr){
            auto splitsrbinname = "bin_"+split_cat_name+"_"+split_bin.bin.binnames.at(icr);
            auto crbinname = "bin_"+control_region+"_"+TString(crMapping.at(split_cat_name)).ReplaceAll("NoDPhi_","_")+"_"+split_bin.bin.binnames.at(icr);
            results[mergedsr_binname]; // touch it: initialize it if not, otherwise should append (5-6j, and >=7j)
            results[mergedsr_binname].push_back("<"+splitsrbinname+">*("+crbinname+")");
          }
        }
      }
    }
  }

  map<std::string, std::string> rltstr;
  for (const auto &s : results){
    rltstr[s.first] = joinString(s.second, " + ").Data();
    cout << s.first << endl << rltstr[s.first] << endl << endl;
  }
  return rltstr;
}

map<std::string, std::string> lepcrBinMap = makeBinMap("lepcr");
map<std::string, std::string> phocrBinMap = makeBinMap("phocr");
map<std::string, std::string> qcdcrBinMap = makeBinMap("qcdcr");
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
#endif /* ESTTOOLS_LMPARAMETERS_HH_ */
