#ifndef ESTTOOLS_HMPARAMETERS_HH_
#define ESTTOOLS_HMPARAMETERS_HH_

#include "../utils/EstHelper.hh"

namespace EstTools{

const TString datadir = ".";
const TString lumistr = "36.2";

TString getLumi(){return lumistr(TRegexp("[0-9]+.[0-9]"));}

// lumi and base weight
const TString wgtvar = lumistr+"*weight*topptWeight";
//const TString wgtvar = lumistr+"*weight*truePUWeight*btagWeight";

// photon trigger eff.
const TString phowgt = wgtvar;
//const TString phowgt = wgtvar + "*qcdRespTailWeight";

// No Lepton SF
const TString lepvetowgt = wgtvar;
const TString lepselwgt  = wgtvar;
const TString vetoes = " && nvetolep==0 && nvetotau==0";

// Tag-and-Probe Lepton SF
//const TString lepvetowgt = wgtvar + "*leptnpweightHM*lepvetoweightHM";
//const TString lepselwgt  = wgtvar + "*leptnpweightHM";
//const TString vetoes = " && nvetolep==0 && (nvetotau==0 || (ismc && npromptgentau>0))";

// 1LCR Lepton SF
//const TString lepvetowgt = wgtvar + "*lepvetoweight";
//const TString lepselwgt  = wgtvar + "*lepselweight";
//const TString vetoes = " && ((nvetolep==0 && nvetotau==0) || (ismc && (ngoodgenele>0 || ngoodgenmu>0 || npromptgentau>0)))";

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
//const TString qcdwgt = wgtvar + "*qcdRespTailWeight";
const TString qcdwgt = wgtvar;
//const TString qcdvetowgt = lepvetowgt + "*qcdRespTailWeight";
const TString qcdvetowgt = lepvetowgt;

// signal weights
//const TString sigwgt = lepvetowgt + "*btagFastSimWeight*isrWeightTight*(1.0*(mtcsv12met<=175)+sdtopFastSimWeight*sdwFastSimWeight*(mtcsv12met>175))";
const TString sigwgt = lepvetowgt;

// triggers
const TString trigSR = " && (passmetmht || ismc)";
const TString trigPhoCR = " && passtrigphoOR && origmet<200";
const TString trigDiLepCR = " && passtrigdilepOR && dileppt>200";
const TString datasel = " && passjson && (passmetfilters || process==10) && j1chEnFrac>0.1 && j1chEnFrac<0.99";
const TString qcdSpikeRemovals = " && (!(run==1 && lumi==46160 && event==331634716)) && (!(run==1 && lumi==91626 && event==208129617))";

// ------------------------------------------------------------------------
// search regions and control regions

const TString baseline = "met>250 && njets>=5 && nbjets>=1 && nlbjets>=2 && dphij1met>0.5 && dphij2met>0.5 && dphij3met>0.5 && dphij4met>0.5";
const TString baseNoDPhi = "met>250 && njets>=5 && nbjets>=1 && nlbjets>=2";
const TString dphi = " && dphij1met>0.5 && dphij2met>0.5 && dphij3met>0.5 && dphij4met>0.5" + qcdSpikeRemovals,
    dphi_invert = " && (dphij1met<0.1 || dphij2met<0.1 || dphij3met < 0.1)" + qcdSpikeRemovals;
const TString baseForNorm = "met>250 && njets>=5 && nlbjets>=2 && nbjets>=1";

const TString nt0 =  "nsdtop==0";
const TString nt1 =  "nsdtop==1";
const TString nt2 =  "nsdtop>=2";
const TString nw0 =  "nsdw==0";
const TString nw1 =  "nsdw==1";
const TString nw2 =  "nsdw>=2";
const TString nrt0 = "nrestop==0";
const TString nrt1 = "nrestop==1";
const TString nrt2 = "nrestop>=2";

std::vector<TString> srbins{
  "nb1_mtb0_nj7_nrt0",
  "nb2_mtb0_nj7_nrt0",
  "nb1_mtb0_nj7_nrt1",
  "nb2_mtb0_nj7_nrt1",

  // 0 taggged
  "nb1_mtb175_nj7_nt0_nrt0_nw0",
  "nb2_mtb175_nj7_nt0_nrt0_nw0",

  // 1 tagged
  "nb1_mtb175_nj5t_nt0_nrt0_nw1",
  "nb2_mtb175_nj5t_nt0_nrt0_nw1",
  "nb1_mtb175_nj5t_nt0_nrt1_nw0",
  "nb2_mtb175_nj5t_nt0_nrt1_nw0",
  "nb1_mtb175_nj5t_nt1_nrt0_nw0",
  "nb2_mtb175_nj5t_nt1_nrt0_nw0",

  // 1+1
  "nb1_mtb175_nj5t_nt1_nrt0_nw1",
  "nb2_mtb175_nj5t_nt1_nrt0_nw1",
  "nb1_mtb175_nj5t_nt0_nrt1_nw1",
  "nb2_mtb175_nj5t_nt0_nrt1_nw1",

  "nb2_mtb175_nj5t_nt1_nrt1_nw0",

  // 2
  "nb2_mtb175_nj5t_nt2_nrt0_nw0",
  "nb2_mtb175_nj5t_nt0_nrt2_nw0",
  "nb2_mtb175_nj5t_nt0_nrt0_nw2",
};

const TString lowmt = "mtcsv12met<=175", highmt = "mtcsv12met>175",
    nb1 = "nbjets==1", nb2 = "nbjets>=2", nbgeq1 = "nbjets>=1",
    mednj = "njets>=5 && njets<7", highnj = "njets>=7", nj5g = "njets>=5",
    nw1g = "1==0",
    nt1g = "1==0";

std::map<TString, TString> normMap{
  {"nb1", nb1},
  {"nb2", nb2}
};

std::map<TString, TString> srcuts{

  {"nb1_mtb0_nj7_nrt0",          addCuts({nrt0, nb1, lowmt, highnj})},
  {"nb2_mtb0_nj7_nrt0",          addCuts({nrt0, nb2, lowmt, highnj})},
  {"nb1_mtb0_nj7_nrt1",          addCuts({nrt1, nb1, lowmt, highnj})},
  {"nb2_mtb0_nj7_nrt1",          addCuts({nrt1, nb2, lowmt, highnj})},


  // 0
  {"nb1_mtb175_nj7_nt0_nrt0_nw0",  addCuts({nt0, nrt0, nw0, nb1, highmt, highnj})},
  {"nb2_mtb175_nj7_nt0_nrt0_nw0",  addCuts({nt0, nrt0, nw0, nb2, highmt, highnj})},

  // 1
  {"nb1_mtb175_nj5t_nt0_nrt0_nw1",  addCuts({nt0, nrt0, nw1, nb1, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt0_nrt0_nw1",  addCuts({nt0, nrt0, nw1, nb2, highmt, nj5g})},
  {"nb1_mtb175_nj5t_nt0_nrt1_nw0",  addCuts({nt0, nrt1, nw0, nb1, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt0_nrt1_nw0",  addCuts({nt0, nrt1, nw0, nb2, highmt, nj5g})},
  {"nb1_mtb175_nj5t_nt1_nrt0_nw0",  addCuts({nt1, nrt0, nw0, nb1, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt1_nrt0_nw0",  addCuts({nt1, nrt0, nw0, nb2, highmt, nj5g})},

  // 1+1
  {"nb1_mtb175_nj5t_nt1_nrt0_nw1",  addCuts({nt1, nrt0, nw1, nb1, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt1_nrt0_nw1",  addCuts({nt1, nrt0, nw1, nb2, highmt, nj5g})},
  {"nb1_mtb175_nj5t_nt0_nrt1_nw1",  addCuts({nt0, nrt1, nw1, nb1, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt0_nrt1_nw1",  addCuts({nt0, nrt1, nw1, nb2, highmt, nj5g})},

  {"nb2_mtb175_nj5t_nt1_nrt1_nw0",  addCuts({nt1, nrt1, nw0, nb2, highmt, nj5g})},


  // 2
  {"nb2_mtb175_nj5t_nt2_nrt0_nw0",  addCuts({nt2, nrt0, nw0, nb2, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt0_nrt2_nw0",  addCuts({nt0, nrt2, nw0, nb2, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt0_nrt0_nw2",  addCuts({nt0, nrt0, nw2, nb2, highmt, nj5g})},
};
std::map<TString, TString> srlabels = srcuts;

std::map<TString, std::vector<int>> srMETbins{

  {"nb1_mtb0_nj7_nrt0",  {250, 300, 400, 500, 1000}},
  {"nb2_mtb0_nj7_nrt0",  {250, 300, 400, 500, 1000}},
  {"nb1_mtb0_nj7_nrt1",  {250, 300, 400, 500, 1000}},
  {"nb2_mtb0_nj7_nrt1",  {250, 300, 400, 500, 1000}},

  // 0
  {"nb1_mtb175_nj7_nt0_nrt0_nw0",  {250, 350, 450, 550, 1000}},
  {"nb2_mtb175_nj7_nt0_nrt0_nw0",  {250, 350, 450, 550, 1000}},

  // 1
  {"nb1_mtb175_nj5t_nt0_nrt0_nw1",  {250, 350, 450, 550, 650, 1000}},
  {"nb2_mtb175_nj5t_nt0_nrt0_nw1",  {250, 350, 450, 550, 650, 1000}},
  {"nb1_mtb175_nj5t_nt0_nrt1_nw0",  {250, 350, 450, 550, 650, 1000}},
  {"nb2_mtb175_nj5t_nt0_nrt1_nw0",  {250, 350, 450, 550, 650, 1000}},
  {"nb1_mtb175_nj5t_nt1_nrt0_nw0",  {250, 350, 450, 550, 650, 1000}},
  {"nb2_mtb175_nj5t_nt1_nrt0_nw0",  {250, 350, 450, 550, 650, 1000}},

  // 1+1
  {"nb1_mtb175_nj5t_nt1_nrt0_nw1",  {250, 300, 400, 500, 1000}},
  {"nb2_mtb175_nj5t_nt1_nrt0_nw1",  {250, 300, 400, 500, 1000}},

  {"nb1_mtb175_nj5t_nt0_nrt1_nw1",  {250, 300, 400, 500, 1000}},
  {"nb2_mtb175_nj5t_nt0_nrt1_nw1",  {250, 300, 400, 500, 1000}},

  {"nb2_mtb175_nj5t_nt1_nrt1_nw0",  {250, 450, 1000}},

  // 2
  {"nb2_mtb175_nj5t_nt2_nrt0_nw0",  {250, 1000}},
  {"nb2_mtb175_nj5t_nt0_nrt2_nw0",  {250, 1000}},
  {"nb2_mtb175_nj5t_nt0_nrt0_nw2",  {250, 1000}},

};

// new cr
std::map<TString, TString> newcrcuts{
  {"nb1_mtb0_nj7_nrt0",          addCuts({nb1, lowmt, highnj})},
  {"nb2_mtb0_nj7_nrt0",          addCuts({nb2, lowmt, highnj})},
  {"nb1_mtb0_nj7_nrt1",          addCuts({nb1, lowmt, highnj})},
  {"nb2_mtb0_nj7_nrt1",          addCuts({nb2, lowmt, highnj})},

  // 0
  {"nb1_mtb175_nj7_nt0_nrt0_nw0",  addCuts({nb1, highmt, highnj})},
  {"nb2_mtb175_nj7_nt0_nrt0_nw0",  addCuts({nb2, highmt, highnj})},

  // 1
  {"nb1_mtb175_nj5t_nt0_nrt0_nw1",  addCuts({nb1, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt0_nrt0_nw1",  addCuts({nb2, highmt, nj5g})},
  {"nb1_mtb175_nj5t_nt0_nrt1_nw0",  addCuts({nb1, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt0_nrt1_nw0",  addCuts({nb2, highmt, nj5g})},
  {"nb1_mtb175_nj5t_nt1_nrt0_nw0",  addCuts({nb1, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt1_nrt0_nw0",  addCuts({nb2, highmt, nj5g})},

  // 1+1
  {"nb1_mtb175_nj5t_nt1_nrt0_nw1",  addCuts({nb1, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt1_nrt0_nw1",  addCuts({nb2, highmt, nj5g})},
  {"nb1_mtb175_nj5t_nt0_nrt1_nw1",  addCuts({nb1, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt0_nrt1_nw1",  addCuts({nb2, highmt, nj5g})},

  {"nb2_mtb175_nj5t_nt1_nrt1_nw0",  addCuts({nb2, highmt, nj5g})},


  // 2
  {"nb2_mtb175_nj5t_nt2_nrt0_nw0",  addCuts({nb2, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt0_nrt2_nw0",  addCuts({nb2, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt0_nrt0_nw2",  addCuts({nb2, highmt, nj5g})},
};
std::map<TString, TString> newcrlabels = newcrcuts;

// ICHEP cr
std::map<TString, TString> crcuts{
  {"nb1_mtb0_nj7_nrt0",          addCuts({nb1, nrt0, lowmt, highnj})},
  {"nb2_mtb0_nj7_nrt0",          addCuts({nb2, nrt0, lowmt, highnj})},
  {"nb1_mtb0_nj7_nrt1",          addCuts({nb1, nrt1, lowmt, highnj})},
  {"nb2_mtb0_nj7_nrt1",          addCuts({nb2, nrt1, lowmt, highnj})},

  // 0
  {"nb1_mtb175_nj7_nt0_nrt0_nw0",  addCuts({nb1, nt0, nrt0, nw0, highmt, highnj})},
  {"nb2_mtb175_nj7_nt0_nrt0_nw0",  addCuts({nb2, nt0, nrt0, nw0, highmt, highnj})},

  // 1
  {"nb1_mtb175_nj5t_nt0_nrt0_nw1",  addCuts({nb1, nt0, nrt0, nw1, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt0_nrt0_nw1",  addCuts({nb2, nt0, nrt0, nw1, highmt, nj5g})},
  {"nb1_mtb175_nj5t_nt0_nrt1_nw0",  addCuts({nb1, nt0, nrt1, nw0, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt0_nrt1_nw0",  addCuts({nb2, nt0, nrt1, nw0, highmt, nj5g})},
  {"nb1_mtb175_nj5t_nt1_nrt0_nw0",  addCuts({nb1, nt1, nrt0, nw0, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt1_nrt0_nw0",  addCuts({nb2, nt1, nrt0, nw0, highmt, nj5g})},

  // 1+1
  {"nb1_mtb175_nj5t_nt1_nrt0_nw1",  addCuts({nb1, nt1, nrt0, nw1, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt1_nrt0_nw1",  addCuts({nb2, nt1, nrt0, nw1, highmt, nj5g})},
  {"nb1_mtb175_nj5t_nt0_nrt1_nw1",  addCuts({nb1, nt0, nrt1, nw1, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt0_nrt1_nw1",  addCuts({nb2, nt0, nrt1, nw1, highmt, nj5g})},

  {"nb2_mtb175_nj5t_nt1_nrt1_nw0",  addCuts({nb2, nt1, nrt1, nw0, highmt, nj5g})},


  // 2
  {"nb2_mtb175_nj5t_nt2_nrt0_nw0",  addCuts({nb2, nt2, nrt0, nw0, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt0_nrt2_nw0",  addCuts({nb2, nt0, nrt2, nw0, highmt, nj5g})},
  {"nb2_mtb175_nj5t_nt0_nrt0_nw2",  addCuts({nb2, nt0, nrt0, nw2, highmt, nj5g})},

};
std::map<TString, TString> crlabels = crcuts;

std::map<TString, std::vector<int>> crMETbins = srMETbins;
std::map<TString, std::vector<int>> qcdcrMETbins = crMETbins;



map<TString, Category> srCatMap(){
  map<TString, Category> cmap;
  for (auto &name : srbins){
    cmap[name] = Category(name, srcuts.at(name), srlabels.at(name), BinInfo("met", "#slash{E}_{T}", srMETbins.at(name), "GeV"));
  }
  return cmap;
}

map<TString, Category> phoCatMap(){
  map<TString, Category> cmap;
  for (auto &name : srbins){
    cmap[name] = Category(name, crcuts.at(name), crlabels.at(name), BinInfo("met", "#slash{E}_{T}^{#gamma}", crMETbins.at(name), "GeV"));
  }
  return cmap;
}

map<TString, Category> lepCatMap(){
  TString varlabel = ADD_LEP_TO_MET ? "p_{T}^{W}" : "#slash{E}_{T}";
  map<TString, Category> cmap;
  for (auto &name : srbins){
    // No integration in Nb now!
    //    cmap[name] = Category(name, crcuts.at(name), crlabels.at(name), BinInfo("met", varlabel, crMETbins.at(name), "GeV"));
    cmap[name] = Category(name, srcuts.at(name), srlabels.at(name), BinInfo("met", varlabel, crMETbins.at(name), "GeV"));
  }
  return cmap;
}

map<TString, Category> qcdCatMap(){
  map<TString, Category> cmap;
  for (auto &name : srbins){
    cmap[name] = Category(name, crcuts.at(name), crlabels.at(name), BinInfo("met", "#slash{E}_{T}", qcdcrMETbins.at(name), "GeV"));
  }
  return cmap;
}

map<TString, Category> zllCatMap{
  {"on-z",  Category("on-z",  "dilepmass > 80 && dilepmass < 100",                      "on Z",   BinInfo("met", "#slash{E}_{T}^{ll}", vector<double>{200, 1000}, "GeV"))},
  {"off-z", Category("off-z", "dilepmass > 20 && (dilepmass < 80 || dilepmass > 100)",  "off Z",  BinInfo("met", "#slash{E}_{T}^{ll}", vector<double>{200, 1000}, "GeV"))}
};

//AP - new cats
map<TString, Category> newlepCatMap(){
  TString varlabel = ADD_LEP_TO_MET ? "p_{T}^{W}" : "#slash{E}_{T}";
  map<TString, Category> cmap;
  for (auto &name : srbins){
    // No integration in Nb now!
    //    cmap[name] = Category(name, crcuts.at(name), crlabels.at(name), BinInfo("met", varlabel, crMETbins.at(name), "GeV"));
    cmap[name] = Category(name, newcrcuts.at(name), newcrlabels.at(name), BinInfo("met", varlabel, crMETbins.at(name), "GeV"));
  }
  return cmap;
}

// ------------------------------------------------------------------------
// samples

BaseConfig phoConfig(){
  BaseConfig     config;

  config.inputdir = "/tmp/trees";
  config.outputdir = "/tmp/plots/HighMass/znunu";
  config.header = "#sqrt{s} = 13 TeV, "+lumistr+" fb^{-1}";

  config.addSample("singlepho",   "Data",           datadir+"/photoncr/singlepho",  "1.0",  datasel + trigPhoCR);

  config.addSample("photon",      "Photon",         "photoncr/photon",     phowgt, datasel + trigPhoCR);
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

  config.inputdir = "/tmp/trees";
  config.outputdir = "/tmp/plots/HighMass/zllcr";
  config.header = "#sqrt{s} = 13 TeV, "+lumistr+" fb^{-1}";

  config.addSample("dyll",      "Z#rightarrowll+jets",    "zllcr/zll",                      lepselwgt, datasel + trigDiLepCR);
  config.addSample("ttbar",     "t#bar{t}",               "zllcr/ttbar",                    lepselwgt, datasel + trigDiLepCR);
  config.addSample("doublelep", "Data",                   datadir+"/zllcr/doublelep",       "1.0",     datasel + trigDiLepCR);

  config.sel = baseForNorm;
  config.catMaps = zllCatMap;
  for (auto &c : zllCatMap) config.categories.push_back(c.first);

  return config;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BaseConfig lepConfig(){
  BaseConfig     config;

  config.inputdir = "/tmp/trees";
  config.outputdir = "/tmp/plots/HighMass/LLB";
  config.header = "#sqrt{s} = 13 TeV, "+lumistr+" fb^{-1}";

  if (ADD_LEP_TO_MET){
    config.addSample("singlelep",   "Data",          datadir+"/lepcr/singlelep",       "1.0",     datasel + trigLepCR + lepcrsel);
    config.addSample("ttbar",       "t#bar{t}",      "lepcr/ttbar",        onelepcrwgt, datasel + trigLepCR + lepcrsel);
    config.addSample("wjets",       "W+jets",        "lepcr/wjets",        onelepcrwgt, datasel + trigLepCR + lepcrsel);
    config.addSample("tW",          "tW",            "lepcr/tW",              onelepcrwgt, datasel + trigLepCR + lepcrsel);
    config.addSample("ttW",         "ttW",           "lepcr/ttW",             onelepcrwgt, datasel + trigLepCR + lepcrsel);
    //    config.addSample("qcd",         "QCD",           "lepcr/qcd",             onelepcrwgt, datasel + trigLepCR + lepcrsel);
  }else{
    config.addSample("singlelep",   "Data",          datadir+"/sr/met",    "1.0",     datasel + trigSR + revert_vetoes);
    config.addSample("ttbar",       "t#bar{t}",      "sr/ttbar",        lepselwgt, datasel + trigSR + revert_vetoes);
    config.addSample("wjets",       "W+jets",        "sr/wjets",        lepselwgt, datasel + trigSR + revert_vetoes);
    config.addSample("tW",          "tW",            "sr/tW",              lepselwgt, datasel + trigSR + revert_vetoes);
    config.addSample("ttW",         "ttW",           "sr/ttW",             lepselwgt, datasel + trigSR + revert_vetoes);
    //    config.addSample("qcd",         "QCD",           "sr/qcd",             lepselwgt, datasel + trigSR + revert_vetoes);
  }

  config.addSample("ttbar-sr",       "t#bar{t}",      "sr/ttbar",        lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("wjets-sr",       "W+jets",        "sr/wjets",        lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("tW-sr",          "tW",            "sr/tW",              lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("ttW-sr",         "ttW",           "sr/ttW",             lepvetowgt, datasel + trigSR + vetoes);
  //  config.addSample("qcd-sr",         "QCD",           "qcd-std/qcd",        lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("rare-sr",        "Rare",          "sr/rare",            lepvetowgt, datasel + trigSR + vetoes);

  config.sel = baseline;
  config.categories = srbins;
  config.catMaps = srCatMap();

  if(ICHEPCR){
    std::cout << "lepConfig: ICHEP2016 CR" << std::endl;
    config.crCatMaps = lepCatMap();
  }else{
    std::cout << "lepConfig: Moriond2017 CR" << std::endl;
    config.crCatMaps = newlepCatMap();
  }

  return config;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BaseConfig srConfig(){
  BaseConfig     config;

  config.inputdir = "/tmp/trees";
  config.outputdir = "/tmp/plots/testSR";
  config.header = "#sqrt{s} = 13 TeV, "+lumistr+" fb^{-1}";

  //  config.addSample("singlelep",   "Data",          datadir+"/sr/met",    "1.0",     datasel + trigSR + revert_vetoes);
  config.addSample("ttbar",       "t#bar{t}",      "ttbar",        lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("wjets",       "W+jets",        "wjets",        lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("znunu",       "Z#rightarrow#nu#nu", "znunu",   lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("qcd",         "QCD",           "qcd",          lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("tW",          "tW",            "tW",           lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("ttW",         "ttW",           "ttW",          lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("ttZ",         "ttZ",           "ttZ",          lepvetowgt, datasel + trigSR + vetoes);
  config.addSample("diboson",     "Diboson",       "diboson",      lepvetowgt, datasel + trigSR + vetoes);

  config.addSample("T2tt_450_250",  "T2tt(450,250)",  "T2tt_450_250",   sigwgt, datasel + trigSR + vetoes);
  config.addSample("T2tt_1000_500", "T2tt(1000,500)", "T2tt_1000_500",  sigwgt, datasel + trigSR + vetoes);
  config.addSample("T2tt_1100_1",   "T2tt(1100,1)",   "T2tt_1000_1",    sigwgt, datasel + trigSR + vetoes);
  config.addSample("T2bW_850_1",    "T2bW(850,1)",    "T2bW_850_1",     sigwgt, datasel + trigSR + vetoes);

  COLOR_MAP["T2tt_450_250"] = kRed;
  COLOR_MAP["T2tt_1000_500"] = kGreen+3;
  COLOR_MAP["T2tt_1100_1"] = kBlack;
  COLOR_MAP["T2bW_850_1"] = kMagenta;

  config.sel = baseline;
  config.categories = srbins;
  config.catMaps = srCatMap();

  return config;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BaseConfig qcdConfig(){
  BaseConfig     config;

  config.inputdir = "/tmp/trees";
  config.outputdir = "/tmp/plots/HighMass/QCD";
  config.header = "#sqrt{s} = 13 TeV, "+lumistr+" fb^{-1}";

  // qcdcr
  config.addSample("data-cr",     "Data",          datadir+"/sr/met",      "1.0",       datasel + trigSR + vetoes + dphi_invert);
  config.addSample("qcd-cr",      "QCD",           "qcd-std/qcd",          qcdwgt,      datasel + trigSR + dphi_invert);

  config.addSample("qcd-withveto-cr",  "QCD",      "qcd-std/qcd",          qcdvetowgt,  datasel + trigSR + vetoes + dphi_invert);
  config.addSample("qcd-withveto-sr",  "QCD",      "qcd-std/qcd",          qcdvetowgt,  datasel + trigSR + vetoes + dphi);

  // qcdcr: other bkg subtraction
  config.addSample("ttbar-cr",       "t#bar{t}",      "sr/ttbar",     lepvetowgt,  datasel + trigSR + vetoes + dphi_invert);
  config.addSample("wjets-cr",       "W+jets",        "sr/wjets",     lepvetowgt,  datasel + trigSR + vetoes + dphi_invert);
  config.addSample("tW-cr",          "tW",            "sr/tW",           lepvetowgt,  datasel + trigSR + vetoes + dphi_invert);
  config.addSample("ttW-cr",         "ttW",           "sr/ttW",          lepvetowgt,  datasel + trigSR + vetoes + dphi_invert);
  config.addSample("znunu-cr",       "Znunu",         "sr/znunu",        lepvetowgt,  datasel + trigSR + vetoes + dphi_invert);

  // onelepcr: norm correction for other bkg subtraction
  config.addSample("data-norm",      "Data",          datadir+"/sr/met",          "1.0",       datasel + trigSR + dphi + revert_vetoes);
  config.addSample("ttbar-norm",     "t#bar{t}",      "sr/ttbar",     lepselwgt,   datasel + trigSR + dphi + revert_vetoes);
  config.addSample("wjets-norm",     "W+jets",        "sr/wjets",     lepselwgt,   datasel + trigSR + dphi + revert_vetoes);
  config.addSample("tW-norm",        "tW",            "sr/tW",           lepselwgt,   datasel + trigSR + dphi + revert_vetoes);
  config.addSample("ttW-norm",       "ttW",           "sr/ttW",          lepselwgt,   datasel + trigSR + dphi + revert_vetoes);
  config.addSample("qcd-norm",       "QCD",           "qcd-std/qcd",     lepselwgt,   datasel + trigSR + dphi + revert_vetoes);

  // qcdsr
  config.addSample("qcd-sr",         "QCD",           "qcd-std/qcd",     qcdwgt,      datasel + trigSR + dphi);

  config.sel = baseNoDPhi;
  config.categories = srbins;
  config.catMaps = srCatMap();
  config.crCatMaps = qcdCatMap();

  return config;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BaseConfig sigConfig(){
  BaseConfig     config;

  config.inputdir = "/tmp/trees";
  config.outputdir = "/tmp/plots/HighMass/sig";
  config.header = "#sqrt{s} = 13 TeV, "+lumistr+" fb^{-1}";

  config.addSample("data-sr",        "Data",             datadir+"/sr/met",                    "1.0",  datasel + trigSR + vetoes);

  //  config.addSample("T2tt_425_325",  "T2tt(425,325)",  "sr/T2tt_425_325",  sigwgt, datasel + trigSR + vetoes);
  config.addSample("T2tt_600_300",  "T2tt(600,300)",  "sig/T2tt_600_300",  sigwgt, datasel + trigSR + vetoes);
  config.addSample("T2tt_850_50",   "T2tt(850,50)",   "sig/T2tt_850_50",   sigwgt, datasel + trigSR + vetoes);
  //  config.addSample("T2tt_950_1",    "T2tt(950,1)",    "sig/T2tt_950_1",  sigwgt, datasel + trigSR + vetoes);

  config.sel = baseline;
  config.categories = srbins;
  config.catMaps = srCatMap();

  return config;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*
const TString nt1 =  "nsdtop==1";
const TString nt2 =  "nsdtop>=2";
const TString nw0 =  "nsdw==0";
const TString nw1 =  "nsdw==1";
const TString nw2 =  "nsdw>=2";
const TString nrt0 = "nrestop==0";
const TString nrt1 = "nrestop==1";
const TString nrt2 = "nrestop>=2";
*/

map<TString, BinInfo> varDict {
  {"norm",      BinInfo("met", "#slash{E}_{T}", vector<int>{0, 1000}, "GeV")},
  {"met",       BinInfo("met", "#slash{E}_{T}", vector<int>{250, 350, 450, 550, 650, 1000}, "GeV")},
  {"metzg",     BinInfo("met", "#slash{E}_{T}^{#gamma/ll}", vector<int>{200, 300, 400, 500, 600, 800}, "GeV")},
  {"origmet",   BinInfo("origmet", "Original #slash{E}_{T}", 16, 0, 800, "GeV")},
  {"njets",     BinInfo("njets", "N_{j}", 12, -0.5, 11.5)},
  {"ncttstd",   BinInfo("ncttstd", "N_{ctt}", 3, -0.5, 2.5)},
  {"nt",        BinInfo("nsdtop", "N_{T}^{T}", 3, -0.5, 2.5)},
  {"nw",        BinInfo("nsdw", "N_{W}^{T}", 3, -0.5, 2.5)},
  {"nrt",       BinInfo("nrestop", "N_{resT}^{M}", 3, -0.5, 2.5)},
  {"nlbjets",   BinInfo("nlbjets", "N_{B}^{loose}", 5, -0.5, 4.5)},
  {"nbjets",    BinInfo("nbjets",  "N_{B}^{medium}", 5, -0.5, 4.5)},
  {"dphij1met", BinInfo("dphij1met", "#Delta#phi(j_{1},#slash{E}_{T})", 30, 0, 3)},
  {"dphij2met", BinInfo("dphij2met", "#Delta#phi(j_{2},#slash{E}_{T})", 30, 0, 3)},
  {"mtcsv12met",BinInfo("mtcsv12met", "min(m_{T}(b_{1},#slash{E}_{T}),m_{T}(b_{2},#slash{E}_{T}))", 6, 0, 300)},
  {"leptonpt",  BinInfo("leptonpt", "p_{T}^{lep} [GeV]", 16, 0, 800)},
  {"leptonptovermet",  BinInfo("leptonpt/met", "p_{T}^{lep}/#slash{E}_{T}", 20, 0, 1.)}
};

}
#endif /* ESTTOOLS_HMPARAMETERS_HH_ */
