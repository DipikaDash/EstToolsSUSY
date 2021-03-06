/*
 * Znunu.C
 *
 *  Created on: Sep 23, 2015
 *      Author: hqu
 */

#include <fstream>

#include "Syst_SR_Parameters.hh"
//#include "Syst_LowMET_Parameters.hh"

#include "../../EstMethods/LLBEstimator.hh"
#include "../../EstMethods/ZnunuEstimator.hh"
#include "../../EstMethods/QCDEstimator.hh"

using namespace EstTools;

vector<Quantity> getZnunuPred(){
  auto phocfg = phoConfig();
  ZnunuEstimator z(phocfg);
  z.zllcr_cfg = zllConfig();
  z.zll_normMap = normMap;
  z.phocr_normMap = phoNormMap;
  z.pred();
  z.printYields();
  return z.yields.at("_TF");
}

vector<Quantity> getQCDPred(){
  auto qcdcfg = qcdConfig();
  QCDEstimator q(qcdcfg);
  q.runBootstrapping = false;
  q.pred();
  q.printYields();
  return q.yields.at("_TF");
}

map<TString, vector<Quantity>> getLLBPred(){
  auto llbcfg = lepConfig();
  LLBEstimator l(llbcfg);
  l.pred();
  l.printYields();
  Quantity::removeNegatives(l.yields.at("ttZ-sr"));
  Quantity::removeNegatives(l.yields.at("diboson-sr"));

  return {
    {"ttbarplusw", l.yields.at("_TF")},
    {"ttZ",        l.yields.at("ttZ-sr")},
    {"diboson",    l.yields.at("diboson-sr")},
  };
}


void SystTopW(std::string outfile_path = "values_unc_wtoptag.conf"){

//  vector<TString> sdmvaSyst = {"sdMVAWeight", "sdMVAWeight_STATS_UP"};
//  vector<TString> resTopSyst = {"resTopWeight", "resTopWeight_STATS_UP"};

  vector<TString> sdmvaSystMistag = {"sdMVAWeight_MISTAG_STATS_W", "sdMVAWeight_MISTAG_STATS_T", "sdMVAWeight_MISTAG_NB", "sdMVAWeight_MISTAG_PS"};
  vector<TString> resTopSystMistag = {"resTopWeight_MISTAG_STATS", "resTopWeight_MISTAG_NB", "resTopWeight_MISTAG_PS"};

  vector<TString> sdmvaSystEff = {"sdMVAWeight_STATS_W", "sdMVAWeight_STATS_T", "sdMVAWeight_PS", "sdMVAWeight_GEN", "sdMVAWeight_MISTAG_UP_W", "sdMVAWeight_MISTAG_UP_T"};
  vector<TString> resTopSystEff = {"resTopWeight_STATS", "resTopWeight_PS", "resTopWeight_GEN", "resTopWeight_MISTAG_UP", "resTopWeight_NMATCH"};


  vector<TString> bkgnames  = {"qcd", "znunu", "diboson", "ttZ", "ttbarplusw"};
  map<TString, map<TString, vector<Quantity>>> proc_syst_pred; // {proc: {syst: yields}}
  for (auto &bkg : bkgnames){
    proc_syst_pred[bkg] = map<TString, vector<Quantity>>();
  }

  inputdir = "/data/hqu/trees/20170221_wtopSyst";

  // nominal
  {
    sys_name = "nominal";
    proc_syst_pred["znunu"][sys_name] = getZnunuPred();
    proc_syst_pred["qcd"][sys_name]   = getQCDPred();
    auto llb = getLLBPred();
    for (auto &p : llb) proc_syst_pred[p.first][sys_name] = p.second;
  }

  // mistag: relevant for all bkgs
  for (auto &sdsyst : sdmvaSystMistag){
    sys_name = sdsyst; sdmvawgt = sdsyst; restopwgt = "resTopWeight";
    cout << "\n\n ====== Using weights " << sdmvawgt << " and " << restopwgt << "======\n\n";
    proc_syst_pred["znunu"][sys_name] = getZnunuPred();
    proc_syst_pred["qcd"][sys_name]   = getQCDPred();
    auto llb = getLLBPred();
    for (auto &p : llb) proc_syst_pred[p.first][sys_name] = p.second;
  }
  for (auto &resyst : resTopSystMistag){
    sys_name = resyst; sdmvawgt = "sdMVAWeight"; restopwgt = resyst;
    cout << "\n\n ====== Using weights " << sdmvawgt << " and " << restopwgt << "======\n\n";
    proc_syst_pred["znunu"][sys_name] = getZnunuPred();
    proc_syst_pred["qcd"][sys_name]   = getQCDPred();
    auto llb = getLLBPred();
    for (auto &p : llb) proc_syst_pred[p.first][sys_name] = p.second;
  }


  // eff: only for ttbarplusw, ttZ and VV
  for (auto &sdsyst : sdmvaSystEff){
    sys_name = sdsyst; sdmvawgt = sdsyst; restopwgt = "resTopWeight";
    cout << "\n\n ====== Using weights " << sdmvawgt << " and " << restopwgt << "======\n\n";
    auto llb = getLLBPred();
    for (auto &p : llb) proc_syst_pred[p.first][sys_name] = p.second;
  }
  for (auto &resyst : resTopSystEff){
    sys_name = resyst; sdmvawgt = "sdMVAWeight"; restopwgt = resyst;
    cout << "\n\n ====== Using weights " << sdmvawgt << " and " << restopwgt << "======\n\n";
    auto llb = getLLBPred();
    for (auto &p : llb) proc_syst_pred[p.first][sys_name] = p.second;
  }

  cout << "\n\n Write unc to " << outfile_path << endl;
  ofstream outfile(outfile_path);
  auto config = lepConfig();

  for (auto &bkg : bkgnames){
    auto nominal_pred = proc_syst_pred[bkg]["nominal"];
    for (auto &sPair : proc_syst_pred[bkg]){
      if(sPair.first=="nominal") continue;
      auto uncs = sPair.second / nominal_pred;

      unsigned ibin = 0;
      for (auto &cat_name : config.categories){
        auto &cat = config.catMaps.at(cat_name);
        for (unsigned ix = 0; ix < cat.bin.nbins; ++ix){
          auto xlow = toString(cat.bin.plotbins.at(ix), 0);
          auto xhigh = (ix==cat.bin.nbins-1) ? "inf" : toString(cat.bin.plotbins.at(ix+1), 0);
          auto binname = "bin_" + cat_name + "_" + cat.bin.var + xlow + "to" + xhigh;
          auto uncType = sPair.first;
//          outfile << binname << "\t" << uncType << "\t" << bkg << "\t" << uncs.at(ibin).value << endl;
          double val = uncs.at(ibin).value;
          if (val>2 || std::isnan(val)) {
            cout << "Invalid unc, set to 100%: " << binname << "\t" << uncType << "\t" << bkg << "\t" << uncs.at(ibin).value << endl;
            val = 2;
          }else if (val<0.5){
            cout << "Invalid unc, set to -100%: " << binname << "\t" << uncType << "\t" << bkg << "\t" << uncs.at(ibin).value << endl;
            val = 0.001;
          }
          outfile << binname << "\t" << uncType << "\t" << bkg << "\t" << val << endl;
          ++ibin;
        }
      }
    }
  }

  outfile.close();

  cout << "Done!" << endl;

}
