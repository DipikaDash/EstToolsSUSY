#include "../EstMethods/QCDEstimator.hh"

#include "SRParameters.hh"

using namespace EstTools;
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

vector<Quantity> QCDPred(){
  auto config = qcdConfig();

  QCDEstimator z(config);

  //  z.runBootstrapping = false;
  z.splitTF = SPLITTF;
  z.pred();
  z.naiveTF();
  z.printYields();

  std::map<TString,int> digits;
  digits["_DATA"] = 0; // indicate it's data for proper formatting
  digits["_QCDTF_CR_to_SR_noextrap"] = -3;
  digits["_QCDTF_SR_extrap"] = -3;

  z.printYieldsTableLatex({"_DATA", "_TF", "_pred"}, labelMap, "yields_qcd_lm.tex","lm", digits); //LM
  if(z.splitTF){
    z.printYieldsTableLatex({"_DATA", "_TF", "_QCDTF_CR_to_SR_noextrap", "_QCDTF_SR_extrap", "_pred"}, labelMap, "/tmp/hqu/yields_qcd_hm.tex","hm", digits);
  }else{
    z.printYieldsTableLatex({"_DATA", "_TF", "_pred"}, labelMap, "yields_qcd_hm.tex","hm", digits);
  }

  return z.yields.at("_pred");
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void plotQCDCR(){
  auto config = qcdConfig();
  config.catMaps = config.crCatMaps;

  BaseEstimator z(config.outputdir);
  z.setConfig(config);

  z.plotDataMC({"ttbar-cr", "wjets-cr", "minor-cr", "rare-cr", "znunu-cr", "qcd-withveto-cr"}, "data-cr", false);

}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void plotQCDInclusive(){
  auto config = qcdConfig();

  config.sel ="";// "njets>=2 && met>250 && metovsqrtht>10";
  //  config.addSample("qcd",  "QCD",      "qcd-4bd/qcd",          qcdvetowgt,  datasel + trigSR + vetoes + dphi_invert);
  config.addSample("qcd",  "QCD",      "qcd_orig_2018",          wgtvar,  datasel + trigSR + vetoes + dphi_invert);

  //  config.addSample("qcd",  "QCD",      "qcd-std/qcd",          qcdvetowgt,  datasel + trigSR + vetoes + dphi_invert);
  BaseEstimator z(config);

  vector<TString> mc_samples = {"ttbar-cr", "wjets-cr","minor-cr","znunu-cr","rare-cr","qcd-cr"};
  TString data_sample = "data-cr";

  map<TString, BinInfo> varDict {
    //  {"met",       BinInfo("MET_pt", "#slash{E}_{T}", 20, 0, 2000, "GeV")},
    //	  {"ht",        BinInfo("Stop0l_HT",  "H_{T}", 60, 0, 3000, "GeV")},
    {"njets",     BinInfo("nJet", "N_{j}", 12, -0.5, 11.5,"")}, 
    {"njets202p2",     BinInfo("njet_j202p2", "N_{j202p2}", 12, -0.5, 11.5,"")},                                                                                                  
      {"njets302p2",     BinInfo("njet_j302p2", "N_{j302p2}", 12, -0.5, 11.5,"")},
      {"njets30",     BinInfo("njet_j30", "N_{j30}", 12, -0.5, 11.5,"")},                                                                                                      
	{"njets30to100",     BinInfo("njet_j30to100", "N_{j30to100}", 12, -0.5, 11.5,"")},                                                                                         
	  {"njets100to200",     BinInfo("njet_j100to200", "N_{j100to200}", 12, -0.5, 11.5,"")},                                                                                      
	    {"njets200to400",     BinInfo("njet_j200to400", "N_{j200to400}", 12, -0.5, 11.5,"")},                                                                                    
	      {"njets400to600",     BinInfo("njet_j400to600", "N_{j400to600}", 12, -0.5, 11.5,"")},
		{"njets600to1k",     BinInfo("njet_j600to1000", "N_{j600to1000}", 12, -0.5, 11.5,"")},
		  {"njets1kto1p5k",     BinInfo("njet_j1kto1p5k", "N_{j1kto1p5k}", 12, -0.5, 11.5,"")},
		    {"njetsj20bar",     BinInfo("njet_j20bar", "N_{j20bar}", 12, -0.5, 11.5,"")},
		      {"njetsj20end",     BinInfo("njet_j20end", "N_{j20end}", 12, -0.5, 11.5,"")},
			{"njetsj30bar",     BinInfo("njet_j30bar", "N_{j30bar}", 12, -0.5, 11.5,"")},
			  {"njetsj30end",     BinInfo("njet_j30end", "N_{j30end}", 12, -0.5, 11.5,"")},

 
	      //	      {"j1pt",        BinInfo("JetPass_pt[0]",  "jet1p_{T}", 50, 0, 2000, "GeV")},
	     //   {"j1eta",        BinInfo("JetPass_eta[0]",  "jet1_{eta}", 10, -5, 5, "")},
	     // 	 {"j1phi",        BinInfo("JetPass_phi[0]",  "jet1_{phi}", 10, -5, 5, "")},
    //       {"njettight",        BinInfo("nJets2p2",  "N_{jtight}", 12, -0.5, 11.5, "")},
          // {"j2pt",        BinInfo("JetPass_pt[1]",  "jet2p_{T}", 50, 0, 2000, "GeV")}
	  // {"j2eta",        BinInfo("JetPass_eta[0]",  "jet1_{eta}", 10, -5, 5, "")}
          // {"j2phi",        BinInfo("JetPass_phi[1]",  "jet2_{phi}", 10, -5, 5, "")}
          // {"j1nhef",    BinInfo("Jet_neHEF[0]",  "j1neHEnFrac", 20, 0, 1, "")},                                           
			    // {"j1nemef",    BinInfo("Jet_neEmEF[0]",  "j1neEmEnFrac", 20, 0, 1, "")},
	       // {"j1chef",    BinInfo("Jet_chHEF[0]",  "j1chHEnFrac", 20, 0, 1, "")},
	       //	{"j1cemef",    BinInfo("Jet_chEmEF[0]",  "j1chEmEnFrac", 20, 0, 1, "")},

	  //    {"nt",        BinInfo("nsdtoploose", "N_{t}", 2, -0.5, 1.5)},
    //    {"nw",        BinInfo("nsdwloose", "N_{W}", 2, -0.5, 1.5)},
    //    {"nlbjets",   BinInfo("nlbjets", "N_{B}^{loose}", 5, -0.5, 4.5)},
    //    {"nbjets",    BinInfo("nbjets",  "N_{B}^{medium}", 5, -0.5, 4.5)},
    //    {"dphij1met", BinInfo("dphij1met", "#Delta#phi(j_{1},#slash{E}_{T})", 32, 0, 3.2)},
    //    {"dphij2met", BinInfo("dphij2met", "#Delta#phi(j_{2},#slash{E}_{T})", 32, 0, 3.2)},
    //    {"dphij3met", BinInfo("dphij3met", "#Delta#phi(j_{2},#slash{E}_{T})", 32, 0, 3.2)},
    //    {"mtcsv12met",BinInfo("mtcsv12met", "min(m_{T}(b_{1},#slash{E}_{T}),m_{T}(b_{2},#slash{E}_{T}))", 6, 0, 300)},
    //    {"metovsqrtht",BinInfo("metovsqrtht", "#slash{E}_{T}/#sqrt{H_{T}}", 10, 0, 20)},
    //    {"dphij1lmet",BinInfo("dphij1lmet", "#Delta#phi(j_{1}^{ISR},#slash{E}_{T})", vector<double>{0, 2, 3})},
    //    {"njl",       BinInfo("njl", "N_{j}^{ISR}", 5, -0.5, 4.5)},
    //    {"j1lpt",     BinInfo("j1lpt", "p_{T}(j_{1}^{ISR}) [GeV]", 20, 0, 1000)},PV_npvsGood
    //    {"csvj1pt",   BinInfo("csvj1pt", "p_{T}(b_{1}) [GeV]", 8, 20, 100)}
    //     {"Jet_ID",   BinInfo("Pass_JetID", "tight JetID", 2, 0, 2)}                                                                           
		  
			    };
  for (auto &var : varDict){
    z.plotDataMC(var.second, mc_samples, data_sample, Category::dummy_category(),false,"",true );
  }



}
