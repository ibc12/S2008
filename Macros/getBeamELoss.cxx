#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"

#include <memory>
#include <stdexcept>
#include <string>

void getBeamELoss()
{
    std::string beam {"20Na"};

    std::string conf {};
    if(beam == "20Mg")
        conf = "../configs/data_20Mg.conf";
    else if(beam == "20Na")
        conf = "../configs/data.conf";
    else
        throw std::runtime_error("Undefined beam particle");

    ActRoot::DataManager dataman {conf, ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chainFilter {dataman.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(chainFilter.get());

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    auto gated {df.Filter("fLightIdx != -1")};

    auto def {gated.Define("ELoss",
                           [](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
                           {
                               auto& cl {tpc.fClusters[mer.fBeamIdx]};
                               double eloss {};
                               for(const auto& v : cl.GetVoxels())
                               {
                                   if(v.GetPosition().X() < 8)
                                       eloss += v.GetCharge();
                               }
                               return eloss;
                           },
                           {"MergerData", "TPCData"})};

    auto hELoss {def.Histo1D({"hELoss", "Beam ELoss;#DeltaE_{beam} [MeV]", 4000, 0, 60000}, "ELoss")};


    // Draw
    auto* c0 {new TCanvas {"c0", "Beam ELoss"}};
    hELoss->SetTitle(beam.c_str());
    hELoss->DrawClone();

    // Save to disk
    auto outfile {std::make_unique<TFile>(TString::Format("./Outputs/beam_eloss_%s.root", beam.c_str()), "recreate")};
    outfile->cd();
    hELoss->SetNameTitle(TString::Format("hELoss%s", beam.c_str()), beam.c_str());
    hELoss->Write();
    outfile->Close();
}
