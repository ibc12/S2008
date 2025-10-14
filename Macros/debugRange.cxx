#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDF/HistoModels.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"

#include "Math/Vector3Dfwd.h"

void debugRange()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    dataman.SetRuns(31, 35);
    auto chain {dataman.GetChain()};
    auto chainFilter {dataman.GetChain(ActRoot::ModeType::EFilter)};

    chain->AddFriend(chainFilter.get());

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    auto gated {
        df.Filter([](ActRoot::MergerData& m) { return m.fHeavyIdx != -1 && !m.fLight.IsL1(); }, {"MergerData"})};

    auto def {gated
                  .Define("End",
                          [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                          {
                              // Get heavy cluster
                              auto heavy {tpc.fClusters[m.fHeavyIdx]};
                              heavy.SetUseExtVoxels(true);
                              // RP
                              auto rp {tpc.fRPs.front()};
                              // Line
                              auto line {heavy.GetLine()};
                              // Align using RP in pad units
                              line.AlignUsingPoint(rp);
                              // Sort alogn direction
                              auto& voxels {heavy.GetRefToVoxels()};
                              // std::sort(voxels.begin(), voxels.end());
                              heavy.SortAlongDir();
                              // Get last
                              auto end {heavy.GetVoxels().back().GetPosition()};
                              end += ROOT::Math::XYZVectorF {0.5, 0.5, 0.5};
                              // Scale point
                              end.SetX(end.X() * 2);
                              end.SetY(end.Y() * 2);
                              end.SetZ(end.Z() * 2.208);
                              // Scale line
                              line.Scale(2, 2.208);
                              // Project
                              auto proj {line.ProjectionPointOnLine(end)};
                              return proj;
                          },
                          {"MergerData", "TPCData"})
                  .Define("EndX", "End.X()")
                  .Define("EndY", "End.Y()")
                  .Define("EndZ", "End.Z()")};

    // Histograms
    ROOT::RDF::TH1DModel mPos {"hPos", "Position;Pos;Counts", 256, 0, 256};
    auto hRPx {def.Histo1D(mPos, "fRP.fCoordinates.fX")};
    auto hRPy {def.Histo1D(mPos, "fRP.fCoordinates.fY")};
    auto hRPz {def.Histo1D(mPos, "fRP.fCoordinates.fZ")};

    auto hEndX {def.Histo1D(mPos, "EndX")};
    auto hEndY {def.Histo1D(mPos, "EndY")};
    auto hEndZ {def.Histo1D(mPos, "EndZ")};

    auto* c0 {new TCanvas {"c0", "Debug range"}};
    c0->DivideSquare(6);
    c0->cd(1);
    hRPx->DrawClone();
    c0->cd(2);
    hRPy->DrawClone();
    c0->cd(3);
    hRPz->DrawClone();
    c0->cd(4);
    hEndX->DrawClone();
    c0->cd(5);
    hEndY->DrawClone();
    c0->cd(6);
    hEndZ->DrawClone();
}
