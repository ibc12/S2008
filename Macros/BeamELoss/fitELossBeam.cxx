#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TString.h"
#include "TCanvas.h"
#include "ActSRIM.h"
#include "TF1.h"
#include "ActLine.h"

void fitELossBeam()
{
    // Get histos from files
    TFile* file = TFile::Open("./Outputs/DE_20Ne_10pads.root");
    TH1D* h = file->Get<TH1D>("hDE_20Ne");

    auto srim {new ActPhysics::SRIM};
    srim->ReadTable("20Ne", "../../Calibrations/SRIM/20Ne_950mbar_95-5.txt.txt");
    double E {5.5 * 20};
    auto eLoss {E - srim->Slow("20Ne", E, 20)};

    // Fit the histograms into gausians
    h->Fit("gaus", "Q");
    // Get functions to plot them
    auto f {h->GetFunction("gaus")};

    // Let's do two lines, starting in 0,0 and going to eLoss from srim, and the center of the gausian fit to calibrate
    auto mean {f->GetParameter(1)};

    ROOT::Math::XYZPointF point_0 {0,0,0};
    ROOT::Math::XYZPointF point {eLoss, mean, 0};
    auto line {new ActRoot::Line(point_0, point)};
    // Get the sigmas and convert it to energy
    auto sigma_raw {f->GetParameter(2)};
    auto sigma_MeV {line->MoveToY(sigma_raw).X()};

    // Do the couts for the mean and the sigmas
    std::cout << "20Ne: E_loss from SRIM: " << eLoss << " MeV, mean from fit: " << mean << " u.a., sigma from fit: " << sigma_raw << " u.a. = " << sigma_MeV << " MeV" << std::endl;
    
    // Plot
    auto c{new TCanvas("c")};
    h->DrawClone();
    f->SetLineColor(kRed);
    f->Draw("same");
}