#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TVirtualPad.h"
void compEcmBefAfter()
{
    // Open file
    auto* file {new TFile {"../PostAnalysis/Outputs/gated_ecm.root"}};
    file->ls();
    // Get histograms
    auto* hBef2d {file->Get<TH2D>("hEpRMg")};
    hBef2d->SetDirectory(nullptr);
    auto* hAft2d {file->Get<TH2D>("hEpRMg")};
    hAft2d->SetDirectory(nullptr);
    // Add them
    hBef2d->Add(hAft2d);

    // And now Ecm
    auto* hBef {file->Get<TH1D>("hECM;1")};
    hBef->SetDirectory(nullptr);
    auto* hAft {file->Get<TH1D>("hECM;2")};
    hAft->SetDirectory(nullptr);
    // Normalize
    // hBef->Scale(2 * hAft->Integral() / hBef->Integral());
    // Set titles and so
    hBef->SetTitle("Before");
    hAft->SetTitle("After");
    hBef->SetLineColor(1);
    hAft->SetLineColor(9);


    // Plot
    auto* c0 {new TCanvas {"c0", "Ecm canvas"}};
    c0->DivideSquare(2);
    c0->cd(1);
    hBef2d->Draw("colz");
    c0->cd(2);
    hBef->Draw("hist");
    hAft->Draw("hist same");
    gPad->BuildLegend();
}
