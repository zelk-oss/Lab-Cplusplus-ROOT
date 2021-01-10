#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TStyle.h"
#include <iostream>

void Analize(int genLoops = 1e5) {
  TFile *file = new TFile("rivelatore.root", "UPDATE");

  // resuming histograms from ROOT file
  TH1F *hPTypes = (TH1F *)file->Get("hPTypes");
  TH1F *hTheta = (TH1F *)file->Get("hTheta");
  TH1F *hPhi = (TH1F *)file->Get("hPhi");
  TH1F *hImpulse = (TH1F *)file->Get("hImpulse");
  TH1F *hTImpulse = (TH1F *)file->Get("hTImpulse");
  TH1F *hEnergy = (TH1F *)file->Get("hEnergy");
  TH1F *hInvMass = (TH1F *)file->Get("hInvMass");
  TH1F *hSameCharge = (TH1F *)file->Get("hSameCharge");
  TH1F *hOddCharge = (TH1F *)file->Get("hOddCharge");
  TH1F *hPKconcordant = (TH1F *)file->Get("hPKconcordant");
  TH1F *hPKdiscordant = (TH1F *)file->Get("hPKdiscordant");
  TH1F *hKDecay = (TH1F *)file->Get("hKDecay");
  TH1F *hDiff1 = (TH1F *)file->Get("hDiff1");
  TH1F *hDiff2 = (TH1F *)file->Get("hDiff2");

  hKDecay->Sumw2();

  // histograms cosmetics
  hPTypes->SetMinimum(0);
  hPTypes->SetFillColor(kBlue - 4);
  hPTypes->GetXaxis()->SetTitle("Particle types");
  hPTypes->GetYaxis()->SetTitle("Number of occurrences");
  hPTypes->GetXaxis()->SetBinLabel(1, "Pions +");
  hPTypes->GetXaxis()->SetBinLabel(2, "Pions -");
  hPTypes->GetXaxis()->SetBinLabel(3, "Kaons +");
  hPTypes->GetXaxis()->SetBinLabel(4, "Kaons +");
  hPTypes->GetXaxis()->SetBinLabel(5, "Protons+");
  hPTypes->GetXaxis()->SetBinLabel(6, "Protons-");
  hPTypes->GetXaxis()->SetBinLabel(7, "K*");

  hTheta->SetFillColor(kMagenta - 9);
  hTheta->GetXaxis()->SetTitle("Theta (radiants)");
  hTheta->GetYaxis()->SetTitle("Number of occurrences");

  hPhi->SetFillColor(kGreen);
  hPhi->GetXaxis()->SetTitle("Phi (radiants)");
  hPhi->GetYaxis()->SetTitle("Number of occurrences");

  hImpulse->SetFillColor(kTeal - 9);
  hImpulse->GetXaxis()->SetTitle("Impulse (GeV)");
  hImpulse->GetYaxis()->SetTitle("Number of occurrences");

  hTImpulse->SetFillColor(kOrange);
  hTImpulse->GetXaxis()->SetTitle("Transverse impulse (GeV)");
  hTImpulse->GetYaxis()->SetTitle("Number of occurrences");
  hTImpulse->Draw();

  hEnergy->SetFillColor(kYellow);
  hEnergy->GetXaxis()->SetTitle("Energy (GeV)");
  hEnergy->GetYaxis()->SetTitle("Number of occurrences");
  hEnergy->Draw();

  TCanvas *massCanvas = new TCanvas("massCanvas", "Invariant mass histograms");
  massCanvas->Divide(2, 3);

  hInvMass->SetFillColor(kBlue);
  hInvMass->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
  hInvMass->GetYaxis()->SetTitle("Number of occurrences");
  massCanvas->cd(1);
  hInvMass->Draw();

  hSameCharge->SetFillColor(kOrange - 3);
  hSameCharge->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
  hSameCharge->GetYaxis()->SetTitle("Number of Occurrences");
  massCanvas->cd(2);
  hSameCharge->Draw();

  hOddCharge->SetFillColor(kRed - 4);
  hOddCharge->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
  hOddCharge->GetYaxis()->SetTitle("Number of Occurrences");
  massCanvas->cd(3);
  hOddCharge->Draw();

  hPKconcordant->SetFillColor(kMagenta - 9);
  hPKconcordant->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
  hPKconcordant->GetYaxis()->SetTitle("Number of Occurrences");
  massCanvas->cd(4);
  hPKconcordant->Draw();

  hPKdiscordant->SetFillColor(kBlue - 3);
  hPKdiscordant->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
  hPKdiscordant->GetYaxis()->SetTitle("Number of Occurrences");
  massCanvas->cd(5);
  hPKdiscordant->Draw();

  hKDecay->SetFillColor(kCyan - 8);
  hKDecay->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
  hKDecay->GetYaxis()->SetTitle("Number of Occurrences");

  hDiff1->SetFillColor(kCyan - 8);
  hDiff1->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
  hDiff1->GetYaxis()->SetTitle("Number of Occurrences");

  hDiff2->SetFillColor(kCyan - 8);
  hDiff2->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
  hDiff2->GetYaxis()->SetTitle("Number of Occurrences");


  /////////////////////////////////////////////////////
  // Analizing particle types histograms
  std::cout << "_____________________________________________" << '\n';
  std::cout << "\033[0;31mParticle Types distribution stats:\033[0m\n";
  for (int i = 0; i < 9;
       i++) // checking particle types histo underflow and overflow as well
  {
    std::cout << "Bin " << i << " Entries fraction = "
              << hPTypes->GetBinContent(i) / (genLoops * 100) << " ± "
              << hPTypes->GetBinError(i) / (genLoops * 100) << '\n';
  }
  std::cout << "_____________________________________________"
            << "\n\n";


  /////////////////////////////////////////////////////
  // analizing angles distributions
  std::cout << "_____________________________________________" << '\n';
  std::cout << "\033[0;31mPolar angle uniform distribution fit:\033[0m\n";
  hTheta->Fit("pol0", "Q"); // uniform distribution fit
  hTheta->GetFunction("pol0")->SetLineColor(kRed);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  hTheta->GetFunction("pol0")->Draw("SAME");
  std::cout << "Chi Square = " << hTheta->GetFunction("pol0")->GetChisquare()
            << '\n';
  std::cout << "DOF = " << hTheta->GetFunction("pol0")->GetNDF() << '\n';
  std::cout << "Reduced Chi Square = "
            << hTheta->GetFunction("pol0")->GetChisquare() /
                   hTheta->GetFunction("pol0")->GetNDF()
            << "\n";

  std::cout << "_____________________________________________"
            << "\n\n";

  ////////////////////////////////////////////////////////

  std::cout << "_____________________________________________" << '\n';
  std::cout << "\033[0;31mAzimuthal angle uniform distribution fit:\033[0m\n";
  hPhi->Fit("pol0", "Q"); // uniform distribution fit
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  hPhi->GetFunction("pol0")->SetLineColor(kRed);
  hPhi->GetFunction("pol0")->Draw("SAME");
  std::cout << "Chi Square = " << hPhi->GetFunction("pol0")->GetChisquare()
            << '\n';
  std::cout << "DOF = " << hPhi->GetFunction("pol0")->GetNDF() << '\n';
  std::cout << "Reduced Chi Square = "
            << hPhi->GetFunction("pol0")->GetChisquare() /
                   hPhi->GetFunction("pol0")->GetNDF()
            << "\n";

  std::cout << "_____________________________________________"
            << "\n\n";

  //////////////////////////////////////////////////////////
  // analizing impulse distribution
  std::cout << "_____________________________________________" << '\n';
  std::cout << "\033[0;31mImpulse exponential fit:\033[0m\n";
  hImpulse->Fit("expo", "Q"); // exponential dist. fit
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  hImpulse->GetFunction("expo")->SetLineColor(kRed);
  hImpulse->GetFunction("expo")->Draw("SAME");
  std::cout << "Chi Square = " << hImpulse->GetFunction("expo")->GetChisquare()
            << '\n';
  std::cout << "DOF = " << hImpulse->GetFunction("expo")->GetNDF() << '\n';
  std::cout << "Reduced Chi Square = "
            << hImpulse->GetFunction("expo")->GetChisquare() /
                   hImpulse->GetFunction("expo")->GetNDF()
            << "\n";
  std::cout << "Mean = " << hImpulse->GetMean() << " ± "
            << hImpulse->GetMeanError() << "(GeV)\n";
  std::cout << "Mean construction value is 1 GeV. Are values compatible?\n";

  std::cout << "_____________________________________________"
            << "\n\n";

  /////////////////////////////////////////////////////////
  // Decay histograms analysis
  std::cout << "_____________________________________________" << '\n';
  std::cout << "\033[0;31mDecay of K* gaussian fit:\033[0m\n";
  hKDecay->Fit("gaus", "", "", 0.72, 1.08);
  hKDecay->GetFunction("gaus")->SetLineColor(kRed);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  hKDecay->GetFunction("gaus")->Draw("SAME");
  std::cout << "Chi Square = " << hKDecay->GetFunction("gaus")->GetChisquare()
            << '\n';
  std::cout << "DOF = " << hKDecay->GetFunction("gaus")->GetNDF() << '\n';
  std::cout << "Reduced Chi Square = "
            << hKDecay->GetFunction("gaus")->GetChisquare() /
                   hKDecay->GetFunction("gaus")->GetNDF()
            << "\n";

  std::cout << "_____________________________________________"
            << "\n\n";

  std::cout << "_____________________________________________" << '\n';
  std::cout << "\033[0;31mOpposite and Same charge particles difference "
               "gaussian fit:\033[0m\n";
  hDiff1->Fit("gaus", "", "", 0.72, 1.08);
  hDiff1->GetFunction("gaus")->SetLineColor(kRed);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  hDiff1->GetFunction("gaus")->Draw("SAME");

  std::cout << "Chi Square = " << hDiff1->GetFunction("gaus")->GetChisquare()
            << '\n';
  std::cout << "DOF = " << hDiff1->GetFunction("gaus")->GetNDF() << '\n';
  std::cout << "Reduced Chi Square = "
            << hDiff1->GetFunction("gaus")->GetChisquare() /
                   hDiff1->GetFunction("gaus")->GetNDF()
            << "\n";

  std::cout << "_____________________________________________"
            << "\n\n";

  std::cout << "_____________________________________________" << '\n';
  std::cout << "\033[0;31mOpposite and Same charge Pions and Kaons difference "
               "gaussian fit:\033[0m\n";
  hDiff2->Fit("gaus", "", "", 0.72, 1.08);
  hDiff2->GetFunction("gaus")->SetLineColor(kRed);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  hDiff2->DrawCopy("Hist");
  hDiff2->DrawCopy("Same", "Func");
  std::cout << "Chi Square = " << hDiff2->GetFunction("gaus")->GetChisquare()
            << '\n';
  std::cout << "DOF = " << hDiff2->GetFunction("gaus")->GetNDF() << '\n';
  std::cout << "Reduced Chi Square = "
            << hDiff2->GetFunction("gaus")->GetChisquare() /
                   hDiff2->GetFunction("gaus")->GetNDF()
            << "\n";
  std::cout << "_____________________________________________"
            << "\n\n";

  // particle types, impulse and angles canvas
  TCanvas *c1 = new TCanvas("c1", "Detector statistics");
  c1->Divide(2, 2);
  c1->cd(1);
  hPTypes->Draw();
  c1->cd(2);
  hImpulse->Draw();
  c1->cd(3);
  hTheta->Draw();
  c1->cd(4);
  hPhi->Draw();

  // k* decay and difference histograms canvas
  TCanvas *c2 = new TCanvas("c2", "K* decay statistics");
  c2->Divide(2, 2);
  c2->cd(1);
  hKDecay->Draw("bar0");
  c2->cd(2);
  hDiff1->Draw();
  c2->cd(3);
  hDiff2->Draw();

  c1->SaveAs("types-impulse-angles.pdf");
  c2->SaveAs("Kstar-stats.pdf");
  massCanvas->SaveAs("invariant-masses.pdf");
  c1->SaveAs("types-impulse-angles.root");
  c2->SaveAs("Kstar-stats.root");
  massCanvas->SaveAs("invariant-masses.root");

  file->Close();
}