//to execute program:
//check runROOT.sh file is current directory, then type
//chmod +x runROOT.sh
//./runROOT.sh

// standard library includes
#include <cmath>
#include <cstdlib>
#include <iostream>

// ROOT includes
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TRandom.h"

// program files
#include "Particle.h"
#include "ParticleType.h"
#include "ResonanceType.h"

int Generate(int genLoops = 1e5) {
  int arrayDim = 130;
  int N = 100;
  // array of particle, dimension takes account of potential decay product
  Particle particle[arrayDim];
  // Particle types added to relative array
  Particle::AddParticleType("Pion+", 0.13957, +1);
  Particle::AddParticleType("Pion-", 0.13957, -1);
  Particle::AddParticleType("Kaon+", 0.49367, +1);
  Particle::AddParticleType("Kaon-", 0.49367, -1);
  Particle::AddParticleType("Proton+", 0.93827, +1);
  Particle::AddParticleType("Proton-", 0.93827, -1);
  Particle::AddParticleType("K*", 0.89166, 0, 0.050);

  // ROOT file
  TFile *file = new TFile("rivelatore.root", "RECREATE");

  // particle types histogram, generation percentages are later defined.
  TH1F *hPTypes = new TH1F("hPTypes", "Particle Types distribution", 7, 0., 7);
  TH1F *hTheta =
      new TH1F("hTheta", "Polar angle distribution", 400, 0, TMath::Pi());
  TH1F *hPhi =
      new TH1F("hPhi", "Azimuthal angle distribution", 400, 0, 2 * TMath::Pi());
  TH1F *hImpulse = new TH1F("hImpulse", "Impulse distribution", 400, 0, 8);
  TH1F *hTImpulse = new TH1F(
      "hTImpulse", "Transverse (X & Y) Impulse distribution", 400, 0, 8);
  TH1F *hEnergy = new TH1F("hEnergy", "Particle Energy", 400, 0, 4.5);

  // invariant mass combinations histograms declarations
  TH1F *hInvMass =
      new TH1F("hInvMass", "Invariant mass: all particle", 400, 0, 5);
  TH1F *hSameCharge = new TH1F(
      "hSameCharge", "Invariant mass: same charge particle", 400, 0, 5);
  TH1F *hOddCharge = new TH1F(
      "hOddCharge", "Invariant mass: opposite charge particle", 400, 0, 5);
  TH1F *hPKconcordant =
      new TH1F("hPKconcordant", "Invariant mass: Pions and Kaons, same charge",
               400, 0, 5);
  TH1F *hPKdiscordant =
      new TH1F("hPKdiscordant",
               "Invariant mass: Pions and Kaons, opposite charge", 400, 0, 5);
  TH1F *hKDecay = new TH1F("hKDecay", "K* from decay", 400, 0, 5);

  gRandom->SetSeed();

  // events loop, number of iterations = genLoops
  for (int i = 0; i < genLoops; i++) {
    double phi;   // azimuthal coordinate
    double theta; // polar coordinate
    double P;     // impulse

    int extraPos = 0; // starting point for extra particle (from K* decay)

    for (int j = 0; j < N; j++) {
      // initialization of angle coordinate and impulse variables.
      phi = gRandom->Uniform(2 * TMath::Pi());
      theta = gRandom->Uniform(TMath::Pi());
      P = gRandom->Exp(1); // medium impulse = 1Gev
      double Px = P * TMath::Sin(theta) * TMath::Cos(phi);
      double Py = P * TMath::Sin(theta) * TMath::Sin(phi);
      double Pz = P * TMath::Cos(theta);

      double transverseImpulse =
          std::sqrt(Px * Px + Py * Py); // impulse on X and Y axis
      particle[j].SetP(Px, Py, Pz);

      // random number, uniformly distributed
      double index = gRandom->Uniform(1);
      // particle types percentages defined here
      if (index < 0.4) {
        particle[j].SetParticle("Pion+");
        // particle[j].GetCharge();

      } else if (index < 0.8) {
        particle[j].SetParticle("Pion-");

      } else if (index < 0.85) {
        particle[j].SetParticle("Kaon+");

      } else if (index < 0.9) {
        particle[j].SetParticle("Kaon-");

      } else if (index < 0.945) {
        particle[j].SetParticle("Proton+");

      } else if (index < 0.99) {
        particle[j].SetParticle("Proton-");
      } else if (index < 0.995)
            { //K* into Pion+ Kaon-
                particle[j].SetParticle("K*");
                particle[N + extraPos].SetParticle("Pion+");
                particle[N + extraPos + 1].SetParticle("Kaon-");
                particle[j].Decay2body(particle[N + extraPos], particle[N + extraPos + 1]);
                extraPos++;
                extraPos++;
            }
            else
            { //K* into Pion- Kaon+
                particle[j].SetParticle("K*");
                particle[N + extraPos].SetParticle("Pion-");
                particle[N + extraPos + 1].SetParticle("Kaon+");
                particle[j].Decay2body(particle[N + extraPos], particle[N + extraPos + 1]);
                extraPos++;
                extraPos++;
            }
      // particle types histogram filled according to percentages distribution
      hPTypes->Fill(particle[j].GetIndex());
hTheta->Fill(theta);
      hPhi->Fill(phi);
      hImpulse->Fill(P);
      hTImpulse->Fill(transverseImpulse);
      hEnergy->Fill(particle[j].GetEnergy());
    }

    // filling invariant mass histogram
   int newArrayDim = N + extraPos;
    for (int h = 0; h < newArrayDim - 1; h++) {
      for (int k = h + 1; k < newArrayDim; k++) {
        int hCharge = particle[h].GetCharge();
        int kCharge = particle[k].GetCharge();
        int hIndex = particle[h].GetIndex();
        int kIndex = particle[k].GetIndex();

        double invMass = particle[h].InvMass(particle[k]);
        hInvMass->Fill(invMass); // filling invariant mass (with no charge
                                 // constraints) histogram

        if ((hCharge > 0 && kCharge > 0) || (hCharge < 0 && kCharge < 0)) {
          hSameCharge->Fill(invMass); // if confronted particle have the same
                                      // charge hSameCharge is filled
        }

        if (((hCharge > 0) && (kCharge < 0)) ||
            ((hCharge < 0) && (kCharge > 0))) {
          hOddCharge->Fill(invMass); // particle with opposite charge
        }

        if (((hIndex == 0) && (kIndex == 2)) ||
            ((hIndex == 1) && (kIndex == 3)) ||
            ((hIndex == 2) && (kIndex == 0)) ||
            ((hIndex == 3) && (kIndex == 1))) {
          hPKconcordant->Fill(invMass); // Pions and Kaons with equal charge
        }

        if (((hIndex == 0) && (kIndex == 3)) ||
            ((hIndex == 1) && (kIndex == 2)) ||
            ((hIndex == 3) && (kIndex == 0)) ||
            ((hIndex == 2) && (kIndex == 1))) {
          hPKdiscordant->Fill(invMass); // Pions and Kaons with opposite charge
        }
      }
    }

    // if any K* particle decayed => filling of relative invariant mass
    // histogram
    if (extraPos != 0) {
      for (int f = 0; f < extraPos; f += 2) {
        hKDecay->Fill(particle[N + f].InvMass(particle[N + f + 1]));
      }
    }
  }

  // definition of difference histograms
  // hDiff1 = hOddCharge - hSameCharge
  TH1F *hDiff1 = new TH1F(
      "hDiff1", "Opposite and Same charge particle difference", 400, 0, 5);
  hDiff1->Sumw2();
  hDiff1->Add(hOddCharge, hSameCharge, 1, -1);
  hDiff1->SetEntries(hDiff1->Integral());

  // hDiff2 = Pions-Kaons opposite charge - Pions-Kaons same charge
  TH1F *hDiff2 =
      new TH1F("hDiff2", "Opposite and Same charge Pions and Kaons difference",
               400, 0, 5);
  hDiff2->Sumw2();
  hDiff2->Add(hPKdiscordant, hPKconcordant, 1, -1);
  hDiff2->SetEntries(hDiff2->Integral());

  /////////////Saving histograms to file/////////////////

  TCanvas *c3 = new TCanvas("c3", "Types, angles, impulse");
  c3->Divide(2, 2);
  c3->cd(1);
  hPTypes->Write();
  c3->cd(2);
  hTheta->Write();
  c3->cd(3);
  hPhi->Write();
  c3->cd(4);
  hImpulse->Write();

  TCanvas *c4 = new TCanvas("c4", "Invariant Mass Decay");
  c4->Divide(3, 1);
  c4->cd(1);
  hKDecay->Write();
  c4->cd(2);
  hDiff1->Write();
  c4->cd(3);
  hDiff2->Write();

  hTImpulse->Write();
  hEnergy->Write();
  hInvMass->Write();
  hSameCharge->Write();
  hOddCharge->Write();
  hPKconcordant->Write();
  hPKdiscordant->Write();

  file->Close();

  return 0;
}
