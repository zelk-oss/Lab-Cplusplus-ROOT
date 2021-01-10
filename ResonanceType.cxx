#include "ResonanceType.h"

ResonanceType::ResonanceType(std::string fName, double const fMass,
                             int const fCharge, double const fWidth)
    : ParticleType(fName, fMass, fCharge), fWidth_(fWidth) {}
double ResonanceType::GetWidth() const { return fWidth_; }
void ResonanceType::Print() const {
  ParticleType::Print();
  std::cout << "Particle width " << fWidth_ << '\n';
}
