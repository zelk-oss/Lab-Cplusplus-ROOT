#ifndef RESONANCETYPE_H
#define RESONANCETYPE_H

#include "ParticleType.h"

class ResonanceType : public ParticleType {
public:
ResonanceType() = default;
ResonanceType(std::string fName, double const fMass, int const fCharge, double const fWidth);
double GetWidth() const override;
void Print() const override;

private:
double const fWidth_;
};
#endif //RESONANCETYPE_H