#ifndef PARTICLE_H
#define PARTICLE_H

#include "TObject.h"
#include "TLorentzVector.h"

class Particle : public TObject {
public:

    Int_t   pdgId, pdgId_Mother1, pdgId_Mother2;
    Int_t   status;

    Float_t Px, Py, Pz, Energy, Mass;
    Float_t PT, Eta, Phi;

    Int_t   Charge;
    Int_t   Mother1, Mother2;
    Int_t   Color1, Color2;

    Float_t Lifetime;
    Int_t   Helicity;

    TLorentzVector p4;

    Particle() :
        pdgId(0), pdgId_Mother1(0), pdgId_Mother2(0), status(0),
        Px(0), Py(0), Pz(0), Energy(0), Mass(0),
        PT(0), Eta(0), Phi(0),
        Charge(0), Mother1(0), Mother2(0),
        Color1(0), Color2(0),
        Lifetime(0), Helicity(0)
    {}

    const TLorentzVector& P4() const { return p4; }

    void SetP4(float px, float py, float pz, float e) {
        Px = px; Py = py; Pz = pz; Energy = e;
        p4.SetPxPyPzE(px, py, pz, e);
        PT  = p4.Pt();
        Eta = p4.Eta();
        Phi = p4.Phi();
    }

    ClassDef(Particle, 1); 
};

#endif