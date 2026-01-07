#ifndef PARTICLE_H
#define PARTICLE_H

#include "TObject.h"
#include "TLorentzVector.h"

class Particle : public TObject {
public:

    Int_t   Status;
    Int_t   pdgID;
    Int_t   pdgID_Mother1;
    Int_t   pdgID_Mother2;
    Int_t   pdgID_Daughter1;   
    Int_t   pdgID_Daughter2;   

    Float_t Px, Py, Pz, Energy, Mass;
    Float_t PT, Eta, Phi;

    Float_t Charge;
    TLorentzVector p4;

    Float_t Lifetime;
    Int_t   Helicity;
    
    Int_t   Color1, Color2;

    Int_t   Index;
    Int_t   Mother1;
    Int_t   Mother2;
    Int_t   Daughter1;         
    Int_t   Daughter2;         


    Particle() :
        pdgID(0), pdgID_Mother1(0), pdgID_Mother2(0), Status(0),
        Px(0), Py(0), Pz(0), Energy(0), Mass(0),
        PT(0), Eta(0), Phi(0),
        Charge(0), Mother1(0), Mother2(0),
        Color1(0), Color2(0),
        Lifetime(0), Helicity(0), Index(0),
        Daughter1(0), Daughter2(0), pdgID_Daughter1(0), pdgID_Daughter2(0)
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
