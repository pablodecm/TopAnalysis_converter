
#pragma once

/**
 * Container class with derived high level features for mva training
 */
class EventHighLevel {
  public:
    float met;
    float dilept_inv_mass;
    float dilepton_MT2;
    float d_phi_l_l;
    float d_phi_min_met_j;
    float d_phi_met_l0;
    float d_phi_met_l1;
    float d_phi_met_metll;

    EventHighLevel() {}
    virtual ~EventHighLevel() {}

    ClassDef(EventHighLevel, 1);

};
