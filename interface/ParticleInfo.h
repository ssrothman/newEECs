#ifndef SRothman_EECs_ParticleInfo_h
#define SRothman_EECs_ParticleInfo_h

#include <fastjet/PseudoJet.hh>

#include <stdexcept>

namespace eec {

  // Attach with p.set_user_info(new EECParticleInfo(q)); fastjet takes ownership.
  class EECParticleInfo : public fastjet::PseudoJet::UserInfoBase {
  public:
    explicit EECParticleInfo(int charge, int pdgId = 0) : charge_(charge), pdgId_(pdgId) {}
    int charge() const { return charge_; }
    int pdgId() const { return pdgId_; }

  private:
    int charge_;
    int pdgId_;
  };

  inline int particleCharge(const fastjet::PseudoJet& p) {
    if (!p.has_user_info<EECParticleInfo>())
      throw std::invalid_argument("eec: particle has no EECParticleInfo user_info, charge unavailable");
    return p.user_info<EECParticleInfo>().charge();
  }

}  // namespace eec

#endif
