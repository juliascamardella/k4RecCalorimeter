#include "SimulateSiPMwithContrib.h"
#include "DD4hep/DD4hepUnits.h"
#include <algorithm>
#include <cmath>
#include <numeric>

DECLARE_COMPONENT(SimulateSiPMwithContrib)

SimulateSiPMwithContrib::SimulateSiPMwithContrib(const std::string& aName, ISvcLocator* aSvcLoc)
    : Gaudi::Algorithm(aName, aSvcLoc) {

  // register DataHandles with hit collection names
  declareProperty("inputHitCollection", m_simHits, "input calo collection name");
  declareProperty("outputHitCollection", m_digiHits, "output calo collection name");
  declareProperty("outputLinkCollection", m_hitLinks, "output links collection name");
}

StatusCode SimulateSiPMwithContrib::initialize() {
  StatusCode sc = Gaudi::Algorithm::initialize();

  if (sc.isFailure())
    return sc;

  // Initialize random services
  m_randSvc = service("RndmGenSvc", false);

  if (!m_randSvc) {
    error() << "Couldn't get RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_rndmExp.initialize(m_randSvc, Rndm::Exponential(1.0)).isFailure() ||
      m_rndmFlat.initialize(m_randSvc, Rndm::Flat(0., 1.)).isFailure()) {
    error() << "Couldn't initialize RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }

  // scintillation time profile: check and build cumulative distribution of the components
  const auto& taus = m_scintDecayTimes.value();
  const auto& fracs = m_scintDecayFractions.value();
  if (taus.empty() || taus.size() != fracs.size()) {
    error() << "scintDecayTimes and scintDecayFractions must be non-empty and have the same size" << endmsg;
    return StatusCode::FAILURE;
  }
  const double sumFrac = std::accumulate(fracs.begin(), fracs.end(), 0.);
  if (sumFrac <= 0.) {
    error() << "scintDecayFractions must sum to a positive number" << endmsg;
    return StatusCode::FAILURE;
  }
  m_cumFractions.clear();
  double acc = 0.;
  for (double f : fracs) {
    acc += f / sumFrac;
    m_cumFractions.push_back(acc);
  }
  m_cumFractions.back() = 1.; // protect against round-off
  if (!m_isCherenkov) {
    for (size_t i = 0; i < taus.size(); ++i)
      info() << "Scintillation component " << i << ": tau = " << taus[i] << " ns, fraction = "
             << fracs[i] / sumFrac << endmsg;
  }

  // initialize SiPM properties
  sipm::SiPMProperties properties;
  properties.setSignalLength(m_sigLength);
  properties.setSize(m_sipmSize);
  properties.setDcr(m_Dcr);
  properties.setXt(m_Xt);
  properties.setSampling(m_sampling);
  properties.setRecoveryTime(m_recovery);
  properties.setPitch(m_cellPitch);
  properties.setAp(m_afterpulse);
  properties.setFallTimeFast(m_falltimeFast);
  properties.setRiseTime(m_risetime);
  properties.setSnr(m_snr);
  properties.setPdeType(sipm::SiPMProperties::PdeType::kSimplePde); // this sets PDE to a single number
  properties.setPde(1.0); // PDE (and geometrical efficiency) are already applied in the SIM
  // properties.setPdeSpectrum(m_wavelen, m_sipmEff); // if you want wavelen dependency

  // set other parameters if provided
  for (const auto& [key, value] : m_params.value()) {
    if (key == "pde" || key == "PDE" || key == "Pde")
      warning() << "'pde' set in params: PDE is already applied in the SIM, it will be applied twice!" << endmsg;
    properties.setProperty(key, value);
  }

  // consistency of the time windows
  const double maxT = (m_maxPeTime.value() > 0.) ? m_maxPeTime.value() : m_gateStart.value() + m_gateL.value();
  if (maxT >= m_sigLength.value())
    warning() << "Latest p.e. time considered (" << maxT << " ns) is >= signalLength (" << m_sigLength.value()
              << " ns): photo-electrons beyond the SimSiPM window will be discarded anyway. Increase signalLength."
              << endmsg;

  m_sensor = std::make_unique<sipm::SiPMSensor>(properties); // must be constructed from SiPMProperties

  info() << "SimulateSiPMwithContrib initialized (isCherenkov = " << m_isCherenkov.value() << ")" << endmsg;
  info() << properties << endmsg; // sipm::SiPMProperties has std::ostream& operator<<

  return StatusCode::SUCCESS;
}

StatusCode SimulateSiPMwithContrib::execute(const EventContext&) const {
  const edm4hep::SimCalorimeterHitCollection* scintHits = m_simHits.get();
  edm4hep::CalorimeterHitCollection* digiHits = m_digiHits.createAndPut();
  auto* links = m_hitLinks.createAndPut();

  const auto& taus = m_scintDecayTimes.value();
  const double maxT =
      std::min((m_maxPeTime.value() > 0.) ? m_maxPeTime.value() : m_gateStart.value() + m_gateL.value(),
               m_sigLength.value());

  // loop through the hits (each hit corresponds to a crystal / fibre)
  for (unsigned int idx = 0; idx < scintHits->size(); idx++) {
    const auto& scintHit = scintHits->at(idx);

    if (!scintHit.isAvailable()) {
      std::cout << "ERROR: Hit not available!" << std::endl;
      continue;
    }

    std::vector<double> vecTimes;
    std::vector<double> vecWavelens;

    unsigned long long nTot = 0, nOut = 0;

    // loop through the hit contributions
    for (auto contrib = scintHit.contributions_begin(); contrib != scintHit.contributions_end(); ++contrib) {

      if (!contrib->isAvailable()) {
        std::cerr << "ERROR: Hit contribution not available!" << std::endl;
        continue;
      }
      // number of photo-electrons: stored as n x GeV in the SIM, converted back to n by the EDM4hep converter
      // (rounded, to avoid losing a p.e. with a floating-point 2.9999)
      const long long npe = std::llround(contrib->getEnergy());
      if (npe <= 0)
        continue;

      // prompt photon arrival time at SiPM (ns), binned in the SIM
      const double time = contrib->getTime();

      vecTimes.reserve(vecTimes.size() + npe);
      vecWavelens.reserve(vecWavelens.size() + npe);

      for (long long pe = 0; pe < npe; pe++) {
        double thisTime = time;
        // scintillation: sample the emission time for EACH p.e. from the multi-exponential profile
        if (!m_isCherenkov) {
          const double u = m_rndmFlat.shoot();
          size_t k = 0;
          while (k + 1 < m_cumFractions.size() && u > m_cumFractions[k])
            ++k;
          thisTime += taus[k] * m_rndmExp.shoot();
        }
        ++nTot;
        if (thisTime >= maxT) { // cannot contribute to the gate, and is outside the SimSiPM window
          ++nOut;
          continue;
        }
        vecTimes.push_back(thisTime);
        vecWavelens.push_back(1.); // photons wavelengths are dummy, we are using photo-electrons
      }
    } // contrib

    m_nPeTotal += nTot;
    m_nPeOutOfWindow += nOut;

    if (vecTimes.empty())
      continue;

    // sort time of arrivals from shortest to longest
    std::sort(vecTimes.begin(), vecTimes.end());

    m_sensor->resetState();
    m_sensor->addPhotons(vecTimes, vecWavelens); // Sets photon times & wavelengths
    m_sensor->runEvent();                        // Runs the simulation

    // Using only analog signal (ADC conversion is still experimental)
    const sipm::SiPMAnalogSignal anaSignal = m_sensor->signal();

    // if the signal never exceeds the threshold, it will return -1.
    const double integral =
        std::max(0., anaSignal.integral(m_gateStart, m_gateL, m_thres));           // (intStart, intGate, threshold)
    const double toa = std::max(0., anaSignal.toa(m_gateStart, m_gateL, m_thres)); // (intStart, intGate, threshold)

    // Apply threshold to save digi hit.
    // By default this is a zero suppression cut, can the changed to apply 1-suppression, 2-suppression, ...
    // Note: check the sipm gain applied in SimSipm, by default it is 1.
    if (integral < m_suppressionThreshold)
      continue;

    // NOTE: I noticed very few digi hits have TOA > ms, clearly a bug to be investigated
    if (toa > 1e6)
      continue;

    auto digiHit = digiHits->create();
    digiHit.setEnergy(integral * m_scaleADC.value()); // conversion from ADC to GeV
    digiHit.setEnergyError(m_scaleADC.value() * std::sqrt(integral));
    digiHit.setPosition(scintHit.getPosition());
    digiHit.setCellID(scintHit.getCellID());
    // Toa and m_gateStart are in ns
    digiHit.setTime(toa + m_gateStart);

    auto link = links->create();
    link.setFrom(digiHit);
    link.setTo(scintHit);
  }

  return StatusCode::SUCCESS;
}

StatusCode SimulateSiPMwithContrib::finalize() {
  const auto tot = m_nPeTotal.load();
  if (tot > 0)
    info() << "Photo-electrons discarded because later than " << (m_maxPeTime.value() > 0. ? m_maxPeTime.value() : m_gateStart.value() + m_gateL.value())
           << " ns: " << m_nPeOutOfWindow.load() << " / " << tot << " ("
           << 100. * double(m_nPeOutOfWindow.load()) / double(tot) << " %)" << endmsg;
  return Gaudi::Algorithm::finalize();
}
