#ifndef SimulateSiPMwithContrib_h
#define SimulateSiPMwithContrib_h 1

#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/TimeSeriesCollection.h"

#include "k4FWCore/DataHandle.h"

#include "Gaudi/Algorithm.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include <atomic>

// Check for SiPMSensor header location (similar to how DigiSiPM handles this)
#if __has_include("SiPMSensor.h")
#include "SiPMSensor.h"
#elif __has_include("sipm/SiPMSensor.h")
#include "sipm/SiPMSensor.h"
#endif

/** @class SimulateSiPMwithContrib
 *
 *  Applies the SiPM response (SimSiPM) on top of hits whose contributions
 *  carry photo-electrons (energy field = number of p.e.) and their prompt
 *  arrival time at the SiPM (time field), as produced by the IDEA_o2 SIM.
 *
 *  Used for:
 *   - HCAL fibres (Lorenzo Pezzotti), scintillation decay ~ 2.8 ns (single component)
 *   - ECAL crystals (SCEPCal): the scintillation decay is sampled PER PHOTO-ELECTRON
 *     from a sum of exponentials (e.g. BSO: 21.778 ns @ 12%, 97.680 ns @ 88%),
 *     because the SIM does not include the scintillation time.
 *
 *  Simulated SiPM effects: DCR, crosstalk, afterpulsing, cell recovery, signal shape.
 *  Any other SimSiPM property can be set through the "params" map (see below).
 *
 *  @author Lorenzo Pezzotti (original), ECAL adaptation: Julia Scamardella
 *  @date   2026
 */

class SimulateSiPMwithContrib : public Gaudi::Algorithm {
public:
  SimulateSiPMwithContrib(const std::string& name, ISvcLocator* svcLoc);
  virtual ~SimulateSiPMwithContrib() {};

  StatusCode initialize();
  StatusCode execute(const EventContext&) const;
  StatusCode finalize();

private:
  // Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  Rndm::Numbers m_rndmExp;  // Exponential(1): multiplied by the chosen tau
  Rndm::Numbers m_rndmFlat; // Flat(0,1): used to choose the decay component

  mutable k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection> m_simHits{"", Gaudi::DataHandle::Reader, this};
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_digiHits{"", Gaudi::DataHandle::Writer, this};
  mutable k4FWCore::DataHandle<edm4hep::CaloHitSimCaloHitLinkCollection> m_hitLinks{"", Gaudi::DataHandle::Writer,
                                                                                    this};

  std::unique_ptr<sipm::SiPMSensor> m_sensor;

  // cumulative probabilities of the scintillation components (built in initialize())
  std::vector<double> m_cumFractions;

  // book-keeping of p.e. falling outside the SiPM signal window (diagnostic)
  mutable std::atomic<unsigned long long> m_nPeTotal{0};
  mutable std::atomic<unsigned long long> m_nPeOutOfWindow{0};

  // ---------------- SiPM properties explicitly exposed (defaults: HCAL, Hamamatsu S14160-1310PS) ----------------
  // Signal properties
  Gaudi::Property<double> m_sigLength{this, "signalLength", 200., "signal length in ns"};
  Gaudi::Property<double> m_sampling{this, "sampling", 0.1, "SiPM sampling rate in ns"};
  Gaudi::Property<double> m_risetime{this, "risetime", 1., "signal rise time in ns"};
  Gaudi::Property<double> m_falltimeFast{this, "falltimeFast", 6.5, "signal fast component decay time in ns"};

  // SiPM physical properties
  Gaudi::Property<double> m_sipmSize{this, "SiPMsize", 6.0, "width of photosensitive area in mm"};
  Gaudi::Property<double> m_cellPitch{this, "cellpitch", 50., "SiPM cell size in um"};
  Gaudi::Property<double> m_recovery{this, "recovery", 10.,
                                     "SiPM cell recovery time in ns"}; // https://arxiv.org/abs/2001.10322

  // Noise parameters
  Gaudi::Property<double> m_Dcr{this, "DCR", 120e3, "SiPM DCR in Hz"}; // dark count rate
  Gaudi::Property<double> m_Xt{this, "Xtalk", 0.01, "SiPM crosstalk"};
  Gaudi::Property<double> m_afterpulse{this, "afterpulse", 0.03, "afterpulse probability"};
  Gaudi::Property<double> m_snr{this, "SNR", 20., "SNR value in dB"}; // signal-to-noise ratio

  // integration parameters
  Gaudi::Property<double> m_gateStart{this, "gateStart", 5., "Integration gate starting time in ns"};
  Gaudi::Property<double> m_gateL{this, "gateLength", 100.,
                                  "Integration gate length in ns"}; // Should be approx 5 times fallTimeFast (see above)
  Gaudi::Property<double> m_thres{this, "threshold", 1.5,
                                  "Integration threshold"}; // Threshold in p.e. (1.5 to suppress DCR)

  // Any other SimSiPM property, key -> value. Keys (case-insensitive) accepted by SiPMProperties::setProperty:
  //   size [mm], pitch [um], sampling, cellrecovery / recoverytime, signallength, risetime,
  //   falltimefast, falltimeslow, slowcomponentfraction, tauapfast, tauapslow,
  //   ccgv, snr, pde, dcr, xt, dxt, ap
  // ATTENTION: applied AFTER the properties above, so they override them.
  // Do NOT set "pde" here unless you want to apply PDE twice (it is already applied in the SIM).
  Gaudi::Property<std::map<std::string, double>> m_params{this, "params", {}, "optional SimSiPM parameters"};

  // switch to select if you are processing scintillation or Cherenkov hits
  Gaudi::Property<bool> m_isCherenkov{this, "isCherenkov", false,
                                      "switch to select if you are processing scintillation or Cherenkov hits"};
  // SiPM efficiency
  Gaudi::Property<std::vector<double>> m_wavelen{
      this, "wavelength", {1000., 100.}, "wavelength vector in nm (decreasing order)"};
  Gaudi::Property<std::vector<double>> m_sipmEff{this, "sipmEfficiency", {1.0, 1.0}, "SiPM efficiency vs wavelength"};

  // Scintillation time profile: sum of exponentials, sampled per photo-electron.
  // Ignored for Cherenkov. Defaults reproduce the HCAL fibres (single 2.8 ns component).
  //   BSO : scintDecayTimes = [21.778, 97.680], scintDecayFractions = [0.12, 0.88]
  //   BGO : scintDecayTimes = [51., 322.]  (fractions from the paper, Tab. 2)
  Gaudi::Property<std::vector<double>> m_scintDecayTimes{
      this, "scintDecayTimes", {2.8}, "scintillation decay time constants in ns"};
  Gaudi::Property<std::vector<double>> m_scintDecayFractions{
      this, "scintDecayFractions", {1.0}, "scintillation component fractions (normalized internally)"};

  // photo-electrons arriving later than this (absolute time, ns) are dropped: they cannot contribute to the gate
  // and would fall outside the SimSiPM signal window. <=0 means: gateStart + gateLength.
  Gaudi::Property<double> m_maxPeTime{this, "maxPeTime", -1.,
                                      "p.e. arriving after this time (ns) are discarded; <=0: gateStart+gateLength"};

  // Conversion of the SiPM integral (p.e.-equivalent units) to the energy stored in the digi hit.
  // Default 1 = no conversion (HCAL behaviour). ECAL test: 1/1965 (S), 1/97.75 (C) as in CreateOpticalCaloCells.
  // NOTE: this is a LINEAR scale; the real response saturates, so it is only a first-look calibration.
  Gaudi::Property<double> m_scaleADC{this, "scaleADC", 1., "multiplicative constant applied to the integral"};

  // Suppression threshold parameters
  Gaudi::Property<double> m_suppressionThreshold{this, "suppressionTreshold", 0.001,
                                                 "Threshold on ADC integral to reject digi hit"};
};

#endif
