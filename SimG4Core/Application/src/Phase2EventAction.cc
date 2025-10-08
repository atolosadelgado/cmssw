#include "SimG4Core/Application/interface/Phase2EventAction.h"
#include "SimG4Core/Application/interface/SimRunInterface.h"
#include "SimG4Core/Notification/interface/TmpSimEvent.h"
#include "SimG4Core/Notification/interface/TmpSimVertex.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/EndOfEvent.h"
#include "SimG4Core/Notification/interface/CMSSteppingVerbose.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Randomize.hh"

Phase2EventAction::Phase2EventAction(const edm::ParameterSet& p,
                                     SimRunInterface* rm,
                                     SimTrackManager* iManager,
                                     CMSSteppingVerbose* sv)
    : m_runInterface(rm),
      m_trackManager(iManager),
      m_SteppingVerbose(sv),
      m_stopFile(p.getParameter<std::string>("StopFile")),
      m_printRandom(p.getParameter<bool>("PrintRandomSeed")),
      m_debug(p.getUntrackedParameter<bool>("debug", false))
      {
        // xmin/max, binsize in mm
        double xmin = 3000;
        double xmax = 6000;
        double binsize = 0.05;
        // bin of 50um
        int nbins( (xmax-xmin)/binsize);
        hEhgcal_average = std::make_unique<TH1D>("hEhgcal_average","hEhgcal_average;z (mm);E (MeV)", nbins,xmin,xmax);
        hEhgcal_1evt    = std::make_unique<TH1D>("hEhgcal_average","hEhgcal_average;z (mm);E (MeV)", nbins,xmin,xmax);

        hNprimaries = std::make_unique<TH1D>("hNprimaries","hNprimaries", 10000,0,10000);
        hNprimaries->SetDirectory(0);

        hNprimariesTime = std::make_unique<TH2F>("hNprimariesTime","hNprimariesTime", 200,0,1000,200,0,100000);

      }

void Phase2EventAction::BeginOfEventAction(const G4Event* anEvent) {
  // Alvaro histograms
  hEhgcal_1evt->Reset();
  nevents++;
  start_time = std::chrono::steady_clock::now();

  BeginOfEvent e(anEvent);
  m_beginOfEventSignal(&e);

  if (m_printRandom) {
    edm::LogVerbatim("SimG4CoreApplication")
        << "BeginOfEvent " << anEvent->GetEventID() << " Random number: " << G4UniformRand();
  }

  if (nullptr != m_SteppingVerbose) {
    m_SteppingVerbose->beginOfEvent(anEvent);
  }

}

void Phase2EventAction::EndOfEventAction(const G4Event* anEvent) {


  // Alvaro histograms,
  int nPrimaries = 0;
  for (int iv = 0; iv < anEvent->GetNumberOfPrimaryVertex(); ++iv) {
      nPrimaries += anEvent->GetPrimaryVertex(iv)->GetNumberOfParticle();
  }
  hNprimaries->Fill(nPrimaries);
  auto end_time = std::chrono::steady_clock::now();
  double duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
  hNprimariesTime->Fill(nPrimaries,duration_ms);


  if (1 == nevents) {
    // just copy the event
    hEhgcal_average->Add(hEhgcal_1evt.get());
  }
  else {
    int nb = hEhgcal_average->GetNbinsX();
    for (int ib = 1; ib <= nb; ++ib) {
        double x = hEhgcal_1evt->GetBinContent(ib);
        double mean_prev = hEhgcal_average->GetBinContent(ib);
        double M2_prev = 0.0;
        double err_prev = hEhgcal_average->GetBinError(ib); // actualmente contiene sqrt(M2_prev) o 0
        if (err_prev > 0.0) M2_prev = err_prev * err_prev;

        double delta = x - mean_prev;
        double mean_new = mean_prev + delta / nevents;
        double M2_new = M2_prev + delta * (x - mean_new);

        hEhgcal_average->SetBinContent(ib, mean_new);
        hEhgcal_average->SetBinError(ib, std::sqrt(M2_new)); // guardamos sqrt(M2) temporalmente
    }
  }

  int ThreadIndex = m_runInterface->getThreadIndex();
  std::string ofilename = "test" + std::to_string(ThreadIndex) + ".root";
  TFile * ofile = new TFile( ofilename.c_str() ,"recreate");
  TH1D * hEhgcal_averaged_fromofile;
  // TH1D * hEhgcal_averaged_fromofile = (TH1D*) ofile->Get("hEhgcal");
  // create histogram if not existing
  // if(nullptr == hEhgcal_averaged_fromofile)
  {
    hEhgcal_averaged_fromofile = (TH1D*) hEhgcal_average->Clone("hEhgcal");
    hEhgcal_averaged_fromofile->SetDirectory(ofile);
    hEhgcal_averaged_fromofile->Reset();
  }
  // finalize Welford incremental variance calculation
  if (nevents > 1) {
    int nb = hEhgcal_average->GetNbinsX();
    for (int ib = 1; ib <= nb; ++ib) {
      // overwrite mean with mean (stays)
      hEhgcal_averaged_fromofile->SetBinContent(ib, hEhgcal_average->GetBinContent(ib));
      // calculate the mean error from M2
      double sqrtM2 = hEhgcal_average->GetBinError(ib); // this is sqrt(M2)
      double M2 = sqrtM2 * sqrtM2;
      double stddev = (nevents > 1) ? std::sqrt(M2 / (nevents - 1)) : 0.0;
      double err_mean = (nevents > 0) ? stddev / std::sqrt(nevents) : 0.0;
      // overwrite error with mean error
      hEhgcal_averaged_fromofile->SetBinError(ib, err_mean);
    }
  }
  // for(int i = 0; i<= hEhgcal_average->GetNbinsX(); ++i)
  // {
  //   // hEhgcal_averaged_fromofile->SetBinContent(i, hEhgcal_average->GetBinContent(i));
  //   // hEhgcal_averaged_fromofile->SetBinError(i, hEhgcal_average->GetBinError(i));
  //   hEhgcal_averaged_fromofile->SetBinContent(i, hEhgcal_1evt->GetBinContent(i));
  //   hEhgcal_averaged_fromofile->SetBinError(i, hEhgcal_1evt->GetBinError(i));
  // }

  hEhgcal_average->Write();
  hNprimaries->Write();
  hNprimariesTime->Write();

  ofile->Close();
  std::cout << "\tdabadaba" << std::endl;



  ///---------------
    if (m_printRandom) {
    edm::LogVerbatim("SimG4CoreApplication")
        << "Phase2EventAction::EndOfEventAction: " << anEvent->GetEventID() << " Random number: " << G4UniformRand();
  }
  if (!m_stopFile.empty() && std::ifstream(m_stopFile.c_str())) {
    edm::LogWarning("SimG4CoreApplication")
        << "Phase2EventAction::EndOfEventAction: termination signal received at event " << anEvent->GetEventID();
    // soft abort run
    m_runInterface->abortRun(true);
  }
  if (anEvent->GetNumberOfPrimaryVertex() == 0) {
    edm::LogWarning("SimG4CoreApplication") << "Phase2EventACtion::EndOfEventAction: event " << anEvent->GetEventID()
                                            << " must have failed (no G4PrimaryVertices found) and will be skipped";
    return;
  }

  m_trackManager->storeTracks();

  // dispatch now end of event
  EndOfEvent e(anEvent);
  m_endOfEventSignal(&e);

  // delete transient objects
  m_trackManager->reset();
}

void Phase2EventAction::abortEvent() { m_runInterface->abortEvent(); }

void Phase2EventAction::Update_HGCaleprofile(double zpos_mm, double edep_MeV)
{
  hEhgcal_1evt->Fill(zpos_mm, edep_MeV);
  // std::cout << "\t zpos " << zpos_mm << "\t edep " << edep_MeV << std::endl;
}

