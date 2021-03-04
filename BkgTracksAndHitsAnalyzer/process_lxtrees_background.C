//
//

#ifndef __RUN_PROC_HITS_TREE__

void process_lxtrees_background(const char *fnlist = 0, const char *commentstr = 0)
{
   gROOT->ProcessLineSync("#define __RUN_PROC_HITS_TREE__ 1");
   gROOT->ProcessLineSync(".L MHists.C+");
   gROOT->ProcessLineSync(".L ProcessLxSim.C+");
   gROOT->ProcessLine("#include \"process_lxtrees_background.C\"");
   gROOT->ProcessLine(Form("process_hits_tree_draw(\"%s\")", fnlist));
   gROOT->ProcessLine("#undef __RUN_PROC_HITS_TREE__");
}

#else

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <math.h>
#include <chrono>
#include <vector>

#include "TChain.h"

#include "MHists.h"


int ProcessList(const std::string &fnamelist, std::vector<std::string> &flist);
void CreateHistograms(MHists *mh);


int process_hits_tree_draw(const char *fnlist = 0, const char *commentstr = 0)
{

  int debugl = 0; //1;

  if (!fnlist) {
    std::cout << "Usage: root -l process_lxtrees.C\'(\"file_with_list_of_mc_files\")\'\n";
    return -1;
  }
  
  std::string fnamelist(fnlist);
  std::vector<std::string>  flist;
  ProcessList(fnamelist, flist);  
  if (debugl) {
    std::cout << "The following files will be processed:\n";
    std::for_each(flist.begin(), flist.end(), [](const std::string ss) {std::cout << ss << std::endl;});
  }

  MHists *lhist = new MHists();
  CreateHistograms(lhist);

  ProcessLxSim *lxsim = new ProcessLxSim(flist);
  
  std::vector<int> nprimary(3,0);
  int evtoread = 50;
  int nev = 0;
  int nrec = 1;
//   while (nrec && nev < evtoread) {
  while (nrec) {
    nrec = lxsim->ReadNextEvent();
    if (!nrec) continue;
    if (!(nev % 10000)) { std::cout << "Event " << nev << std::endl; }
    
//     lxsim->DumpEvent();

    const auto &tracksv = lxsim->GetTracks();
    const auto voltrck = tracksv.at(0);
    const std::vector<int> &detidv = voltrck.detid;
    if(nev > 1640000 && false)cout << "After detidv " << std::endl;
//     const auto itr = std::find (detidv.begin(), detidv.end(), -1);
//     int primndx = -1;
//     if (itr == detidv.end()) {
//       std::cout << "Warning! Primary track is not found!\n";
//       throw std::logic_error("Primary track is not found!");
//     } else {
//       primndx = itr - detidv.begin();
//     }

//     std::cout << "Event " << voltrck.eventid << std::endl;
//     std::cout << "Primary: detid: " << detidv.at(primndx) << "  index: " << primndx << std::endl;
    int primpdg = 11;
//     int primpdg = voltrck.pdg.at(primndx);
    double primary_energy = 16.5;
//     double primary_energy = voltrck.E.at(primndx);
//     std::cout << "Primary: particle: " << primpdg << "Energ:y " << primary_energy << std::endl;
    int phid = 0;
    switch (primpdg) {
      case -11 : phid = 2; break; 
      case  22 : phid = 1; break;
    }
//     lhist->FillHistW("primary_e", phid, voltrck.E.at(primndx), voltrck.weight);
    nprimary[phid] += 1;

    const auto &hitv = lxsim->GetHits();
    const auto &hittrckv = lxsim->GetHitTracks();
    const auto &htrcks = hittrckv.at(0);

    const int ntracklayers = 16;  
    std::vector<int> ntracker_hits_signal(ntracklayers, 0);
    std::vector<int> ntracker_hits_backgr(ntracklayers, 0);

    const int necallayers = 21;  
    std::vector<int> necal_hits_signal(necallayers, 0);
    std::vector<int> necal_hits_backgr(necallayers, 0);
    double showeren = 0.0;
    double showerx = 0.0;
    std::map<int, std::vector<int> > track_hit_map;
    
    if(nev > 1640000 && false)cout << "After track hit map : nev " << nev << std::endl;
    
    // loop on hits
    for (const auto &hit : hitv) {
      const auto &tridlist = hit.track_list;
      const auto &trxv = hit.trackx;
      const auto &tryv = hit.tracky;
      const auto &trzv = hit.trackz;
      const auto &trtv = hit.trackt;
      const auto &trackedepv = hit.trackedep;
      
      if(nev > 1640000 && false)cout << "Inside hit loop : nev " << nev << std::endl;
      int det_id = hit.detid;
      double ev_weight = hit.weight;
      
      if (det_id >= 1000 && det_id <= 1008) {
        int det = det_id - 1000;
        int layer_id = hit.layerid;
        int detlayer = det * ntracklayers + layer_id;
        lhist->FillHistW("tracking_planes_hits_x", detlayer, hit.cellx, ev_weight);
        lhist->FillHistW("tracking_planes_hits_y", detlayer, hit.celly, ev_weight);
        lhist->FillHistW("tracking_planes_hits_xy", detlayer, hit.cellx, hit.celly, ev_weight);
        lhist->FillHistW("tracking_planes_hits_xy_edep", detlayer, hit.cellx, hit.celly, hit.edep*ev_weight);
        lhist->FillHistW("tracking_planes_hits_edep_log", detlayer, hit.edep, ev_weight); // new histogram
        
        lhist->FillHistW("tracking_planes_n_tracks_per_hit", layer_id, tridlist.size(), ev_weight);
        if (primpdg==-11) {
          lhist->FillHistW("tracking_planes_n_tracks_per_hit_signal", layer_id, tridlist.size(), ev_weight);
          ntracker_hits_signal[layer_id] += 1;
        } else {
          lhist->FillHistW("tracking_planes_n_tracks_per_hit_background", layer_id, tridlist.size(), ev_weight);
          ntracker_hits_backgr[layer_id] += 1;
        }
        if(nev > 1640000 && false)cout << "Inside hit and detid loop : nev " << nev << std::endl;
        std::vector<int> ntracks_per_hits_signal(ntracklayers, 0);
        std::vector<int> ntracks_per_hits_backgr(ntracklayers, 0);
        std::vector<int> ntracks_per_hits_signal_ecut(ntracklayers, 0);
        std::vector<int> ntracks_per_hits_backgr_ecut(ntracklayers, 0);
        std::vector<int> ntracks_per_hits_signal_zcut(ntracklayers, 0);
        std::vector<int> ntracks_per_hits_backgr_zcut(ntracklayers, 0);
        /// count the tracks after the vertex cut
        int trackNumber = 0;
        //loop on hit tracks
        for (int trkid : tridlist) {
            
          if(tridlist.size()!=trackedepv.size())continue;
          if(tridlist.size()!=trtv.size())continue;
          if(tridlist.size()!=htrcks.pproc.size())continue;
          
          
              
          const auto titr = std::find(htrcks.trackid.begin(), htrcks.trackid.end(), trkid);
          if (titr == htrcks.trackid.end()) {
            throw std::logic_error("Track not found in HitTracks tree!"); 
          }
          
          
          if(nev > 1640000 && false)cout << "Inside hit and detid and tracks loop : nev " << nev << std::endl;
          int vndx = titr - htrcks.trackid.begin();
//           std::cout << "track edep: " << trackedepv.at(vndx) << std::endl;
//           std::cout << "track time: " << trtv.at(vndx) << std::endl;
          lhist->FillHistW("tracking_planes_hit_track_e", layer_id, htrcks.E.at(vndx), ev_weight);
          lhist->FillHistW("tracking_planes_hit_track_proc", layer_id, htrcks.pproc.at(vndx), ev_weight);
          lhist->FillHistW("tracking_planes_hit_track_pdg", layer_id, htrcks.pdg.at(vndx), ev_weight);
          lhist->FillHistW("tracking_planes_hit_track_edep_all_log", layer_id, trackedepv.at(vndx), ev_weight);/// new histogram
          if(nev > 1640000 && false)cout << "Inside hit and detid and tracks loop after vector : nev " << nev << std::endl;
          
          if(htrcks.pdg.at(vndx) == 2112){
            lhist->FillHistW("tracking_planes_hit_track_edep_neutron_log", layer_id, trackedepv.at(vndx), ev_weight);/// new histogram
//             lhist->FillHistW("tracking_planes_hit_track_time_neutron", layer_id, trtv.at(vndx), ev_weight);/// new histogram
//             lhist->FillHistW("tracking_planes_hit_track_proc_neutron", layer_id, htrcks.pproc.at(vndx), ev_weight);
            lhist->FillHistW("tracking_planes_hits_xy_edep_neutron", detlayer, hit.cellx, hit.celly, trackedepv.at(vndx)*ev_weight);
          }
          
          if(nev > 1640000 && false)cout << "Inside hit and detid and tracks loop after neutron : nev " << nev << std::endl;
          
          if(htrcks.pdg.at(vndx) == 2212){
            lhist->FillHistW("tracking_planes_hit_track_edep_proton_log", layer_id, trackedepv.at(vndx), ev_weight);/// new histogram
//             lhist->FillHistW("tracking_planes_hit_track_time_proton", layer_id, trtv.at(vndx), ev_weight);/// new histogram
//             lhist->FillHistW("tracking_planes_hit_track_proc_proton", layer_id, htrcks.pproc.at(vndx), ev_weight);
            lhist->FillHistW("tracking_planes_hits_xy_edep_proton", detlayer, hit.cellx, hit.celly, trackedepv.at(vndx)*ev_weight);
          }
          
          if(nev > 1640000 && false)cout << "Inside hit and detid and tracks loop after proton : nev " << nev << std::endl;
          
          if(htrcks.pdg.at(vndx) == 111){
            lhist->FillHistW("tracking_planes_hit_track_edep_pi0_log", layer_id, trackedepv.at(vndx), ev_weight);/// new histogram
//             lhist->FillHistW("tracking_planes_hit_track_time_pi0", layer_id, trtv.at(vndx), ev_weight);/// new histogram
//             lhist->FillHistW("tracking_planes_hit_track_proc_pi0", layer_id, htrcks.pproc.at(vndx), ev_weight);
            lhist->FillHistW("tracking_planes_hits_xy_edep_pi0", detlayer, hit.cellx, hit.celly, trackedepv.at(vndx)*ev_weight);
          }
          
          if(nev > 1640000 && false)cout << "Inside hit and detid and tracks loop after pion : nev " << nev << std::endl;
          
          if(abs(htrcks.pdg.at(vndx)) == 211){
            lhist->FillHistW("tracking_planes_hit_track_edep_pi_log", layer_id, trackedepv.at(vndx), ev_weight);/// new histogram
//             lhist->FillHistW("tracking_planes_hit_track_time_pi", layer_id, trtv.at(vndx), ev_weight);/// new histogram
//             lhist->FillHistW("tracking_planes_hit_track_proc_pi", layer_id, htrcks.pproc.at(vndx), ev_weight);
            lhist->FillHistW("tracking_planes_hits_xy_edep_pi", detlayer, hit.cellx, hit.celly, trackedepv.at(vndx)*ev_weight);
          }
          
          if(nev > 1640000 && false)cout << "Inside hit and detid and tracks loop after pi+/- : nev " << nev << std::endl;
          
          if(abs(htrcks.pdg.at(vndx)) == 13){
            lhist->FillHistW("tracking_planes_hit_track_edep_muon_log", layer_id, trackedepv.at(vndx), ev_weight);/// new histogram
//             lhist->FillHistW("tracking_planes_hit_track_time_muon", layer_id, trtv.at(vndx), ev_weight);/// new histogram
//             lhist->FillHistW("tracking_planes_hit_track_proc_muon", layer_id, htrcks.pproc.at(vndx), ev_weight);
            lhist->FillHistW("tracking_planes_hits_xy_edep_muon", detlayer, hit.cellx, hit.celly, trackedepv.at(vndx)*ev_weight);
          }
          
          if(nev > 1640000 && false)cout << "Inside hit and detid and tracks loop after muon : nev " << nev << std::endl;
          
//           /// remove the hotspot region
//           //if(!((htrcks.vtxz.at(vndx) > 3800 && htrcks.vtxz.at(vndx) < 6500) && htrcks.vtxx.at(vndx) < -25))
//           trackNumber++;
//           
//           if ( abs(htrcks.vtxz.at(vndx)-trzv[0]) > 0.025 ) {
//             lhist->FillHistW("tracking_planes_hit_track_e_zcut", layer_id, htrcks.E.at(vndx), ev_weight);
//             lhist->FillHistW("tracking_planes_hit_track_pdg_zcut", layer_id, htrcks.pdg.at(vndx), ev_weight);
//             track_hit_map[htrcks.trackid.at(vndx)].push_back(hit.hitid);
//             if (primpdg==-11) {
//               lhist->FillHistW("tracking_planes_hit_track_e_signal_zcut", layer_id, htrcks.E.at(vndx), ev_weight);
//               lhist->FillHistW("tracking_planes_hit_track_pdg_signal_zcut", layer_id, htrcks.pdg.at(vndx), ev_weight);
//               lhist->FillHistW("tracking_planes_hit_track_id_signal_zcut", layer_id, htrcks.trackid.at(vndx), ev_weight);
//             } else {
//               lhist->FillHistW("tracking_planes_hit_track_e_background_zcut", layer_id, htrcks.E.at(vndx), ev_weight);
//               lhist->FillHistW("tracking_planes_hit_track_pdg_background_zcut", layer_id, htrcks.pdg.at(vndx), ev_weight);
//               lhist->FillHistW("tracking_planes_hit_track_id_background_zcut", layer_id, htrcks.trackid.at(vndx), ev_weight);
//             }
//           }
          
        }
//         /// fill the histogram for number of tracks after vertex cut
//         lhist->FillHistW("tracking_planes_n_tracks_per_hit_aftervtxcut", detlayer, trackNumber, ev_weight);
//         if(trackNumber!=0){
//             lhist->FillHistW("tracking_planes_hits_x_aftervtxcut", detlayer, hit.cellx, ev_weight);
//             lhist->FillHistW("tracking_planes_hits_y_aftervtxcut", detlayer, hit.celly, ev_weight);
//             lhist->FillHistW("tracking_planes_hits_xy_aftervtxcut", detlayer, hit.cellx, hit.celly, ev_weight);
//             lhist->FillHistW("tracking_planes_hits_xy_edep_aftervtxcut", detlayer, hit.cellx, hit.celly, hit.edep*ev_weight);
//         }
      }
      
//       //ECAL
//       if (det_id >= 2000 && det_id <= 2001) {
//         const int nlayers = necallayers;  
//         int det = det_id - 2000;
//         int layer_id = hit.layerid;
//         int detlayer = det * nlayers + layer_id;
//         double ecalw = 1.0;
// 
//         if (primpdg==-11) {
//           lhist->FillHistW("ecal_hit_edep_signal", layer_id, hit.edep, ecalw);
//           lhist->FillHistW("ecal_hit_z_signal", det, layer_id, ecalw);
//           lhist->FillHistW("ecal_hit_z_signal_edep", det, layer_id, hit.edep*ecalw);
//           lhist->FillHistW("ecal_hit_zx_signal_edep", det, layer_id, hit.cellx, hit.edep*ecalw);
//           lhist->FillHistW("ecal_hit_xy_signal_edep", det, hit.cellx, hit.celly, hit.edep*ecalw);
//           lhist->FillHistW("ecal_hit_xy_layer_signal_edep", detlayer, hit.cellx, hit.celly, hit.edep*ecalw);
//           lhist->FillHist("ecal_hit_signal_weight", 0, ev_weight);
//           showeren += hit.edep;
//           showerx += hit.edep * hit.cellx;
//           necal_hits_signal[layer_id] += 1;
//         } else {
//           lhist->FillHistW("ecal_hit_edep_background", layer_id, hit.edep, ev_weight);
//           lhist->FillHistW("ecal_hit_z_background", det, layer_id, ev_weight);
//           lhist->FillHistW("ecal_hit_z_background_edep", det, layer_id, hit.edep*ev_weight);
//           lhist->FillHistW("ecal_hit_zx_background_edep", det, layer_id, hit.cellx, hit.edep*ev_weight);
//           lhist->FillHistW("ecal_hit_xy_background_edep", det, hit.cellx, hit.celly, hit.edep*ev_weight);
//           lhist->FillHistW("ecal_hit_xy_layer_background_edep", detlayer, hit.cellx, hit.celly, hit.edep*ev_weight);
//           necal_hits_backgr[layer_id] += 1;
//           if ((primpdg==11 && primary_energy > 16.4)) {
//             lhist->FillHistW("ecal_hit_edep_background_dump", layer_id, hit.edep, ev_weight);
//             lhist->FillHistW("ecal_hit_z_background_dump", det, layer_id, ev_weight);
//             lhist->FillHistW("ecal_hit_z_background_dump_edep", det, layer_id, hit.edep*ev_weight);
//             lhist->FillHistW("ecal_hit_zx_background_dump_edep", det, layer_id, hit.cellx, hit.edep*ev_weight);
//             lhist->FillHistW("ecal_hit_xy_background_dump_edep", det, hit.cellx, hit.celly, hit.edep*ev_weight);
//             lhist->FillHistW("ecal_hit_xy_layer_background_dump_edep", detlayer, hit.cellx, hit.celly, hit.edep*ev_weight);
//           }
//           
//           if ((primpdg==11 && primary_energy <= 16.4)) {
//             lhist->FillHistW("ecal_hit_edep_background_el_low", layer_id, hit.edep, ev_weight);
//             lhist->FillHistW("ecal_hit_z_background_el_low", det, layer_id, ev_weight);
//             lhist->FillHistW("ecal_hit_z_background_el_low_edep", det, layer_id, hit.edep*ev_weight);
//             lhist->FillHistW("ecal_hit_zx_background_el_low_edep", det, layer_id, hit.cellx, hit.edep*ev_weight);
//             lhist->FillHistW("ecal_hit_xy_background_el_low_edep", det, hit.cellx, hit.celly, hit.edep*ev_weight);
//             lhist->FillHistW("ecal_hit_xy_layer_background_el_low_edep", detlayer, hit.cellx, hit.celly, hit.edep*ev_weight);
//           }
// 
//           if (primpdg==22) {
//             lhist->FillHistW("ecal_hit_edep_background_gamma", layer_id, hit.edep, ev_weight);
//             lhist->FillHistW("ecal_hit_z_background_gamma", det, layer_id, ev_weight);
//             lhist->FillHistW("ecal_hit_z_background_gamma_edep", det, layer_id, hit.edep*ev_weight);
//             lhist->FillHistW("ecal_hit_zx_background_gamma_edep", det, layer_id, hit.cellx, hit.edep*ev_weight);
//             lhist->FillHistW("ecal_hit_xy_background_gamma_edep", det, hit.cellx, hit.celly, hit.edep*ev_weight);
//             lhist->FillHistW("ecal_hit_xy_layer_background_gamma_edep", detlayer, hit.cellx, hit.celly, hit.edep*ev_weight);
//           }
//         }
//       }
//       
//       if (det_id == 5000) { // Compton electron scintillator.
//         int det = det_id - 5000;
//         if (primpdg==11) {
//           const auto tlitr = std::find(tridlist.begin(), tridlist.end(), 1);
//           if (tlitr != tridlist.end()) {
//             lhist->FillHistW("ce_screen_hit_edep_signal_prime", det, hit.edep, ev_weight);
//             lhist->FillHistW("ce_screen_hit_xy_signal_prime_edep", det, hit.cellx, hit.celly, hit.edep*ev_weight);
//             lhist->FillHistW("ce_screen_hit_n_tracks_signal_prime", det, tridlist.size(), ev_weight);
//           } else {
//             lhist->FillHistW("ce_screen_hit_edep_signal_other", det, hit.edep, ev_weight);
//             lhist->FillHistW("ce_screen_hit_xy_signal_other_edep", det, hit.cellx, hit.celly, hit.edep*ev_weight);
//             lhist->FillHistW("ce_screen_hit_n_tracks_signal_other", det, tridlist.size(), ev_weight);
//           }
//         } else {
//           lhist->FillHistW("ce_screen_hit_edep_background", det, hit.edep, ev_weight);
//           lhist->FillHistW("ce_screen_hit_xy_background_edep", det, hit.cellx, hit.celly, hit.edep*ev_weight);
//           lhist->FillHistW("ce_screen_hit_n_tracks_background", det, tridlist.size(), ev_weight);
//         }
//       }
      
    } //hits
    
//     for (int il = 0; il < ntracklayers; ++il) {
//       if (primpdg==-11) lhist->FillHistW("tracker_n_hits_signal", il, ntracker_hits_signal[il], voltrck.weight);
//       else lhist->FillHistW("tracker_n_hits_background", il, ntracker_hits_backgr[il], voltrck.weight);
//     }
//     for (int il = 0; il < necallayers; ++il) {
//       if (primpdg==-11) lhist->FillHistW("ecal_n_hits_signal", il, necal_hits_signal[il], voltrck.weight);
//       else lhist->FillHistW("ecal_n_hits_background", il, necal_hits_backgr[il], voltrck.weight);
//     }
//     
//     if (showeren>0.0) {
//       lhist->FillHistW("ecal_shower_e_signal", 0, showeren, voltrck.weight);
//       lhist->FillHistW("ecal_shower_x_signal", 0, showerx/showeren, showeren*voltrck.weight);
//       lhist->FillHistW("ecal_shower_xe_signal", 0, showerx/showeren, showeren, voltrck.weight);
//     }
//     
//     for (auto itr = track_hit_map.cbegin(); itr != track_hit_map.cend(); ++itr) {
//       if (primpdg==-11) {  
//         if (itr->first == 1) {
//           lhist->FillHistW("tracking_tracks_nhits_signal", 0, itr->second.size(), voltrck.weight);
//         } else {
//           lhist->FillHistW("tracking_tracks_nhits_other", 0, itr->second.size(), voltrck.weight);
//         }
//       } else {
//         lhist->FillHistW("tracking_tracks_nhits_background", 0, itr->second.size(), voltrck.weight);
//       }
//     }
    
    track_hit_map.clear();
    
    ++nev;
  }

  std::string suffix("_hits");
  std::string foutname = fnamelist.substr(fnamelist.find_last_of("/")+1);
  foutname = foutname.substr(0, foutname.find_last_of("."));
  foutname += suffix + std::string(".root");
  lhist->SaveHists(foutname);

//  delete hitstree;
  return 0;  
}



int ProcessList(const std::string &fnamelist, std::vector<std::string> &flist)
{
  std::fstream  fdata;
  fdata.open(fnamelist, std::ios::in);
  if (!fdata.is_open()) {
    throw std::runtime_error(std::string("Error reding data from the file ") + fnamelist);
  }
  
  unsigned long lid = 0;
  while (!fdata.eof()) {
    std::string  ffname;
    double fweight;
    fdata >> ffname;
    if (!fdata.fail()) { 
//       std::cout << "File name " << ffname << " is read from the list file" << std::endl;
      flist.push_back(ffname);
    }
    else if (fdata.eof()) { break; }
    else {
      std::cout << "ProcessList(..)  :  Error reading data from the file " << fnamelist 
                << ",  line: " << lid << ". Exit." << std::endl;
      fdata.close();          
      return -2;
    }
    ++lid;
  }
  
  fdata.close();

  return 0;
}  
  


void CreateHistograms(MHists *mh) 
{
  std::cout << "Creating histograms\n";

  mh->AddHists("primary_e", 3, 20000, 0.0, 20.0);
  
  //// variable binned in X axis histograms
  Double_t logmin           = -6;
  Double_t logmax           = 1;
  const Int_t nbins         = 450;
  Double_t logbinwidth      = (logmax-logmin)/nbins;
  Double_t xpoints[nbins+2] = {0.0};
  for(Int_t i=0 ; i<=nbins ; i++) xpoints[i] = TMath::Power( 10,(logmin + i*logbinwidth) );
  xpoints[451] = 2*TMath::Power( 10,1);

  int nlayers = 16;
  int nsensor = 9;
  int ndet = nsensor * nlayers;
  int npixx = 1024;
  int npixy = 512;
  mh->AddHists("tracking_planes_hits_x", ndet, npixx, 0.0, npixx);
  mh->AddHists("tracking_planes_hits_y", ndet, npixy, 0.0, npixy);
  mh->AddHists("tracking_planes_hits_xy", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
  mh->AddHists("tracking_planes_hits_xy_edep", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
  mh->AddHists("tracking_planes_hits_edep", ndet, 20000, 0.0, 20.0);
//   mh->AddHists("tracking_planes_hits_x_aftervtxcut", ndet, npixx, 0.0, npixx);
//   mh->AddHists("tracking_planes_hits_y_aftervtxcut", ndet, npixy, 0.0, npixy);
//   mh->AddHists("tracking_planes_hits_xy_aftervtxcut", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
//   mh->AddHists("tracking_planes_hits_xy_edep_aftervtxcut", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
//   
//   
//   mh->AddHists("tracking_planes_n_tracks_per_hit_aftervtxcut", ndet, 100, 0.0, 100.0);
  
  mh->AddHists("tracking_planes_n_tracks_per_hit", nlayers, 100, 0.0, 100.0);
  mh->AddHists("tracking_planes_n_tracks_per_hit_signal", nlayers, 100, 0.0, 100.0);
  mh->AddHists("tracking_planes_n_tracks_per_hit_background", nlayers, 100, 0.0, 100.0);
        
  mh->AddHists("tracking_planes_hit_track_e", nlayers, 20000, 0.0, 20.0);
  mh->AddHists("tracking_planes_hit_track_proc", nlayers, 2100, 0.0, 2100.0);
  mh->AddHists("tracking_planes_hit_track_pdg", nlayers, 100, -50.0, 50.0);
  
  
  /// new hostogram for QGSP physics list
  mh->AddHistsLogBin("tracking_planes_hits_edep_log", ndet, nbins+1, xpoints);
  mh->AddHistsLogBin("tracking_planes_hit_track_edep_neutron_log", nlayers, nbins+1, xpoints);/// new histogram
//   mh->AddHists("tracking_planes_hit_track_time_neutron", nlayers, 1000, -50, 50);/// new histogram
//   mh->AddHists("tracking_planes_hit_track_proc_neutron", nlayers, 2100, 0.0, 2100.0);
  mh->AddHists("tracking_planes_hits_xy_edep_neutron", ndet, npixx, 0.0, npixx, npixy, 0.0, npixx);
  
  mh->AddHistsLogBin("tracking_planes_hit_track_edep_proton_log", nlayers, nbins+1, xpoints);/// new histogram
//   mh->AddHists("tracking_planes_hit_track_time_proton", nlayers, 1000, -50, 50);/// new histogram
//   mh->AddHists("tracking_planes_hit_track_proc_proton", nlayers, 2100, 0.0, 2100.0);
  mh->AddHists("tracking_planes_hits_xy_edep_proton", ndet, npixx, 0.0, npixx, npixy, 0.0, npixx);
  
  mh->AddHistsLogBin("tracking_planes_hit_track_edep_pi0_log", nlayers, nbins+1, xpoints);/// new histogram
//   mh->AddHists("tracking_planes_hit_track_time_pi0", nlayers, 1000, -50, 50);/// new histogram
//   mh->AddHists("tracking_planes_hit_track_proc_pi0", nlayers, 2100, 0.0, 2100.0);
  mh->AddHists("tracking_planes_hits_xy_edep_pi0", ndet, npixx, 0.0, npixx, npixy, 0.0, npixx);
  
  
  mh->AddHistsLogBin("tracking_planes_hit_track_edep_pi_log", nlayers, nbins+1, xpoints);/// new histogram
//   mh->AddHists("tracking_planes_hit_track_time_pi", nlayers, 1000, -50, 50);/// new histogram
//   mh->AddHists("tracking_planes_hit_track_proc_pi", nlayers, 2100, 0.0, 2100.0);
  mh->AddHists("tracking_planes_hits_xy_edep_pi", ndet, npixx, 0.0, npixx, npixy, 0.0, npixx);
  
  mh->AddHistsLogBin("tracking_planes_hit_track_edep_muon_log", nlayers, nbins+1, xpoints);/// new histogram
//   mh->AddHists("tracking_planes_hit_track_time_muon", nlayers, 1000, -50, 50);/// new histogram
//   mh->AddHists("tracking_planes_hit_track_proc_muon", nlayers, 2100, 0.0, 2100.0);
  mh->AddHists("tracking_planes_hits_xy_edep_muon", ndet, npixx, 0.0, npixx, npixy, 0.0, npixx);


          
//   mh->AddHists("tracking_planes_hit_track_e_zcut", nlayers, 20000, 0.0, 20.0);
//   mh->AddHists("tracking_planes_hit_track_pdg_zcut", nlayers, 100, -50.0, 50.0);
//   mh->AddHists("tracking_planes_hit_track_e_signal_zcut", nlayers, 20000, 0.0, 20.0);
//   mh->AddHists("tracking_planes_hit_track_pdg_signal_zcut", nlayers, 100, -50.0, 50.0);
//   mh->AddHists("tracking_planes_hit_track_e_background_zcut", nlayers, 20000, 0.0, 20.0);
//   mh->AddHists("tracking_planes_hit_track_pdg_background_zcut", nlayers, 100, -50.0, 50.0);
//   
//   mh->AddHists("tracking_planes_hit_track_id_signal_zcut", nlayers, 100, 0.0, 100.0);
//   mh->AddHists("tracking_planes_hit_track_id_background_zcut", nlayers, 100, 0.0, 100.0);
  
  mh->AddHists("tracker_n_hits_signal", nlayers, 1000, 0.0, 1000.0);
  mh->AddHists("tracker_n_hits_background", nlayers, 1000, 0.0, 1000.0);

  mh->AddHists("tracking_tracks_nhits_signal", 1, 10, 0.0, 10.0);
  mh->AddHists("tracking_tracks_nhits_other", 1, 10, 0.0, 10.0);
  mh->AddHists("tracking_tracks_nhits_background", 1, 10, 0.0, 10.0);
  
  
  
  
  
  mh->AddHistsLogBin("tracking_planes_hits_edep_log", ndet, nbins+1, xpoints);
  mh->AddHistsLogBin("tracking_planes_hit_track_edep_all_log", ndet, nbins+1, xpoints);
  
//   nlayers = 21;
//   nsensor = 2;
//   ndet = nsensor * nlayers;
//   npixx = 110;
//   npixy = 11;
//   mh->AddHists("ecal_hits_x", ndet, npixx, 0.0, npixx);
//   mh->AddHists("ecal_hits_edep", ndet, 10000, 0.0, 10.0);
//   mh->AddHists("ecal_hits_xy", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
//   mh->AddHists("ecal_hits_ntrck", ndet, 1000, 0.0, 1000);
//   mh->AddHists("ecal_hits_xy_edep", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
//   mh->AddHists("ecal_hits_xy_ntrck", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
// 
//   mh->AddHists("ecal_hit_edep_signal", nlayers, 10000, 0.0, 1.0e-2);
//   mh->AddHists("ecal_hit_z_signal", nsensor, 30, 0.0, 30.0);
//   mh->AddHists("ecal_hit_z_signal_edep", nsensor, 30, 0.0, 30.0);
//   mh->AddHists("ecal_hit_zx_signal_edep", nsensor, 30, 0.0, 30.0, npixx, 0.0, npixx);
//   mh->AddHists("ecal_hit_edep_background", nlayers, 10000, 0.0, 1.0e-2);
//   mh->AddHists("ecal_hit_z_background", nsensor, 30, 0.0, 30.0);
//   mh->AddHists("ecal_hit_z_background_edep", nsensor, 30, 0.0, 30.0);
//   mh->AddHists("ecal_hit_zx_background_edep", nsensor, 30, 0.0, 30.0, npixx, 0.0, npixx);
//   mh->AddHists("ecal_n_hits_signal", nlayers, 1000, 0.0, 1000.0);
//   mh->AddHists("ecal_n_hits_background", nlayers, 1000, 0.0, 1000.0);
//   
//   mh->AddHists("ecal_hit_xy_signal_edep", nsensor, npixx, 0.0, npixx, npixy, 0.0, npixy);
//   mh->AddHists("ecal_hit_xy_layer_signal_edep", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
//   mh->AddHists("ecal_hit_signal_weight", 1, 300000, 0.0, 300000);
//   mh->AddHists("ecal_hit_xy_background_edep", nsensor, npixx, 0.0, npixx, npixy, 0.0, npixy);
//   mh->AddHists("ecal_hit_xy_layer_background_edep", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
// 
//   mh->AddHists("ecal_hit_edep_background_dump", nlayers, 10000, 0.0, 1.0e-2);
//   mh->AddHists("ecal_hit_z_background_dump", nsensor, 30, 0.0, 30.0);
//   mh->AddHists("ecal_hit_z_background_dump_edep", nsensor, 30, 0.0, 30.0);
//   mh->AddHists("ecal_hit_zx_background_dump_edep", nsensor, 30, 0.0, 30.0, npixx, 0.0, npixx);
//   mh->AddHists("ecal_hit_xy_background_dump_edep", nsensor, npixx, 0.0, npixx, npixy, 0.0, npixy);
//   mh->AddHists("ecal_hit_xy_layer_background_dump_edep", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
// 
//   mh->AddHists("ecal_hit_edep_background_el_low", nlayers, 10000, 0.0, 1.0e-2);
//   mh->AddHists("ecal_hit_z_background_el_low", nsensor, 30, 0.0, 30.0);
//   mh->AddHists("ecal_hit_z_background_el_low_edep", nsensor, 30, 0.0, 30.0);
//   mh->AddHists("ecal_hit_zx_background_el_low_edep", nsensor, 30, 0.0, 30.0, npixx, 0.0, npixx);
//   mh->AddHists("ecal_hit_xy_background_el_low_edep", nsensor, npixx, 0.0, npixx, npixy, 0.0, npixy);
//   mh->AddHists("ecal_hit_xy_layer_background_el_low_edep", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
//           
//   mh->AddHists("ecal_hit_edep_background_gamma", nlayers, 10000, 0.0, 1.0e-2);
//   mh->AddHists("ecal_hit_z_background_gamma", nsensor, 30, 0.0, 30.0);
//   mh->AddHists("ecal_hit_z_background_gamma_edep", nsensor, 30, 0.0, 30.0);
//   mh->AddHists("ecal_hit_zx_background_gamma_edep", nsensor, 30, 0.0, 30.0, npixx, 0.0, npixx);
//   mh->AddHists("ecal_hit_xy_background_gamma_edep", nsensor, npixx, 0.0, npixx, npixy, 0.0, npixy);
//   mh->AddHists("ecal_hit_xy_layer_background_gamma_edep", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
//   
//   mh->AddHists("ecal_shower_e_signal", 1, 1000, 0.0, 1.0);
//   mh->AddHists("ecal_shower_x_signal", 1, 110, 0.0, 110.0);
//   mh->AddHists("ecal_shower_xe_signal", 1, 110, 0.0, 110.0, 1000, 0.0, 1.0);
//   
//   ndet = 2;
//   npixx = 150;
//   npixy = 25;
//   mh->AddHists("lyso_hits_x", ndet, npixx, 0.0, npixx);
//   mh->AddHists("lyso_hits_xy", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
//   mh->AddHists("lyso_hits_xy_edep", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
//   mh->AddHists("lyso_hits_edep", ndet, 10000, 0.0, 10.0);
// 
//   ndet = 8;
//   mh->AddHists("gammamon_hit_edep", ndet, 10000, 0.0, 10000.0);
//   
//   ndet = 1;
//   npixx = 1000;
//   npixy = 200;
//   mh->AddHists("ce_screen_hit_xy_signal_prime_edep", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
//   mh->AddHists("ce_screen_hit_edep_signal_prime", ndet, 1000.0, 0.0, 1.0e-2);
//   mh->AddHists("ce_screen_hit_n_tracks_signal_prime", ndet, 50.0, 0.0, 50.0);
//   mh->AddHists("ce_screen_hit_xy_signal_other_edep", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
//   mh->AddHists("ce_screen_hit_edep_signal_other", ndet, 1000.0, 0.0, 1.0e-2);
//   mh->AddHists("ce_screen_hit_n_tracks_signal_other", ndet, 50.0, 0.0, 50.0);
//   mh->AddHists("ce_screen_hit_xy_background_edep", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
//   mh->AddHists("ce_screen_hit_edep_background", ndet, 1000.0, 0.0, 1.0e-2);
//   mh->AddHists("ce_screen_hit_n_tracks_background", ndet, 50.0, 0.0, 50.0);

}


#endif

