void cscSelector::SaveCSCWithMuon()
{


   endcapL.clear(); stationL.clear(); ringL.clear(); chamberL.clear();
   vector<int> ZmuIndex; ZmuIndex.clear();

   for (int i = 0; i < *muons_nMuons; i++) {

       if (!(muons_isPFMuon[i] || muons_isGlobalMuon[i] || muons_isTrackerMuon[i]) ) continue;
       if (!(abs(muons_eta[i]) < 2.4 && abs(muons_eta[i]) < 2.4)) continue;
       if (!(abs(muons_dz[i]) < 1 && abs(muons_dxy[i]) < 0.05)) continue;
       if (muons_pt[i] < 5) continue;

       ZmuIndex.push_back(i);

   }

   for (int i = 0; i < int(ZmuIndex.size()); i++) {
       int tmpIndex = ZmuIndex[i];
       for (int j = 0; j < int(muons_cscSegmentRecord_endcap[tmpIndex].size()); j++) {
           int endcap = muons_cscSegmentRecord_endcap[tmpIndex][j];
           int station = muons_cscSegmentRecord_station[tmpIndex][j];
           int ring = muons_cscSegmentRecord_ring[tmpIndex][j];
           int chamber = muons_cscSegmentRecord_chamber[tmpIndex][j];

           if ( (find(endcapL.begin(), endcapL.end(), endcap) != endcapL.end() ) &&
                (find(stationL.begin(), stationL.end(), station) != stationL.end() ) &&
                (find(ringL.begin(), ringL.end(), ring) != ringL.end() ) &&
                (find(chamberL.begin(), chamberL.end(), chamber) != chamberL.end() )
              ) continue;

           endcapL.push_back(endcap);
           stationL.push_back(station);
           ringL.push_back(ring);
           chamberL.push_back(chamber);

           }
       }


}

