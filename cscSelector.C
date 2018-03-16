#define cscSelector_cxx
// The class definition in cscSelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("cscSelector.C")
// root> T->Process("cscSelector.C","some options")
// root> T->Process("cscSelector.C+")
//


#include "cscSelector.h"
#include <TH2.h>
#include <TStyle.h>
#include <algorithm>  

#include "TMatrixDSparse.h"
#include "clib/fillMatrix.h"
#include "clib/CSCWithMuon.h"
#include "clib/scanPattern.h"

#include "clib/CSCWireSegmentPattern.h"
#include "clib/CSCStripSegmentPattern.h"

/*
void cscSelector::Initialize()
{

   endcapL.clear(); stationL.clear(); ringL.clear(); chamberL.clear();

}
*/

void cscSelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

//   nSegPerChamber = new TH1F("nSegPerChamebr","",10,0,10);
//   nRHPerSeg = new TH1F("nRHPerSeg","",7,0,7);

   TString option = GetOption();
}

void cscSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t cscSelector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetEntry(entry);
  
//   if (entry > 20000) return kTRUE; 
   //cout << "Event: " << *Event << endl;
   //if (*Event == 59284848/*479434223*/) cout << *cscSegments_nSegments << endl;

   SaveCSCWithMuon();

   // if no CSC contains a muon segment, skip this event
   if (int(endcapL.size()) == 0) return kTRUE;

   CountObjectsInCSCs(true,true,true,true,true,false,false,false); //RH,Seg,Wire,Strip,Comparator,ALCT,CLCT,CLT

// SOMEHOW CANNOT PROCESS ALL DIGIS TOGHETHER !!!

   for (int i = 0; i < int(endcapL.size()); i++) {

       int endcap_s = endcapL[i];
       int station_s = stationL[i];
       int ring_s = ringL[i];
       int chamber_s = chamberL[i];
//if (!(endcap_s ==2 && station_s ==4 && ring_s==1 && chamber_s==15)) continue;
       // use chamber Only with n seg
//       if (nSegmentL[i].size() != 1) continue;
       if (nSegmentL[i].size() != 4) continue;
       int segIndex = nSegmentL[i][0];
//       if (cscSegments_nRecHits[segIndex] != nRH[i].size()) continue;
//       if (!(cscSegments_nRecHits[segIndex] == 6 && nRH[i].size() == 6 ) ) continue;

/*
       if (!(cscSegments_nRecHits[nSegmentL[i][0]] == 6 && 
             cscSegments_nRecHits[nSegmentL[i][1]] == 6 && 
             cscSegments_nRecHits[nSegmentL[i][2]] == 6 && 
             cscSegments_nRecHits[nSegmentL[i][3]] == 6 && nRH[i].size() < 30) ) continue;
*/

// fill matrix

       if (nWire[i].size() <= 0 || nComparator[i].size() <=0) continue;    

       int nWireGroups = -1;
       if (station_s == 1) nWireGroups = 48;
       if (station_s >= 1 && ring_s == 2) nWireGroups = 64;
       if (station_s == 2 && ring_s == 1) nWireGroups = 112;
       if (station_s > 2 && ring_s == 1) nWireGroups = 96;
//cout << "nWireGroup " <<  nWireGroups << endl;
       TMatrixDSparse wireMatrix(6,nWireGroups);
       FillWireMatrix(nWire[i], wireMatrix);

       int nStrips = 80;
       if (station_s == 1 && ring_s == 4) nStrips = 112;//48;
       if (station_s == 1 && ring_s == 1) nStrips = 112;//64;
       if (station_s == 1 && ring_s == 3) nStrips = 64;

       bool doStagger = true;
       if (station_s == 1 && (ring_s == 1 || ring_s == 4) ) doStagger = false;

       TMatrixDSparse comparatorMatrix(6,2*nStrips+2);
       FillComparatorMatrix(nComparator[i], comparatorMatrix, doStagger);


       bool reverseRowIndex = false;
       if (station_s == 3 || station_s == 4) reverseRowIndex = true;

       vector<CSC1DSeg> allWireSegs = MakeScans(wireMatrix, reverseRowIndex, w_rows, w_cols, nWGsInPatterns, patternRanks_w, 4);
       vector<CSC1DSeg> allComparatorSegs = MakeScans(comparatorMatrix, reverseRowIndex, s_rows, s_cols, nHalfStrips, patternRanks_s, 9);

       vector<CSC1DSeg> allWireSegs_old = MakeScans(wireMatrix, reverseRowIndex, w_rows_old, w_cols_old, nWGsInPatterns_old, patternRanks_w_old, 1);
       vector<CSC1DSeg> allComparatorSegs_old = MakeScans(comparatorMatrix, reverseRowIndex, s_rows_old, s_cols_old, nHalfStrips_old, patternRanks_s_old, 9);

if (allComparatorSegs.size()==2 && allComparatorSegs_old.size()==1) 
   {
// cout << "nWireSeg: " << allWireSegs.size() << ", nStripSeg: " << allComparatorSegs.size() << ", nWireSeg_old: " << allWireSegs_old.size() << ", nStripSeg_old: " << allComparatorSegs_old.size() << endl;

   cout << "wide pattern ranks: " << allComparatorSegs[0].patternRank << " and " << allComparatorSegs[1].patternRank << endl;
   cout << "narrow pattern ranks: " << allComparatorSegs_old[0].patternRank << endl;

   PrintSparseMatrix(comparatorMatrix);

   }

       // one wire and com seg per chamber

/*
       if (int(allWireSegs.size()) != 1 || int(allComparatorSegs.size()) != 1) continue;
       CSC1DSeg wireSeg = allWireSegs[0];
       CSC1DSeg comSeg = allComparatorSegs[0];
       if (wireSeg.nHits != 6 || comSeg.nHits != 6) continue;

       wireTime_comTime->Fill(wireSeg.MeanTime(), comSeg.MeanTime());
 
       wireTime_comTime_l1->Fill((wireSeg.hitsTime)[0], (comSeg.hitsTime)[0]);
       wireTime_comTime_l2->Fill((wireSeg.hitsTime)[1], (comSeg.hitsTime)[1]);
       wireTime_comTime_l3->Fill((wireSeg.hitsTime)[2], (comSeg.hitsTime)[2]);
       wireTime_comTime_l4->Fill((wireSeg.hitsTime)[3], (comSeg.hitsTime)[3]);
       wireTime_comTime_l5->Fill((wireSeg.hitsTime)[4], (comSeg.hitsTime)[4]);
       wireTime_comTime_l6->Fill((wireSeg.hitsTime)[5], (comSeg.hitsTime)[5]);

       vector<int> sortedWireTime = wireSeg.SortedTime();
       vector<int> sortedComTime = comSeg.SortedTime();

       wireTime_comTime_r1->Fill(sortedWireTime[0],sortedComTime[0]); 
       wireTime_comTime_r2->Fill(sortedWireTime[1],sortedComTime[1]);
       wireTime_comTime_r3->Fill(sortedWireTime[2],sortedComTime[2]);
       wireTime_comTime_r4->Fill(sortedWireTime[3],sortedComTime[3]);
       wireTime_comTime_r5->Fill(sortedWireTime[4],sortedComTime[4]);
       wireTime_comTime_r6->Fill(sortedWireTime[5],sortedComTime[5]);


       // if having 2 wire and strip segs
       if (int(allWireSegs.size()) != 2 || int(allComparatorSegs.size()) != 2) continue; 
       CSC1DSeg wireSeg1 = allWireSegs[0];
       CSC1DSeg wireSeg2 = allWireSegs[1];
       CSC1DSeg comSeg1 = allComparatorSegs[0];
       CSC1DSeg comSeg2 = allComparatorSegs[1];

       if (wireSeg1.nHits != 6 || comSeg1.nHits != 6 || wireSeg2.nHits != 6 || comSeg2.nHits != 6) continue;

       vector<int> sortedWireTime1 = wireSeg1.SortedTime();
       vector<int> sortedComTime1 = comSeg1.SortedTime();
       vector<int> sortedWireTime2 = wireSeg2.SortedTime();
       vector<int> sortedComTime2 = comSeg2.SortedTime();

cout << sortedWireTime1[2] << ", " << sortedComTime1[2] << ", " << sortedWireTime2[2] << ", " << sortedComTime2[2] << endl;

     if ( ! (  sortedWireTime1[2] >=7 && sortedWireTime1[2] <= 11 && sortedComTime1[2] >= 6 && sortedComTime1[2] <=10 &&
               sortedWireTime2[2] >=7 && sortedWireTime2[2] <= 11 && sortedComTime2[2] >= 6 && sortedComTime2[2] <=10) ) continue; 
cout << "here" << endl;            
       deltaTime->Fill(abs(wireSeg1.MeanTime()-wireSeg2.MeanTime()), abs(comSeg1.MeanTime()-comSeg2.MeanTime()));

       deltaTime_l1->Fill(abs((wireSeg1.hitsTime)[0]-(wireSeg2.hitsTime)[0]), abs((comSeg1.hitsTime)[0]-(comSeg2.hitsTime)[0]));
       deltaTime_l2->Fill(abs((wireSeg1.hitsTime)[1]-(wireSeg2.hitsTime)[1]), abs((comSeg1.hitsTime)[1]-(comSeg2.hitsTime)[1]));
       deltaTime_l3->Fill(abs((wireSeg1.hitsTime)[2]-(wireSeg2.hitsTime)[2]), abs((comSeg1.hitsTime)[2]-(comSeg2.hitsTime)[2]));
       deltaTime_l4->Fill(abs((wireSeg1.hitsTime)[3]-(wireSeg2.hitsTime)[3]), abs((comSeg1.hitsTime)[3]-(comSeg2.hitsTime)[3]));
       deltaTime_l5->Fill(abs((wireSeg1.hitsTime)[4]-(wireSeg2.hitsTime)[4]), abs((comSeg1.hitsTime)[4]-(comSeg2.hitsTime)[4]));
       deltaTime_l6->Fill(abs((wireSeg1.hitsTime)[5]-(wireSeg2.hitsTime)[5]), abs((comSeg1.hitsTime)[5]-(comSeg2.hitsTime)[5]));

       deltaTime_r1->Fill(abs(sortedWireTime1[0]-sortedWireTime2[0]), abs(sortedComTime1[0]-sortedComTime2[0]));       
       deltaTime_r2->Fill(abs(sortedWireTime1[1]-sortedWireTime2[1]), abs(sortedComTime1[1]-sortedComTime2[1]));
       deltaTime_r3->Fill(abs(sortedWireTime1[2]-sortedWireTime2[2]), abs(sortedComTime1[2]-sortedComTime2[2]));
       deltaTime_r4->Fill(abs(sortedWireTime1[3]-sortedWireTime2[3]), abs(sortedComTime1[3]-sortedComTime2[3]));
       deltaTime_r5->Fill(abs(sortedWireTime1[4]-sortedWireTime2[4]), abs(sortedComTime1[4]-sortedComTime2[4]));
       deltaTime_r6->Fill(abs(sortedWireTime1[5]-sortedWireTime2[5]), abs(sortedComTime1[5]-sortedComTime2[5]));
*/
//cout << "endcap: " << endcap_s << ", station: " << station_s << ", ring: " << ring_s << ", chamber: " << chamber_s << endl;
//cout << "nWireSeg: " << allWireSegs.size() << endl;
//cout << "nStripSeg: " << allComparatorSegs.size() << endl;

    }//loop over chambers with muon

    return kTRUE;


}

void cscSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void cscSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   TCanvas *c1 = new TCanvas("c1", "", 800,800);//200,10,400,400);
   outputRootFile = new TFile("tmpRootPlots/CSCresults_" + tag + ".root","RECREATE");
   outputRootFile->cd();
//c1->SetLogy();
//chi2PerDOF->SetMinimum(0.1);
//nSegPerChamber->SetMaximum(9000);
//nRHPerSeg->SetMaximum(8000);
//   nSegPerChamber->Write();
   nRHPerSeg->Write();
//   chi2PerDOF->Write();
   nWireDigi_Layer->Write();
   nStripDigi_Layer->Write();
   nComparatorDigi_Layer->Write();

   wireTime_layer->Write();
   comparatorTime_layer->Write();
   wireTime_comparatorTime->Write();
   wireTime_comparatorTime_out->Write();

   wireTime_comTime->Write();

   wireTime_comTime_l1->Write();
   wireTime_comTime_l2->Write();
   wireTime_comTime_l3->Write();
   wireTime_comTime_l4->Write();
   wireTime_comTime_l5->Write();
   wireTime_comTime_l6->Write();

   wireTime_comTime_r1->Write();
   wireTime_comTime_r2->Write();
   wireTime_comTime_r3->Write();
   wireTime_comTime_r4->Write();
   wireTime_comTime_r5->Write();
   wireTime_comTime_r6->Write();

   deltaTime->Write();

   deltaTime_l1->Write();
   deltaTime_l2->Write();
   deltaTime_l3->Write();
   deltaTime_l4->Write();
   deltaTime_l5->Write();
   deltaTime_l6->Write();

   deltaTime_r1->Write();
   deltaTime_r2->Write();
   deltaTime_r3->Write();
   deltaTime_r4->Write();
   deltaTime_r5->Write();
   deltaTime_r6->Write();

   nRH_lumi_ME11a->Write();
   nRH_lumi_ME11b->Write();
   nRH_lumi_MEx1->Write();
   nRH_lumi_MEx2->Write();
   nRH_lumi_ME13->Write();

   nSeg_lumi_ME11a->Write();
   nSeg_lumi_ME11b->Write();
   nSeg_lumi_MEx1->Write();
   nSeg_lumi_MEx2->Write();
   nSeg_lumi_ME13->Write();

   nALCT_lumi_ME11a->Write();
   nALCT_lumi_ME11b->Write();
   nALCT_lumi_MEx1->Write();
   nALCT_lumi_MEx2->Write();
   nALCT_lumi_ME13->Write();

   nCLCT_lumi_ME11a->Write();
   nCLCT_lumi_ME11b->Write();
   nCLCT_lumi_MEx1->Write();
   nCLCT_lumi_MEx2->Write();
   nCLCT_lumi_ME13->Write();

   nLCT_lumi_ME11a->Write();
   nLCT_lumi_ME11b->Write();
   nLCT_lumi_MEx1->Write();
   nLCT_lumi_MEx2->Write();
   nLCT_lumi_ME13->Write();

   outputRootFile->Close();

   cout << "root file " << outputRootFile->GetName() << " made " << endl;
//   nSegPerChamber->Draw("hist");
//   nRHPerSeg->Draw("hist");
//   chi2PerDOF->Draw("hist");
//   c1->SaveAs(savedir + "chi2PerDOF_" + tag + ".png");
//   c1->SaveAs(savedir + "chi2PerDOF_" + tag + ".pdf");



//   c1->SaveAs(savedir + "nSegsPerChamber_" + tag + ".png");
//   c1->SaveAs(savedir + "nSegsPerChamber_" + tag + ".pdf");
//   c1->SaveAs(savedir + "nRHPerSeg_" + tag + ".png");
//   c1->SaveAs(savedir + "nRHPerSeg_" + tag + ".pdf");

}



// user defined function

void cscSelector::SetInputs(TString tag_, TString savedir_, bool doME11_)
{

     tag = tag_;
     savedir = savedir_;
     doME11 = doME11_;

}



void cscSelector::FillSegs(int segIndex, int chamberIndex, vector<SegsInChamber>& segs)
{

       if (int(segs.size()) == 0) {

          vector<int> tmpIndex; tmpIndex.push_back(segIndex);
          segs.push_back(make_pair(chamberIndex,tmpIndex));

          } else {

                 int savedIndex = -1;
                 for (int j = 0; j < int(segs.size()); j++) {
                     if (segs[j].first == chamberIndex) savedIndex = j;
                     }

                 if (savedIndex >= 0) {

                    segs[savedIndex].second.push_back(segIndex);

                    } else {

                           vector<int> tmpIndex2; tmpIndex2.push_back(segIndex);
                           segs.push_back(make_pair(chamberIndex,tmpIndex2));

                           }

                 }


}




/*
void cscSelector::MakeWireMatrix()
{


}
*/
void cscSelector::CountObjectsInCSCs(bool doRH, bool doSeg, 
                                     bool doWire, bool doStrip, bool doComparator, 
                                     bool doALCT, bool doCLCT, bool doLCT)
{

// count number of seg/alct/clct/lct in CSCs with muon segment

   for (int i = 0; i < int(endcapL.size()); i++) {

       int endcap = endcapL[i];
       int station = stationL[i];
       int ring = ringL[i];
       int chamber = chamberL[i];

       if (doRH) {
          nRH[i].clear();
          for (int j = 0; j < *recHits2D_nRecHits2D; j++) {
              int endcap_rh = recHits2D_ID_endcap[j]; 
              int station_rh = recHits2D_ID_station[j];
              int ring_rh = recHits2D_ID_ring[j];
              int chamber_rh = recHits2D_ID_chamber[j];

              if (endcap==endcap_rh && station == station_rh && ring == ring_rh && chamber == chamber_rh) nRH[i].push_back(j);// += 1;
          }
       }

       if (doSeg) {
          nSegmentL[i].clear();
          for (int j = 0; j < *cscSegments_nSegments; j++) {
              int endcap_s = cscSegments_ID_endcap[j]; 
              int station_s = cscSegments_ID_station[j];
              int ring_s = cscSegments_ID_ring[j];
              int chamber_s = cscSegments_ID_chamber[j];

              if (endcap==endcap_s && station == station_s && ring == ring_s && chamber == chamber_s) nSegmentL[i].push_back(j); // += 1;
          }
       }

       if (doWire) {
          nWire[i].clear();
          for (int j = 0; j < *firedWireDigis_nWireDigis; j++) {

              int endcap_w = firedWireDigis_ID_endcap[j];
              int station_w = firedWireDigis_ID_station[j];
              int ring_w = firedWireDigis_ID_ring[j];
              int chamber_w = firedWireDigis_ID_chamber[j];

              if (endcap==endcap_w && station == station_w && ring == ring_w && chamber == chamber_w) nWire[i].push_back(j); // += 1;
          }
       }

       if (doStrip) {
          nStrip[i].clear();
          for (int j = 0; j < *firedStripDigis_nStripDigis; j++) {
              int endcap_sp = firedStripDigis_ID_endcap[j];
              int station_sp = firedStripDigis_ID_station[j];
              int ring_sp = firedStripDigis_ID_ring[j];
              int chamber_sp = firedStripDigis_ID_chamber[j];

              if (endcap==endcap_sp && station == station_sp && ring == ring_sp && chamber == chamber_sp) nStrip[i].push_back(j); // += 1;
          }
       }

       if (doComparator) {
          nComparator[i].clear();
          for (int j = 0; j < *comparatorDigis_nDigis; j++) {
              int endcap_c = comparatorDigis_ID_endcap[j];
              int station_c = comparatorDigis_ID_station[j];
              int ring_c = comparatorDigis_ID_ring[j];
              int chamber_c = comparatorDigis_ID_chamber[j];

              if (endcap==endcap_c && station == station_c && ring == ring_c && chamber == chamber_c) nComparator[i].push_back(j); // += 1;
          }
       }
/*
       if (doALCT) {
          nALCT[i].clear();
          for (int j = 0; j < *alct_nAlcts; j++) {
              int endcap_alct = alct_ID_endcap[j];
              int station_alct = alct_ID_station[j];
              int ring_alct = alct_ID_ring[j];
              int chamber_alct = alct_ID_chamber[j];

              if (endcap==endcap_alct && station == station_alct && ring == ring_alct && chamber == chamber_alct) nALCT[i].push_back(j);// += 1;
          }
       }

       if (doCLCT) {
          nCLCT[i].clear();
          for (int j = 0; j < *clct_nClcts; j++) {
              int endcap_clct = clct_ID_endcap[j];
              int station_clct = clct_ID_station[j];
              int ring_clct = clct_ID_ring[j];
              int chamber_clct = clct_ID_chamber[j];

              if (endcap==endcap_clct && station == station_clct && ring == ring_clct && chamber == chamber_clct) nCLCT[i].push_back(j);// += 1;
          }
       }

       if (doLCT) {
          nLCT[i].clear();
          for (int j = 0; j < *correlatedLct_nLcts; j++) {
              int endcap_lct = correlatedLct_ID_endcap[j];
              int station_lct = correlatedLct_ID_station[j];
              int ring_lct = correlatedLct_ID_ring[j];
              int chamber_lct = correlatedLct_ID_chamber[j];

              if (endcap==endcap_lct && station == station_lct && ring == ring_lct && chamber == chamber_lct) nLCT[i].push_back(j);// += 1;
          }
       }
*/
   }


/*   
   for (int i = 0; i < int(endcapL.size()); i++) {

       // select chamber

       int endcap_s = endcapL[i];
       int station_s = stationL[i];
       int ring_s = ringL[i];
       int chamber_s = chamberL[i];

       if (station_s == 1 && ring_s == 4) {
          nRH_lumi_ME11a->Fill(*instLumi, nRH[i].size() );
          nSeg_lumi_ME11a->Fill(*instLumi, nSegmentL[i].size() );
          nALCT_lumi_ME11a->Fill(*instLumi, nALCT[i].size() );
          nCLCT_lumi_ME11a->Fill(*instLumi, nCLCT[i].size() );
          nLCT_lumi_ME11a->Fill(*instLumi, nLCT[i].size() );
          }

       if (station_s == 1 && ring_s == 1) {
          nRH_lumi_ME11b->Fill(*instLumi, nRH[i].size() );
          nSeg_lumi_ME11b->Fill(*instLumi, nSegmentL[i].size() );
          nALCT_lumi_ME11b->Fill(*instLumi, nALCT[i].size() );
          nCLCT_lumi_ME11b->Fill(*instLumi, nCLCT[i].size() );
          nLCT_lumi_ME11b->Fill(*instLumi, nLCT[i].size() );
          }

       if (station_s != 1 && ring_s == 1) {
          nRH_lumi_MEx1->Fill(*instLumi, nRH[i].size() );
          nSeg_lumi_MEx1->Fill(*instLumi, nSegmentL[i].size() );
          nALCT_lumi_MEx1->Fill(*instLumi, nALCT[i].size() );
          nCLCT_lumi_MEx1->Fill(*instLumi, nCLCT[i].size() );
          nLCT_lumi_MEx1->Fill(*instLumi, nLCT[i].size() );
          }

       if (ring_s == 2) {
          nRH_lumi_MEx2->Fill(*instLumi, nRH[i].size() );
          nSeg_lumi_MEx2->Fill(*instLumi, nSegmentL[i].size() );
          nALCT_lumi_MEx2->Fill(*instLumi, nALCT[i].size() );
          nCLCT_lumi_MEx2->Fill(*instLumi, nCLCT[i].size() );
          nLCT_lumi_MEx2->Fill(*instLumi, nLCT[i].size() );
          }

       if (ring_s == 3) {
          nRH_lumi_ME13->Fill(*instLumi, nRH[i].size() );
          nSeg_lumi_ME13->Fill(*instLumi, nSegmentL[i].size() );
          nALCT_lumi_ME13->Fill(*instLumi, nALCT[i].size() );
          nCLCT_lumi_ME13->Fill(*instLumi, nCLCT[i].size() );
          nLCT_lumi_ME13->Fill(*instLumi, nLCT[i].size() );
          }

       }
*/
}


void cscSelector::PrintSparseMatrix(TMatrixDSparse inputMatrix) 
{

     int * rowIndex = inputMatrix.GetRowIndexArray();
     int * colIndex = inputMatrix.GetColIndexArray();
     double * pData = inputMatrix.GetMatrixArray();

     int nRow = inputMatrix.GetNrows();
     int nCol = inputMatrix.GetNcols();

     TH2F* h_matrix = new TH2F("h_matrix","",nCol,0,nCol,nRow,0,nRow);

     for (int i = 0; i < nRow; i++) {
         int sIndex = rowIndex[i];
         int eIndex = rowIndex[i+1];
         for (int j = sIndex; j < eIndex; j++) {
             int icol = colIndex[j];
             double data = pData[j];
//             cout << i+inputMatrix.GetRowLwb()+1 << ", " << icol+inputMatrix.GetColLwb() << ", " << data << endl;
             h_matrix->SetBinContent(icol+inputMatrix.GetColLwb(),i+inputMatrix.GetRowLwb()+1,data);
             }
         }
//cout << h_matrix->Integral() << endl;
     WriteTH2F(h_matrix);
     delete h_matrix;      
}

void cscSelector::WriteTH2F(TH2F* hist) {
//cout << hist->Integral() << endl;

     for (int i = 1; i < hist->GetNbinsY()+1; i++) {
         for (int j = 1; j < hist->GetNbinsX()+1; j++) {
             if (hist->GetBinContent(j,i)==0)
                {std::cout << "-";} else {std::cout << hist->GetBinContent(j,i)-1;}

             }
         std::cout << std::endl;
         }

}
