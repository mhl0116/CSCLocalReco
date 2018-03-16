#include <vector>
#include <TROOT.h>
#include "TMatrixDSparse.h"
#include "TArrayI.h"
#include "TArrayD.h"
#include <vector>
#include <utility>      // std::pair, std::make_pair
#include <algorithm>    // std::find

vector<CSC1DSeg> cscSelector::MakeScans(TMatrixDSparse inputMatrix, bool reverseRowIndex, int** rows, int** cols, int* nHits, int* ranks, int nPatterns) {
  
    vector<CSC1DSeg> allSegs1D; 

    for (int i = 6; i >= 3; i--) { // nHits/Seg loop
     for (int j = 0; j < nPatterns; j++) { // patterns loop

         int* rowIndex = rows[j];
         int* colIndex = cols[j];
         int nhits = nHits[j];
         int rank = ranks[j];

         vector<CSC1DSeg> new1DSegs = ScanPattern(inputMatrix, i, reverseRowIndex, rowIndex, colIndex, nhits, rank);
         // add to allSegs1D
         for (int k = 0; k < int(new1DSegs.size()); k++) allSegs1D.push_back(new1DSegs[k]);
//cout << "allSegs1D size: " << allSegs1D.size() << endl;
         }
        } 

    return allSegs1D;
}


vector<CSC1DSeg> cscSelector::ScanPattern(TMatrixDSparse& inputMatrix, int nHitsPerSeg, bool reverseRowIndex, //for ME 3 and ME 4
                                                         int* r_index, int* c_index, int nHitsInPattern, int patternRank) {

     vector<CSC1DSeg> cscSegs1D; // container saving 1D segment
     vector<int> usedCols; // save used cols (wire group or half strip number), avoid repeat scan
     int* rowArray = inputMatrix.GetRowIndexArray(); // https://root.cern.ch/root/html534/guides/users-guide/LinearAlgebra.html
     int* colArray = inputMatrix.GetColIndexArray();
     int nElement = rowArray[inputMatrix.GetNrows()]; // only scan cols with non-zero element
     int nCols = inputMatrix.GetNcols();

     for (int i = 0; i < nElement; i++) {
         int keyPos = colArray[i];
         if (find(usedCols.begin(), usedCols.end(), keyPos) != usedCols.end() ) continue;  

//         int c_index_tmp[nHitsInPattern] = {};  

         int timeHits[6] = {0,0,0,0,0,0};
         int nNonZeroHits[6] = {0,0,0,0,0,0};
         vector<pair<int, int> > nonZeroIndexs;

         for (int j = 0; j < nHitsInPattern; j++) {
             int c_index_tmp = c_index[j] + keyPos;
//             c_index_tmp[j] += keyPos; // shift pattern to corresponding key wire group / half strip
           
             if (c_index_tmp < 0 || c_index_tmp >= nCols) continue;
             
             int rowIndex = r_index[j];
             int colIndex = c_index_tmp;
             if (reverseRowIndex) rowIndex = 5-rowIndex;
             int thisElement = inputMatrix[rowIndex][colIndex]; // if wire group or comparator/strip fired, this value should be positive (by construction), other wise 0

             timeHits[rowIndex] += thisElement;
             if (thisElement > 0) {
                nNonZeroHits[rowIndex]+=1;
                nonZeroIndexs.push_back(make_pair(rowIndex,colIndex));
                }
             }

         int nHits = CountNonZero(nNonZeroHits,6);         
         if (nHits < nHitsPerSeg) continue; // a segment is found only if number of hits is bigger than nHitsPerSeg (6,5,4,3)

         CSC1DSeg tmp1DSeg;
         tmp1DSeg.keyPos = keyPos; usedCols.push_back(keyPos);
         tmp1DSeg.nHits = nHits;
         tmp1DSeg.patternRank = patternRank; 

         for (int k = 0; k < 6; k++) {
             if (nNonZeroHits[k]>0) (tmp1DSeg.hitsTime)[k] = (timeHits[k]/nNonZeroHits[k]);
             else (tmp1DSeg.hitsTime)[k] = 0;
             } // note, this is int, may need to think to change it to double      
         // clean up used wire/strip hits
         for (int k = 0; k < int(nonZeroIndexs.size()); k++) inputMatrix[nonZeroIndexs[k].first][nonZeroIndexs[k].second] = 0;
//         Print1DSeg(tmp1DSeg);
         cscSegs1D.push_back(tmp1DSeg);

         }
 
    return cscSegs1D;  
}


int cscSelector::CountNonZero(int * array, int size) {

    int count = 0;

    for (int i = 0; i < size; i++) {

     if (array[i] > 0) count+=1;

    }

    return count;
}


void cscSelector::Print1DSeg(CSC1DSeg seg1D) {

     cout << "keyPos: " << seg1D.keyPos << endl;
     cout << "nHits: "  << seg1D.nHits << endl;
     cout << "time: ";
     for (int i = 0; i < 6; i++) cout << (seg1D.hitsTime)[i] << " ";

     cout << endl;
     cout << seg1D.MeanTime() << endl;
}

