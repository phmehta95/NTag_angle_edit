#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <limits>

#include <TRandom3.h>
#include <TCollection.h>
#include <TList.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include "TMath.h"
#include <geotnkC.h>

#include "Calculator.hh"
<<<<<<< HEAD
#include "skparmC.h"
#include "sktqC.h"
#include "skheadC.h"
#include "apmringC.h"
#include "apmueC.h"
#include "apmsfitC.h"
#include "appatspC.h"
#include "apringspC.h"
#include "skroot_loweC.h"
#include "spliTChanOutC.h"
#include "fitqunoutC.h"
#include "geotnkC.h"
#include "neworkC.h"
#include "nbnkC.h"
#include "skonl/softtrg_cond.h"

#include "SKLibs.hh"
#include "SKIO.hh"
#include "GetStopMuVertex.hh"
#include "NTagBankIO.hh"
#include "Calculator.hh"
#include "NoiseManager.hh"
#include "EventNTagManager.hh"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TreeManager.h"
=======
>>>>>>> 030e87beaa8c96090057e292df276de5d7db0b56

std::default_random_engine c_ranGen;
TRandom3 ranGen;


float Dot(const float a[3], const float b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

float Norm(const float vec[3])
{
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

float Norm(float x, float y, float z)
{
    return sqrt(x*x + y*y + z*z);
}

float GetDistance(const float vec1[3], const float vec2[3])
{
    float tmpVec[3];

    for (int i = 0; i < 3; i++)
        tmpVec[i] = vec1[i] - vec2[i];

    return Norm(tmpVec);
}

float GetRMS(const std::vector<float>& vec)
{
    float N  = static_cast<float>(vec.size());
    float mean = 0.;
    float var  = 0.;

    for (auto const& value: vec)
        mean += value / N;
    for (auto const& value: vec)
        var += (value-mean)*(value-mean) / (N-1);

    return sqrt(var);
}

float GetLegendreP(int i, float& x)
{
    float result = 0.;

    switch (i) {
        case 1:
            result = x; break;
        case 2:
            result = (3*x*x-1)/2.; break;
        case 3:
            result = (5*x*x*x-3*x)/2; break;
        case 4:
            result = (35*x*x*x*x-30*x*x+3)/8.; break;
        case 5:
            result = (63*x*x*x*x*x-70*x*x*x+15*x)/8.; break;
    }

    return result;
}

float GetOpeningAngle(TVector3 uA, TVector3 uB, TVector3 uC)
{
    // make sure the inputs are unit vectors
    // uA = uA.Unit(); uB = uB.Unit(); uC = uC.Unit();

    // sides of the triangle formed by the three unit vectors
    double a = (uA-uB).Mag();
    double b = (uC-uA).Mag();
    double c = (uB-uC).Mag();

    if (a*b*c == 0) {
        //double angleAB = (180./M_PI) * uA.Angle(uB)/2.;
        //double angleAC = (180./M_PI) * uA.Angle(uC)/2.;
        return uA.Angle(uB) == 0 ? (uA.Angle(uC) == 0 ? 0 : uA.Angle(uC)) : uA.Angle(uB);
    }

    else {
        // circumradius of the triangle
        double r = a*b*c / sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c));

        if (r>=1)
            return 90.; // prevents NaN
        else
            return (180./M_PI) * asin(r);
    }
}

float GetDWall(TVector3 vtx)
{
    float vertex[3] = {(float)vtx.x(), (float)vtx.y(), (float)vtx.z()};
    return wallsk_(vertex);
}

float GetDWallInDirection(TVector3 vtx, TVector3 dir)
{
    dir = dir.Unit();

    float dot = vtx.Dot(dir) - vtx.z()*dir.z();
    float dirSq = dir.Perp2(); float vtxSq = vtx.Perp2();

    // Calculate distance to barrel and distance to top/bottom
    float distR = (-dot + sqrt(dot*dot + dirSq*(RINTK*RINTK - vtxSq))) / dirSq;
    float distZ = dir.z() > 0 ? (ZPINTK-vtx.z())/dir.z() : (ZMINTK-vtx.z())/dir.z();

    // Return the smaller
    return distR < distZ ? distR : distZ;
}

unsigned int GetMinIndex(std::vector<float>& vec)
{
    std::vector<float>::iterator iter = std::min_element(vec.begin(), vec.end());
    size_t index = std::distance(vec.begin(), iter);
    return index;
}

unsigned int GetMaxIndex(std::vector<float>& vec)
{
    std::vector<float>::iterator iter = std::max_element(vec.begin(), vec.end());
    size_t index = std::distance(vec.begin(), iter);
    return index;
}

void SetSeed(int seed)
{
<<<<<<< HEAD
  c_ranGen.seed(seed);
  ranGen.SetSeed(seed);
=======
    c_ranGen.seed(seed);
    ranGen.SetSeed(seed);
>>>>>>> 030e87beaa8c96090057e292df276de5d7db0b56
}

TString PickFile(TString dirPath, const char* extension)
{
    std::vector<TString> list = GetListOfFiles(dirPath, extension);
    return PickRandom(list);
}

TString PickSubdirectory(TString dirPath)
{
    std::vector<TString> list = GetListOfSubdirectories(dirPath);
    return PickRandom(list);
}

std::vector<TString> GetListOfFiles(TString dirPath, const char* extension, bool recursive)
{
    std::vector<TString> list;

    TSystemDirectory dir(dirPath, dirPath);
    TList *files = dir.GetListOfFiles();

    if (files) {
        TSystemFile *file;
        TString fileName, filePath;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fileName = file->GetName();
            filePath = file->GetTitle();
            if (recursive && file->IsDirectory() && fileName != "." && fileName != "..") {
                std::cout << "Checking files in directory: " << filePath << std::endl;
                std::vector<TString> subdirList = GetListOfFiles(filePath, extension);
                list.insert(list.end(), subdirList.begin(), subdirList.end());
            }
            if (!file->IsDirectory() && fileName.EndsWith(extension)) {
                list.push_back(filePath + "/" + fileName);
            }
        }
    }

    return list;
}

std::vector<TString> GetListOfSubdirectories(TString dirPath)
{
    std::vector<TString> list;

    TSystemDirectory dir(dirPath, dirPath);
    TList *files = dir.GetListOfFiles();

    if (files) {
        TSystemFile *file;
        TString fileName, filePath;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fileName = file->GetName();
            filePath = file->GetTitle();
            if (file->IsDirectory() && fileName != "." && fileName != "..") {
                list.push_back(filePath);
            }
        }
    }

    return list;
}

std::vector<unsigned int> GetRangeIndex(const std::vector<double>& sortedVec, double low, double high)
{
    std::vector<unsigned int> vIndex;
    //std::cout << "low: " << low << " high: " << high << std::endl;

    auto start = std::lower_bound(sortedVec.begin(), sortedVec.end(), low);
    auto end   = std::upper_bound(sortedVec.begin(), sortedVec.end(), high);

    //std::cout << "start_i: " << start-sortedVec.begin() << " end_i: " << end-sortedVec.begin() << std::endl;

    for (auto it=start; it!=end; ++it) {
        vIndex.push_back(it-sortedVec.begin());
    }

    return vIndex;
}

std::vector<std::pair<float, int>> Histogram(std::vector<float> vec, int nBins, float min, float max)
{
    std::vector<std::pair<float, int>> hist(nBins, std::pair<float, int>(0, 0));
    const float bWidth = (max - min) / float(nBins);

    std::sort(vec.begin(), vec.end());
    
    int index = 0;
    auto binStartItr = vec.begin();
    for (auto& bin: hist) {
        const float bMin = min + index*bWidth;
        const float bMax = min + (index+1)*bWidth;

        auto start = std::lower_bound(binStartItr, vec.end(), bMin);
        auto end = std::lower_bound(start, vec.end(), bMax*0.99999);
        binStartItr = end;
        
        bin.first = (bMin+bMax)/2.; // bin center
        bin.second = end - start;   // bin value
        index++;
    }
    
    return hist;
}

std::vector<std::string> Split(std::string target, std::string delim)
{
    std::vector<std::string> v;

    if (!target.empty()) {
        while (target.size() ){
            int index = target.find(delim);
            if (index != std::string::npos) {
                v.push_back(target.substr(0,index));
                target = target.substr(index+delim.size());
                if (target.size()==0) v.push_back(target);
            }
            else {
                v.push_back(target);
                target = "";
            }
        }         
    }
    return v;
}



  Bool_t checkMiss (const Int_t cab)
{
  for (Int_t i=0; i<NMIS; i++) {
    if ( MISCH[i] == cab ) return kTRUE;
  }
  return kFALSE;
}

int getNhits(Float_t *v, int start_index, Float_t width, int nhits)
// Get nhits from v[0] in time window = width
// v[] - hit timing, sorted in time
// width - time window width
// start_index - index of start pmt hit
// nhits - total number  of v[] hits
{
    int i = start_index;
    while (1) {
        i++;
        if((i > nhits-1 ) || (TMath::Abs((v[i]-v[start_index])) > width)) break;
    }
    return TMath::Abs(i - start_index);
}





void MaxHitsInWindow(Float_t *t, int ALLHITS, float width, int &maxnhits, int &maxindex)
{
    maxnhits = 0;
    maxindex = 0;

    //Finding max N15 and its index, limiting potential number of hits
    for (int in = 0; in < ALLHITS; in++) {

        // Limit to -500 to 1500 ns time window
        //if ( t[in] < tchelow)  continue;
        //if ( t[in] > tchehigh) continue;

        // Calculate hits in 15 ns window
        int N = getNhits(t, in, width, ALLHITS);
        if ( N > maxnhits) {
            maxnhits = N;
            maxindex = in;
        }
    }//ALLHITS
}


Double_t lfrsqrt(Double_t a2, Double_t b2, Double_t c2)
{
  Double_t rsquare=0.0;
  Double_t tmpr2=0.0;

  //++allCount;

  tmpr2= a2*(2.*b2-a2) + b2*(2.*c2-b2) + c2*(2.*a2-c2);
  if (tmpr2>0){
    rsquare= a2*b2*c2/tmpr2;
  }else {
      //++negativeCount;
      //printf("lfrsqrd: Problem: %f, %f, %f, %f.\n", tmpr2, a2, b2, c2);
    return 2;
  }

  return rsquare;
}




void MinimumWidthWindow(Float_t *t, int ALLHITS, float width, int &maxnhits, int &maxindex)
{
    // ----------------------------------------------------------------------
    // ---- get time ordered tof subtracted (vertex)
    // ----------------------------------------------------------------------

    // v is already sorted


    // ----------------------------------------------------------------------
    // ---- scan for the requested number of hits
    // ----------------------------------------------------------------------

    double dt   = 100000.;
    for (int i = 0; i < ALLHITS-maxnhits; i++) {
      //if ( t[i] < tchelow)  continue;
      //if ( t[i] > tchehigh) continue;

        if (t[i+maxnhits] - t[i] < dt) {
            dt = t[i+maxnhits] - t[i];
        }
    }


    // ----------------------------------------------------------------------
    // ---- if dt too big, find maximal number of hits fitting into twindow
    // ----------------------------------------------------------------------

    if (dt > width) {
    //    if (verbose)
    //  std::cout << "dt too wide: "std::endl;
        MaxHitsInWindow(t, ALLHITS, width, maxnhits, maxindex);
    }
    //else {
    //   if (verbose)
    //std::cout << "MinimumWidthWindow, maxnhits=" << maxnhits << ", maxindex=" << maxindex << std::endl;
    //}

}

float GetCherenkovAngle(float *bonsaivertex, bool limithits)
{

    int nstep = 100;
    angledown = 0;
    angleup   = 90;
    tchelow   = -500;
    tchehigh  = 1500;
    hcheren = new TH1D("hcheren", "Opening Cherenkov Angle from 3-PMT combinations", nstep, angledown, angleup);

  
  Int_t goodcount=0, ipmt;
  for (int i=0;i<sktqz_.nqiskz;i++){
        if (sktqz_.ihtiflz[i]&0x02 && sktqz_.icabiz[i] <= MAXPM && !checkMiss(sktqz_.icabiz[i]) && sktqz_.ihtiflz[i]&0x01){ //check if hits are in gate and not on bad/missing channels and hits that are main event
            goodcount+=1;
  
}
  }
   const Int_t MAXN10 = 200;

    const int ALLHITS = goodcount;
    const Float_t C_WATER = 21.5833;

    Int_t   cabiz[ALLHITS], cabiz2[ALLHITS];
    Float_t tiskz[ALLHITS], tiskz2[ALLHITS];
    Float_t qiskz[ALLHITS], qiskz2[ALLHITS];
    Int_t   index[ALLHITS], nindex[MAXN10];
    Double_t dhit[3], dtmp;
    Float_t hitv_x[ALLHITS];
    Float_t hitv_y[ALLHITS];
    Float_t hitv_z[ALLHITS];
    Double_t rsqrd, opang, rnosq;
    Int_t icount=0;

    // Copy the TQ arrays
    
      for (int i=0; i<sktqz_.nqiskz; i++)
        {
          if (sktqz_.ihtiflz[i]&0x02 && sktqz_.icabiz[i] <= MAXPM  && !checkMiss(sktqz_.icabiz[i]\
) && sktqz_.ihtiflz[i]&0x01) //check if hits are in gate and not on bad/missing channels and hits that are main event
            {
              cabiz2[icount] = sktqz_.icabiz[i];
              tiskz2[icount] = sktqz_.tiskz[i];
              qiskz2[icount] = sktqz_.qiskz[i];
              icount+=1;
            }


	}
	
// TOF subtraction for all hits in copied arrays
    for (int i=0; i<ALLHITS; i++) {
        //if (i % 100 == 0)     cout << "position : " << cabiz2[i]-1<<  geopmt_.xyzpm[cabiz2[i]-1][0]<<endl;
        Float_t tof;
        tof = TMath::Sqrt((bonsaivertex[0] - geopmt_.xyzpm[cabiz2[i]-1][0]) * (bonsaivertex[0] - geopmt_.xyzpm[cabiz2[i]-1][0])
                          +(bonsaivertex[1] - geopmt_.xyzpm[cabiz2[i]-1][1]) * (bonsaivertex[1] - geopmt_.xyzpm[cabiz2[i]-1][1])
                          +(bonsaivertex[2] - geopmt_.xyzpm[cabiz2[i]-1][2]) * (bonsaivertex[2] - geopmt_.xyzpm[cabiz2[i]-1][2])) / C_WATER;
        tiskz2[i] -= tof;

	//std::cout << "TOF CORRECTED HIT TIMES: " << tiskz2[i] << std::endl;
    }


    // Sort hits by TOF-corrected time
    TMath::Sort(ALLHITS, tiskz2, index, kFALSE); // In increasing order
    for (int i=0; i<ALLHITS; i++){
        cabiz[i] = cabiz2[ index[i] ];
        tiskz[i] = tiskz2[ index[i] ];
        qiskz[i] = qiskz2[ index[i] ];
    }

    int N15      = 0;
    int N15index = 0;
    //std::cout << "DBG1" << std::endl;
    if (limithits)   {
      //std::cout << "DBG2" << std::endl;

       N15 = 183;
       //std::cout << "DBG3" << std::endl;

       MinimumWidthWindow(tiskz, ALLHITS, 15., N15, N15index);
       //std::cout << "DBG30" <<std::endl;
    }
    //std::cout << "DBG4" << std::endl;
    else {
    //  std::cout << "DBG5" << std::endl;
    MaxHitsInWindow(tiskz, ALLHITS, 15., N15, N15index);
    }



    
    int maxncmb = N15*N15*N15 / 6;
    std::vector<double> abc2(N15*N15, 0.0);
    //std::vector<double> abc2(ALLHITS*ALLHITS, 0.0);
    std::vector<double> maxcmb(maxncmb, 0.0);

    hcheren->Reset();

    for (int ii=N15index; ii < N15index+N15; ii++){
    //for (int ii=N15index; ii < ; ii++){
      //std::cout << "ii" << ii << std::endl;
      //std::cout << "N15index" << N15index << std::endl;
      //std::cout << "N15" << N15 << std::endl;
        for (int jj=0; jj<3; jj++){
	  dhit[jj] = geopmt_.xyzpm[cabiz[ii]-1][jj]- bonsaivertex[jj];
	  //std::cout << "dhit" << dhit[jj] << std::endl;
	  //  std::cout << "dhit0" << dhit[0] << std::endl;
	  //  std::cout << "dhit1" << dhit[1] << std::endl;
	  //  std::cout << "dhit2" << dhit[2] << std::endl;
        }
        //dtmp= 1.0/TMath::Sqrt(dhit[0]*dhit[0]+dhit[1]*dhit[1]+dhit[2]*dhit[2]);
	dtmp=1.0/(TMath::Sqrt(dhit[0]*dhit[0]+dhit[1]*dhit[1]+dhit[2]*dhit[2]));
	//std::cout << "DTMP " << dtmp << std::endl;
	int intmp= ii-N15index;
	//int intmp=ii;
	hitv_x[intmp]= dhit[0]*dtmp;
        hitv_y[intmp]= dhit[1]*dtmp;
        hitv_z[intmp]= dhit[2]*dtmp;
    }

    int tmpindex_a=0, tmpindex_b=0, tmpindex_c=0;
    //Calculate lengths of difference vectors between directions
    for (int aa = 0; aa < N15-1; aa++) {
      //for (int aa = 0; aa < ALLHITS-1; aa++) {
      for (int bb=aa+1; bb<N15; bb++) {
	//for (int bb=aa+1; bb<ALLHITS; bb++) {
	tmpindex_a=aa*N15+bb;
	// tmpindex_a=aa*ALLHITS+bb;
            abc2[tmpindex_a] = (hitv_x[aa]-hitv_x[bb])*(hitv_x[aa]-hitv_x[bb])+(hitv_y[aa]-hitv_y[bb])*(hitv_y[aa]-hitv_y[bb])+(hitv_\
z[aa]-hitv_z[bb])*(hitv_z[aa]-hitv_z[bb]);
	    //std::cout << "abc2: " << abc2[tmpindex_a] << std::endl;
	}
    }
       //Do all direction triangles and fill histogram
          for (int aa=0;aa<N15-2;aa++){
      //for (int aa=0;aa<ALLHITS-2;aa++){
	  for (int bb=aa+1;bb<N15-1;bb++){
	//for (int bb=aa+1;bb<ALLHITS-1;bb++){
            tmpindex_a=aa*N15+bb;
            //dtmp= abc2[aa][bb];
            dtmp= abc2[tmpindex_a];
	    //std::cout << "dtmp " << dtmp << std::endl;
            for (int cc=bb+1;cc<N15;cc++){
                tmpindex_b=aa*N15+cc;
                tmpindex_c=bb*N15+cc;
		//std::cout << "tmpindex_b" << tmpindex_b << std::endl;
		//std::cout << "tmpindex_c" << tmpindex_c << std::endl;
		rsqrd = lfrsqrt(dtmp, abc2[tmpindex_b], abc2[tmpindex_c]);
		//std::cout << "rsqrd " << rsqrd << std::endl;
                rnosq = TMath::Sqrt(rsqrd);
		//std::cout << "RNOSQ " << rnosq << std::endl;
                opang = TMath::ASin(rnosq)*180.0/TMath::Pi();
		//return opang;
                // if rnosq>1-->return pi/2
		//std::cout << "OPANG " << opang << std::endl;
                hcheren->Fill(opang);
		//hcheren->Draw();
            }
        }
    }


      
    Double_t unpack[nstep], b2deg=0.0, lsum=0.0;
    int nwindow=7, midpoint=3, look=100;
    int ilowbnd=41, ihibnd=67;
    int nrange, llook, ij=0, locale=0;
    Double_t lheight=0, peakangle=0.0;

    //std::cout << "NSTEP" << nstep << std::endl;
    //std::cout << "UNPACK" << unpack[nstep] << std::endl;

    nrange = nstep-nwindow+1;
    //std::cout << "NRANGE" << nrange << std::endl;
    b2deg  = 90.0/100.0;

    anglebinnum= nstep;
    for(int hh=0;hh<nstep;hh++){
        unpack[hh]=hcheren->GetBinContent(hh+1);
        plotcontent[hh] = hcheren->GetBinContent(hh);
	//std::cout << "UNPACKnew" << unpack[nstep] << std::endl;
    }
    for(int ii=0;ii<nrange;ii++){
        lsum=0.0;
        for (int jj=0;jj<nwindow;jj++){
            ij=ii+jj;
            lsum+=unpack[ij];
	    //std::cout << "lsum" << lsum << std::endl;
        }
        if (lheight<lsum){
            lheight=lsum;
            locale=ii+1+midpoint;}
	//std::cout << "locale" << locale << std::endl;
    }
    peakangle=locale*b2deg;
    //hcheren->Draw();
    //std::cout << "peakangle" << peakangle << std::endl;
    //fEventVariables.Set("CherenkovAngle", peakangle);
    delete hcheren;
    return peakangle;
    //    printf("Cherenkov angle is %f\n",  peakangle);
    
    

}

