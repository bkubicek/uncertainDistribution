#ifndef _BKRANDOMCPP
#define _BKRANDOMCPP


#include <vector>
#include <random>
#include <tuple>
//#include <math.h>
#include <cmath>
#include <utility>

const double pi = 3.14159265358979323846;

typedef std::vector<double> VecD;
typedef std::tuple< VecD, VecD> TupleVecD;
typedef std::pair< double, double > PairDouble;

inline void normalizeVec(VecD &v)
{
    double wsum=0;
    for(unsigned i=0;i<v.size();i++)
        if (v[i]<0)
            v[i]=0;
    for(unsigned i=0;i<v.size();i++)
        wsum+=v[i];
    for(unsigned i=0;i<v.size();i++)
        v[i]*=1/wsum;
}



inline TupleVecD getNormal(double mean, double stddev, double sigmaMax, unsigned N)
{    
    double leftVal= mean-sigmaMax;
    double rightVal= mean+sigmaMax;
    double delta = (rightVal-leftVal)/N;
    VecD centers;
    centers.resize(N);
    VecD weights(N,0);
    
    
    centers[0]=leftVal+delta/2;
    for(unsigned i=1;i<N;i++)
        centers[i]=centers[i-1]+delta;
    double preFac=1/std::sqrt(2*pi)/stddev;
    double expFac=0.5/stddev/stddev;
    for(unsigned i=0;i<N;i++)
        weights[i]=preFac*std::exp(-std::pow(centers[i]-mean,2) *expFac);
    
    normalizeVec(weights);
    
    return {centers, weights};
}

inline TupleVecD getTri(double leftVal, double centerVal, double rightVal, unsigned N)
{
    double delta = (rightVal-leftVal)/N;
    VecD centers;
    centers.resize(N);
    VecD weights(N,0);
    
    double centerPos=(centerVal-leftVal)/delta;
    
    double stepL = 1/centerPos;
    double stepR = 1/(N-centerPos);
    
    
    centers[0]=leftVal+delta/2.;
    for(unsigned i=1;i<N;i++)
        centers[i]=centers[i-1]+delta;
    
    weights[0] = stepL*delta/2.;
    for(unsigned i=1;i<N;i++)
        if (i<centerPos)
            weights[i]=weights[i-1]+stepL;
        else
            weights[i]=weights[i-1]-stepR;
    
    normalizeVec(weights);
    
    return {centers, weights};
    
}


inline TupleVecD getRect(double leftVal, double rightVal, unsigned N)
{    
    double delta = (rightVal-leftVal)/N;
    VecD centers;
    centers.resize(N+2);
    VecD weights(N+2,1/double(N));
    
    
    centers[0]=leftVal-delta/2;
    weights[0]=0;
    for(unsigned i=1;i<N+1;i++)
        centers[i]=centers[i-1]+delta;
    centers[N+1]=centers[N]+delta;
    weights[N+1]=0;
    
    return {centers, weights};
}



inline bool comparatorPair ( const PairDouble &l, const PairDouble &r)
   { return l.first < r.first; }
  
 
inline TupleVecD resampleDistribution(const VecD &centers, const VecD &weights, unsigned newN)
{
    VecD newCenters;
    VecD newWeights;
    newCenters.resize(newN);
    newWeights.resize(newN,0);
    
    
    double newLeft= *std::min_element(centers.begin(), centers.end());
    double newRight= *std::max_element(centers.begin(), centers.end());
    double delta = (newRight-newLeft)/newN;
    
    newCenters[0]=newLeft+delta/2;
    for(unsigned i=1;i<newN;i++)
        newCenters[i]=newCenters[i-1]+delta;
    
    for(unsigned i=0;i<centers.size();i++)
    {
        double indexPosF=(centers[i]-newLeft)/delta;
        int indexPos= std::round(indexPosF);
        if(indexPos<0)
            indexPos=0;
        if(indexPos>=newN)
            indexPos=newN-1;
        newWeights[indexPos]+=weights[i];
    }
    
    return {newCenters, newWeights};
}

inline TupleVecD addDistributions(
    const VecD &centers1, const VecD &weights1, 
    const VecD &centers2, const VecD &weights2, unsigned newN)
{
    VecD resCenters;
    VecD resWeights;
    resCenters.reserve( centers1.size()*centers2.size() );
    resWeights.reserve( centers1.size()*centers2.size() );
    for (unsigned i =0;i<centers1.size();i++)
    {
        double c1=centers1[i];
        double w1=weights1[i];
        for (unsigned j =0;j<centers2.size();j++)
        {
                double c2=centers2[j];
                double w2=weights2[j];
                resCenters.push_back(c1+c2);
                resWeights.push_back(w1*w2);
                
        }
    }
    normalizeVec(resWeights);
    
    return  resampleDistribution(resCenters,resWeights,newN);

}


inline TupleVecD multiplyDistributions(
     const VecD &centers1, const VecD &weights1, 
     const VecD &centers2, const VecD &weights2, unsigned newN)
{
    VecD resCenters;
    VecD resWeights;
    resCenters.reserve( centers1.size()*centers2.size() );
    resWeights.reserve( centers1.size()*centers2.size() );
    for (unsigned i =0;i<centers1.size();i++)
    {
        double c1=centers1[i];
        double w1=weights1[i];
        for (unsigned j =0;j<centers2.size();j++)
        {
                double c2=centers2[j];
                double w2=weights2[j];
                resCenters.push_back(c1*c2);
                resWeights.push_back(w1*w2);
                
        }
    }
    normalizeVec(resWeights);
    
    return  resampleDistribution(resCenters,resWeights,newN);

}

inline TupleVecD divideDistributions(
   const VecD &centers1, const VecD &weights1, 
   const VecD &centers2, const VecD &weights2, unsigned newN)
{
    VecD resCenters;
    VecD resWeights;
    resCenters.reserve( centers1.size()*centers2.size() );
    resWeights.reserve( centers1.size()*centers2.size() );
    for (unsigned i =0;i<centers1.size();i++)
    {
        double c1=centers1[i];
        double w1=weights1[i];
        for (unsigned j =0;j<centers2.size();j++)
        {
                double c2=centers2[j];
                double w2=weights2[j];
                resCenters.push_back(c1/c2);
                resWeights.push_back(w1*w2);
                
        }
    }
    normalizeVec(resWeights);

    return  resampleDistribution(resCenters,resWeights,newN);
}




#endif