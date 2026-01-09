/**
 * Some C++ routines for handling covulation of distributions
 *
 * The distributions are given by two arrays, each beeing a std::vector.
 * The "centers" array is the center points of the intervals in the x-axis
 * The "weights" array the likelyhood of the distribution at the corresponding center point
 * The weights sum up 1. This is different from "they integrate to 1". 
 * By having this nomenclature, the weights actually descripe the likelyhood of beeing in that interval
 *
 * The purpose is to have a toolkit for error propagation calculations
 * Might be usefull for people who perform accredited testing, and need to proove accuracy.
 * 
 * all definitions are inline, so that they can "legally" recide in a header file only
 * this simplifies compilation, as no corresponding c++ file needs to be defined
 *
 * initial written by Bernhard Kubicek, AIT Austrian institute of technology
 *  originally for the SolidPV EU Project, EURAMET
 * 
 * License: LGPL-3.0-or-later
 */

#ifndef _BKRANDOMCPP
#define _BKRANDOMCPP


#include <vector>
#include <random>
#include <tuple>
#include <cmath>
#include <utility>

const double pi = 3.14159265358979323846;

//! an Enum to allow to define if one wants to add, multiply or divide two distributions
enum DistOps {OpAdd, OpMul, OpDiv};

//! simpler type names
typedef std::vector<double> VecD;
typedef std::tuple< VecD, VecD> TupleVecD;
typedef std::pair< double, double > PairDouble;


//! take a weights vector, and scale it so that the sum is 1
inline void normalizeVec(VecD &v)
{
    // override negative weights, that are not possible in distributions
    for(unsigned i=0;i<v.size();i++)
        if (v[i]<0)
            v[i]=0;
    // find sum
    double wsum=0;
    for(unsigned i=0;i<v.size();i++)
        wsum+=v[i];
    // scale to one
    for(unsigned i=0;i<v.size();i++)
        v[i]*=1/wsum;
}


//! calculate a normal distribution represented in the system of interval centers and weights
//! @param mean self explenatory
//! @param stddev  self explenatory
//! @param sigmaMax is the maximum sigma, until the distribution should be evaluated. "3-4" sigmas is reasonable
//! @param N is the number of samples used.
inline TupleVecD getNormal(double mean, double stddev, double sigmaMax, unsigned N)
{    
    // find the interval
    double leftVal= mean-sigmaMax*stddev;
    double rightVal= mean+sigmaMax*stddev;
    double delta = (rightVal-leftVal)/N;
    
    //allocate the result vectors
    VecD centers(N,0);
    VecD weights(N,0);
    
    // write the centers
    centers[0]=leftVal+delta/2;
    for(unsigned i=1;i<N;i++)
        centers[i]=centers[i-1]+delta;
        
    // optimized factors for the noraml distribution function
    double preFac=1/std::sqrt(2*pi)/stddev;
    double expFac=0.5/stddev/stddev;
    
    for(unsigned i=0;i<N;i++)
        weights[i]=preFac*std::exp(-std::pow(centers[i]-mean,2) *expFac);
    
    //normalize so that not the integral but the sum is one.
    normalizeVec(weights);
    
    return {centers, weights};
}

//! calculate a triangular distribution represented in the system of interval centers and weights
//! @param leftVal the smallest pos; at which the raising flank start
//! @param centerVal the pos; at which the peak is reached
//! @param rightVal the largest pos; at which the falling flank stops
//! @param N is the number of samples used. 
inline TupleVecD getTri(double leftVal, double centerVal, double rightVal, unsigned N)
{
    double delta = (rightVal-leftVal)/N;
    VecD centers(N,0);
    VecD weights(N,0);
    
    // index of center Position
    double centerPos=(centerVal-leftVal)/delta; 
    
    // steepness of rise and fall
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

//! calculate a triangular distribution represented in the system of interval centers and weights
//! @param leftVal the smallest pos; at which the raising flank start
//! @param centerVal the pos; at which the peak is reached
//! @param rightVal the largest pos; at which the falling flank stops
//! @param N is the number of samples used. 
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


//! helper function while sorting a vector of pairs of doubles
//! so that the first element is used for sorting, and the second is moved correspondingly
inline bool comparatorPair ( const PairDouble &l, const PairDouble &r)
{ 
    return l.first < r.first; 
}

//! find the min/max of a vector
//! @param v the vector
//! @param dmin the double that the minimum is written into
//! @param dmax the double that the maximum is written into
inline void getLimits(const VecD &v, double &dmin, double &dmax)
{
    dmin= *std::min_element(v.begin(), v.end());
    dmax= *std::max_element(v.begin(), v.end());
}
  
//! resample a distributino in the system of interval centers and weights to a smaller size
//! @param centers the original interval centers
//! @param weights the likelyhood of the corresponding interval centers
//! @param newN the new size should be smaller!
inline TupleVecD resampleDistribution(const VecD &centers, const VecD &weights, unsigned newN)
{
    // allocate results vectors
    VecD newCenters;
    VecD newWeights;
    newCenters.resize(newN);
    newWeights.resize(newN,0);
    
    //find range
    double newLeft, newRight; //init. not necessary as always overwritten in getLimits
    getLimits(centers, newLeft, newRight);
    double delta = (newRight-newLeft)/newN;
    
    //write resampled center positions
    newCenters[0]=newLeft+delta/2;
    for(unsigned i=1;i<newN;i++)
        newCenters[i]=newCenters[i-1]+delta;
    
    // go through all weights and add them in the new interval
    // can have aliasing effects
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

//! interact one distribution with another
//! One can add/multiply/divide two distributions.
//! adding means, that one e.g. has two sticks with the likely lengths are the distribution centers
//! and one wants to get the total length of two sticks combined
//! multiply means, that e.g. one has a inaccuaratly measured current, and an inaccurate Resistor
//! and wants to know the likely distribution of the voltage= current*Resistance
//! divide means, that e.g. one has a inaccuaratly measured voltage, and an inaccurate Resistor
//! and wants to know the likely distribution of the current=voltage/Resistance
//!
//! this mulipurpose routine avoids code duplication
//!
//! all the possible pair choices of distribution 1 and 2 are evaluated, e.g. d1_sample*d2_sample
//! and put into an distribution of a reasonable size newN. This avoids that there are O^2 choices
//! and the dramatic slowdown/memory increase if the full resolution would be outputted
//! however this allows for aliasing effects
//! 
//! @param op: One out of "OpAdd, OpMul, OpDiv" to define what should be done
//! @param weights the likelyhood of the corresponding interval centers
//! @param newN the new size should be smaller!
inline TupleVecD operateDistributionsAndResample( const DistOps op,
    const VecD &centers1, const VecD &weights1, 
    const VecD &centers2, const VecD &weights2, unsigned newN=0)
{
    // have a default reasonable choice for newN
    if (newN==0)
        newN = std::min(centers1.size(), centers2.size()) +1;
    
    // allocate the centers/weights for the results
    VecD resCenters( newN, 0);
    VecD resWeights( newN, 0);
    
    //get limits for the operated centers
    double min1,max1,min2,max2;
    getLimits(centers1, min1, max1);
    getLimits(centers2, min2, max2);
    
    //find min/max of the resulting distribution
    double newLeft, newRight;
    switch(op)
    {
    case OpAdd: 
            newLeft  =min1+min2;
            newRight =max1+max2;
        break;
    case OpMul:
            newLeft  =min1*min2;
            newRight =max1*max2;
        break;
    case OpDiv: //this only works if there is no zero crossing in distribution 2
            newLeft  =min1/max2;
            newRight =max1/min2;
        break;
    }
    
    // define the new centers, weights are more compliacted
    double delta = (newRight-newLeft)/newN;
    
    resCenters[0]=newLeft+delta/2;
    for(unsigned i=1;i<newN;i++)
        resCenters[i]=resCenters[i-1]+delta; 
        
    // iterate through all combinations of points
    // to find the likelyhood of the operated new value 
    // and add it to the weights array
    for (unsigned i =0;i<centers1.size();i++)
    {
        const double &c1=centers1[i];
        const double &w1=weights1[i];
        for (unsigned j =0;j<centers2.size();j++)
        {
                const double &c2=centers2[j];
                const double &w2=weights2[j];
                
                double newC; //the new center pos
                switch(op)
                {
                case OpAdd: newC=c1+c2; break;
                case OpMul: newC=c1*c2; break;
                case OpDiv: newC=c1/c2; break;
                }
                
                //find the index at which this new Pos is in the results
                double indexPosF=(newC-newLeft)/delta;
                int indexPos= std::floor(indexPosF);
                double aliasingWeight=(indexPosF-indexPos); //how closely is it hitting at the integer pos
                if(indexPos<0)
                    indexPos=0;
                if(indexPos>=newN)
                    indexPos=newN-1;
                
                if ( indexPos==newN-1 )
                    resWeights[indexPos]+=w1*w2;
                else 
                {
                    //anti-aliasing
                    resWeights[indexPos]+=w1*w2 * (1-aliasingWeight);
                    resWeights[indexPos+1]+=w1*w2 * (aliasingWeight);
                }
                
        }
    }
    // polish the results weight so that they are mostly summed to 1
    normalizeVec(resWeights);
    
    return  {resCenters,resWeights};
}



#endif