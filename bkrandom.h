#ifndef _BKRANDOMCPP
#define _BKRANDOMCPP


#include <vector>
#include <random>
#include <tuple>
//#include <math.h>
#include <cmath>
#include <utility>

constexpr double pi = 3.14159265358979323846;

int testab(int a)
{
    return a+1;
}

void normalizeVec(std::vector<double> &v)
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

std::vector<double>  getNormal(double mean, double stddev, unsigned N)
{
    
    std::random_device r;
    std::default_random_engine e1(r());
    //std::mt19937 e2(seed2);
    
    std::normal_distribution d{mean, stddev};
    
    std::vector<double> result;
    result.resize(N);
    
    for(unsigned i=0;i<N;i++)
        result[i]=d(e1);
    
    return result;
}

std::tuple< std::vector<double>, std::vector<double> >
  getNormalW(double mean, double stddev, double sigmaMax, unsigned N)
{    
    double leftVal= mean-sigmaMax;
    double rightVal= mean+sigmaMax;
    double delta = (rightVal-leftVal)/N;
    std::vector<double> centers;
    centers.resize(N);
    std::vector<double> weights(N,0);
    
    
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

std::tuple< std::vector<double>, std::vector<double> >
    getTriW(double leftVal, double centerVal, double rightVal, unsigned N)
{
    double delta = (rightVal-leftVal)/N;
    std::vector<double> centers;
    centers.resize(N);
    std::vector<double> weights(N,0);
    
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

std::vector<double>  getRect(double leftVal, double rightVal, unsigned N)
{
    
    std::random_device r;
    std::default_random_engine e1(r());
    //std::mt19937 e2(seed2);
    
    std::uniform_real_distribution d(leftVal, rightVal);
    
    std::vector<double> result;
    result.resize(N);
    
    for(unsigned i=0;i<N;i++)
        result[i]=d(e1);
    
    return result;
}

std::tuple< std::vector<double>, std::vector<double> >
  getRectW(double leftVal, double rightVal, unsigned N)
{    
    double delta = (rightVal-leftVal)/N;
    std::vector<double> centers;
    centers.resize(N+2);
    std::vector<double> weights(N+2,1/double(N));
    
    
    centers[0]=leftVal-delta/2;
    weights[0]=0;
    for(unsigned i=1;i<N+1;i++)
        centers[i]=centers[i-1]+delta;
    centers[N+1]=centers[N]+delta;
    weights[N+1]=0;
    
    return {centers, weights};
}

bool comparatorPair ( const std::pair<double,double> &l, const std::pair<double,double> &r)
   { return l.first < r.first; }
  
 
std::tuple< std::vector<double>, std::vector<double> >
    resampleDistribution(std::vector<double> centers, std::vector<double> weights, unsigned newN)
{
    std::vector<double> newCenters;
    std::vector<double> newWeights;
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
std::tuple< std::vector<double>, std::vector<double> >
    addDistributions(std::vector<double> centers1, std::vector<double> weights1, 
                     std::vector<double> centers2, std::vector<double> weights2, unsigned newN)
{
    std::vector<double> resCenters;
    std::vector<double> resWeights;
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


std::tuple< std::vector<double>, std::vector<double> >
    multiplyDistributions(std::vector<double> centers1, std::vector<double> weights1, 
                          std::vector<double> centers2, std::vector<double> weights2, unsigned newN)
{
    std::vector<double> resCenters;
    std::vector<double> resWeights;
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

std::tuple< std::vector<double>, std::vector<double> >
    divideDistributions(std::vector<double> centers1, std::vector<double> weights1, 
                          std::vector<double> centers2, std::vector<double> weights2, unsigned newN)
{
    std::vector<double> resCenters;
    std::vector<double> resWeights;
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