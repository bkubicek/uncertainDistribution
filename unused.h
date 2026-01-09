int testab(int a)
{
    return a+1;
}

vDouble  getNormal(double mean, double stddev, unsigned N)
{
    
    std::random_device r;
    std::default_random_engine e1(r());
    //std::mt19937 e2(seed2);
    
    std::normal_distribution d{mean, stddev};
    
    vDouble result;
    result.resize(N);
    
    for(unsigned i=0;i<N;i++)
        result[i]=d(e1);
    
    return result;
}


vDouble getRect(double leftVal, double rightVal, unsigned N)
{
    
    std::random_device r;
    std::default_random_engine e1(r());
    //std::mt19937 e2(seed2);
    
    std::uniform_real_distribution d(leftVal, rightVal);
    
    vDouble result;
    result.resize(N);
    
    for(unsigned i=0;i<N;i++)
        result[i]=d(e1);
    
    return result;
}