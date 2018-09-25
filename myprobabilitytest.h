#ifndef MYPROBABILITYTEST_H
#define MYPROBABILITYTEST_H

#include "myoptional.h"
#include "mytests.h"
#include <random>

class Probability_test
{
public:
    template<class DataFunctionGenerator, class CumulativeProbabilityFunction>
    Op_void Probabity_test(const DataFunctionGenerator& f, const CumulativeProbabilityFunction& Prob, std::mt19937_64& mt, std::size_t sample_size, double p_accepted, double p_rejected, double sample_increase_factor,std::size_t max_samples)


    {

        std::size_t total_samples=0;
        double p=0;
        std::stringstream ss;
        while(total_samples<max_samples)
        {
            auto sample=f(mt,sample_size);
            total_samples+=sample_size;
            double p=Prob(sample,ss);
            if (p>p_accepted)
            {
                return Op_void(true,
                               "Accepted with sample size="+std::to_string(sample_size)+
                               ", p("+std::to_string(p)+")>p_accepted ("+std::to_string(p_accepted)+"); "+ss.str());
            }
            else if (p<p_rejected)
            {
                return Op_void(false,
                               "Rejected with sample size="+std::to_string(sample_size)+
                               ", p("+std::to_string(p)+")<p_rejected ("+std::to_string(p_rejected)+"); "+ss.str());

            }
            else sample_size*=sample_increase_factor;
        }
        return Op_void(false,
                       "Undecided after "+std::to_string(total_samples)+" total samples. Last with sample size="+
                       std::to_string(sample_size)+"; p_accepted ("+std::to_string(p_accepted)+")>p("+std::to_string(p)+
                       +")<p_rejected ("+std::to_string(p_rejected)+"); "+ss.str());
    }


    template<class DataFunctionGenerator, class CumulativeProbabilityFunction>
    Op_void Probabity_test(const DataFunctionGenerator& f, const CumulativeProbabilityFunction& Prob, std::mt19937_64& mt)
    {
         return Probabity_test(f,Prob,mt,sample_size_,p_accepted_,p_rejected_,sample_increase_factor_,max_samples_);
    }

    Probability_test(std::size_t sample_size, double p_accepted, double p_rejected, double sample_increase_factor,std::size_t max_samples)
        :sample_size_{sample_size},p_accepted_{p_accepted},p_rejected_{p_rejected}, sample_increase_factor_{sample_increase_factor},max_samples_{max_samples}
    {}

private:
     std::size_t sample_size_;
     double p_accepted_;
     double p_rejected_;
     double sample_increase_factor_;
     std::size_t max_samples_;
};





#endif // MYPROBABILITYTEST_H
