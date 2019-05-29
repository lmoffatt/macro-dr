#ifndef MATRIXCALCULATIONS_H
#define MATRIXCALCULATIONS_H

#include "Matrix.h"
#include "matrixderivative.h"
#include "matrixerror.h"

struct calculate_eigenvalues
{
    auto operator()(const M_Matrix<double>& x)
    {
        assert(x.nrows()==x.nrows());
        typedef myOptional_t<Matrix_Decompositions::eigensystem_type> Op;

        Error<M_Matrix<double>,::norm_1,true> Ex(x);
        auto eig=Matrix_Decompositions::EigenSystem_full_real_eigenvalue(Ex);
        if (eig.has_value())
        {
            auto [VR, landa, VL, CL,CV]=std::move(eig).value();
            assert(VR*VL==Matrix_Generators::eye<double>(VR.ncols()));
            assert(VR*landa*VL==x);
            return Op(Matrix_Decompositions::eigensystem_type(std::move(VR.center()),std::move(landa.center()),std::move(VL.center()),
                                                              std::move(CL), std::move(CV)));

        }
        else {
            return Op(false,eig.error());
        }
    }
    template<class Norm, bool diff>
    auto operator()(const Error<M_Matrix<double>,Norm,diff>& x)
    {
        assert(x.nrows()==x.ncols());
       return Matrix_Decompositions::EigenSystem_full_real_eigenvalue(x);
    }

};






#endif // MATRIXCALCULATIONS_H
