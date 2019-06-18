#ifndef MATRIXERROR_H
#define MATRIXERROR_H
#include "myError.h"

#include "Matrix.h"

template<typename T, class Norm, bool diff>
class Error<M_Matrix<T>,Norm,diff>
{
    M_Matrix<T> center_;
    M_Matrix<T> norm_;
public:
    auto& center()const {return center_;}
    auto& center() {return center_;}
    auto& norm()const {return norm_;}
    auto& norm() {return norm_;}

    auto& error()const {return Norm::pow_inv(norm_);}

    auto nrows()const {return center().nrows();}
    auto ncols()const {return center().ncols();}
    auto size()const {return center().size();}
    auto type()const {return center().type();}

    Error<double,Norm,diff> getvalue() const {
        assert(size() == 1);
        return Error<double,Norm,diff>(center()[0],norm()[0]);
    }


    Error<T,Norm,diff> operator[](std::size_t i)const {return Error<T,Norm,diff>(center()[i],norm()[i]);}

    Error<T,Norm,diff> operator()(std::size_t i, std::size_t j)const {return Error<T,Norm,diff>(center()(i,j),norm()(i,j));}

    template <class F> Error apply(const F &f) const& {
        M_Matrix<T> c(nrows(), ncols(), type());
        M_Matrix<T> n(nrows(), ncols(), type());
        for (std::size_t i = 0; i < size(); ++i)
        {
            auto v=std::invoke(f, Error<T,Norm,diff>(center()[i],norm()[i]));
            c[i]=v.center();
            n[i]=v.norm();
        }
        return Error(std::move(c),std::move(n));
    }



    Error()=default;

    Error(M_Matrix<T>&& val, M_Matrix<T>&& n): center_{std::move(val)},norm_{std::move(n)} {
        assert(norm()>=0.0);
    }
    Error(const M_Matrix<T>& val): center_{val},norm_{val.apply(&Norm::pow)*Norm::pow(std::numeric_limits<T>::epsilon()*10)} {}

    Error(M_Matrix<T>&& val): center_{std::move(val)},norm_{center().apply(&Norm::pow)*Norm::pow(std::numeric_limits<T>::epsilon()*10)} {}

    Error(std::size_t nrows_, std::size_t ncols_, typename M_Matrix<T>::TYPE t, T x)
        : center_(nrows_,ncols_,t,x),norm_(nrows_,ncols_,t,Norm::pow(x*std::numeric_limits<double>::epsilon())){}

    Error(std::size_t nrows_, std::size_t ncols_, typename M_Matrix<T>::TYPE t)
        : center_(nrows_,ncols_,t),norm_(nrows_,ncols_,t){}

    Error(std::size_t nrows_, std::size_t ncols_, const T& x)
        : center_(nrows_,ncols_,x),norm_(nrows_,ncols_,Norm::pow(x*std::numeric_limits<double>::epsilon())){}


    void set(std::size_t i, std::size_t j, const Error<T,Norm,diff>& value)
    {
        center_(i,j)=value.center();
        norm_(i,j)=value.norm();
    }

    void add(std::size_t i, std::size_t j, const Error<T,Norm,diff>& value)
    {
        center_(i,j)+=value.center();
        norm_(i,j)+=value.norm();
    }



    void set(std::size_t i, const Error<T,Norm,diff>& value)
    {
        center_[i]=value.center();
        norm_[i]=value.norm();
    }


    Error& operator +=(const Error& other){
        center_+=other.center();
        norm_+=other.norm();
        return *this;
    }
    Error& operator -=(const Error& other){
        center_-=other.center();
        norm_+=other.norm();
        return *this;
    }



};
template<typename T, class Norm, bool diff>
std::ostream& operator<<(std::ostream& os, const Error<M_Matrix<T>,Norm,diff>& x)
{
    os<<x.center()<<"+/-"<<x.norm()<<"\n";
    return os;
}

template<typename T, class Norm, bool diff>
bool operator==(const Error<M_Matrix<T>,Norm,diff>& one, const Error<M_Matrix<T>,Norm,diff>& other)
{
    if (one.center().nrows()!=other.center().nrows()) return false;
    if( one.center().size()!=other.center().size()) return false;
    for (std::size_t i=0; i<one.center().size(); ++i)
        if(!(Error<T,Norm,diff>(one.center()[i],one.norm()[i])==Error<T,Norm,diff>(other.center()[i],other.norm()[i])))
            return false;
    return true;
}
template<typename T, class Norm, bool diff>
bool operator==(const Error<M_Matrix<T>,Norm,diff>& one, const M_Matrix<T>& other)
{
    if (one.center().nrows()!=other.nrows()) return false;
    if( one.center().ncols()!=other.ncols()) return false;
    for (std::size_t i=0; i<one.center().nrows(); ++i)
        for (std::size_t j=0; j<one.center().ncols(); ++j)
            if(!(Error<T,Norm,diff>(one.center()(i,j),one.norm()(i,j))==other(i,j)))
            return false;
    return true;
}


template <class T, class Norm, bool diff,typename...Ts>
auto are_Equal_v(const Error<M_Matrix<T>,Norm,diff> &one, const Error<M_Matrix<T>,Norm,diff>& other, double /*eps*/, [[maybe_unused]] double epsf,std::ostream &os, Ts...context)->std::enable_if_t<!std::is_pointer_v<T>,bool>
{
    if (!(one==other))
    {
        using namespace io;
        (os<<...<<context);
        return false;
    }
    else return true;
}



template<typename T, class Norm, bool diff>
Error<M_Matrix<T>,Norm,diff> operator+(const Error<M_Matrix<T>,Norm,diff>& one, const Error<M_Matrix<T>,Norm,diff>& other)
{
    Error<M_Matrix<T>,Norm,diff> out(one);
    out+=other;
    return out;
}

template<typename T, class Norm, bool diff>
Error<M_Matrix<T>,Norm,diff> operator-(const Error<M_Matrix<T>,Norm,diff>& one, const Error<M_Matrix<T>,Norm,diff>& other)
{
    Error<M_Matrix<T>,Norm,diff> out(one);
    out-=other;
    return out;
}

template<typename T, class Norm, bool diff>
Error<M_Matrix<T>,Norm,diff> operator*(const Error<M_Matrix<T>,Norm,diff>& one, const Error<M_Matrix<T>,Norm,diff>& other)
{
    auto val=one.center()*other.center();
    auto norm=one.center().apply(&Norm::pow)*other.norm()+one.norm()*other.center().apply(&Norm::pow);

    return Error<M_Matrix<T>,Norm,diff>(std::move(val),std::move(norm));
   }

   template<typename T, class Norm, bool diff>
   Error<M_Matrix<T>,Norm,diff> operator*(const Error<M_Matrix<T>,Norm,diff>& one, const Error<double,Norm,diff>& other)
   {
       auto val=one.center()*other.center();
       auto norm=one.center().apply(&Norm::pow)*other.norm()+one.norm()*Norm::pow(other.center());

       return Error<M_Matrix<T>,Norm,diff>(std::move(val),std::move(norm));
   }

   template<typename T, class Norm, bool diff>
   Error<M_Matrix<T>,Norm,diff> operator*( const Error<double,Norm,diff>& first,const Error<M_Matrix<T>,Norm,diff>& second)
   {
       auto val=first.center()*second.center();
       auto norm=first.norm()*second.center().apply(&Norm::pow)+Norm::pow(first.center())*second.norm();

       return Error<M_Matrix<T>,Norm,diff>(std::move(val),std::move(norm));
   }





   template<typename T, class Norm, bool diff>
   Error<M_Matrix<T>,Norm,diff> operator*(const Error<M_Matrix<T>,Norm,diff>& one, const M_Matrix<T>& other)
   {
       auto val=one.center()*other;
       auto norm=one.norm()*other.apply(&Norm::pow);

       return Error<M_Matrix<T>,Norm,diff>(std::move(val),std::move(norm));
   }

   template<typename T, class Norm, bool diff>
   Error<M_Matrix<T>,Norm,diff> operator*(const M_Matrix<T>& one, const Error<M_Matrix<T>,Norm,diff>& other)
   {
       auto val=one*other.center();
       auto norm=one.apply(&Norm::pow)*other.norm();

       return Error<M_Matrix<T>,Norm,diff>(std::move(val),std::move(norm));
   }

   template<typename T, class Norm, bool diff>
   Error<M_Matrix<T>,Norm,diff> operator*(const Error<M_Matrix<T>,Norm,diff>& one, const T& other)
   {
       auto val=one.center()*other;
       auto norm=one.norm()*Norm::pow(other);

       return Error<M_Matrix<T>,Norm,diff>(std::move(val),std::move(norm));
   }

   template<typename T, class Norm, bool diff>
   Error<M_Matrix<T>,Norm,diff> operator*(const T& one, const Error<M_Matrix<T>,Norm,diff>& other)
   {
       auto val=one*other.center();
       auto norm=one.apply(&Norm::pow)*other.norm();

       return Error<M_Matrix<T>,Norm,diff>(std::move(val),std::move(norm));
   }




 template <class T, class Norm, bool diff>
   myOptional_t<Error<M_Matrix<T>,Norm,diff>> inv(const Error<M_Matrix<T>,Norm,diff> &x) {
       typedef myOptional_t<Error<M_Matrix<T>,Norm,diff>> Op;

       auto invx = matrix_inverse::Matrix_inverse(x.center());
       if (!invx.has_value())
           return Op(false, "cannot invert x " + invx.error());
       else {
           auto inv_x_n=invx.value().first.apply(&Norm::pow);
           auto norm = inv_x_n* x.norm() * inv_x_n;
           return Op(
               Error<M_Matrix<T>,Norm,diff>(std::move(invx).value().first, std::move(norm)));
       }
   }

   template <> struct E<Matrix_Decompositions::eigensystem_type> {
       typedef std::tuple<Error<M_Matrix<double>,::norm_1,true>, Error<M_Matrix<double>,::norm_1,true>,
                          Error<M_Matrix<double>,::norm_1,true>,M_Matrix<double>,M_Matrix<double>>
           type;
   };

   namespace Matrix_Decompositions {



   template <class Norm, bool diff>
    auto EigenSystem_full_real_eigenvalue(const Error<M_Matrix<double>,Norm,diff> &Ex) {
       typedef myOptional_t<
           std::tuple<Error<M_Matrix<double>,Norm,diff>, Error<M_Matrix<double>,Norm,diff>,
                      Error<M_Matrix<double>,Norm,diff>, M_Matrix<double>, M_Matrix<double>>>
           Op;
       auto res = Matrix_Decompositions::EigenSystem_full_real_eigenvalue(Ex.center());
       if (!res)
           return Op(false,
                     "Derivative error cannot calculete function value" + res.error());
       else {
           auto [VR, landa, VL, CL,CV] = std::move(res).value();

           Error<M_Matrix<double>,Norm,diff> Elanda(std::move(landa),CL.apply(&Norm::pow)*Norm::pow(std::numeric_limits<double>::epsilon()));


           M_Matrix<double> VRN(VR.nrows(), VR.ncols(), VR.type());
           for (std::size_t i = 0; i < VR.nrows(); ++i) {
               for (std::size_t j = 0; j < VR.ncols(); ++j) {
                   VRN(i,j)=Norm::pow(CV[j]*std::numeric_limits<double>::epsilon());
               }

               }
           Error<M_Matrix<double>,Norm,diff> DVR(std::move(VR), std::move(VRN));
           auto DVL = inv(DVR);
           if (!DVL)
               return Op(false, " fails to invert the left eigenvector");

           return Op(std::tuple(std::move(DVR), std::move(Elanda),
                                std::move(DVL).value(), std::move(CL), std::move(CV)));
       }
   }
   } // namespace Matrix_Decompositions


   namespace Matrix_Unary_Transformations
   {
   template <typename T, class Norm,bool diff> Error<M_Matrix<T>,Norm,diff> diag(const Error<M_Matrix<T>,Norm,diff> &x) {
       return Error<M_Matrix<T>,Norm,diff>(diag(x.center()),diag(x.norm()));
      }
   }
   template <typename T, class Norm,bool diff> Error<M_Matrix<T>,Norm,diff> quadraticForm_XTX(const Error<M_Matrix<T>,Norm,diff> &x) {
       M_Matrix<T> out_norm(x.center().ncols(), x.center().ncols(), Matrix_TYPE::SYMMETRIC, T{});
       for (std::size_t i = 0; i < x.ncols(); ++i)
           for (std::size_t j = 0; j < i + 1; ++j)
               for (std::size_t k = 0; k < x.nrows(); ++k)
                   out_norm(i, j) += Norm::pow(x.center()(k, i)) * x.norm()(k, j)+
                       x.norm()(k, i) * Norm::pow(x.center()(k, j));
       return Error<M_Matrix<T>,Norm,diff>(quadraticForm_XTX(x.center()),std::move(out_norm));
   }

   template <class Norm,bool diff>
   auto TranspMult
       (const Error<M_Matrix<double>,Norm,diff> &one, const Error<M_Matrix<double>,Norm,diff> &other)
   {
       M_Matrix<double> one_norm=one.center().apply(&Norm::pow);
       M_Matrix<double> other_norm=other.center().apply(&Norm::pow);
       return Error<M_Matrix<double>,Norm,diff>(Matrix_Binary_Transformations::TranspMult(one.center(),other.center()),
Matrix_Binary_Transformations::TranspMult(one.norm(),other_norm)+Matrix_Binary_Transformations::TranspMult(one_norm,other.norm()));

   }

   template <typename T, typename S,class Norm,bool diff>
   auto quadraticForm_BT_A_B
       (const Error<M_Matrix<T>,Norm,diff> &one, const Error<M_Matrix<S>,Norm,diff> &other)
           -> Error<M_Matrix<decltype(std::declval<T>() * std::declval<S>())>,Norm,diff> {

           typedef decltype (std::declval<T>()*std::declval<S>()) R;
           auto one_norm=one.center().apply(&Norm::pow);
           auto other_norm=other.center().apply(&Norm::pow);
           return Error<M_Matrix<R>,Norm,diff>(Matrix_Binary_Transformations::quadraticForm_BT_A_B(one.center(),other.center()),
                                                 Matrix_Binary_Transformations::quadraticForm_BT_A_B(one.norm(),other_norm)
                                                     +Matrix_Binary_Transformations::quadraticForm_BT_A_B(one_norm,other.norm()));

       }

   template <class F, class Norm,bool diff>
   Error<double,Norm,diff> accumulate(const Error<M_Matrix<double>,Norm,diff> &x, const F &f, double start) {
       Error<double,Norm,diff> out(start);
       for (std::size_t i=0; i<x.size(); ++i)
           out = f(out,x.center()[i], x.norm()[i]);
       return out;
   }



   template <typename T, class Norm,bool diff> Error<double,Norm,diff> max(const Error<M_Matrix<T>,Norm,diff>  &x) {
       using std::max;
       return accumulate(x, [](Error<double,Norm,diff> one, double val, double norm) {
               if (one.center()>val) return one; else return Error<double, Norm, diff>(val,norm); },
           -std::numeric_limits<double>::infinity());
   }


   template <class Norm,bool diff>
   Error<double,Norm,diff> min(const Error<M_Matrix<double>,Norm,diff>  &x) {
       using std::min;
       return accumulate(x, [](Error<double,Norm,diff> one, double val, double norm) {
               if (one.center()<val) return one; else return Error<double, Norm, diff>(val,norm); },
           std::numeric_limits<double>::infinity());
   }

   template <typename T, class Norm,bool diff,std::enable_if_t<!is_Matrix_v<std::decay_t<T>>, int> = 0>
   Error<double,Norm,diff> elemDivSafe(const Error<T,Norm,diff>& x, const Error<T,Norm,diff>& y,
                             double eps = std::numeric_limits<double>::epsilon()) {
       if (std::abs(y.center()) > eps)
           return x / y;
       else
           return Error<double,Norm,diff>(0);
   }

   template <typename T, typename S, class Norm,bool diff>
   auto elemDivSafe(const Error<M_Matrix<T>,Norm,diff>   &x, const Error<M_Matrix<S>,Norm,diff>   &y,
                    [[maybe_unused]] double eps = std::numeric_limits<double>::epsilon()  ) {
       typedef decltype(std::declval<T>() * std::declval<S>()) R;
       //double norm = Matrix_Unary_Functions::norm_1(y);
       Error<M_Matrix<R>,Norm,diff> out(x.nrows(),x.ncols(),x.type());
       for (std::size_t i=0; i<x.size(); ++i)
           out.set(i,elemDivSafe(x[i],y[i]));

       return out;
   }

   template <typename T, typename S, class Norm,bool diff>
   auto elemMult(const Error<M_Matrix<T>,Norm,diff>   &x, const Error<M_Matrix<S>,Norm,diff>   &y) {
       typedef decltype(std::declval<T>() * std::declval<S>()) R;
       assert(x.nrows()==y.nrows());
       assert(x.ncols()==y.ncols());
       //double norm = Matrix_Unary_Functions::norm_1(y);
       Error<M_Matrix<R>,Norm,diff> out(x.nrows(),x.ncols(),x.type());
       for (std::size_t i=0; i<x.size(); ++i)
           out.set(i,x[i]*y[i]);

       return out;
   }

   template <typename T, class Norm,bool diff>
   auto TransposeSum(const Error<M_Matrix<T>,Norm,diff>   &x) {
       return Error<M_Matrix<T>,Norm,diff>(TransposeSum(x.center()), TransposeSum(x.norm()));


   }




   template <typename T, class Norm,bool diff>
   Error<M_Matrix<double>,Norm,diff> expm_taylor(const Error<M_Matrix<double>,Norm,diff> & x, std::size_t order=6)
   {
       Error<M_Matrix<double>,Norm,diff> out=x+Matrix_Generators::eye<double>(x.ncols());
       Error<M_Matrix<double>,Norm,diff> xr=x;
       double a=1.0;
       for (std::size_t n=2; n<order; ++n)
       {
           a/=n;
           xr=xr*x;
           out=out+xr*a;
       }
       double error=expm_taylor_error(center(x),order);
       auto error_norm=M_Matrix<double>(x.nrows(),x.ncols(),Matrix_TYPE::FULL,Norm::pow(error));
       out.norm()+=error_norm;
       return out;

   }

   template <typename T, class Norm,bool diff>
   Error<M_Matrix<T>,Norm,diff>  calc_Commutator_Taylor_Series(const Error<M_Matrix<T>,Norm,diff> & GP, const Error<M_Matrix<T>,Norm,diff>  & Qx,std::int8_t order)
   {
       Error<M_Matrix<T>,Norm,diff> dGn = GP;
       auto Gtot = dGn;
       double a = 1.0;
       for (std::size_t n = 2; n + 1 < order; ++n) {
           a /= n;
           dGn = Qx * dGn - dGn * Qx;
           Gtot += dGn * a;
       }
       double error=calc_Commutator_Taylor_Series_error(center(Qx),center(GP),order);
       auto error_norm=M_Matrix<double>(Qx.nrows(),Qx.ncols(),Matrix_TYPE::FULL,Norm::pow(error));
       Gtot.norm()+=error_norm;
       return Gtot;

   }

   template <typename T, class Norm,bool diff>
   void next_derivative(std::vector<Error<M_Matrix<T>,Norm,diff>> &dGvar_n, Error<M_Matrix<T>,Norm,diff> &Tau_n,
                        const Error<M_Matrix<T>,Norm,diff> &Q, const Error<M_Matrix<T>,Norm,diff> &GP) {
       Error<M_Matrix<T>,Norm,diff> neg_2Q(Q * (-2.0));
       for (auto &e : dGvar_n)
           e = e * neg_2Q;
       Tau_n = Q * Tau_n + Tau_n * Q;
       dGvar_n.push_back(Tau_n * GP);
   }

   template <typename T, class Norm,bool diff>
   Error<M_Matrix<T>,Norm,diff>  calc_square_Commutator_Taylor_Series(const Error<M_Matrix<T>,Norm,diff> &Qx, const Error<M_Matrix<T>,Norm,diff> &P,
                                                    const Error<M_Matrix<T>,Norm,diff> &G, std::size_t order) {
       auto Tau_n = G;
       auto GP = G * P;
       std::vector<double> coeff(1, 1.0);
       std::vector<Error<M_Matrix<T>,Norm,diff> > dGvar_n(1, G * GP);
       auto Gvar_tot = dGvar_n[0];
       double a = 2.0;
       for (std::size_t n = 3; n + 1 < order; ++n) {
           a /= n;
           pascal_triangle(coeff);
           next_derivative(dGvar_n, Tau_n, Qx, GP);
           for (std::size_t i = 0; i < coeff.size(); ++i)
               Gvar_tot += dGvar_n[i] * (coeff[i] * a);
       }
       double error=calc_square_Commutator_Taylor_Series_error(center(Qx),center(G),order);
       auto error_norm=M_Matrix<double>(Qx.nrows(),Qx.ncols(),Matrix_TYPE::FULL,Norm::pow(error));
       Gvar_tot.norm()+=error_norm;
       return Gvar_tot;
   }



   template <typename T, class Norm,bool diff>
   auto mean(const Error<M_Matrix<T>,Norm,diff> &A, const Error<M_Matrix<T>,Norm,diff> &B)
        {
       assert(A.nrows() == B.nrows());
       assert(B.ncols() == A.ncols());
       auto out=(A+B)*0.5;
       out.norm()+=(A-B).apply(&Norm::pow);
       return out;

   }

   template <typename T, class Norm,bool diff>
   auto variance(const Error<M_Matrix<T>,Norm,diff> &A, const Error<M_Matrix<T>,Norm,diff> &B)
  {
       assert(A.nrows() == B.nrows());
       assert(B.ncols() == A.ncols());
       auto out= (A-B).apply([](auto x){return sqr(x)/2.0;});
       out.norm()+=out.center().apply(Norm::pow);
       return out;

   }


#endif // MATRIXERROR_H
