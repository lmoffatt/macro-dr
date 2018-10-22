#ifndef MATRIXDERIVATIVE_H
#define MATRIXDERIVATIVE_H
#include "Matrix.h"




template<>
class Derivative<double>
{
    double f_;
    M_Matrix<double> const * x_;
    M_Matrix<double> dfdx_;
public:
    typedef double primitive_type;
    constexpr static auto const className=my_static_string("Derivative_")+my_trait<primitive_type>::className;
    typedef   Derivative<primitive_type> self_type ;

    static auto get_constructor_fields()
    {
        return std::make_tuple(
                    grammar::field(C<self_type>{},"f",&self_type::f),
                    grammar::field(C<self_type>{},"elogL",&self_type::x),
                    grammar::field(C<self_type>{},"vlogL",&self_type::dfdx)
                    );
    }
    double f()const { return f_;}
    M_Matrix<double> const & x() const { return *x_;}
    M_Matrix<double> const& dfdx()const { return dfdx_;}
    Derivative(double fx,  const M_Matrix<double>& myx, M_Matrix<double>&& der): f_{fx},x_{&myx}, dfdx_{std::move(der)}{}
    Derivative()=default;
    Derivative(const M_Matrix<double>& myx): f_{0},x_{&myx}, dfdx_{myx.nrows(),myx.ncols(),myx.type(),0.0}{}
    Derivative& operator+=(const Derivative& other)
    {
        assert(x()==other.x());
        f_+=other.f();
        dfdx_+=other.dfdx();
        return *this;
    }

};






template<>
class Derivative<M_Matrix<double>>
{
    M_Matrix<double> f_;
    M_Matrix<double>const * x_;
    M_Matrix<M_Matrix<double>> dfdx_;
public:
    Derivative& set_by_f_index(std::size_t i, std::size_t j, const Derivative<double>& dfij)
    {
        assert(x()==dfij.x());
        f()(i,j)=dfij.f();
        for (std::size_t i=0; i<dfdx().size(); ++i)
            dfdx()[i](i,j)=dfij.dfdx()[i];
        return *this;


    }

    M_Matrix<double> const & f()const { return f_;}
    M_Matrix<double>  & f(){ return f_;}
    M_Matrix<double> const & x() const { return *x_;}
    M_Matrix<M_Matrix<double>> const& dfdx()const { return dfdx_;}
    M_Matrix<M_Matrix<double>> & dfdx() { return dfdx_;}
    Derivative(M_Matrix<double>&& fx, const  M_Matrix<double>& myx, M_Matrix<M_Matrix<double>>&& der): f_{std::move(fx)},x_{&myx}, dfdx_{std::move(der)}{}

    template<class M_f, class M_d>
    Derivative(M_f fx, const  M_Matrix<double>& myx, M_d der): f_{std::forward<M_f>(fx)},x_{&myx}, dfdx_{std::forward<M_d>(der)}{}



    Derivative(std::size_t nrows, std::size_t ncols, M_Matrix<double>::TYPE t, const M_Matrix<double>& x)
        :f_{nrows,ncols,t},x_{&x},dfdx_{x.nrows(),x.ncols(),x.type(),M_Matrix<double>(nrows,ncols,t)}{}
    Derivative()=default;

    M_Matrix<Derivative<double>> to_Matrix()const
    {
        M_Matrix<Derivative<double>> out(f().nrows(),f().ncols(),f().type());
        for (std::size_t i=0; i<f().size(); ++i)
        {
            M_Matrix<double> dfdxi(x());
            for (std::size_t j=0; j<x().size(); ++j)
                dfdxi[j]=dfdx()[j][i];
            out[i]=Derivative<double>(f()[i],x(),std::move(dfdxi));
        }
        return out;
    }
    Derivative (const M_Matrix<Derivative<double>>& x):f_{x.nrows(),x.ncols(),x.type()},x_{&x[0].x()},
    dfdx_{x[0].x().nrows(),x[0].x().ncols(),x[0].x().type()}
    {
        for (std::size_t i=0; i<x.size(); ++i)
        {
            f_[i]=x[i].f();
            for (std::size_t j=0; j<x[0].x().size(); ++j)
                dfdx_[j][i]=x[i].dfdx()[j];
        }

    }


    template<class F>
    Derivative apply(const F& f)
    {
        auto M=to_Matrix();
        return Derivative(M.apply(f));
    }
    Derivative& operator+=(const Derivative& other)
    {
        assert(x()==other.x());
        f_+=other.f();
        dfdx_+=other.dfdx();
        return *this;
    }
};

template<class T>
Derivative<T> operator +(const Derivative<T>& x,  Constant<T>&& c)
{
    return Derivative<T>(x.f()+c.value,x.x(),x.dfdx());
}



auto exp(D,double x)
{
    return std::exp(x);
}


auto exp(const Derivative<double>& x)
{
    return Derivative<double>(std::exp(x.f()),x.x(),exp(D(),x.f())*x.dfdx());
}

auto exp(const Derivative<M_Matrix<double>>& x)
{
    return Derivative<M_Matrix<double>>(exp(x.f()),x.x(),elemMult_a(x.dfdx(),exp(x.f())));
}

auto abs(D,double x)
{
    if (x>=0) return  1.0;
    else return -1.0;
}



auto abs(const Derivative<double>& x)
{
    return Derivative<double>(std::abs(x.f()),x.x(),abs(D(),x.f())*x.dfdx());
}




template <class> struct is_Derivative: public std::false_type{};

template <typename T>
struct is_Derivative<Derivative<T>>: public std::true_type{};

template<class D>
constexpr static bool is_Derivative_v=is_Derivative<D>::value;

template<typename T>
auto elemMult(const Derivative<M_Matrix<T>>& x,const Derivative<M_Matrix<T>>& y)
{
    assert(x.x()==y.y());
    return Derivative<M_Matrix<T>
            >(elemMult(x.f(),y.f()),x.x(),
              elemMult_a(x.dfdx(),y.f())+elemMult_a(x.f(),y.dfdx()));
}

template<typename T>
auto elemDiv(const Derivative<M_Matrix<T>>& x,const Derivative<M_Matrix<T>>& y)
{
    assert(x.x()==y.y());
    auto f=elemDiv(x.f(),y.f());
    auto dfdx=elemDiv_a(x.dfdx(),y.f());
    auto dfdy=elemMult_a(f,elemDiv_a(y.dfdx(),y.f()));

    return Derivative<M_Matrix<T>>(std::move(f),x.x(),dfdx+dfdy);
}

template<typename T>
auto elemDivSafe(const Derivative<M_Matrix<T>>& x,const Derivative<M_Matrix<T>>& y, double eps=std::numeric_limits<double>::epsilon())
{
    assert(x.x()==y.y());
    auto f=elemDivSafe(x.f(),y.f(),eps);
    auto dfdx=elemDivSafe_a(x.dfdx(),y.f(),eps);
    auto dfdy=elemMult_a(elemDivSafe_a(y.dfdx(),y.f(),eps),f);

    return Derivative<M_Matrix<T>>(std::move(f),x.x(),dfdx+dfdy);
}



template< class T>
auto compose (const Derivative<T>& dfdx, const Derivative<M_Matrix<double>>& dxdy)
{
    assert(dfdx.x()==dxdy.f());
    return Derivative<T>(dfdx.f(), dxdy.x(), dfdx.dfdx()*dxdy.dfdx());
}


template<typename T, typename S>
auto operator*(const Derivative<T>& x,S t)
->Derivative<std::enable_if_t<!is_Derivative_v<S>,T>
>
{
    return  Derivative<T>(x.f()*t,x.x(),x.dfdx()*t);
}

template<typename T, typename S>
auto operator*(S t,const Derivative<T>& x)
->Derivative<std::enable_if_t<!is_Derivative_v<S>,T>
>
{
    return  Derivative<T>(t*x.f(),x.x(),t*x.dfdx());
}





template<typename T>
Derivative<T> operator+(const Derivative<T>& x,const Derivative<T>& y)

{
    assert(x.x()==y.x());
    return  Derivative<T>(x.f()+y.f(),x.x(),x.dfdx()+y.dfdx());
}

template<class T>
Derivative<M_Matrix<T>>  TransposeSum(const Derivative<M_Matrix<T>>& x)
{
    assert(x.f().ncols()==x.f().nrows());
    auto f=TransposeSum(x.f());
    auto df=x.dfdx().apply([](auto& x){return TransposeSum(x);});
    return Derivative<M_Matrix<T>>(std::move(f),x.x(),std::move(df));
}

template<class T>
Derivative<M_Matrix<T>>  Transpose(const Derivative<M_Matrix<T>>& x)
{
    assert(x.f().ncols()==x.f().nrows());
    auto f=Transpose(x.f());
    auto df=x.dfdx().apply([](auto& x){return Transpose(x);});
    return Derivative<M_Matrix<T>>(std::move(f),x.x(),std::move(df));
}


template<typename T>
Derivative<M_Matrix<T>> operator-(const Derivative<M_Matrix<T>>& x,const Derivative<M_Matrix<T>>& y)
{
    assert(x.x()==y.x());
    return  Derivative<M_Matrix<T>>(x.f()-y.f(),x.x(),x.dfdx()-y.dfdx());
}



template<typename T>
Derivative<T>
operator *
(const Derivative<T>& one,
 const Derivative<T>& other)
{
    assert(one.x()==other.x());
    return  Derivative<T>(one.f()*other.f(),one.x(),one.dfdx()*other.f()+one.f()*other.dfdx());
}


template <class T>
myOptional_t<Derivative<M_Matrix<T>>>
inv(const Derivative<M_Matrix<T>>& x)
{
    typedef  myOptional_t<Derivative<M_Matrix<T>>> Op;

    auto invx=matrix_inverse::Matrix_inverse(x.f());
    if (!invx.has_value()) return Op(false,"cannot invert x "+invx.error());
    else
    {
        //        auto fxxf=Permute<T,1,3,0,2>(x.dfdx());
        auto df=-invx.value()*x.dfdx()*invx.value();
        //      df=Permute<T,2,0,3,1>(df);
        return Op(Derivative<M_Matrix<T>>(std::move(invx).value(),x.x(),std::move(df)));
    }

}
namespace Matrix_Decompositions {

auto EigenSystem_full_real_eigenvalue(const Derivative<M_Matrix<double>>& Dx)
{
    typedef myOptional_t<std::tuple<Derivative<M_Matrix<double>>,Derivative<M_Matrix<double>>,Derivative<M_Matrix<double>>>>
            Op;
    auto res=Matrix_Decompositions::EigenSystem_full_real_eigenvalue(Dx.f());
    if (!res)
        return Op(false,"Derivative error cannot calculete function value"+res.error());
    else
    {
        auto [W, landa, V]=std::move(res).value();
                // auto fxxf=Permute<double,1,3,0,2>(Dx.dfdx());
                auto fxxf=Dx.dfdx();
                M_Matrix<M_Matrix<double>> dlanda(landa.nrows(),Dx.x().nrows(),
                M_Matrix<double>(landa.ncols(),Dx.x().ncols(),0));
                M_Matrix<M_Matrix<double>> dV(V.nrows(),Dx.x().nrows(),
                M_Matrix<double>(V.ncols(),Dx.x().ncols(),0));
                for (std::size_t i=0; i<landa.nrows(); ++i)
        {
            auto v=V(i,":");
            auto vT=Transpose(v);
            for (std::size_t j1=0; j1<Dx.x().nrows(); ++j1)
                for (std::size_t j2=0; j2<Dx.x().ncols(); ++j2)
                    dlanda(i,j1)(i,j2)=(vT*fxxf(j1,j2)*v).getvalue();

            auto vinv=pinv(landa(i,i)*Matrix_Generators::eye<double>(landa.nrows())-Dx.f());
            if (!vinv) return Op(false, " Error calculating Derivative, problem with pseudoinverse "+vinv.error());
            for (std::size_t j1=0; j1<Dx.x().nrows(); ++j1)
                for (std::size_t j2=0; j2<Dx.x().ncols(); ++j2)
                {
                    auto dvi=vinv.value()*fxxf(j1,j2)*v;
                    for (std::size_t i2=0; i2<V.ncols(); ++i2)
                        dV(i,j1)(i2,j2)=dvi(i2,0);
                }

        }
        Derivative<M_Matrix<double>> Dlanda(std::move(landa),Dx.x(),std::move(dlanda));
        Derivative<M_Matrix<double>> DV(std::move(V),Dx.x(),std::move(dV));
        auto DW=inv(DV);
        if (!DW) return Op(false, " fails to invert the left eigenvector");

        return Op(std::tuple(std::move(DW).value(),Dlanda,DV));
    }


}
} // namespace Matrix_Decompositions

template<>
struct myDerivative<Matrix_Decompositions::eigensystem_type>
{
    typedef std::tuple<Derivative<M_Matrix<double>>, Derivative<M_Matrix<double>>, Derivative<M_Matrix<double>>> type;
};


template<typename T>
Derivative<M_Matrix<T>> diag(const Derivative<M_Matrix<T>>& x)
{
    return Derivative<M_Matrix<T>>(Matrix_Unary_Transformations::diag(x.f()),x.x(),
                                   x.dfdx().apply([](auto& df){return diag(df);})
                                   );
}





#endif // MATRIXDERIVATIVE_H
