#ifndef MAXIMUM_ENTROPY_H
#define MAXIMUM_ENTROPY_H
#include "Matrix.h"


struct maxent
{
    class step
    {
    private:
        std::size_t n_iter_;
        M_Matrix<double> pi_;
        M_Matrix<double> fi_;
        double S_;
        M_Matrix<double> landa_;
        M_Matrix<double> G_;
        M_Matrix<double> H_;
    public:
        std::size_t n_iter()const { return  n_iter_;}
        M_Matrix<double> const& pi()const { return  pi_;}
        M_Matrix<double> const& fi()const { return  fi_;}
        double S()const { return  S_;}
        M_Matrix<double> const& landa()const { return  landa_;}
        M_Matrix<double> const& G()const { return  G_;}
        M_Matrix<double> const& H()const { return  H_;}


        step()=default;

        step(std::size_t n_iter__,
             M_Matrix<double> const& pi__,
             M_Matrix<double> const&fi__,
             double S__,
             M_Matrix<double>const& landa__,
             M_Matrix<double>const& G__,
             M_Matrix<double> const& H__):
                                            n_iter_{n_iter__},pi_{pi__},fi_{fi__},S_{S__},landa_{landa__},G_{G__},H_{H__}
        {}
        step(std::size_t n_iter__,
             M_Matrix<double>&& pi__,
             M_Matrix<double>&& fi__,
             double S__,
             M_Matrix<double>&& landa__,
             M_Matrix<double>&& G__,
             M_Matrix<double>&& H__):
                                       n_iter_{n_iter__},pi_{std::move(pi__)},fi_{std::move(fi__)},S_{S__},landa_{std::move(landa__)},G_{std::move(G__)},H_{std::move(H__)}
        {}
        typedef step self_type;
        static auto get_constructor_fields()
        {
            return std::make_tuple(
                grammar::field(C<self_type>{},"n_iter",&self_type::n_iter),
                grammar::field(C<self_type>{},"P_i",&self_type::pi),
                grammar::field(C<self_type>{},"F_i",&self_type::fi),
                grammar::field(C<self_type>{},"S",&self_type::S),
                grammar::field(C<self_type>{},"landa",&self_type::landa),
                grammar::field(C<self_type>{},"G",&self_type::G),
                grammar::field(C<self_type>{},"H",&self_type::H));

        }
    };

static auto z_i(const M_Matrix<double>& f_ij, const M_Matrix<double>& landa)
{
    typedef myOptional_t<M_Matrix<double>> Op;
    assert(f_ij.ncols()==landa.ncols());
    assert(landa.nrows()==1);
    auto out= multTransp(f_ij, landa).apply([](auto&x){ return std::exp(-x);});
    if (isfinite(out))
        return Op(std::move(out));
    else {
        return Op(false,"some elements non finite: "+ToString(out));
    }

}


static auto calc_step(std::size_t nstep,const M_Matrix<double>& landa,const M_Matrix<double>& F, const M_Matrix<double>& f_ij)
{
    typedef myOptional_t<step> Op;
    assert(F.ncols()==f_ij.ncols());
    assert(F.nrows()==1);
    assert(F.ncols()==landa.ncols());
    assert(landa.nrows()==1);
    auto zi=z_i(f_ij,landa);
    if (!zi) return Op(false, "landa values to big: "+ToString(landa)+"z_i"+zi.error());
    double Zv=sum(zi.value());
    double S=std::log(Zv)+multTransp(landa,F).getvalue();
    auto pi=zi.value()/Zv;;
    auto fi=TranspMult(pi,f_ij);
    auto ffij=quadraticForm_BT_A_B(diag(pi),f_ij);
    return Op(step(nstep+1,std::move(pi), std::move(fi),S,landa,F-fi,ffij-quadraticForm_XTX(fi)));
}

class end_criteria
{
public:
    bool operator()(step& current,step& previous )const
    {
        if (maxAbs(current.landa())>max_landa_)
        {
            current=std::move(previous);
            return true;
        }
        else if (current.n_iter()>maxiter_) return true;
        else if (maxAbs(current.G())<maxG_) return true;
        else return maxAbs(current.G()-previous.G())<maxG_;
    }
    end_criteria(std::size_t maxiter, double maxGradient, double max_landa)
        : maxiter_{maxiter},maxG_{maxGradient}, max_landa_{max_landa}{}


private:
    std::size_t maxiter_;
    double maxG_;
    double max_landa_;
};


static auto next(const step& current,const M_Matrix<double>& F, const M_Matrix<double>& f_ij, double lm, double maxstep)
{
    typedef myOptional_t<step> Op;
  //  auto Hinv=inv(current.H()+diag(current.H()+eye<double>(current.H().ncols()))*lm);
    auto Hinv=inv(current.H()+diag(current.H())*lm);
    if (!Hinv.has_value())
    {
        return Op(false,"fails on inverse of Hessian "+Hinv.error());
    }
    else
    {
        auto step=current.G()*Hinv.value().first;
        step=step-mean(step);
        double m=maxAbs(step);
        if ( m>maxstep)
            step=step*maxstep/m;
        auto next_landa=current.landa()-step;
        return Op(calc_step(current.n_iter(),next_landa,F,f_ij));
    }
}

static
    std::vector<std::size_t> non_zero_index(const M_Matrix<double>& cov, double eps=std::sqrt(std::numeric_limits<double>::epsilon()))
{
    assert(cov.nrows()==cov.ncols());
    assert(cov.isSymmetric());

    std::vector<std::size_t> index;
    for (std::size_t i=0; i<cov.nrows(); ++i)
    {
        assert(cov(i,i)>=0);
        if (cov(i,i)>eps) index.push_back(i);
    }
    return index;
}

static M_Matrix<double> sub_matrix(const M_Matrix<double>& m, const std::vector<std::size_t>& index)
{
    assert(m.ncols()>=index.size());
    assert(m.ncols()==m.nrows());

        M_Matrix<double> out(index.size(),index.size(),m.type());
        for (std::size_t i=0; i<index.size(); ++i)
            for (std::size_t j=0; j<index.size(); ++j)
                out(i,j)=m(index[i],index[j]);
        return out;

}


static M_Matrix<double> sub_vector_column(const M_Matrix<double>& m, const std::vector<std::size_t>& index)
{
    assert(m.ncols()>=index.size());


        M_Matrix<double> out(m.nrows(),index.size());
        for (std::size_t i=0; i<out.nrows(); ++i)
        for (std::size_t j=0; j<index.size(); ++j)
            out(i,j)=m(i,index[j]);
        return out;

 }

 static step sub_step(const step& s, const std::vector<std::size_t>& index)
 {
     return step(s.n_iter(),
                 s.pi(),sub_vector_column(s.fi(),index),
                 s.S(),sub_vector_column(s.landa(),index),sub_vector_column(s.G(),index),sub_matrix(s.H(),index));
 }



 static M_Matrix<double> super_matrix(const M_Matrix<double>& Super,const M_Matrix<double>& sub, const std::vector<std::size_t>& index)
 {
     assert(sub.ncols()==index.size());
     auto out=Super;
     assert (sub.ncols()==sub.nrows());
         for (std::size_t i=0; i<index.size(); ++i)
             for (std::size_t j=0; j<index.size(); ++j)
                 out(index[i],index[j])=sub(i,j);
         return out;

 }

 static M_Matrix<double> super_vector_column(const M_Matrix<double>& Super,const M_Matrix<double>& sub, const std::vector<std::size_t>& index)
 {
     assert(sub.ncols()==index.size());
     auto out=Super;

         assert(sub.nrows()==Super.nrows());
         for (std::size_t i=0; i<out.nrows(); ++i)
             for (std::size_t j=0; j<index.size(); ++j)
                 out(i,index[j])=sub(i,j);
         return out;

 }


 static step super_step(const step& super, const step& sub,  const std::vector<std::size_t>& index)
 {
     return step(sub.n_iter(),
                 sub.pi(),super_vector_column(super.fi(),sub.fi(),index),
                 sub.S(),super_vector_column(super.landa(),sub.landa(),index),super_vector_column(super.G(),sub.G(),index),super_matrix(super.H(),sub.H(),index));
 }


 template<class fun>
 static auto optimize_it(const step& initial, const M_Matrix<double>& F, const M_Matrix<double>& f_ij,const fun& criteria , double maxstep)
 {
     double LM=1;
     typedef myOptional_t<step> Op;
     assert(F.ncols()==f_ij.ncols());
     assert(F.nrows()==1);
     auto previous=initial;
     auto res=next(previous,F,f_ij,LM,maxstep);
     if (!res) return Op(false,"failed because "+res.error()+" current step: "+ToString(previous));
     auto current=std::move(res).value();
     while(!criteria(current,previous))
     {
         auto res=next(current,F,f_ij,LM,maxstep);
         if (!res) return Op(false,"failed because "+res.error()+" current step: "+ToString(current));
         else {
                 previous=std::move(current);
                 current=std::move(res).value();

         }
     }
     return Op(current);


 }


template<class fun>
 static auto optimize(const M_Matrix<double>& F, const M_Matrix<double>& f_ij,const fun& criteria , double maxstep)
{
    typedef myOptional_t<step> Op;
    assert(F.ncols()==f_ij.ncols());
    assert(F.nrows()==1);
    auto landa=M_Matrix<double>(1,F.ncols(),0.0);
    auto current_o=calc_step(0,landa,F,f_ij);
    if (!current_o) return Op(false,"fails on start: "+current_o.error());
    auto current=std::move(current_o).value();
    auto index=non_zero_index(current.H());
    if (index.size()==current.G().size())
    {
        return optimize_it(current,F,f_ij,criteria,maxstep);
    }
    else
    {
        auto sub_current=sub_step(current,index);
        auto sub_F=sub_vector_column(F,index);
        auto sub_fij=sub_vector_column(f_ij,index);
        auto sub_res=optimize_it(sub_current,sub_F,sub_fij,criteria,maxstep);
        if (!sub_res) return Op(false,sub_res.error());
        else  return Op(super_step(current,sub_res.value(),index));
    }


}





};




#endif // MAXIMUM_ENTROPY_H
