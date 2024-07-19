// Mikael Jagan: from https://github.com/jaganmn/misc/blob/master/tmb_distributions/distributions.h
namespace glmmtmb
{

template<class Type>
Type mvlgamma(Type x, int p = 1)
{
    Type res = lgamma(x);
    if (p == 1)
    {
        return res;
    }
    for (int i = 1; i < p; ++i)
    {
        res += lgamma(x - Type(0.5 * i));
    }
    res += Type(0.25 * p * (p - 1) * log(M_PI));
    return res;
}
  
template<class Type>
Type dlkj(const vector<Type> &x, Type eta, int give_log = 0)
{
    int len = x.size();
    if (len == 0)
    {
        return ( give_log ? Type(0.0) : Type(1.0) );
    }
    int n = 0.5 * (1.0 + sqrt(1.0 + 8.0 * len));
    matrix<Type> L(n, n);
    L.setIdentity();
    for (int i = 0, k = 0; i < n; ++i)
    {
        for (int j = 0; j < i; ++j, ++k)
	{
	    L(i, j) = x(k);
	}
    }
    Type log_det_X = -(L.array() * L.array()).rowwise().sum().log().sum();
    Type log_res = (eta - Type(1.0)) * log_det_X;
    return ( give_log ? log_res : exp(log_res) );
}

template<class Type>
Type dwishart(const vector<Type> &x, Type df, const vector<Type> &scale, int give_log = 0)
{
    int len = x.size();
    int n = 0.5 * (-1.0 + sqrt(1.0 + 8.0 * len));

    matrix<Type> L_X(n, n);
    L_X.setIdentity();
    for (int i = 0, k = n; i < n; ++i)
    {
        for (int j = 0; j < i; ++j, ++k)
	{
	    L_X(i, j) = x(k);
	}
    }

    matrix<Type> L_S = L_X;
    for (int i = 0, k = n; i < n; ++i)
    {
        for (int j = 0; j < i; ++j, ++k)
	{
	    L_S(i, j) = scale(k);
	}
    }
    
    vector<Type> log_diag_LLT_X = (L_X.array() * L_X.array()).rowwise().sum().log();
    vector<Type> log_diag_LLT_S = (L_S.array() * L_S.array()).rowwise().sum().log();
    
    Type log_det_X = Type(2.0) * x.head(n).sum()     - log_diag_LLT_X.sum();
    Type log_det_S = Type(2.0) * scale.head(n).sum() - log_diag_LLT_S.sum();
    
    matrix<Type> invL_S = atomic::matinv(L_S);
    matrix<Type> A = (invL_S.transpose() * invL_S).array() * (L_X * L_X.transpose()).array();
    vector<Type> log_diag_D = x.head(n) - scale.head(n) - Type(0.5) * (log_diag_LLT_X - log_diag_LLT_S);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
	{
	    A(i, j) *= exp(log_diag_D(i) + log_diag_D(j));
	}
    }

    Type log_res = Type(-0.5) *
        (
	 df * log_det_S +
	 (-df + Type(n + 1)) * log_det_X +
	 df * Type(n * M_LN2) +
	 Type(2.0) * mvlgamma(Type(0.5) * df, n) +
	 A.sum() /* this term is 'tr(invS * X)' */
	);
    return ( give_log ? log_res : exp(log_res) ); 
}

template<class Type>
Type dinvwishart(const vector<Type> &x, Type df, const vector<Type> &scale, int give_log = 0)
{
    int len = x.size();
    int n = 0.5 * (-1.0 + sqrt(1.0 + 8.0 * len));

    matrix<Type> L_X(n, n);
    L_X.setIdentity();
    for (int i = 0, k = n; i < n; ++i)
    {
        for (int j = 0; j < i; ++j, ++k)
	{
	    L_X(i, j) = x(k);
	}
    }

    matrix<Type> L_S = L_X;
    for (int i = 0, k = n; i < n; ++i)
    {
        for (int j = 0; j < i; ++j, ++k)
	{
	    L_S(i, j) = scale(k);
	}
    }

    vector<Type> log_diag_LLT_X = (L_X.array() * L_X.array()).rowwise().sum().log();
    vector<Type> log_diag_LLT_S = (L_S.array() * L_S.array()).rowwise().sum().log();
    
    Type log_det_X = Type(2.0) * x.head(n).sum()     - log_diag_LLT_X.sum();
    Type log_det_S = Type(2.0) * scale.head(n).sum() - log_diag_LLT_S.sum();
    
    matrix<Type> invL_X = atomic::matinv(L_X);
    matrix<Type> A = (L_S * L_S.transpose()).array() * (invL_X.transpose() * invL_X).array();
    vector<Type> log_diag_D = scale.head(n) - x.head(n) - Type(0.5) * (log_diag_LLT_S - log_diag_LLT_X);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
	{
	    A(i, j) *= exp(log_diag_D(i) + log_diag_D(j));
	}
    }

    Type log_res = Type(-0.5) *
        (
	 -df * log_det_S +
	 (df + Type(n + 1)) * log_det_X +
	 df * Type(n * M_LN2) +
	 Type(2.0) * mvlgamma(Type(0.5) * df, n) +
	 A.sum() /* this term is 'tr(S * invX)' */
	);
    return ( give_log ? log_res : exp(log_res) ); 
}

} /* namespace distributions */

