// Generated by rstantools.  Do not edit by hand.

/*
    MEDdist is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MEDdist is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MEDdist.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_expectile_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_expectile");
    reader.add_event(53, 51, "end", "model_expectile");
    return reader;
}
template <typename T0__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__>::type, Eigen::Dynamic, Eigen::Dynamic>
square_root(const Eigen::Matrix<T0__, Eigen::Dynamic, Eigen::Dynamic>& A,
                const int& R, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 4;
        validate_non_negative_index("eigvecs", "R", R);
        validate_non_negative_index("eigvecs", "R", R);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> eigvecs(R, R);
        stan::math::initialize(eigvecs, DUMMY_VAR__);
        stan::math::fill(eigvecs, DUMMY_VAR__);
        stan::math::assign(eigvecs,eigenvectors_sym(A));
        current_statement_begin__ = 5;
        validate_non_negative_index("eigvals", "R", R);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> eigvals(R);
        stan::math::initialize(eigvals, DUMMY_VAR__);
        stan::math::fill(eigvals, DUMMY_VAR__);
        stan::math::assign(eigvals,eigenvalues_sym(A));
        current_statement_begin__ = 6;
        return stan::math::promote_scalar<fun_return_scalar_t__>(mdivide_right(multiply(eigvecs, diag_matrix(stan::math::sqrt(eigvals))), eigvecs));
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct square_root_functor__ {
    template <typename T0__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__>::type, Eigen::Dynamic, Eigen::Dynamic>
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, Eigen::Dynamic>& A,
                const int& R, std::ostream* pstream__) const {
        return square_root(A, R, pstream__);
    }
};
template <typename T1__>
typename boost::math::tools::promote_args<T1__>::type
log_norm_const(const int& d,
                   const T1__& rho, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T1__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 10;
        local_scalar_t__ num(DUMMY_VAR__);
        (void) num;  // dummy to suppress unused var warning
        stan::math::initialize(num, DUMMY_VAR__);
        stan::math::fill(num, DUMMY_VAR__);
        stan::math::assign(num,(((d - 1) * stan::math::log(2)) + stan::math::log((1 + inv_sqrt((1 - pow(rho, 2)))))));
        current_statement_begin__ = 11;
        local_scalar_t__ den(DUMMY_VAR__);
        (void) den;  // dummy to suppress unused var warning
        stan::math::initialize(den, DUMMY_VAR__);
        stan::math::fill(den, DUMMY_VAR__);
        stan::math::assign(den,(d * stan::math::log((stan::math::sqrt((1 - rho)) + stan::math::sqrt((1 + rho))))));
        current_statement_begin__ = 12;
        return stan::math::promote_scalar<fun_return_scalar_t__>((num - den));
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct log_norm_const_functor__ {
    template <typename T1__>
        typename boost::math::tools::promote_args<T1__>::type
    operator()(const int& d,
                   const T1__& rho, std::ostream* pstream__) const {
        return log_norm_const(d, rho, pstream__);
    }
};
template <bool propto, typename T0__, typename T1__, typename T2__, typename T3__, typename T4__, typename T5__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__>::type>::type
mymed_lpdf(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& x,
               const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& m,
               const Eigen::Matrix<T2__, Eigen::Dynamic, Eigen::Dynamic>& sigma_inv,
               const Eigen::Matrix<T3__, Eigen::Dynamic, Eigen::Dynamic>& sigma_inv_sqrt,
               const Eigen::Matrix<T4__, Eigen::Dynamic, 1>& nu,
               const T5__& r,
               const int& d, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__>::type>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 16;
        local_scalar_t__ out(DUMMY_VAR__);
        (void) out;  // dummy to suppress unused var warning
        stan::math::initialize(out, DUMMY_VAR__);
        stan::math::fill(out, DUMMY_VAR__);
        stan::math::assign(out,multi_normal_prec_log(x, m, sigma_inv));
        current_statement_begin__ = 17;
        stan::math::assign(out, (out - log_norm_const(d, r, pstream__)));
        current_statement_begin__ = 18;
        stan::math::assign(out, (out - (multiply(multiply(multiply((r / 2), transpose(subtract(x, m))), sigma_inv_sqrt), nu) * stan::math::sqrt(quad_form(sigma_inv, subtract(x, m))))));
        current_statement_begin__ = 19;
        return stan::math::promote_scalar<fun_return_scalar_t__>(out);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
template <typename T0__, typename T1__, typename T2__, typename T3__, typename T4__, typename T5__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__>::type>::type
mymed_lpdf(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& x,
               const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& m,
               const Eigen::Matrix<T2__, Eigen::Dynamic, Eigen::Dynamic>& sigma_inv,
               const Eigen::Matrix<T3__, Eigen::Dynamic, Eigen::Dynamic>& sigma_inv_sqrt,
               const Eigen::Matrix<T4__, Eigen::Dynamic, 1>& nu,
               const T5__& r,
               const int& d, std::ostream* pstream__) {
    return mymed_lpdf<false>(x,m,sigma_inv,sigma_inv_sqrt,nu,r,d, pstream__);
}
struct mymed_lpdf_functor__ {
    template <bool propto, typename T0__, typename T1__, typename T2__, typename T3__, typename T4__, typename T5__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__>::type>::type
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& x,
               const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& m,
               const Eigen::Matrix<T2__, Eigen::Dynamic, Eigen::Dynamic>& sigma_inv,
               const Eigen::Matrix<T3__, Eigen::Dynamic, Eigen::Dynamic>& sigma_inv_sqrt,
               const Eigen::Matrix<T4__, Eigen::Dynamic, 1>& nu,
               const T5__& r,
               const int& d, std::ostream* pstream__) const {
        return mymed_lpdf(x, m, sigma_inv, sigma_inv_sqrt, nu, r, d, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_expectile
  : public stan::model::model_base_crtp<model_expectile> {
private:
        int d;
        int n;
        double rho;
        vector_d nu_tilde;
        std::vector<vector_d> Y;
public:
    model_expectile(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_expectile(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_expectile_namespace::model_expectile";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 24;
            context__.validate_dims("data initialization", "d", "int", context__.to_vec());
            d = int(0);
            vals_i__ = context__.vals_i("d");
            pos__ = 0;
            d = vals_i__[pos__++];
            check_greater_or_equal(function__, "d", d, 1);
            current_statement_begin__ = 25;
            context__.validate_dims("data initialization", "n", "int", context__.to_vec());
            n = int(0);
            vals_i__ = context__.vals_i("n");
            pos__ = 0;
            n = vals_i__[pos__++];
            check_greater_or_equal(function__, "n", n, 0);
            current_statement_begin__ = 26;
            context__.validate_dims("data initialization", "rho", "double", context__.to_vec());
            rho = double(0);
            vals_r__ = context__.vals_r("rho");
            pos__ = 0;
            rho = vals_r__[pos__++];
            check_greater_or_equal(function__, "rho", rho, 0);
            check_less_or_equal(function__, "rho", rho, 1);
            current_statement_begin__ = 27;
            validate_non_negative_index("nu_tilde", "d", d);
            context__.validate_dims("data initialization", "nu_tilde", "vector_d", context__.to_vec(d));
            nu_tilde = Eigen::Matrix<double, Eigen::Dynamic, 1>(d);
            vals_r__ = context__.vals_r("nu_tilde");
            pos__ = 0;
            size_t nu_tilde_j_1_max__ = d;
            for (size_t j_1__ = 0; j_1__ < nu_tilde_j_1_max__; ++j_1__) {
                nu_tilde(j_1__) = vals_r__[pos__++];
            }
            stan::math::check_unit_vector(function__, "nu_tilde", nu_tilde);
            current_statement_begin__ = 28;
            validate_non_negative_index("Y", "d", d);
            validate_non_negative_index("Y", "n", n);
            context__.validate_dims("data initialization", "Y", "vector_d", context__.to_vec(n,d));
            Y = std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >(n, Eigen::Matrix<double, Eigen::Dynamic, 1>(d));
            vals_r__ = context__.vals_r("Y");
            pos__ = 0;
            size_t Y_j_1_max__ = d;
            size_t Y_k_0_max__ = n;
            for (size_t j_1__ = 0; j_1__ < Y_j_1_max__; ++j_1__) {
                for (size_t k_0__ = 0; k_0__ < Y_k_0_max__; ++k_0__) {
                    Y[k_0__](j_1__) = vals_r__[pos__++];
                }
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 32;
            validate_non_negative_index("mu", "d", d);
            num_params_r__ += d;
            current_statement_begin__ = 33;
            validate_non_negative_index("Sigma_inv", "d", d);
            validate_non_negative_index("Sigma_inv", "d", d);
            num_params_r__ += (d + ((d * (d - 1)) / 2));
            current_statement_begin__ = 34;
            validate_non_negative_index("diag_sigma_sq", "d", d);
            num_params_r__ += d;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_expectile() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 32;
        if (!(context__.contains_r("mu")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable mu missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("mu");
        pos__ = 0U;
        validate_non_negative_index("mu", "d", d);
        context__.validate_dims("parameter initialization", "mu", "vector_d", context__.to_vec(d));
        Eigen::Matrix<double, Eigen::Dynamic, 1> mu(d);
        size_t mu_j_1_max__ = d;
        for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
            mu(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(mu);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable mu: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 33;
        if (!(context__.contains_r("Sigma_inv")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable Sigma_inv missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("Sigma_inv");
        pos__ = 0U;
        validate_non_negative_index("Sigma_inv", "d", d);
        validate_non_negative_index("Sigma_inv", "d", d);
        context__.validate_dims("parameter initialization", "Sigma_inv", "matrix_d", context__.to_vec(d,d));
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Sigma_inv(d, d);
        size_t Sigma_inv_j_2_max__ = d;
        size_t Sigma_inv_j_1_max__ = d;
        for (size_t j_2__ = 0; j_2__ < Sigma_inv_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < Sigma_inv_j_1_max__; ++j_1__) {
                Sigma_inv(j_1__, j_2__) = vals_r__[pos__++];
            }
        }
        try {
            writer__.cov_matrix_unconstrain(Sigma_inv);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable Sigma_inv: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 34;
        if (!(context__.contains_r("diag_sigma_sq")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable diag_sigma_sq missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("diag_sigma_sq");
        pos__ = 0U;
        validate_non_negative_index("diag_sigma_sq", "d", d);
        context__.validate_dims("parameter initialization", "diag_sigma_sq", "vector_d", context__.to_vec(d));
        Eigen::Matrix<double, Eigen::Dynamic, 1> diag_sigma_sq(d);
        size_t diag_sigma_sq_j_1_max__ = d;
        for (size_t j_1__ = 0; j_1__ < diag_sigma_sq_j_1_max__; ++j_1__) {
            diag_sigma_sq(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_lb_unconstrain(0, diag_sigma_sq);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable diag_sigma_sq: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 32;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> mu;
            (void) mu;  // dummy to suppress unused var warning
            if (jacobian__)
                mu = in__.vector_constrain(d, lp__);
            else
                mu = in__.vector_constrain(d);
            current_statement_begin__ = 33;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> Sigma_inv;
            (void) Sigma_inv;  // dummy to suppress unused var warning
            if (jacobian__)
                Sigma_inv = in__.cov_matrix_constrain(d, lp__);
            else
                Sigma_inv = in__.cov_matrix_constrain(d);
            current_statement_begin__ = 34;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> diag_sigma_sq;
            (void) diag_sigma_sq;  // dummy to suppress unused var warning
            if (jacobian__)
                diag_sigma_sq = in__.vector_lb_constrain(0, d, lp__);
            else
                diag_sigma_sq = in__.vector_lb_constrain(0, d);
            // transformed parameters
            current_statement_begin__ = 38;
            validate_non_negative_index("sigma_inv_sqrt", "d", d);
            validate_non_negative_index("sigma_inv_sqrt", "d", d);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> sigma_inv_sqrt(d, d);
            stan::math::initialize(sigma_inv_sqrt, DUMMY_VAR__);
            stan::math::fill(sigma_inv_sqrt, DUMMY_VAR__);
            stan::math::assign(sigma_inv_sqrt,square_root(Sigma_inv, d, pstream__));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 38;
            size_t sigma_inv_sqrt_j_1_max__ = d;
            size_t sigma_inv_sqrt_j_2_max__ = d;
            for (size_t j_1__ = 0; j_1__ < sigma_inv_sqrt_j_1_max__; ++j_1__) {
                for (size_t j_2__ = 0; j_2__ < sigma_inv_sqrt_j_2_max__; ++j_2__) {
                    if (stan::math::is_uninitialized(sigma_inv_sqrt(j_1__, j_2__))) {
                        std::stringstream msg__;
                        msg__ << "Undefined transformed parameter: sigma_inv_sqrt" << "(" << j_1__ << ", " << j_2__ << ")";
                        stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable sigma_inv_sqrt: ") + msg__.str()), current_statement_begin__, prog_reader__());
                    }
                }
            }
            // model body
            current_statement_begin__ = 43;
            lp_accum__.add(normal_log<propto__>(mu, rep_vector(0, d), diag_sigma_sq));
            current_statement_begin__ = 44;
            for (int j = 1; j <= d; ++j) {
                current_statement_begin__ = 45;
                lp_accum__.add(inv_gamma_log<propto__>(get_base1(diag_sigma_sq, j, "diag_sigma_sq", 1), 1, 1));
            }
            current_statement_begin__ = 47;
            lp_accum__.add(wishart_log<propto__>(Sigma_inv, d, diag_matrix(rep_vector(1, d))));
            current_statement_begin__ = 48;
            for (int i = 1; i <= n; ++i) {
                current_statement_begin__ = 49;
                lp_accum__.add(mymed_lpdf<propto__>(get_base1(Y, i, "Y", 1), mu, Sigma_inv, sigma_inv_sqrt, nu_tilde, rho, d, pstream__));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("mu");
        names__.push_back("Sigma_inv");
        names__.push_back("diag_sigma_sq");
        names__.push_back("sigma_inv_sqrt");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(d);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(d);
        dims__.push_back(d);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(d);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(d);
        dims__.push_back(d);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_expectile_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> mu = in__.vector_constrain(d);
        size_t mu_j_1_max__ = d;
        for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
            vars__.push_back(mu(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Sigma_inv = in__.cov_matrix_constrain(d);
        size_t Sigma_inv_j_2_max__ = d;
        size_t Sigma_inv_j_1_max__ = d;
        for (size_t j_2__ = 0; j_2__ < Sigma_inv_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < Sigma_inv_j_1_max__; ++j_1__) {
                vars__.push_back(Sigma_inv(j_1__, j_2__));
            }
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> diag_sigma_sq = in__.vector_lb_constrain(0, d);
        size_t diag_sigma_sq_j_1_max__ = d;
        for (size_t j_1__ = 0; j_1__ < diag_sigma_sq_j_1_max__; ++j_1__) {
            vars__.push_back(diag_sigma_sq(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 38;
            validate_non_negative_index("sigma_inv_sqrt", "d", d);
            validate_non_negative_index("sigma_inv_sqrt", "d", d);
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> sigma_inv_sqrt(d, d);
            stan::math::initialize(sigma_inv_sqrt, DUMMY_VAR__);
            stan::math::fill(sigma_inv_sqrt, DUMMY_VAR__);
            stan::math::assign(sigma_inv_sqrt,square_root(Sigma_inv, d, pstream__));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t sigma_inv_sqrt_j_2_max__ = d;
                size_t sigma_inv_sqrt_j_1_max__ = d;
                for (size_t j_2__ = 0; j_2__ < sigma_inv_sqrt_j_2_max__; ++j_2__) {
                    for (size_t j_1__ = 0; j_1__ < sigma_inv_sqrt_j_1_max__; ++j_1__) {
                        vars__.push_back(sigma_inv_sqrt(j_1__, j_2__));
                    }
                }
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_expectile";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t mu_j_1_max__ = d;
        for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t Sigma_inv_j_2_max__ = d;
        size_t Sigma_inv_j_1_max__ = d;
        for (size_t j_2__ = 0; j_2__ < Sigma_inv_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < Sigma_inv_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "Sigma_inv" << '.' << j_1__ + 1 << '.' << j_2__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        size_t diag_sigma_sq_j_1_max__ = d;
        for (size_t j_1__ = 0; j_1__ < diag_sigma_sq_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "diag_sigma_sq" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t sigma_inv_sqrt_j_2_max__ = d;
            size_t sigma_inv_sqrt_j_1_max__ = d;
            for (size_t j_2__ = 0; j_2__ < sigma_inv_sqrt_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < sigma_inv_sqrt_j_1_max__; ++j_1__) {
                    param_name_stream__.str(std::string());
                    param_name_stream__ << "sigma_inv_sqrt" << '.' << j_1__ + 1 << '.' << j_2__ + 1;
                    param_names__.push_back(param_name_stream__.str());
                }
            }
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t mu_j_1_max__ = d;
        for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t Sigma_inv_j_1_max__ = (d + ((d * (d - 1)) / 2));
        for (size_t j_1__ = 0; j_1__ < Sigma_inv_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "Sigma_inv" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t diag_sigma_sq_j_1_max__ = d;
        for (size_t j_1__ = 0; j_1__ < diag_sigma_sq_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "diag_sigma_sq" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t sigma_inv_sqrt_j_2_max__ = d;
            size_t sigma_inv_sqrt_j_1_max__ = d;
            for (size_t j_2__ = 0; j_2__ < sigma_inv_sqrt_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < sigma_inv_sqrt_j_1_max__; ++j_1__) {
                    param_name_stream__.str(std::string());
                    param_name_stream__ << "sigma_inv_sqrt" << '.' << j_1__ + 1 << '.' << j_2__ + 1;
                    param_names__.push_back(param_name_stream__.str());
                }
            }
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_expectile_namespace::model_expectile stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
