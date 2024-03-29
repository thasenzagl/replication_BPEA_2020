// Code generated by Stan version 2.21.0

#include <stan/model/model_header.hpp>

namespace condhet_horseshoe_model_namespace {

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
    reader.add_event(0, 0, "start", "D:/Dropbox/research/Projects/SP20_BPEA/Code/estim/condhet_horseshoe.stan");
    reader.add_event(43, 41, "end", "D:/Dropbox/research/Projects/SP20_BPEA/Code/estim/condhet_horseshoe.stan");
    return reader;
}

class condhet_horseshoe_model
  : public stan::model::model_base_crtp<condhet_horseshoe_model> {
private:
        int T;
        int q;
        int p;
        vector_d Y;
        matrix_d W;
        matrix_d X;
        double min_scale;
public:
    condhet_horseshoe_model(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }

    condhet_horseshoe_model(stan::io::var_context& context__,
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

        static const char* function__ = "condhet_horseshoe_model_namespace::condhet_horseshoe_model";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        try {
            // initialize data block variables from context__
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "T", "int", context__.to_vec());
            T = int(0);
            vals_i__ = context__.vals_i("T");
            pos__ = 0;
            T = vals_i__[pos__++];
            check_greater_or_equal(function__, "T", T, 1);

            current_statement_begin__ = 3;
            context__.validate_dims("data initialization", "q", "int", context__.to_vec());
            q = int(0);
            vals_i__ = context__.vals_i("q");
            pos__ = 0;
            q = vals_i__[pos__++];
            check_greater_or_equal(function__, "q", q, 0);

            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "p", "int", context__.to_vec());
            p = int(0);
            vals_i__ = context__.vals_i("p");
            pos__ = 0;
            p = vals_i__[pos__++];
            check_greater_or_equal(function__, "p", p, 0);

            current_statement_begin__ = 5;
            validate_non_negative_index("Y", "T", T);
            context__.validate_dims("data initialization", "Y", "vector_d", context__.to_vec(T));
            Y = Eigen::Matrix<double, Eigen::Dynamic, 1>(T);
            vals_r__ = context__.vals_r("Y");
            pos__ = 0;
            size_t Y_j_1_max__ = T;
            for (size_t j_1__ = 0; j_1__ < Y_j_1_max__; ++j_1__) {
                Y(j_1__) = vals_r__[pos__++];
            }

            current_statement_begin__ = 6;
            validate_non_negative_index("W", "T", T);
            validate_non_negative_index("W", "q", q);
            context__.validate_dims("data initialization", "W", "matrix_d", context__.to_vec(T,q));
            W = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(T, q);
            vals_r__ = context__.vals_r("W");
            pos__ = 0;
            size_t W_j_2_max__ = q;
            size_t W_j_1_max__ = T;
            for (size_t j_2__ = 0; j_2__ < W_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < W_j_1_max__; ++j_1__) {
                    W(j_1__, j_2__) = vals_r__[pos__++];
                }
            }

            current_statement_begin__ = 7;
            validate_non_negative_index("X", "T", T);
            validate_non_negative_index("X", "p", p);
            context__.validate_dims("data initialization", "X", "matrix_d", context__.to_vec(T,p));
            X = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(T, p);
            vals_r__ = context__.vals_r("X");
            pos__ = 0;
            size_t X_j_2_max__ = p;
            size_t X_j_1_max__ = T;
            for (size_t j_2__ = 0; j_2__ < X_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_j_1_max__; ++j_1__) {
                    X(j_1__, j_2__) = vals_r__[pos__++];
                }
            }

            current_statement_begin__ = 8;
            context__.validate_dims("data initialization", "min_scale", "double", context__.to_vec());
            min_scale = double(0);
            vals_r__ = context__.vals_r("min_scale");
            pos__ = 0;
            min_scale = vals_r__[pos__++];
            check_greater_or_equal(function__, "min_scale", min_scale, 0);


            // initialize transformed data variables
            // execute transformed data statements

            // validate transformed data

            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 12;
            validate_non_negative_index("gamma_mu", "q", q);
            num_params_r__ += q;
            current_statement_begin__ = 13;
            validate_non_negative_index("gamma_sigma", "q", q);
            num_params_r__ += q;
            current_statement_begin__ = 14;
            validate_non_negative_index("beta_mu_std", "p", p);
            num_params_r__ += p;
            current_statement_begin__ = 15;
            validate_non_negative_index("beta_sigma_std", "p", p);
            num_params_r__ += p;
            current_statement_begin__ = 16;
            validate_non_negative_index("lambda_mu", "p", p);
            num_params_r__ += p;
            current_statement_begin__ = 17;
            validate_non_negative_index("lambda_sigma", "p", p);
            num_params_r__ += p;
            current_statement_begin__ = 18;
            num_params_r__ += 1;
            current_statement_begin__ = 19;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    ~condhet_horseshoe_model() { }


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

        current_statement_begin__ = 12;
        if (!(context__.contains_r("gamma_mu")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable gamma_mu missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("gamma_mu");
        pos__ = 0U;
        validate_non_negative_index("gamma_mu", "q", q);
        context__.validate_dims("parameter initialization", "gamma_mu", "vector_d", context__.to_vec(q));
        Eigen::Matrix<double, Eigen::Dynamic, 1> gamma_mu(q);
        size_t gamma_mu_j_1_max__ = q;
        for (size_t j_1__ = 0; j_1__ < gamma_mu_j_1_max__; ++j_1__) {
            gamma_mu(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(gamma_mu);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable gamma_mu: ") + e.what()), current_statement_begin__, prog_reader__());
        }

        current_statement_begin__ = 13;
        if (!(context__.contains_r("gamma_sigma")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable gamma_sigma missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("gamma_sigma");
        pos__ = 0U;
        validate_non_negative_index("gamma_sigma", "q", q);
        context__.validate_dims("parameter initialization", "gamma_sigma", "vector_d", context__.to_vec(q));
        Eigen::Matrix<double, Eigen::Dynamic, 1> gamma_sigma(q);
        size_t gamma_sigma_j_1_max__ = q;
        for (size_t j_1__ = 0; j_1__ < gamma_sigma_j_1_max__; ++j_1__) {
            gamma_sigma(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(gamma_sigma);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable gamma_sigma: ") + e.what()), current_statement_begin__, prog_reader__());
        }

        current_statement_begin__ = 14;
        if (!(context__.contains_r("beta_mu_std")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta_mu_std missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta_mu_std");
        pos__ = 0U;
        validate_non_negative_index("beta_mu_std", "p", p);
        context__.validate_dims("parameter initialization", "beta_mu_std", "vector_d", context__.to_vec(p));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_mu_std(p);
        size_t beta_mu_std_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_mu_std_j_1_max__; ++j_1__) {
            beta_mu_std(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta_mu_std);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta_mu_std: ") + e.what()), current_statement_begin__, prog_reader__());
        }

        current_statement_begin__ = 15;
        if (!(context__.contains_r("beta_sigma_std")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta_sigma_std missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta_sigma_std");
        pos__ = 0U;
        validate_non_negative_index("beta_sigma_std", "p", p);
        context__.validate_dims("parameter initialization", "beta_sigma_std", "vector_d", context__.to_vec(p));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_sigma_std(p);
        size_t beta_sigma_std_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_sigma_std_j_1_max__; ++j_1__) {
            beta_sigma_std(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta_sigma_std);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta_sigma_std: ") + e.what()), current_statement_begin__, prog_reader__());
        }

        current_statement_begin__ = 16;
        if (!(context__.contains_r("lambda_mu")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable lambda_mu missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("lambda_mu");
        pos__ = 0U;
        validate_non_negative_index("lambda_mu", "p", p);
        context__.validate_dims("parameter initialization", "lambda_mu", "vector_d", context__.to_vec(p));
        Eigen::Matrix<double, Eigen::Dynamic, 1> lambda_mu(p);
        size_t lambda_mu_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < lambda_mu_j_1_max__; ++j_1__) {
            lambda_mu(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(lambda_mu);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable lambda_mu: ") + e.what()), current_statement_begin__, prog_reader__());
        }

        current_statement_begin__ = 17;
        if (!(context__.contains_r("lambda_sigma")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable lambda_sigma missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("lambda_sigma");
        pos__ = 0U;
        validate_non_negative_index("lambda_sigma", "p", p);
        context__.validate_dims("parameter initialization", "lambda_sigma", "vector_d", context__.to_vec(p));
        Eigen::Matrix<double, Eigen::Dynamic, 1> lambda_sigma(p);
        size_t lambda_sigma_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < lambda_sigma_j_1_max__; ++j_1__) {
            lambda_sigma(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(lambda_sigma);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable lambda_sigma: ") + e.what()), current_statement_begin__, prog_reader__());
        }

        current_statement_begin__ = 18;
        if (!(context__.contains_r("tau_mu")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable tau_mu missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("tau_mu");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "tau_mu", "double", context__.to_vec());
        double tau_mu(0);
        tau_mu = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(min_scale, tau_mu);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable tau_mu: ") + e.what()), current_statement_begin__, prog_reader__());
        }

        current_statement_begin__ = 19;
        if (!(context__.contains_r("tau_sigma")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable tau_sigma missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("tau_sigma");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "tau_sigma", "double", context__.to_vec());
        double tau_sigma(0);
        tau_sigma = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(min_scale, tau_sigma);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable tau_sigma: ") + e.what()), current_statement_begin__, prog_reader__());
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
            current_statement_begin__ = 12;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> gamma_mu;
            (void) gamma_mu;  // dummy to suppress unused var warning
            if (jacobian__)
                gamma_mu = in__.vector_constrain(q, lp__);
            else
                gamma_mu = in__.vector_constrain(q);

            current_statement_begin__ = 13;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> gamma_sigma;
            (void) gamma_sigma;  // dummy to suppress unused var warning
            if (jacobian__)
                gamma_sigma = in__.vector_constrain(q, lp__);
            else
                gamma_sigma = in__.vector_constrain(q);

            current_statement_begin__ = 14;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta_mu_std;
            (void) beta_mu_std;  // dummy to suppress unused var warning
            if (jacobian__)
                beta_mu_std = in__.vector_constrain(p, lp__);
            else
                beta_mu_std = in__.vector_constrain(p);

            current_statement_begin__ = 15;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta_sigma_std;
            (void) beta_sigma_std;  // dummy to suppress unused var warning
            if (jacobian__)
                beta_sigma_std = in__.vector_constrain(p, lp__);
            else
                beta_sigma_std = in__.vector_constrain(p);

            current_statement_begin__ = 16;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> lambda_mu;
            (void) lambda_mu;  // dummy to suppress unused var warning
            if (jacobian__)
                lambda_mu = in__.vector_constrain(p, lp__);
            else
                lambda_mu = in__.vector_constrain(p);

            current_statement_begin__ = 17;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> lambda_sigma;
            (void) lambda_sigma;  // dummy to suppress unused var warning
            if (jacobian__)
                lambda_sigma = in__.vector_constrain(p, lp__);
            else
                lambda_sigma = in__.vector_constrain(p);

            current_statement_begin__ = 18;
            local_scalar_t__ tau_mu;
            (void) tau_mu;  // dummy to suppress unused var warning
            if (jacobian__)
                tau_mu = in__.scalar_lb_constrain(min_scale, lp__);
            else
                tau_mu = in__.scalar_lb_constrain(min_scale);

            current_statement_begin__ = 19;
            local_scalar_t__ tau_sigma;
            (void) tau_sigma;  // dummy to suppress unused var warning
            if (jacobian__)
                tau_sigma = in__.scalar_lb_constrain(min_scale, lp__);
            else
                tau_sigma = in__.scalar_lb_constrain(min_scale);

            // transformed parameters
            current_statement_begin__ = 23;
            validate_non_negative_index("beta_mu", "p", p);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta_mu(p);
            stan::math::initialize(beta_mu, DUMMY_VAR__);
            stan::math::fill(beta_mu, DUMMY_VAR__);
            stan::math::assign(beta_mu,multiply(tau_mu, elt_multiply(lambda_mu, beta_mu_std)));

            current_statement_begin__ = 24;
            validate_non_negative_index("beta_sigma", "p", p);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta_sigma(p);
            stan::math::initialize(beta_sigma, DUMMY_VAR__);
            stan::math::fill(beta_sigma, DUMMY_VAR__);
            stan::math::assign(beta_sigma,multiply(tau_sigma, elt_multiply(lambda_sigma, beta_sigma_std)));

            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning

            current_statement_begin__ = 23;
            size_t beta_mu_j_1_max__ = p;
            for (size_t j_1__ = 0; j_1__ < beta_mu_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(beta_mu(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: beta_mu" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable beta_mu: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            current_statement_begin__ = 24;
            size_t beta_sigma_j_1_max__ = p;
            for (size_t j_1__ = 0; j_1__ < beta_sigma_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(beta_sigma(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: beta_sigma" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable beta_sigma: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }

            // model body

            current_statement_begin__ = 29;
            lp_accum__.add(cauchy_log<propto__>(gamma_mu, 0, 5));
            current_statement_begin__ = 30;
            lp_accum__.add(cauchy_log<propto__>(gamma_sigma, 0, 5));
            current_statement_begin__ = 31;
            lp_accum__.add(normal_log<propto__>(beta_mu_std, 0, 1));
            current_statement_begin__ = 32;
            lp_accum__.add(normal_log<propto__>(beta_sigma_std, 0, 1));
            current_statement_begin__ = 33;
            lp_accum__.add(cauchy_log<propto__>(lambda_mu, 0, 1));
            current_statement_begin__ = 34;
            lp_accum__.add(cauchy_log<propto__>(lambda_sigma, 0, 1));
            current_statement_begin__ = 35;
            lp_accum__.add(cauchy_log<propto__>(tau_mu, 0, 1));
            current_statement_begin__ = 36;
            lp_accum__.add(cauchy_log<propto__>(tau_sigma, 0, 1));
            current_statement_begin__ = 39;
            lp_accum__.add(normal_log<propto__>(Y, add(multiply(W, gamma_mu), multiply(X, beta_mu)), stan::math::exp(add(multiply(W, gamma_sigma), multiply(X, beta_sigma)))));

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
        names__.push_back("gamma_mu");
        names__.push_back("gamma_sigma");
        names__.push_back("beta_mu_std");
        names__.push_back("beta_sigma_std");
        names__.push_back("lambda_mu");
        names__.push_back("lambda_sigma");
        names__.push_back("tau_mu");
        names__.push_back("tau_sigma");
        names__.push_back("beta_mu");
        names__.push_back("beta_sigma");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(q);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(q);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(p);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(p);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(p);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(p);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(p);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(p);
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
        static const char* function__ = "condhet_horseshoe_model_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning

        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> gamma_mu = in__.vector_constrain(q);
        size_t gamma_mu_j_1_max__ = q;
        for (size_t j_1__ = 0; j_1__ < gamma_mu_j_1_max__; ++j_1__) {
            vars__.push_back(gamma_mu(j_1__));
        }

        Eigen::Matrix<double, Eigen::Dynamic, 1> gamma_sigma = in__.vector_constrain(q);
        size_t gamma_sigma_j_1_max__ = q;
        for (size_t j_1__ = 0; j_1__ < gamma_sigma_j_1_max__; ++j_1__) {
            vars__.push_back(gamma_sigma(j_1__));
        }

        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_mu_std = in__.vector_constrain(p);
        size_t beta_mu_std_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_mu_std_j_1_max__; ++j_1__) {
            vars__.push_back(beta_mu_std(j_1__));
        }

        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_sigma_std = in__.vector_constrain(p);
        size_t beta_sigma_std_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_sigma_std_j_1_max__; ++j_1__) {
            vars__.push_back(beta_sigma_std(j_1__));
        }

        Eigen::Matrix<double, Eigen::Dynamic, 1> lambda_mu = in__.vector_constrain(p);
        size_t lambda_mu_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < lambda_mu_j_1_max__; ++j_1__) {
            vars__.push_back(lambda_mu(j_1__));
        }

        Eigen::Matrix<double, Eigen::Dynamic, 1> lambda_sigma = in__.vector_constrain(p);
        size_t lambda_sigma_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < lambda_sigma_j_1_max__; ++j_1__) {
            vars__.push_back(lambda_sigma(j_1__));
        }

        double tau_mu = in__.scalar_lb_constrain(min_scale);
        vars__.push_back(tau_mu);

        double tau_sigma = in__.scalar_lb_constrain(min_scale);
        vars__.push_back(tau_sigma);

        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        if (!include_tparams__ && !include_gqs__) return;

        try {
            // declare and define transformed parameters
            current_statement_begin__ = 23;
            validate_non_negative_index("beta_mu", "p", p);
            Eigen::Matrix<double, Eigen::Dynamic, 1> beta_mu(p);
            stan::math::initialize(beta_mu, DUMMY_VAR__);
            stan::math::fill(beta_mu, DUMMY_VAR__);
            stan::math::assign(beta_mu,multiply(tau_mu, elt_multiply(lambda_mu, beta_mu_std)));

            current_statement_begin__ = 24;
            validate_non_negative_index("beta_sigma", "p", p);
            Eigen::Matrix<double, Eigen::Dynamic, 1> beta_sigma(p);
            stan::math::initialize(beta_sigma, DUMMY_VAR__);
            stan::math::fill(beta_sigma, DUMMY_VAR__);
            stan::math::assign(beta_sigma,multiply(tau_sigma, elt_multiply(lambda_sigma, beta_sigma_std)));

            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning

            // write transformed parameters
            if (include_tparams__) {
                size_t beta_mu_j_1_max__ = p;
                for (size_t j_1__ = 0; j_1__ < beta_mu_j_1_max__; ++j_1__) {
                    vars__.push_back(beta_mu(j_1__));
                }
                size_t beta_sigma_j_1_max__ = p;
                for (size_t j_1__ = 0; j_1__ < beta_sigma_j_1_max__; ++j_1__) {
                    vars__.push_back(beta_sigma(j_1__));
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
        return "condhet_horseshoe_model";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t gamma_mu_j_1_max__ = q;
        for (size_t j_1__ = 0; j_1__ < gamma_mu_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "gamma_mu" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t gamma_sigma_j_1_max__ = q;
        for (size_t j_1__ = 0; j_1__ < gamma_sigma_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "gamma_sigma" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t beta_mu_std_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_mu_std_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_mu_std" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t beta_sigma_std_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_sigma_std_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_sigma_std" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t lambda_mu_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < lambda_mu_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lambda_mu" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t lambda_sigma_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < lambda_sigma_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lambda_sigma" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau_mu";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau_sigma";
        param_names__.push_back(param_name_stream__.str());

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
            size_t beta_mu_j_1_max__ = p;
            for (size_t j_1__ = 0; j_1__ < beta_mu_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "beta_mu" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t beta_sigma_j_1_max__ = p;
            for (size_t j_1__ = 0; j_1__ < beta_sigma_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "beta_sigma" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }

        if (!include_gqs__) return;
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t gamma_mu_j_1_max__ = q;
        for (size_t j_1__ = 0; j_1__ < gamma_mu_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "gamma_mu" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t gamma_sigma_j_1_max__ = q;
        for (size_t j_1__ = 0; j_1__ < gamma_sigma_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "gamma_sigma" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t beta_mu_std_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_mu_std_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_mu_std" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t beta_sigma_std_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_sigma_std_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_sigma_std" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t lambda_mu_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < lambda_mu_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lambda_mu" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t lambda_sigma_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < lambda_sigma_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lambda_sigma" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau_mu";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau_sigma";
        param_names__.push_back(param_name_stream__.str());

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
            size_t beta_mu_j_1_max__ = p;
            for (size_t j_1__ = 0; j_1__ < beta_mu_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "beta_mu" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t beta_sigma_j_1_max__ = p;
            for (size_t j_1__ = 0; j_1__ < beta_sigma_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "beta_sigma" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }

        if (!include_gqs__) return;
    }

}; // model

}  // namespace

typedef condhet_horseshoe_model_namespace::condhet_horseshoe_model stan_model;

#ifndef USING_R

stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}

#endif

