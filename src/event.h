#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include <iostream>
#include "stats.hpp"
class Event {
    public:
        //事件定义
        arma::colvec Time;      //事件时间
        arma::ivec   Outcome;   //事件类型
        arma::uvec   Group;     //发生组别, 注意下标从0开始
        arma::colvec Gaus;      //高斯结果
        
        //参数定义
        arma::colvec beta;          
        double theta;
        double ptheta;
        arma::colvec alpha;
        arma::colvec gamma;
        arma::colvec pgamma;
        arma::colvec tau;
        arma::colvec eta;
        arma::colvec omega;
        double Imin;
        double loglikelihood;
        
        //高斯参数定义
        arma::colvec kernal_scale;
        arma::colvec kernal_scale_prior;
        double kernal_variance;
        arma::uword MPA_Level;

        //先验定义
        arma::mat beta_priors;
        arma::colvec theta_priors;
        arma::colvec ptheta_priors;
        arma::mat alpha_priors;
        arma::mat gamma_priors;
        arma::mat pgamma_priors;
        arma::mat tau_priors;
        arma::mat eta_priors;
        arma::mat omega_priors;
        double Imin_priors;

        //特殊事件定义
        double Time_Last;
        arma::ivec N;
        arma::uword Imin_group;
        arma::uword num_group;

        //debug
        arma::uword num_insert;
        arma::uword num_delete;
        // 成员函数声明
        //Event(arma::colvec Time, arma::ivec Outcome, arma::uvec Group); 

        void sort(); 
        void gibbs_sampling();
        void Gaus_sampling();
        void log_likelihood();
        arma::mat Gaus_Output(arma::uword group);

        Event(arma::colvec time, arma::ivec outcome, arma::uvec group, arma::ivec n, arma::colvec kernel_prior) 
        {
            arma::uvec idx;

            idx = arma::sort_index(time);
            Time = time.elem(idx);
            Outcome = outcome.elem(idx);
            Group = group.elem(idx);
            N = n;

            // set derivative variables
            num_group = Group.max() + 1;
            Imin_group = Group(0);
  
            beta_priors.set_size(num_group, 2);
            theta_priors.set_size(2);
            ptheta_priors.set_size(2);
            alpha_priors.set_size(num_group, 2);
            gamma_priors.set_size(num_group, 2);
            pgamma_priors.set_size(num_group, 2);
            tau_priors.set_size(num_group, 2);
            eta_priors.set_size(num_group, 2);
            omega_priors.set_size(num_group, 2);

            beta_priors.fill(1);
            theta_priors.fill(1);
            ptheta_priors.fill(1);
            alpha_priors.fill(1);
            gamma_priors.fill(1);
            pgamma_priors.fill(1);
            tau_priors.fill(1);
            eta_priors.fill(1);
            omega_priors.fill(1);
            Imin_priors = 0;


            beta.set_size(num_group);
            alpha.set_size(num_group);
            gamma.set_size(num_group);
            pgamma.set_size(num_group);
            tau.set_size(num_group);
            eta.set_size(num_group);
            omega.set_size(num_group);

            Time_Last = time.max();
            Imin = time.min();

            kernal_scale_prior = kernel_prior;
            kernal_scale.set_size(num_group);
            for(arma::uword i = 0; i < num_group; i++)
            {
                kernal_scale(i) = stats::rexp(kernal_scale_prior(i));
            }
            
            kernal_variance = 2;
            Gaus.set_size(time.size());
            Gaus.fill(arma::datum::nan);
            MPA_Level = 1000;
            
        }

};

void Event :: sort()
{
    arma::uvec idx;

    idx = arma::sort_index(Time);
    Time = Time.elem(idx);
    Outcome = Outcome.elem(idx);
    Group = Group.elem(idx);

}



void Event :: gibbs_sampling()
{   
    
    // declare variables 
    arma::ivec x(num_group);                        
    // Initialize to be zero
    arma::ivec y(num_group, arma::fill::zeros);
    arma::ivec z(num_group, arma::fill::zeros);
    arma::ivec u(num_group, arma::fill::zeros); 
    arma::ivec v(num_group, arma::fill::zeros); //住院人数
    arma::ivec w(num_group, arma::fill::zeros); //康复人数

    arma::colvec beta_int(num_group, arma::fill::zeros);
    double theta_int = 0;
    arma::colvec gamma_int(num_group, arma::fill::zeros);
    arma::colvec eta_int(num_group, arma::fill::zeros);
    arma::colvec alpha_int(num_group, arma::fill::zeros);
    arma::colvec tau_int(num_group, arma::fill::zeros);
    arma::colvec omega_int(num_group, arma::fill::zeros);

    arma::ivec m(num_group, arma::fill::zeros);
    arma::ivec M(num_group, arma::fill::zeros);
    arma::sword q = 0;
    arma::sword Q = 0;
    arma::ivec e(num_group, arma::fill::zeros);
    arma::ivec E(num_group, arma::fill::zeros);
    arma::ivec f(num_group, arma::fill::zeros);
    arma::ivec b(num_group, arma::fill::zeros);
    arma::ivec a(num_group, arma::fill::zeros);
    arma::ivec g(num_group, arma::fill::zeros);
    
    double dt;

    // Start of the disease
    x = N;
    
    x(Imin_group) -=1;
    y(Imin_group) +=1;
   

    for (arma::uword i = 1; i < Time.n_elem; i++)
    {    
        dt = Time(i) - Time(i - 1);
        switch (Outcome(i))
        {
        case -1: //Thinned
            M(Group(i)) += 1;
            //Covert ivec to colvec
            beta_int    += dt * arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            theta_int   += arma::sum(dt * arma::conv_to<arma::colvec>::from(y));
            gamma_int   += dt * arma::conv_to<arma::colvec>::from(z);
            eta_int     += dt * arma::conv_to<arma::colvec>::from(u);
            alpha_int   += dt * arma::conv_to<arma::colvec>::from(w);
            tau_int     += dt * arma::conv_to<arma::colvec>::from(v);
            omega_int   += dt * arma::conv_to<arma::colvec>::from(y);

            break;
        case 1: // Infection
            m(Group(i)) += 1;
            //Covert ivec to colvec
            beta_int    += dt * arma::conv_to<arma::colvec>::from((x * arma::as_scalar(arma::sum(y + z + u)))); 
            theta_int   += arma::sum(dt * arma::conv_to<arma::colvec>::from(y));
            gamma_int   += dt * arma::conv_to<arma::colvec>::from(z);
            eta_int     += dt * arma::conv_to<arma::colvec>::from(u);
            alpha_int   += dt * arma::conv_to<arma::colvec>::from(w);
            tau_int     += dt * arma::conv_to<arma::colvec>::from(v);
            omega_int   += dt * arma::conv_to<arma::colvec>::from(y);

            x(Group(i)) -= 1;
            y(Group(i)) += 1;


            break;
        case 2: // Diagnosis
            q += 1;
            //Covert ivec to colvec
            beta_int    += dt * arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            theta_int   += arma::sum(dt * arma::conv_to<arma::colvec>::from(y));
            gamma_int   += dt * arma::conv_to<arma::colvec>::from(z);
            eta_int     += dt * arma::conv_to<arma::colvec>::from(u);
            alpha_int   += dt * arma::conv_to<arma::colvec>::from(w);
            tau_int     += dt * arma::conv_to<arma::colvec>::from(v);
            omega_int   += dt * arma::conv_to<arma::colvec>::from(y);

            y(Group(i)) -= 1;
            z(Group(i)) += 1;
            break;
        case 3: // Undetected
            Q += 1;
            //Covert ivec to colvec
            beta_int    += dt * arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            theta_int   += arma::sum(dt * arma::conv_to<arma::colvec>::from(y));
            gamma_int   += dt * arma::conv_to<arma::colvec>::from(z);
            eta_int     += dt * arma::conv_to<arma::colvec>::from(u);
            alpha_int   += dt * arma::conv_to<arma::colvec>::from(w);
            tau_int     += dt * arma::conv_to<arma::colvec>::from(v);
            omega_int   += dt * arma::conv_to<arma::colvec>::from(y);

            y(Group(i)) -= 1;
            u(Group(i)) += 1;
            break;
        case 4: // Diagnosis Hospitalized
            e(Group(i)) += 1;
            //Covert ivec to colvec
            beta_int    += dt * arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            theta_int   += arma::sum(dt * arma::conv_to<arma::colvec>::from(y));
            gamma_int   += dt * arma::conv_to<arma::colvec>::from(z);
            eta_int     += dt * arma::conv_to<arma::colvec>::from(u);
            alpha_int   += dt * arma::conv_to<arma::colvec>::from(w);
            tau_int     += dt * arma::conv_to<arma::colvec>::from(v);
            omega_int   += dt * arma::conv_to<arma::colvec>::from(y);

            z(Group(i)) -= 1;
            v(Group(i)) += 1;
            break;
        case 5: // Diagnosis Recovered
            E(Group(i)) += 1;
            //Covert ivec to colvec
            beta_int    += dt * arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            theta_int   += arma::sum(dt * arma::conv_to<arma::colvec>::from(y));
            gamma_int   += dt * arma::conv_to<arma::colvec>::from(z);
            eta_int     += dt * arma::conv_to<arma::colvec>::from(u);
            alpha_int   += dt * arma::conv_to<arma::colvec>::from(w);
            tau_int     += dt * arma::conv_to<arma::colvec>::from(v);
            omega_int   += dt * arma::conv_to<arma::colvec>::from(y);

            z(Group(i)) -= 1;
            w(Group(i)) += 1;
            break;
        case 6: // Undetected Recovered
            f(Group(i)) += 1;
            //Covert ivec to colvec
            beta_int    += dt * arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            theta_int   += arma::sum(dt * arma::conv_to<arma::colvec>::from(y));
            gamma_int   += dt * arma::conv_to<arma::colvec>::from(z);
            eta_int     += dt * arma::conv_to<arma::colvec>::from(u);
            alpha_int   += dt * arma::conv_to<arma::colvec>::from(w);
            tau_int     += dt * arma::conv_to<arma::colvec>::from(v);
            omega_int   += dt * arma::conv_to<arma::colvec>::from(y);

            u(Group(i)) -= 1;
            w(Group(i)) += 1;
            break;
        case 7: //  Hospitalized Suspected
            a(Group(i)) += 1;
            //Covert ivec to colvec
            beta_int    += dt * arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            theta_int   += arma::sum(dt * arma::conv_to<arma::colvec>::from(y));
            gamma_int   += dt * arma::conv_to<arma::colvec>::from(z);
            eta_int     += dt * arma::conv_to<arma::colvec>::from(u);
            alpha_int   += dt * arma::conv_to<arma::colvec>::from(w);
            tau_int     += dt * arma::conv_to<arma::colvec>::from(v);
            omega_int   += dt * arma::conv_to<arma::colvec>::from(y);

            v(Group(i)) -= 1;
            x(Group(i)) += 1;
            break;
        case 8: // Recovered Suspected
            b(Group(i)) += 1;
            //Covert ivec to colvec
            beta_int    += dt * arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            theta_int   += arma::sum(dt * arma::conv_to<arma::colvec>::from(y));
            gamma_int   += dt * arma::conv_to<arma::colvec>::from(z);
            eta_int     += dt * arma::conv_to<arma::colvec>::from(u);
            alpha_int   += dt * arma::conv_to<arma::colvec>::from(w);
            tau_int     += dt * arma::conv_to<arma::colvec>::from(v);
            omega_int   += dt * arma::conv_to<arma::colvec>::from(y);

            w(Group(i)) -= 1;
            x(Group(i)) += 1;
            break;
        case 9: // Death
            g(Group(i)) += 1;
            //Covert ivec to colvec
            beta_int    += dt * arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            theta_int   += arma::sum(dt * arma::conv_to<arma::colvec>::from(y));
            gamma_int   += dt * arma::conv_to<arma::colvec>::from(z);
            eta_int     += dt * arma::conv_to<arma::colvec>::from(u);
            alpha_int   += dt * arma::conv_to<arma::colvec>::from(w);
            tau_int     += dt * arma::conv_to<arma::colvec>::from(v);
            omega_int   += dt * arma::conv_to<arma::colvec>::from(y);

            z(Group(i)) -= 1;
            break;
        case 10: //Vaccination
            //Covert ivec to colvec
            beta_int    += dt * arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            theta_int   += arma::sum(dt * arma::conv_to<arma::colvec>::from(y));
            gamma_int   += dt * arma::conv_to<arma::colvec>::from(z);
            eta_int     += dt * arma::conv_to<arma::colvec>::from(u);
            alpha_int   += dt * arma::conv_to<arma::colvec>::from(w);
            tau_int     += dt * arma::conv_to<arma::colvec>::from(v);
            omega_int   += dt * arma::conv_to<arma::colvec>::from(y);

            x(Group(i)) -= 1;
            x(Group(i) + 1) += 1;
            break;
        case 11: //Time Last
            //Covert ivec to colvec
            beta_int    += dt * arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            theta_int   += arma::sum(dt * arma::conv_to<arma::colvec>::from(y));
            gamma_int   += dt * arma::conv_to<arma::colvec>::from(z);
            eta_int     += dt * arma::conv_to<arma::colvec>::from(u);
            alpha_int   += dt * arma::conv_to<arma::colvec>::from(w);
            tau_int     += dt * arma::conv_to<arma::colvec>::from(v);
            omega_int   += dt * arma::conv_to<arma::colvec>::from(y);
            break;
        default:
            std::cerr << "Error message : Undefined Outcome Type" << std::endl;
            break;
        }
        bool mask = arma::any(x < 0 or y < 0 or z < 0 or u < 0 or v < 0 or w < 0);
        if(mask)
        {
            std::cerr << "Negetive Num of people"<< y << i << std::endl;
            std::abort();
        }

    }

    // sampling
    bool delta;
    double Imin_diff;
    double Imin_propose;
    arma::uvec I12_idx;
    arma::uvec Imin_fail;
    #pragma omp parallel for
    for (arma::uword i = 0; i < num_group; i++)
    {   
        delta = i == Imin_group;
        beta(i) = stats::rgamma(m(i) + M(i) + beta_priors(i, 0) - delta, 1/(beta_int(i) + beta_priors(i, 1)));
        
        pgamma(i) = stats::rbeta(pgamma_priors(i, 0) + e(i), pgamma_priors(i, 1) + E(i));
        gamma(i) = stats::rgamma(gamma_priors(i, 0) + e(i) + E(i), 1/(gamma_priors(i, 1) + gamma_int(i)));
        eta(i) = stats::rgamma(eta_priors(i,0) + f(i), 1/(eta_priors(i,1) + eta_int(i)));
        alpha(i) = stats::rgamma(alpha_priors(i, 0) + b(i), 1/(alpha_priors(i, 1) + alpha_int(i)));
        tau(i) = stats::rgamma(tau_priors(i, 0) + a(i), 1/(tau_priors(i, 1) + tau_int(i)));
        omega(i) = stats::rgamma(omega_priors(i, 0) + g(i), 1/(omega_priors(i, 1) + omega_int(i)));


    }

    ptheta = stats::rbeta(ptheta_priors(0) + q, ptheta_priors(1) + Q);
    theta = stats::rgamma(q + Q + theta_priors(0), 1/(theta_int + theta_priors(1)));
    
    Imin_diff = stats::rexp(theta + Imin_priors + beta(Imin_group) * (sum(N) - 1));
    I12_idx = arma::find(Outcome == 1, 2);
    Imin_propose = Time(I12_idx(1)) - Imin_diff;
    double Imin_old = Time.min();
    Time(0) = Imin_propose;
    if(!Time.is_sorted())
    {
        Time(0) = Imin_old;
    }
    else
    {
        Imin =  Imin_propose;
    }
    



}

arma::mat rdist(arma::colvec x, arma::colvec y) 
{
    x = arma::sort(x);
    y = arma::sort(y);

    arma::mat distmat(x.size(), y.size(), arma::fill::none);
    #pragma omp parallel for schedule(static)
    for (arma::uword i = 0; i < y.size(); i++)
    {
        // 遍历y
        distmat.col(i) = x - y(i);
    }
    
    return arma::abs(distmat);
}

void Event :: Gaus_sampling()
{
    // declare variables
    arma::uvec gaus_idx;
    arma::mat dist_mat;
    arma::mat cov_mat;
    arma::colvec mean;
    arma::mat adj_dist_mat;
    arma::mat adj_mat;
    arma::colvec pseudo_gaus;
    arma::colvec MPA_Time;

     for(arma::uword i = 0; i < num_group; i ++)
     { 
        gaus_idx = arma::find((Outcome == 1 or Outcome == -1) and Group == i and Time != Imin);
        
        if (0 < gaus_idx.size()  and gaus_idx.size() < MPA_Level)
        {   // Normal Sampling
            dist_mat = rdist(Time.elem(gaus_idx), Time.elem(gaus_idx));
            cov_mat = pow(kernal_variance, 2) * arma::exp(-arma::pow(dist_mat, 2) / (2 * pow(kernal_scale(i), 2)));
            Gaus.elem(gaus_idx) = arma::mvnrnd(mean.zeros(gaus_idx.size()), cov_mat);

        } 
        else if ( gaus_idx.size() >= MPA_Level)
        {   // Mean Projection Approximation
            // Sample appropriate  pseudo data
            MPA_Time = Time(arma::sort(gaus_idx(arma::randperm(gaus_idx.size(), MPA_Level))));
            dist_mat = rdist(MPA_Time, MPA_Time);
            cov_mat = pow(kernal_variance, 2) * arma::exp(-arma::pow(dist_mat, 2) / (2 * pow(kernal_scale(i), 2)));
            cov_mat.diag() += arma::randu();    // Add some noise to make sure cov_mat is sympd;
            pseudo_gaus = arma::mvnrnd(mean.zeros(MPA_Level), cov_mat);

            // Make MPA adjustments
            adj_dist_mat = rdist(Time(gaus_idx), MPA_Time);
            adj_mat = pow(kernal_variance, 2) * arma::exp(-arma::pow(adj_dist_mat, 2) / (2 * pow(kernal_scale(i), 2)));
            Gaus.elem(gaus_idx) = adj_mat * arma::inv_sympd(cov_mat, arma::inv_opts::allow_approx) * pseudo_gaus;
        }
        else std::cerr << "Group"<< i << "Doesn't have Infection or Thinned" << std::endl;
        
     }
}

double logit(double x)
{
    return pow(1 + exp(-x), -1);
}

void Event :: log_likelihood()
{
    // declare variables
    arma::ivec x(num_group);                        
    // Initialize to be zero
    arma::ivec y(num_group, arma::fill::zeros);
    arma::ivec z(num_group, arma::fill::zeros);
    arma::ivec u(num_group, arma::fill::zeros); 
    arma::ivec v(num_group, arma::fill::zeros); //住院人数
    arma::ivec w(num_group, arma::fill::zeros); //康复人数

    arma::colvec Integral1(num_group, arma::fill::zeros);
    arma::colvec Integral2(num_group, arma::fill::zeros);
    arma::colvec Integral3(num_group, arma::fill::zeros);
    arma::colvec Integral4(num_group, arma::fill::zeros);
    arma::colvec Integral5(num_group, arma::fill::zeros);
    arma::colvec Integral6(num_group, arma::fill::zeros);
    arma::colvec Integral7(num_group, arma::fill::zeros);

    double dt;
    double res = 0;
    bool mask;

    // Start of the disease
    x = N;
    x(Imin_group) -=1;
    y(Imin_group) +=1;

    for (arma::uword i = 1; i < Time.n_elem; i++)
    {   
        dt = Time(i) - Time(i - 1);
        
        switch (Outcome(i))
        {
        case -1: //Thinned
            res += log(beta(Group(i))) + log(x(Group(i))) + log(arma::sum(y + z + u)) + log(logit(-Gaus(i)));

            //Covert ivec to colvec
            Integral1   = beta * dt % arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            Integral2   = theta * dt * arma::conv_to<arma::colvec>::from(y);
            Integral3   = gamma * dt % arma::conv_to<arma::colvec>::from(z);
            Integral4   = eta * dt % arma::conv_to<arma::colvec>::from(u);
            Integral5   = tau * dt % arma::conv_to<arma::colvec>::from(v);
            Integral6   = alpha * dt % arma::conv_to<arma::colvec>::from(w);
            Integral7   = omega * dt % arma::conv_to<arma::colvec>::from(z);

            res = res - arma::sum(Integral1 + Integral2 + Integral3 + Integral4 + Integral5 + Integral6 + Integral7);


            break;
        case 1: // Infection
            res += log(beta(Group(i))) + log(x(Group(i))) + log(arma::sum(y + z + u)) + log(logit(Gaus(i)));

            //Covert ivec to colvec
            Integral1   = beta * dt % arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            Integral2   = theta * dt * arma::conv_to<arma::colvec>::from(y);
            Integral3   = gamma * dt % arma::conv_to<arma::colvec>::from(z);
            Integral4   = eta * dt % arma::conv_to<arma::colvec>::from(u);
            Integral5   = tau * dt % arma::conv_to<arma::colvec>::from(v);
            Integral6   = alpha * dt % arma::conv_to<arma::colvec>::from(w);
            Integral7   = omega * dt % arma::conv_to<arma::colvec>::from(z);

            res = res - arma::sum(Integral1 + Integral2 + Integral3 + Integral4 + Integral5 + Integral6 + Integral7);
            x(Group(i)) -= 1;
            y(Group(i)) += 1;

            break;
        case 2: // Diagnosis
            res += log(ptheta) + log(theta) + log(y(Group(i)));

            //Covert ivec to colvec
            Integral1   = beta * dt % arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            Integral2   = theta * dt * arma::conv_to<arma::colvec>::from(y);
            Integral3   = gamma * dt % arma::conv_to<arma::colvec>::from(z);
            Integral4   = eta * dt % arma::conv_to<arma::colvec>::from(u);
            Integral5   = tau * dt % arma::conv_to<arma::colvec>::from(v);
            Integral6   = alpha * dt % arma::conv_to<arma::colvec>::from(w);
            Integral7   = omega * dt % arma::conv_to<arma::colvec>::from(z);

            res = res - arma::sum(Integral1 + Integral2 + Integral3 + Integral4 + Integral5 + Integral6 + Integral7);

            y(Group(i)) -= 1;
            z(Group(i)) += 1;

            break;
        case 3: // Undetected
            res += log(1 - ptheta) + log(theta) + log(y(Group(i)));

            //Covert ivec to colvec
            Integral1   = beta * dt % arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            Integral2   = theta * dt * arma::conv_to<arma::colvec>::from(y);
            Integral3   = gamma * dt % arma::conv_to<arma::colvec>::from(z);
            Integral4   = eta * dt % arma::conv_to<arma::colvec>::from(u);
            Integral5   = tau * dt % arma::conv_to<arma::colvec>::from(v);
            Integral6   = alpha * dt % arma::conv_to<arma::colvec>::from(w);
            Integral7   = omega * dt % arma::conv_to<arma::colvec>::from(z);

            res = res - arma::sum(Integral1 + Integral2 + Integral3 + Integral4 + Integral5 + Integral6 + Integral7);

            y(Group(i)) -= 1;
            u(Group(i)) += 1;


            break;
        case 4: // Diagnosis Hospitalized
            res += log(pgamma(Group(i))) + log(gamma(Group(i))) + log(z(Group(i)));

            //Covert ivec to colvec
            Integral1   = beta * dt % arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            Integral2   = theta * dt * arma::conv_to<arma::colvec>::from(y);
            Integral3   = gamma * dt % arma::conv_to<arma::colvec>::from(z);
            Integral4   = eta * dt % arma::conv_to<arma::colvec>::from(u);
            Integral5   = tau * dt % arma::conv_to<arma::colvec>::from(v);
            Integral6   = alpha * dt % arma::conv_to<arma::colvec>::from(w);
            Integral7   = omega * dt % arma::conv_to<arma::colvec>::from(z);

            res = res - arma::sum(Integral1 + Integral2 + Integral3 + Integral4 + Integral5 + Integral6 + Integral7);
            z(Group(i)) -= 1;
            v(Group(i)) += 1;

            break;
        case 5: // Diagnosis Recovered
            res += log(1 - pgamma(Group(i))) + log(gamma(Group(i))) + log(z(Group(i)));

            //Covert ivec to colvec
            Integral1   = beta * dt % arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            Integral2   = theta * dt * arma::conv_to<arma::colvec>::from(y);
            Integral3   = gamma * dt % arma::conv_to<arma::colvec>::from(z);
            Integral4   = eta * dt % arma::conv_to<arma::colvec>::from(u);
            Integral5   = tau * dt % arma::conv_to<arma::colvec>::from(v);
            Integral6   = alpha * dt % arma::conv_to<arma::colvec>::from(w);
            Integral7   = omega * dt % arma::conv_to<arma::colvec>::from(z);

            res = res - arma::sum(Integral1 + Integral2 + Integral3 + Integral4 + Integral5 + Integral6 + Integral7);

            z(Group(i)) -= 1;
            w(Group(i)) += 1;

            break;
        case 6: // Undetected Recovered
            res += log(eta(Group(i))) + log(u(Group(i)));

            //Covert ivec to colvec
            Integral1   = beta * dt % arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            Integral2   = theta * dt * arma::conv_to<arma::colvec>::from(y);
            Integral3   = gamma * dt % arma::conv_to<arma::colvec>::from(z);
            Integral4   = eta * dt % arma::conv_to<arma::colvec>::from(u);
            Integral5   = tau * dt % arma::conv_to<arma::colvec>::from(v);
            Integral6   = alpha * dt % arma::conv_to<arma::colvec>::from(w);
            Integral7   = omega * dt % arma::conv_to<arma::colvec>::from(z);

            res = res - arma::sum(Integral1 + Integral2 + Integral3 + Integral4 + Integral5 + Integral6 + Integral7);

            u(Group(i)) -= 1;
            w(Group(i)) += 1;

            break;
        case 7: //  Hospitalized Suspected
            res += log(tau(Group(i))) + log(v(Group(i)));

            //Covert ivec to colvec
            Integral1   = beta * dt % arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            Integral2   = theta * dt * arma::conv_to<arma::colvec>::from(y);
            Integral3   = gamma * dt % arma::conv_to<arma::colvec>::from(z);
            Integral4   = eta * dt % arma::conv_to<arma::colvec>::from(u);
            Integral5   = tau * dt % arma::conv_to<arma::colvec>::from(v);
            Integral6   = alpha * dt % arma::conv_to<arma::colvec>::from(w);
            Integral7   = omega * dt % arma::conv_to<arma::colvec>::from(z);

            res = res - arma::sum(Integral1 + Integral2 + Integral3 + Integral4 + Integral5 + Integral6 + Integral7);

            v(Group(i)) -= 1;
            x(Group(i)) += 1;

            break;
        case 8: // Recovered Suspected
            res += log((Group(i))) + log(w(Group(i)));

            //Covert ivec to colvec
            Integral1   = beta * dt % arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            Integral2   = theta * dt * arma::conv_to<arma::colvec>::from(y);
            Integral3   = gamma * dt % arma::conv_to<arma::colvec>::from(z);
            Integral4   = eta * dt % arma::conv_to<arma::colvec>::from(u);
            Integral5   = tau * dt % arma::conv_to<arma::colvec>::from(v);
            Integral6   = alpha * dt % arma::conv_to<arma::colvec>::from(w);
            Integral7   = omega * dt % arma::conv_to<arma::colvec>::from(z);

            res = res - arma::sum(Integral1 + Integral2 + Integral3 + Integral4 + Integral5 + Integral6 + Integral7);

            w(Group(i)) -= 1;
            x(Group(i)) += 1;


            break;
        case 9: // Death
            res += log(omega(Group(i))) + log(z(Group(i)));

            //Covert ivec to colvec
            Integral1   = beta * dt % arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            Integral2   = theta * dt * arma::conv_to<arma::colvec>::from(y);
            Integral3   = gamma * dt % arma::conv_to<arma::colvec>::from(z);
            Integral4   = eta * dt % arma::conv_to<arma::colvec>::from(u);
            Integral5   = tau * dt % arma::conv_to<arma::colvec>::from(v);
            Integral6   = alpha * dt % arma::conv_to<arma::colvec>::from(w);
            Integral7   = omega * dt % arma::conv_to<arma::colvec>::from(z);

            res = res - arma::sum(Integral1 + Integral2 + Integral3 + Integral4 + Integral5 + Integral6 + Integral7);

            z(Group(i)) -= 1;

            break;
        case 10: //Vaccination
            //Covert ivec to colvec
            Integral1   = beta * dt % arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            Integral2   = theta * dt * arma::conv_to<arma::colvec>::from(y);
            Integral3   = gamma * dt % arma::conv_to<arma::colvec>::from(z);
            Integral4   = eta * dt % arma::conv_to<arma::colvec>::from(u);
            Integral5   = tau * dt % arma::conv_to<arma::colvec>::from(v);
            Integral6   = alpha * dt % arma::conv_to<arma::colvec>::from(w);
            Integral7   = omega * dt % arma::conv_to<arma::colvec>::from(z);

            res = res - arma::sum(Integral1 + Integral2 + Integral3 + Integral4 + Integral5 + Integral6 + Integral7);

            x(Group(i)) -= 1;
            x(Group(i) + 1) += 1;

            break;
        case 11: //Time Last
            //Covert ivec to colvec
            Integral1   = beta * dt % arma::conv_to<arma::colvec>::from((x * arma::sum(y + z + u))); 
            Integral2   = theta * dt * arma::conv_to<arma::colvec>::from(y);
            Integral3   = gamma * dt % arma::conv_to<arma::colvec>::from(z);
            Integral4   = eta * dt % arma::conv_to<arma::colvec>::from(u);
            Integral5   = tau * dt % arma::conv_to<arma::colvec>::from(v);
            Integral6   = alpha * dt % arma::conv_to<arma::colvec>::from(w);
            Integral7   = omega * dt % arma::conv_to<arma::colvec>::from(z);

            res = res - arma::sum(Integral1 + Integral2 + Integral3 + Integral4 + Integral5 + Integral6 + Integral7);

            break;
        default:
            std::cerr << "Error message : Undefined Outcome Type" << Outcome(i) << std::endl;
            break;
        }

        mask = arma::any(x < 0 or y < 0 or z < 0 or u < 0 or v < 0 or w < 0);
        

        if(mask)
        {
            res = -arma::datum::inf;
            loglikelihood = res;
            break;
        }
        else
        {
            continue;
        }
    }

    loglikelihood = res;

}

Event Insert(Event &target,  arma::uword group, arma::sword outcome)
{
    // declare variables
    arma::colvec new_t(1);
    arma::colvec new_gaus(1); new_gaus.fill(arma::datum::nan); //Initialized to be NA;
    arma::uvec idx;
    arma::mat dist_mat;
    arma::mat cov_mat;
    arma::mat cov_inv_mat;
    arma::mat adj_dist_mat;
    arma::mat adj_mat;
    arma::colvec mean;
    arma::mat var;
    arma::uvec insert_idx;
    arma::ivec insert_Outcome(1, arma::fill::value(outcome));
    arma::uvec insert_Group(1, arma::fill::value(group));
    double unif = arma::randu<double>();
    double acc;
    double old_log_likelihood = target.loglikelihood;
    double new_log_likelihood;
    arma::uvec num_idx; // used for calculating num;
    arma::sword num;    // How many corresponding outcomes intotal


    new_t = arma::randu(1, arma::distr_param(target.Imin, target.Time_Last));
    insert_idx = arma::find(target.Time < arma::as_scalar(new_t), 1, "last") + 1;
    num_idx = arma::find((target.Outcome == outcome) and target.Group == group and target.Time != target.Imin);
    num = num_idx.size();


    if(outcome == 1 or outcome == -1)
    {   // Make Gaus Predictions
        
        idx = arma::find((target.Outcome == 1 or target.Outcome == -1) and target.Group == group and target.Time != target.Imin);

        if(0 < idx.size() and idx.size() < target.MPA_Level)
        {   
            dist_mat = rdist(target.Time(idx), target.Time(idx));
            cov_mat = pow(target.kernal_variance, 2) * arma::exp(-arma::pow(dist_mat, 2) / (2 * pow(target.kernal_scale(group), 2)));
            adj_dist_mat = rdist(new_t, target.Time(idx));
            adj_mat = pow(target.kernal_variance, 2) * arma::exp(-arma::pow(adj_dist_mat, 2) / (2 * pow(target.kernal_scale(group), 2)));

            cov_inv_mat = arma::inv_sympd(cov_mat, arma::inv_opts::allow_approx);
            mean = adj_mat * cov_inv_mat * target.Gaus(idx);
            var = pow(target.kernal_variance, 2) - adj_mat * cov_inv_mat * adj_mat.t();
            //var.diag() += arma::randu();
            if(var(0, 0) < 0) var(0, 0) = 0;
            new_gaus = arma::mvnrnd(mean, var);

        }   
        else if(idx.size() >= target.MPA_Level)
        {   // MPA
            arma::uvec MPA_idx = arma::sort(idx(arma::randperm(idx.size(), target.MPA_Level)));
            arma::colvec MPA_Time = target.Time(MPA_idx);

            dist_mat = rdist(MPA_Time, MPA_Time);
            cov_mat = pow(target.kernal_variance, 2) * arma::exp(-arma::pow(dist_mat, 2) / (2 * pow(target.kernal_scale(group), 2)));
            cov_mat.diag() += arma::randu(); //add some noise
            adj_dist_mat = rdist(new_t, MPA_Time);
            adj_mat = pow(target.kernal_variance, 2) * arma::exp(-arma::pow(adj_dist_mat, 2) / (2 * pow(target.kernal_scale(group), 2)));
            if(!cov_mat.is_sympd()) return target;
            cov_inv_mat = arma::inv_sympd(cov_mat, arma::inv_opts::allow_approx);
            mean = adj_mat * cov_inv_mat * target.Gaus(MPA_idx);
            var = (pow(target.kernal_variance, 2) - adj_mat * cov_inv_mat * adj_mat.t()) + arma::randu(); //add some noise;
            if(var(0, 0) < 0) var(0, 0) = 0;
            new_gaus = arma::mvnrnd(mean, var);

        }
        else
        {
            std::cerr << "Group" << group <<" Need At Least One Infection or Thinned to Make Gaus Prediction" << std::endl;
        }

    }


    // insert it and evaluate  likelihood
    target.Time.insert_rows(arma::as_scalar(insert_idx), new_t);
    target.Outcome.insert_rows(arma::as_scalar(insert_idx), insert_Outcome);
    target.Group.insert_rows(arma::as_scalar(insert_idx), insert_Group); 
    target.Gaus.insert_rows(arma::as_scalar(insert_idx), new_gaus);
    target.log_likelihood();
    acc = log(target.Time_Last - target.Imin) + (target.loglikelihood - old_log_likelihood) - log(double(num + 1));
    
    if(acc > log(unif))
    {   
        return target;
    }
    else
    {
        target.Time.shed_rows(insert_idx);
        target.Outcome.shed_rows(insert_idx);
        target.Group.shed_rows(insert_idx);
        target.Gaus.shed_rows(insert_idx);
        target.loglikelihood = old_log_likelihood;
        return target;
    }
}

Event Delete(Event &target, arma::uword group, arma::sword outcome)
{
    // declare variables
    arma::colvec old_time(1);
    arma::ivec old_Outcome(1);
    arma::uvec old_Group(1);
    arma::colvec old_Gaus(1);
    arma::uvec old_idx;
    arma::uvec idx;
    arma::uword num_idx;
    double old_loglikelihood;
    double new_loglikelihood;
    double acc;
    double unif = arma::randu<double>();

    idx = arma::find(target.Outcome == outcome and target.Group == group and target.Time != target.Imin);
    num_idx = idx.size();
    if(num_idx == 0) return target;
    if(num_idx == 1 and (outcome == 1 or outcome == -1)) return target;

    old_idx = idx(arma::randperm(num_idx, 1));

    // preserve old stuff
    old_loglikelihood = target.loglikelihood;
    old_time = target.Time(old_idx);
    old_Outcome = target.Outcome(old_idx);
    old_Group = target.Group(old_idx);
    old_Gaus = target.Gaus(old_idx);

    // delete old stuff
    target.Time.shed_rows(old_idx);   
    target.Outcome.shed_rows(old_idx);
    target.Group.shed_rows(old_idx);
    target.Gaus.shed_rows(old_idx);
    target.log_likelihood();
    new_loglikelihood = target.loglikelihood;
    acc = log(double(num_idx)) + (new_loglikelihood - old_loglikelihood) - log(target.Time_Last - target.Imin);
    
   
    if (acc > log(unif))
    {   
        return target;
    }
    else
    {
        target.loglikelihood = old_loglikelihood;
        target.Time.insert_rows(arma::as_scalar(old_idx), old_time);
        target.Outcome.insert_rows(arma::as_scalar(old_idx), old_Outcome);
        target.Group.insert_rows(arma::as_scalar(old_idx), old_Group);
        target.Gaus.insert_rows(arma::as_scalar(old_idx), old_Gaus);
        return target;
    }    
}


Event Update(Event &target, arma::uword group, arma::sword outcome) {
    // declare variables
    arma::uvec idx;
    arma::colvec new_t(1);
    arma::uvec new_idx(1);
    arma::ivec new_Outcome(1, arma::fill::value(outcome));
    arma::uvec new_Group(1, arma::fill::value(group));
    arma::colvec new_Gaus(1); new_Gaus.fill(arma::datum::nan);
    arma::mat dist_mat;
    arma::mat cov_mat;
    arma::mat cov_inv_mat;
    arma::mat adj_dist_mat;
    arma::mat adj_mat;
    arma::colvec mean;
    arma::mat var;
    double new_loglikelihood;


    arma::colvec old_t(1);
    arma::uvec old_idx(1);
    arma::ivec old_Outcome(1);
    arma::uvec old_Group(1);
    arma::colvec old_Gaus(1);
    
    double old_log_likelihood = target.loglikelihood;

    double acc;
    double unif = arma::randu<double>();


    // new stuff
    new_t = arma::randu(1, arma::distr_param(target.Imin, target.Time_Last));
    new_idx = arma::find(target.Time < arma::as_scalar(new_t), 1, "last") + 1;

    idx = arma::find(target.Outcome == outcome && target.Group == group && target.Time != target.Imin && target.Time != arma::as_scalar(new_t));
    
    if(idx.size() <= 0){return target;} 
    
    if(outcome == 1 or outcome == -1)
    {   // Make Gaussian Predictions
        
        arma::uvec GP_idx = arma::find((target.Outcome == 1 or target.Outcome == -1) and target.Group == group and target.Time != target.Imin);
        if(0 < GP_idx.size() and GP_idx.size() < target.MPA_Level)
        {   // Normal

            dist_mat = rdist(target.Time(GP_idx), target.Time(GP_idx));
            cov_mat = pow(target.kernal_variance, 2) * arma::exp(-arma::pow(dist_mat, 2) / (2 * pow(target.kernal_scale(group), 2)));
            adj_dist_mat = rdist(new_t, target.Time(GP_idx));
            adj_mat = pow(target.kernal_variance, 2) * arma::exp(-arma::pow(adj_dist_mat, 2) / (2 * pow(target.kernal_scale(group), 2)));
            
            cov_inv_mat = arma::inv_sympd(cov_mat, arma::inv_opts::allow_approx);
            mean = adj_mat * cov_inv_mat * target.Gaus(GP_idx);
            var = pow(target.kernal_variance, 2) - adj_mat * cov_inv_mat * adj_mat.t();
            //var.diag() += arma::randu();
            if(var(0, 0) < 0) var(0, 0) = 0;
            new_Gaus = arma::mvnrnd(mean, var);
            
        }   
        else if(GP_idx.size() >= target.MPA_Level)
        {   // MPA
            arma::uvec MPA_idx = arma::sort(GP_idx(arma::randperm(GP_idx.size(), target.MPA_Level)));
            arma::colvec MPA_Time = target.Time(MPA_idx);

            dist_mat = rdist(MPA_Time, MPA_Time);
            cov_mat = pow(target.kernal_variance, 2) * arma::exp(-arma::pow(dist_mat, 2) / (2 * pow(target.kernal_scale(group), 2)));
            cov_mat.diag() += arma::randu(); //add some noise
            adj_dist_mat = rdist(new_t, MPA_Time);
            adj_mat = pow(target.kernal_variance, 2) * arma::exp(-arma::pow(adj_dist_mat, 2) / (2 * pow(target.kernal_scale(group), 2)));
            
            if(!cov_mat.is_sympd()) return target;
            
            cov_inv_mat = arma::inv_sympd(cov_mat, arma::inv_opts::allow_approx);
            mean = adj_mat * cov_inv_mat * target.Gaus(MPA_idx);
            var = (pow(target.kernal_variance, 2) - adj_mat * cov_inv_mat * adj_mat.t()) + arma::randu(); //add some noise;
            if(var(0, 0) < 0) var(0, 0) = 0;
            new_Gaus = arma::mvnrnd(mean, var);
            
        } else{
            std::cerr << "Group" << group << "Doesn't Have Any Gauss Events" << std::endl;
        }

    }
   
     // Update it and evaluate  likelihood(insert first then delete)
    target.Time.insert_rows(arma::as_scalar(new_idx), new_t);
    target.Outcome.insert_rows(arma::as_scalar(new_idx), new_Outcome);
    target.Group.insert_rows(arma::as_scalar(new_idx), new_Group); 
    target.Gaus.insert_rows(arma::as_scalar(new_idx), new_Gaus);

    // old stuff

    idx = arma::find(target.Outcome == outcome && target.Group == group && target.Time != target.Imin && target.Time != arma::as_scalar(new_t));

    old_idx = idx(arma::randperm(idx.size(), 1));
    old_t = target.Time(old_idx);
    old_Group = target.Group(old_idx);
    old_Outcome = target.Outcome(old_idx);
    old_Gaus = target.Gaus(old_idx);

    target.Time.shed_rows(old_idx);
    target.Outcome.shed_rows(old_idx);
    target.Group.shed_rows(old_idx);
    target.Gaus.shed_rows(old_idx);

    target.log_likelihood();
    new_loglikelihood = target.loglikelihood;

    acc = new_loglikelihood - old_log_likelihood;
    if(acc > log(unif))
    {   
        //if(outcome == 1 or outcome == -1) std::cout << "RUA" << std::endl;
        return target;

    }
    else
    {   
        // Revert log likelihood
        target.loglikelihood = old_log_likelihood;
        // Insert what deleted
        target.Time.insert_rows(arma::as_scalar(old_idx), old_t);
        target.Outcome.insert_rows(arma::as_scalar(old_idx), old_Outcome);
        target.Group.insert_rows(arma::as_scalar(old_idx), old_Group);
        target.Gaus.insert_rows(arma::as_scalar(old_idx), old_Gaus);
        // deleted what inserted
        new_idx = arma::find(target.Time == arma::as_scalar(new_t), 1);

        target.Time.shed_rows(new_idx);
        target.Outcome.shed_rows(new_idx);
        target.Group.shed_rows(new_idx);
        target.Gaus.shed_rows(new_idx);


        return target;
    }

}


Event Update_Gaus(Event &target, arma::uword group, double nu)
{
    // declare variables
    arma::colvec new_Gaus(target.Gaus); //Initialized to be the same as old Gaus
    arma::uvec idx;

    arma::mat dist_mat;
    arma::mat cov_mat;
    arma::colvec mean;
    arma::mat adj_dist_mat;
    arma::mat adj_mat;

    arma::colvec old_Gaus = target.Gaus;
    double old_loglikelihood = target.loglikelihood;
    double new_loglikelihood;
    double unif = arma::randu<double>();
    double acc;

    idx = arma::find((target.Outcome == 1 or target.Outcome == -1) and target.Group == group and target.Time != target.Imin);
    
    if(0 < idx.size() and idx.size() < target.MPA_Level)
    {   // Normal Sampling
        dist_mat = rdist(target.Time.elem(idx), target.Time.elem(idx));
        cov_mat = pow(target.kernal_variance, 2) * arma::exp(-arma::pow(dist_mat, 2) / (2 * pow(target.kernal_scale(group), 2)));
        new_Gaus(idx) = sqrt(1 - pow(nu, 2)) * (arma::mvnrnd(mean.zeros(idx.size()), cov_mat)) + nu * target.Gaus(idx);

    }
    else if(idx.size() >= target.MPA_Level)
    {   // MPA Sampling
        arma::colvec pseudo_gaus;
        arma::uvec MPA_idx = arma::sort(idx(arma::randperm(idx.size(), target.MPA_Level)));
        arma::colvec MPA_Time = target.Time(MPA_idx);

        dist_mat = rdist(MPA_Time, MPA_Time);
        cov_mat = pow(target.kernal_variance, 2) * arma::exp(-arma::pow(dist_mat, 2) / (2 * pow(target.kernal_scale(group), 2)));
        cov_mat.diag() += arma::randu();
        if(!cov_mat.is_sympd()) return target;
        pseudo_gaus = arma::mvnrnd(mean.zeros(target.MPA_Level), cov_mat);       
   
        // Make MPA adjustments
        adj_dist_mat = rdist(target.Time(idx), MPA_Time);
        adj_mat = pow(target.kernal_variance, 2) * arma::exp(-arma::pow(adj_dist_mat, 2) / (2 * pow(target.kernal_scale(group), 2)));
        cov_mat.diag() += arma::randu(); // Add some noise to make sure cov_mat is sympd;
        new_Gaus(idx) = adj_mat * arma::inv_sympd(cov_mat, arma::inv_opts::allow_approx) * pseudo_gaus;

    }
    else
    {    
        std::cerr << "Group" << group << "Doesn't Have Any Gauss Events" << std::endl;
    }

    // Evaluate Log likelihood
    target.Gaus = new_Gaus;
    target.log_likelihood();
    new_loglikelihood = target.loglikelihood;
    acc = new_loglikelihood - old_loglikelihood;
    if(acc > log(unif))
    {   
        return target;
    }
    else
    {
        //revert Gaus and likelihood;
        target.Gaus = old_Gaus;
        target.loglikelihood = old_loglikelihood;
        return target;
    }

}


Event Update_Kernel(Event &target, arma::uword group)
{
    //delclare variables
    double mean = stats::rexp(target.kernal_scale_prior(group));
    double var = 1.0;
    double new_Kernel_scale = arma::randn(arma::distr_param(mean, var));
    double old_Kernel_scale = target.kernal_scale(group);

    arma::uvec idx = arma::find((target.Outcome == 1 or target.Outcome == -1) and target.Group == group and target.Time != target.Imin);
    arma::mat dist_mat;
    arma::mat old_cov_mat;
    arma::mat new_cov_mat;
    arma::colvec gaus;

    double acc;
    double acc1;
    double acc2;
    double unif = arma::randu<double>();
    
    if( 0 < idx.size() and idx.size()  < target.MPA_Level)
    {
        dist_mat = rdist(target.Time(idx), target.Time(idx));
        old_cov_mat = pow(target.kernal_variance, 2) * arma::exp(-arma::pow(dist_mat, 2) / (2 * pow(target.kernal_scale(group), 2)));
        new_cov_mat = pow(target.kernal_variance, 2) * arma::exp(-arma::pow(dist_mat, 2) / (2 * pow(new_Kernel_scale, 2)));
        gaus = target.Gaus(idx);
    }
    else if(idx.size() >= target.MPA_Level)
    {   //MPA
        arma::uvec MPA_idx = arma::sort(idx(arma::randperm(idx.size(), target.MPA_Level)));
        arma::colvec MPA_Time = target.Time(MPA_idx);
        
        dist_mat = rdist(MPA_Time, MPA_Time);
        old_cov_mat = pow(target.kernal_variance, 2) * arma::exp(-arma::pow(dist_mat, 2) / (2 * pow(target.kernal_scale(group), 2)));
        old_cov_mat.diag() += arma::randu(); // add some noise
        new_cov_mat = pow(target.kernal_variance, 2) * arma::exp(-arma::pow(dist_mat, 2) / (2 * pow(new_Kernel_scale, 2)));
        new_cov_mat.diag() += arma::randu(); // add some noise
        gaus = target.Gaus(MPA_idx);

    }
    else
    {
        std::cerr << "Group" << group << "Doesn't Have Any Gauss Events" << std::endl;
    }

    if(!new_cov_mat.is_sympd()) return target;

    acc1 = -0.5 * log(arma::det(new_cov_mat)) + 
            (-0.5 * arma::as_scalar(gaus.t() * arma::inv_sympd(new_cov_mat, arma::inv_opts::allow_approx) * gaus) - 
                new_Kernel_scale * target.kernal_scale_prior(group));

    acc2 = -0.5 * log(arma::det(old_cov_mat)) + 
            (-0.5 * arma::as_scalar(gaus.t() * arma::inv_sympd(old_cov_mat, arma::inv_opts::allow_approx) * gaus) - 
                 target.kernal_scale(group) * target.kernal_scale_prior(group));

    acc = acc1 - acc2;
    if(acc > log(unif))
    {
        target.kernal_scale(group) = new_Kernel_scale;
        return target;
    }
    else
    {
        return target;
    }

}


void Data_Readin(std::string filename, arma::colvec &Time, arma::ivec &Outcome, arma::uvec &Group)
{

    arma::mat readin_mat;
    
    bool status = readin_mat.load(filename, arma::csv_ascii);
    if(status)
    {   
        Time = readin_mat.col(0);
        Outcome = arma::conv_to<arma::ivec>::from(readin_mat.col(1));
        Group = arma::conv_to<arma::uvec>::from(readin_mat.col(2));


    }
    else
    {
        std::cerr << "Read in Failed" << std::endl;
    }
    

}


void Data_Output(std::string filename, arma::mat target)
{
    


    bool status = target.save(filename, arma::csv_ascii);
    if(!status)
    {
        std::cerr << "Write out Failed" << std::endl;
    }

}

arma::mat Event::Gaus_Output(arma::uword group){
    arma::uvec idx = arma::find(Outcome == 1 and Group == group);
    arma::mat output(idx.size(), 2, arma::fill::none);

    output.col(0) = Time.elem(idx);
    output.col(1) = Gaus.elem(idx) * beta(group);
        
    return output;
    

}


