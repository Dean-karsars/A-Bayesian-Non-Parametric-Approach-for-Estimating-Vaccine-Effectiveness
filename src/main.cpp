#include "event.h"
#include "indicators.hpp"

int main(int argc, char** argv)
{    
    arma::colvec Time;
    arma::ivec Outcome;
    arma::uvec Group;
    /*arma::ivec N = {5546419, 368594     , 5453  , 9386706  , 3802429   , 440450  , 4456411  , 6521686   , 
                    707373  , 8167974 , 15038367  , 1524279, 9260197 , 25772332  , 4382122  , 4520946 , 
                    19591096  , 5197352  , 1398336 , 12318440  , 6912468   , 589395  , 3694941  , 2074299 };*/
    arma::ivec N = {22494, 115, 0};                
    arma::colvec Kernel_prior(arma::size(N), arma::fill::value(1E-4));
    arma::ivec Outcome_vec = {1, -1, 3, 5, 6, 7, 8};
    // initialize
    
    Data_Readin("/home/yao/C++/Covid_Research/NPBCOVID/data/KingsData.csv", Time, Outcome, Group);
    Event sim(Time, Outcome, Group, N, Kernel_prior); 
    sim.MPA_Level = 200;
    sim.sort();
    sim.gibbs_sampling();
    sim.Gaus_sampling();
    sim.log_likelihood();
    sim.kernal_variance = 1.011;
    std::cout << "log Likelihood" << sim.loglikelihood << std::endl;

    

    arma::uword n_iter = 100;
  
    arma::mat beta_output(n_iter, sim.num_group);
    arma::mat pgamma_output(n_iter, sim.num_group);
    arma::mat gamma_output(n_iter, sim.num_group);
    arma::mat alpha_output(n_iter, sim.num_group);
    arma::mat eta_output(n_iter, sim.num_group);
    arma::mat tau_output(n_iter, sim.num_group);
    arma::mat omega_output(n_iter, sim.num_group);
    arma::mat theta_output(n_iter, 1);
    arma::mat ptheta_output(n_iter, 1);

    indicators::ProgressBar bar{
        indicators::option::BarWidth{50},
        indicators::option::Start{"["},
        indicators::option::Fill{"■"},
        indicators::option::Lead{"■"},
        indicators::option::Remainder{"-"},
        indicators::option::End{" ]"},
        indicators::option::ShowElapsedTime{true},
        indicators::option::ShowRemainingTime{true},
        indicators::option::ShowPercentage{true},
        indicators::option::MaxProgress{n_iter}, 
        indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}
    };
    
    for (arma::uword i = 0; i < n_iter; i++)
    {
        arma::uword group = arma::randi<arma::uword>(arma::distr_param(0, sim.num_group - 1));
        for(arma::uword j = 0; j < Outcome_vec.size(); j++)
        {   
            arma::sword current_outcome = Outcome_vec(j);
            arma::uword Insert_or_delete = arma::randi<arma::uword>(arma::distr_param(0, 1));
            
            if(Insert_or_delete == 0)
            {   //Insert
                
                sim = Insert(sim, group, current_outcome);    
                
            }
            else 
            {   //Delete
  
                sim = Delete(sim, group, current_outcome);
                //std::cout << sim.Gaus(sim.MPA_idx(group)).has_nan()<< "Delete" << std::endl;
                
            }
             


            //Update
   
            sim = Update(sim, group, current_outcome);
            //std::cout << sim.Gaus(sim.MPA_idx(group)).has_nan()<< "Update" << std::endl;

            
        }

        sim = Update_Gaus(sim, group, 0.5);
        sim = Update_Kernel(sim, group);
        sim.gibbs_sampling();
        sim.log_likelihood();
        //std::cout << sim.loglikelihood << std::endl;
        
        beta_output.row(i) = sim.beta.t();
        pgamma_output.row(i) = sim.pgamma.t();
        gamma_output.row(i) = sim.gamma.t();
        alpha_output.row(i) = sim.alpha.t();
        eta_output.row(i) = sim.eta.t();
        tau_output.row(i) = sim.tau.t();
        omega_output.row(i) = sim.omega.t();
        theta_output.row(i) = sim.theta;
        ptheta_output.row(i) = sim.ptheta;
        //bar.set_progress((double(i)/double(n_iter)) * 100);
        if((i % 10000 == 0) & (i != 0))
        {
            Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/beta.csv", beta_output);
            Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/pgamma.csv", pgamma_output);
            Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/gamma.csv", gamma_output);
            Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/alpha.csv", alpha_output);
            Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/eta.csv", eta_output);
            Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/tau.csv", tau_output);
            Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/omega.csv", omega_output);
            Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/theta.csv", theta_output);
            Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/ptheta.csv", ptheta_output);
            Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/TimeTable.csv", 
                        arma::join_horiz(   sim.Time, 
                                            arma::conv_to<arma::colvec>::from(sim.Outcome), 
                                            arma::conv_to<arma::colvec>::from(sim.Group),
                                            arma::conv_to<arma::colvec>::from(sim.Gaus)));
        }


        bar.tick();
    }
        
    sim.Time.brief_print("Time");
    sim.Outcome.brief_print("Outcome");
    sim.Group.brief_print("Group");
    sim.Gaus.brief_print("Gaus");
    std::cout << "log Likelihood" << sim.loglikelihood << std::endl;

    Data_Output("/home/yao/C++/Covid_Research/NPBCOVID/data/beta.csv", beta_output);
    Data_Output("/home/yao/C++/Covid_Research/NPBCOVID/data/pgamma.csv", pgamma_output);
    Data_Output("/home/yao/C++/Covid_Research/NPBCOVID/data/gamma.csv", gamma_output);
    Data_Output("/home/yao/C++/Covid_Research/NPBCOVID/data/alpha.csv", alpha_output);
    Data_Output("/home/yao/C++/Covid_Research/NPBCOVID/data/eta.csv", eta_output);
    Data_Output("/home/yao/C++/Covid_Research/NPBCOVID/data/tau.csv", tau_output);
    Data_Output("/home/yao/C++/Covid_Research/NPBCOVID/data/omega.csv", omega_output);
    Data_Output("/home/yao/C++/Covid_Research/NPBCOVID/data/theta.csv", theta_output);
    Data_Output("/home/yao/C++/Covid_Research/NPBCOVID/data/ptheta.csv", ptheta_output);
    Data_Output("/home/yao/C++/Covid_Research/NPBCOVID/data/TimeTable.csv", 
                arma::join_horiz(   sim.Time, 
                                    arma::conv_to<arma::colvec>::from(sim.Outcome), 
                                    arma::conv_to<arma::colvec>::from(sim.Group),
                                    arma::conv_to<arma::colvec>::from(sim.Gaus)));
    
    Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/beta.csv", beta_output);
    Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/pgamma.csv", pgamma_output);
    Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/gamma.csv", gamma_output);
    Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/alpha.csv", alpha_output);
    Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/eta.csv", eta_output);
    Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/tau.csv", tau_output);
    Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/omega.csv", omega_output);
    Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/theta.csv", theta_output);
    Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/ptheta.csv", ptheta_output);
    Data_Output("/mnt/c/Users/27526/My Documents/YaoGongzheng/Paper/Paper_Code2/Data/TimeTable.csv", 
                arma::join_horiz(   sim.Time, 
                                    arma::conv_to<arma::colvec>::from(sim.Outcome), 
                                    arma::conv_to<arma::colvec>::from(sim.Group),
                                    arma::conv_to<arma::colvec>::from(sim.Gaus)));

    bar.mark_as_completed();
    return 0;
}

