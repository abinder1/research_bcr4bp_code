function sim_config = create_sim_config(sigma, B0, TH0, N_m, N_e, tau_e, options)
    arguments
        sigma double
        B0 double
        TH0 double

        N_m double
        N_e double

        tau_e double
        
        options.INC = deg2rad(5.145);

        options.N = 7

        options.mu_S = 1.32712440018 * 10^11
        options.mu_E = 3.986004418 * 10^5
        options.mu_M = 4.9048695 * 10^3
        
        options.a_S = 149598023
        options.a_M = 384400

        options.use_real_atil = false;
    end
    
    mutil_S = options.mu_S / (options.mu_M + options.mu_E);
    
    atil_S_fake = ((mutil_S + 1) * (N_m / N_e)^2)^(1/3);
    atil_S_real = options.a_S / options.a_M;
    
    sim_config.mu = options.mu_M / (options.mu_M + options.mu_E);
    sim_config.mutil_S = mutil_S;
    sim_config.atil_S_real = atil_S_real;

    if options.use_real_atil == true
        sim_config.atil_S = atil_S_real;
    elseif options.use_real_atil == false
        sim_config.atil_S = atil_S_fake;
    end

    sim_config.tau_e = tau_e;
    sim_config.sigma = sigma;
    sim_config.B0 = B0; 
    sim_config.TH0 = TH0;
    sim_config.INC = options.INC;
    sim_config.N = options.N;
end

