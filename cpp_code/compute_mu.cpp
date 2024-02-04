/*
    // NEW CODE -------------------------------------------------------------
    const amrex::Real Tc = 304.128;
    mu = mu_o * 
        ( 
            (A[0] * A[5] + A[1] + A[2]) / 
            (
                A[0] * (1 - std::exp(-A[5] * rholoc * Vcm / 6)) / (rholoc * Vcm / 6)
                + A[1] * (1 - 0.5 * rholoc * Vcm / 6) /  
                ( 
                    (1 - rholoc * Vcm / 6) * (1 - rholoc * Vcm / 6) * (1 - rholoc * Vcm / 6)
                ) * std::exp(A[4] * rholoc * Vcm / 6)
                + A[2] * (1 - 0.5 * rholoc * Vcm / 6) / ((1 - rholoc * Vcm / 6) * (1 - rholoc * Vcm / 6) * (1 - rholoc * Vcm / 6))       
            )           
            + A[5] * rholoc * Vcm / 6        
        )
        +
        3.6344e-5 * std::sqrt(MW_m * Tcm) / std::pow(Vcm, 2.0 / 3.0) 
        * std::exp(A[7] + A[8] / Tstar + A[9] / Tstar / Tstar) 
        * A[6] * (rholoc * Vcm / 6) * (rholoc * Vcm / 6)
        * (   
            A[0] * (1 - std::exp(-A[5] * rholoc * Vcm / 6)) / (rholoc * Vcm / 6)
            + A[1] * (1 - 0.5 * rholoc * Vcm / 6) / 
            ( 
                (1 - rholoc * Vcm / 6) * (1 - rholoc * Vcm / 6) * (1 - rholoc * Vcm / 6)
            ) * std::exp(A[4] * rholoc * Vcm / 6)
            + A[2] * (1 - 0.5 * rholoc * Vcm / 6) / ((1 - rholoc * Vcm / 6) * (1 - rholoc * Vcm / 6) * (1 - rholoc * Vcm / 6))      
        ) / (A[0] * A[5] + A[1] + A[2])
        ;     
    //----------------------------------------------------------------------
*/ 
