class Parameters {
    Integer N = 1000;                       // Population
    Double p = 1.0;
//    Double p = 5.0 * 1 / (N - 1);      // Prob of link between 2 given nodes in Erdös-Renyi network
                                            // mean nbr of connections: (N-1) * p
    Double R0 = 2.5;                    // Basic Rep Number
    Double recov_scale = 5.0;                 // Recovery gamma dist scale
    Double recov_shape = 1.0;                 // Recovery gamma dist shape
    Double beta = R0 * recov_shape / recov_scale / (p * (N - 1));    // Transmission rate per S contact
}
