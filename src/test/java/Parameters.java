class Parameters {
    Integer N = 30;                       // Population
    Double p = 1.0;
//    Double p = 15.0 * 1 / (N - 1);      // Prob of link between 2 given nodes in Erdös-Renyi network
                                            // mean nbr of connections: (N-1) * p
    Double R0 = 2.5;                    // Basic Rep Number
    Double gamma = 0.2;                 // Recovery rate
    Double beta = R0 * gamma / (p * (N - 1));    // Transmission rate per S contact
}
