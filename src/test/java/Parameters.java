class Parameters {
    int N = 100;                       // Population
    double p = 5.0 / (N - 1);                    // Prob of link between 2 given nodes (Erdös-Renyi network)
    double R0 = 2.5;                    // Basic Rep Number
    double gamma = 0.2;                 // Recovery rate
    double beta = R0 * gamma / (p * (N - 1));    // Transmission rate per S contact
}
