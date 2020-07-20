import java.util.*;

class IC {
    Random random = new Random(1234567890);
    private Integer initialInfecteds = 5;
    NetworkState initialState;

    IC(Double simulationStartTime, Integer N){

        NavigableMap<Integer, Compartments> state = new TreeMap<>();

        // Making everyone S at the start
        for (int i = 1; i <= N; i++){
            state.put(i, Compartments.S);
        }

        // Introducing some Is
        for (int i = 1; i <= initialInfecteds; i++){
            state.put(1 + random.nextInt(N-1), Compartments.I);
        }

        initialState = new NetworkState(simulationStartTime, state);
    }
}
