import java.util.*;

class IC {
    Integer initialInfecteds = 20;
    NavigableMap<Integer, Compartments> state = new TreeMap<>();

    IC(Integer N){
        Random random = new Random(1234567890);

        // Making everyone S at the start
        for (int i = 1; i <= N; i++){
            state.put(i, Compartments.S);
        }

        // Introducing some Is
        for (int i = 1; i <= initialInfecteds; i++){
            state.put(1 + random.nextInt(N-1), Compartments.I);
        }

    }
}
