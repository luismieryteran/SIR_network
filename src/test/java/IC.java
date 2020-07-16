import java.util.*;

class IC {

    SortedMap<Integer, Compartments> state = new TreeMap<>();

    IC(Integer N){
        Random random = new Random(1234567890);

        // Making everyone S at the start
        for (int i = 1; i <= N; i++){
            state.put(i, Compartments.S);
        }

        // Introducing two Is
        state.put(1 + random.nextInt(N-1), Compartments.I);
        state.put(1 + random.nextInt(N-1), Compartments.I);
    }
}
