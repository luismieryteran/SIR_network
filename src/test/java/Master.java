import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class Master {

    public static void main (String args[]) {

        // Defining parameters
        Parameters parameters = new Parameters();

        // IC
        IC ic = new IC(parameters.N);

        Path currentRelativePath = Paths.get("");

        System.out.println(currentRelativePath.toAbsolutePath().toString());

        // Defining network
        NetworkBuilder networkBuilder = new NetworkBuilder();
        try {
            networkBuilder.buildNetwork(parameters.N, parameters.p);
        } catch (IOException e) {
            e.printStackTrace();
        }
        Map<Integer, List<Integer>> networkNeighbors = networkBuilder.networkNeighbors;

        // Simulation
        double t0 = 0.0;   // initial time
        NavigableMap<Double, Map<Integer, Object>> dynamicState = new TreeMap<>();   // state of network in time
        dynamicState.put(t0, ic.state);     // setting ic
        NavigableMap<Double, ReactionSpecification> reactionHistory = new TreeMap<>();   // state of network in time

        Simulation simulation = new Simulation();
        simulation.reactionStepping(dynamicState, reactionHistory, networkNeighbors, parameters);
        simulation.printDynamicState(dynamicState);
        simulation.printSummarizedDynamicState(dynamicState);
        simulation.printReactionHistory(reactionHistory);

    }
}
