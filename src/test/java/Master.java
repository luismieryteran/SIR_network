import org.apache.commons.lang3.tuple.MutablePair;

import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class Master {

    public static void main (String args[]) {

        // Defining parameters
        SimParameters simParameters = new SimParameters();
        Parameters parameters = new Parameters();

        // IC
        IC ic = new IC(parameters.N);

        // Defining network
        NetworkBuilder networkBuilder = new NetworkBuilder();
        networkBuilder.buildNetwork(parameters.N, parameters.p);
        try {
            networkBuilder.printNetwork(simParameters.outputPath, networkBuilder.network);
        } catch (IOException e) {
            e.printStackTrace();
        }
        Map<Integer, List<Integer>> networkNeighbors = networkBuilder.networkNeighbors;

        // Simulation
        Double t0 = 0.0;   // initial time
        Simulation simulation = new Simulation(t0, parameters);
        simulation.seedRNG();
        simulation.openOutputFiles(simParameters.outputPath, parameters.N);

        for (Iterator<Integer> experiment =
             IntStream.rangeClosed(1, simParameters.numberExperiments).iterator();
             experiment.hasNext(); ){

            // Current experiment
            Integer exp = experiment.next();

            // list of reactions and when they occur
            NavigableMap<Double, ReactionSpecification> reactionHistory = new TreeMap<>();

            // state of network in time
            NetworkState networkState = new NetworkState(t0, ic.state);

            // Sim setup
            simulation.simulationSetUp(networkState, networkNeighbors, parameters);

            // Simulation
            simulation.reactionStepping(exp, networkState,
                    reactionHistory, networkNeighbors, parameters);

            // Printing to File
            simulation.printnetworkStateMinimal(exp, ic.state, reactionHistory);
            simulation.printReactionHistory(exp, reactionHistory);
        }

        simulation.closeOutputFiles();

    }
}
