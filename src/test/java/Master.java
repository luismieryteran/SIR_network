import javafx.util.Pair;
import org.apache.commons.lang3.tuple.MutablePair;

import java.io.IOException;
import java.util.*;

import static java.util.stream.IntStream.range;

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
        Simulation simulation = new Simulation(t0);
        simulation.openOutputFiles(simParameters.outputPath, parameters.N);


        for (int experiment = 1; experiment <= simParameters.numberExperiments; experiment++) {
            NavigableMap<Double, ReactionSpecification> reactionHistory = new TreeMap<>();   // list of reactions and when they occur
            MutablePair<Double, SortedMap<Integer, Compartments>> dynamicState = new MutablePair(t0, ic.state);   // state of network in time

            // Simulation
            simulation.reactionStepping(experiment, dynamicState, reactionHistory, networkNeighbors, parameters);

            // Printing to File
            simulation.printDynamicStateMinimal(experiment, ic.state, reactionHistory);
            simulation.printReactionHistory(experiment, reactionHistory);
        }
        simulation.closeOutputFiles();

    }
}
