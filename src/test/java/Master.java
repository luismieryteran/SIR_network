import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class Master {

    public static void main (String args[]) {

        // Defining parameters
        SimParameters simParameters = new SimParameters();
        Parameters parameters = new Parameters();

        // IC
        Double simulationStartTime = 0.0;   // initial time
        IC ic = new IC(simulationStartTime, parameters.N);

        // Defining network
        NetworkBuilder networkBuilder = new NetworkBuilder();
        networkBuilder.buildNetwork(parameters.N, parameters.p);
        try {
            networkBuilder.printNetwork(simParameters.outputPath, networkBuilder.network);
        } catch (IOException e) {
            e.printStackTrace();
        }
        Map<Integer, List<Integer>> networkNeighbors = networkBuilder.networkNeighbors;

        // Simulation setup
        Simulation simulation = new Simulation(ic.initialState, parameters);
        simulation.seedRNG();

        // Opening output files
        PrintOutput.openOutputFiles(simParameters.outputPath, parameters.N);

        for (Iterator<Integer> experiment =
             IntStream.rangeClosed(1, simParameters.numberExperiments).iterator();
             experiment.hasNext(); ){

            // Current experiment
            Integer exp = experiment.next();

            // list of reactions and when they occur
            NavigableMap<Double, ReactionSpecification> reactionHistory = new TreeMap<>();

            // state of network in time
            NetworkState networkState = new NetworkState(ic.initialState.getTime(), ic.initialState.getState());

            // Sim setup
            simulation.simulationSetUp(networkState, networkNeighbors, parameters);

            // Simulation
            simulation.reactionStepping(exp, networkState,
                    reactionHistory, networkNeighbors, parameters);

        }

        PrintOutput.closeOutputFiles();

    }
}
