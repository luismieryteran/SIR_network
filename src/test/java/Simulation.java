import org.apache.commons.math3.distribution.GammaDistribution;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;


class Simulation {
    private Integer randomSeed = 1234567890;
    private Random random;
    private GammaDistribution infectiousPeriodDistribution;

    private Double simulationStartTime;

    private NavigableMap<Double, ReactionSpecification> reactionsToCome = new TreeMap<>();   // list of potential reactions to come
    private Map<Integer, Double> recoveryTimesByNode = new HashMap<>();   // Future recovery Times by node

    Simulation(Double time, Parameters parameters){
        this.simulationStartTime = time;
        this.infectiousPeriodDistribution = new GammaDistribution(parameters.recov_scale, parameters.recov_shape);
    }

    void seedRNG(){
        this.infectiousPeriodDistribution.reseedRandomGenerator(randomSeed);
        this.random = new Random(randomSeed);
    }

    void simulationSetUp(NetworkState initialState,
                         Map<Integer, List<Integer>> networkNeighbors,
                         Parameters parameters){
        Map<Integer, Compartments> ic = initialState.getState();

        recoveryTimesByNode.clear();
        reactionsToCome.clear();

        // Adding recovery and transmission times to initially infected nodes
        for ( Integer node : ic.keySet() ){
            if ( ic.get(node) == Compartments.I ){
                assignRecoveryTimeToInfectedNode(node, initialState, parameters);

                assignTransmissionTimesToSourceNode(node, initialState, networkNeighbors, parameters);
            }
        }
    }

    private Double transmissionTime(Double rate){
        if ( rate == 0 ) {
            return Double.POSITIVE_INFINITY;
        } else {
            return Math.log(1 - random.nextDouble()) / ( - rate );
        }
    }

    private Double recoveryTime(Double rate) {
        return Math.log(1 - random.nextDouble()) / ( - rate );

//        return infectiousPeriodDistribution.sample();
    }

    private void assignRecoveryTimeToInfectedNode(Integer node, NetworkState networkState, Parameters parameters){

        Double currentTime = networkState.getTime();

        Double recoveryTime = currentTime + recoveryTime(1.0 / parameters.recov_scale);
        recoveryTimesByNode.put(node, recoveryTime);

        reactionsToCome.put(recoveryTime,
                new ReactionSpecification(ReactionType.Recovery, Arrays.asList(node)));
    }

    private void assignTransmissionTimesToSourceNode(Integer sourceNode,
                                                     NetworkState networkState,
                                                     Map<Integer, List<Integer>> networkNeighbors, Parameters parameters){

        Double currentTime = networkState.getTime();

        // Cycling over contacts of node to consider transmission
        for (Integer contactNode : networkNeighbors.get(sourceNode)){
            if ( networkState.getState().get(contactNode) == Compartments.S ) {
                Double potentialInfectionTime = currentTime + transmissionTime(parameters.beta);

                // If transmission happens before recovery of source, proceed
                if ( potentialInfectionTime <= recoveryTimesByNode.get(sourceNode) ){
                    reactionsToCome.put(potentialInfectionTime,
                            new ReactionSpecification(ReactionType.Infection, Arrays.asList(sourceNode, contactNode)));
                }
            }
        }
    }

    private Map<Compartments, Long> numberOfNodesByCompartment(NetworkState networkState){
      return networkState.getState()
              .values()
              .stream()
              .collect(Collectors.groupingBy(
                      Function.identity(),
                      Collectors.counting()));
    }

    private void reactionStep(NetworkState networkState,
                              NavigableMap<Double, ReactionSpecification> reactionHistory,
                              Map<Integer, List<Integer>> networkNeighbors, Parameters parameters){

        // next reaction time
        Double reactionTime = reactionsToCome.firstKey();

        // Modifying oldState into newState because of reaction
        SortedMap<Integer, Compartments> newState = new TreeMap<>();
        newState.putAll(networkState.getState());

        if ( reactionsToCome.firstEntry().getValue().reactionType == ReactionType.Infection ){

            Integer sourceNode = reactionsToCome.firstEntry().getValue().reactionNodes.get(0);
            Integer targetNode = reactionsToCome.firstEntry().getValue().reactionNodes.get(1);
            reactionsToCome.remove(reactionTime);

            // Ensuring that target node is still susceptible
            if ( newState.get(targetNode) == Compartments.S ){
                newState.put(targetNode, Compartments.I);

                // Updating state
                networkState.setTime(reactionTime);
                networkState.setState(newState);

                // Updating reaction history
                reactionHistory.put(reactionTime,
                        new ReactionSpecification(ReactionType.Infection, Arrays.asList(sourceNode, targetNode)));

                // Assigning recovery and transmission times for newly infected node
                assignRecoveryTimeToInfectedNode(targetNode, networkState, parameters);
                assignTransmissionTimesToSourceNode(targetNode, networkState, networkNeighbors, parameters);
            }

        } else if ( reactionsToCome.firstEntry().getValue().reactionType == ReactionType.Recovery ) {
            Integer recoveringNode = reactionsToCome.firstEntry().getValue().reactionNodes.get(0);
            reactionsToCome.remove(reactionTime);

            newState.put(recoveringNode, Compartments.R);

            // Updating state
            networkState.setTime(reactionTime);
            networkState.setState(newState);

            // Updating reaction history
            reactionHistory.put(reactionTime,
                    new ReactionSpecification(ReactionType.Recovery, Arrays.asList(recoveringNode)));
        }
    }

    void reactionStepping(Integer experiment, NetworkState networkState,
                          NavigableMap<Double, ReactionSpecification> reactionHistory,
                          Map<Integer, List<Integer>> networkNeighbors, Parameters parameters) {

        // Printing IC
        printnetworkState(experiment, networkState);
        printSummarizednetworkState(experiment, networkState);

        while ( reactionsToCome.size() > 0 ) {
//        for (int i = 1; i <= 10; i++){
            // Single step
            reactionStep(networkState, reactionHistory, networkNeighbors, parameters);

            // Printing new state to files
            if ( parameters.N <= 200 & experiment <= 10){
                printnetworkState(experiment, networkState);
            }
            printSummarizednetworkState(experiment, networkState);
        }

    }


    private FileWriter writerSummarizednetworkState;
    private FileWriter writernetworkState;
    private FileWriter writernetworkStateMinimal;
    private FileWriter writerReactionHistory;
    void openOutputFiles(String outputPath, Integer networkSize) {
        try {
            //------------- Summarized Dynamic State
            writerSummarizednetworkState =
                    new FileWriter(outputPath + "summarizedDynamicState.csv");

            // File headers
            writerSummarizednetworkState.append("iter, t, ");
            for (Compartments compartments : Compartments.values()) {
                writerSummarizednetworkState.append(String.valueOf(compartments));

                if ( compartments == Compartments.values()[Compartments.values().length-1] ) {
                    writerSummarizednetworkState.append("\n");
                } else {
                    writerSummarizednetworkState.append(", ");
                }
            }

            //------------- Detailed Dynamic State
            writernetworkState =
                    new FileWriter(outputPath + "dynamicState.csv");

            // File headers
            writernetworkState.append("iter, t, ");
            for (Integer node = 1; node <= networkSize; node++) {
                writernetworkState.append(node.toString());

                if ( node.equals(networkSize) ) {
                    writernetworkState.append("\n");
                } else {
                    writernetworkState.append(", ");
                }
            }

            //------------- Detailed Dynamic State Minimal
            writernetworkStateMinimal =
                    new FileWriter(outputPath + "dynamicStateMinimal.csv");

            // File headers
            writernetworkStateMinimal.append("iter, t, node, compartment\n");

            //------------- Reaction History
            writerReactionHistory =
                    new FileWriter(outputPath + "reactionHistory.csv");
            // File header
            writerReactionHistory.append("iter, t, ReactionType, ReactionNodes\n");

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    void closeOutputFiles(){
        try {
            writerSummarizednetworkState.close();
            writernetworkState.close();
            writernetworkStateMinimal.close();
            writerReactionHistory.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void printnetworkState(Integer experiment, NetworkState networkState){
        // Output file
        try  {
            // Writing to file
            writernetworkState.append(experiment.toString());
            writernetworkState.append(", ");
            writernetworkState.append(String.valueOf(networkState.getTime()));
            writernetworkState.append(", ");

            for (Iterator<Integer> it = networkState.getState().keySet().iterator(); it.hasNext(); ) {
                writernetworkState.append(String.valueOf(networkState.getState().get(it.next())));

                if ( it.hasNext() ) {
                    writernetworkState.append(", ");
                } else {
                    writernetworkState.append("\n");
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    private void printSummarizednetworkState(Integer experiment,
                                             NetworkState networkState) {

        // Counting number of nodes in each compartment
        Map<Compartments, Long> numbersByCompartment = numberOfNodesByCompartment(networkState);
        Arrays.stream(Compartments.values())
                .forEach(c -> numbersByCompartment.putIfAbsent(c, (long) 0));

        // Output file
        try {
            // Writing to file
            writerSummarizednetworkState.append(experiment.toString());
            writerSummarizednetworkState.append(", ");
            writerSummarizednetworkState.append(String.valueOf(networkState.getTime()));
            writerSummarizednetworkState.append(", ");

            for ( Compartments compartments : Compartments.values() ){
                    writerSummarizednetworkState.append(String.valueOf(numbersByCompartment.get(compartments)));

                    if ( compartments == Compartments.values()[Compartments.values().length-1] ) {
                        writerSummarizednetworkState.append("\n");
                    } else {
                        writerSummarizednetworkState.append(", ");
                    }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    void printnetworkStateMinimal(Integer experiment,
                                  NavigableMap<Integer, Compartments> initialState,
                                  NavigableMap<Double, ReactionSpecification> reactionHistory){
        // Output file
        try {
            // Initial state
            for (Map.Entry<Integer, Compartments> nodeState : initialState.entrySet()) {
                writernetworkStateMinimal.append(experiment.toString());
                writernetworkStateMinimal.append(", ");
                writernetworkStateMinimal.append(simulationStartTime.toString());
                writernetworkStateMinimal.append(", ");
                writernetworkStateMinimal.append(nodeState.getKey().toString());
                writernetworkStateMinimal.append(", ");
                writernetworkStateMinimal.append(nodeState.getValue().toString());
                writernetworkStateMinimal.append("\n");
            }

            // Writing to file
            for (Map.Entry<Double, ReactionSpecification> entry : reactionHistory.entrySet()) {
                writernetworkStateMinimal.append(experiment.toString());
                writernetworkStateMinimal.append(", ");
                writernetworkStateMinimal.append(entry.getKey().toString());
                writernetworkStateMinimal.append(", ");
                if (entry.getValue().reactionType == ReactionType.Infection) {
                    writernetworkStateMinimal.append(entry.getValue().reactionNodes.get(1).toString());
                    writernetworkStateMinimal.append(", ");
                    writernetworkStateMinimal.append(Compartments.I.toString());
                } else if (entry.getValue().reactionType == ReactionType.Recovery) {
                    writernetworkStateMinimal.append(entry.getValue().reactionNodes.get(0).toString());
                    writernetworkStateMinimal.append(", ");
                    writernetworkStateMinimal.append(Compartments.R.toString());
                }
                writernetworkStateMinimal.append("\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    void printReactionHistory(Integer experiment, Map<Double, ReactionSpecification> reactionHistory){
        // Output file
        try {
            // Writing to file
            for (Map.Entry<Double, ReactionSpecification> entry : reactionHistory.entrySet()) {
                writerReactionHistory.append(experiment.toString());
                writerReactionHistory.append(", ");
                writerReactionHistory.append(entry.getKey().toString());
                writerReactionHistory.append(", ");
                writerReactionHistory.append(entry.getValue().reactionType.toString());
                writerReactionHistory.append(", ");
                writerReactionHistory.append(entry.getValue().reactionNodes.toString().replace(",", ""));
                writerReactionHistory.append("\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
