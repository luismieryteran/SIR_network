import org.apache.commons.lang3.tuple.MutablePair;
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

    void simulationSetUp(MutablePair<Double, SortedMap<Integer, Compartments>> initialState,
                         Map<Integer, List<Integer>> networkNeighbors,
                         Parameters parameters){
        Map<Integer, Compartments> ic = initialState.getValue();

        recoveryTimesByNode.clear();
        reactionsToCome.clear();

        // Adding recovery and transmission times to initially infected nodes
        for ( Integer node : ic.keySet() ){
            if ( ic.get(node) == Compartments.I ){
                assignRecoveryTimeToInfectedNode(node, initialState);

                assignTransmissionTimesToSourceNode(node, initialState, networkNeighbors, parameters);
            }
        }
    }

    private Double transmissionTime(Double lambda){
        if ( lambda == 0 ) {
            return Double.POSITIVE_INFINITY;
        } else {
            return Math.log(1 - random.nextDouble()) / ( - lambda );
        }
    }

    private Double recoveryTime() {
        return infectiousPeriodDistribution.sample();
    }

    private void assignRecoveryTimeToInfectedNode(Integer node, MutablePair<Double, SortedMap<Integer, Compartments>> dynamicState){

        Double currentTime = dynamicState.getKey();

        Double recoveryTime = currentTime + recoveryTime();
        recoveryTimesByNode.put(node, recoveryTime);

        reactionsToCome.put(recoveryTime,
                new ReactionSpecification(ReactionType.Recovery, Arrays.asList(node)));
    }

    private void assignTransmissionTimesToSourceNode(Integer sourceNode,
                                                     MutablePair<Double, SortedMap<Integer, Compartments>> dynamicState,
                                                     Map<Integer, List<Integer>> networkNeighbors, Parameters parameters){

        Double currentTime = dynamicState.getKey();

        // Cycling over contacts of node to consider transmission
        for (Integer contactNode : networkNeighbors.get(sourceNode)){
            if ( dynamicState.getValue().get(contactNode) == Compartments.S ) {
                Double potentialInfectionTime = currentTime + transmissionTime(parameters.beta);

                // If transmission happens before recovery of source, proceed
                if ( potentialInfectionTime <= recoveryTimesByNode.get(sourceNode) ){
                    reactionsToCome.put(potentialInfectionTime,
                            new ReactionSpecification(ReactionType.Infection, Arrays.asList(sourceNode, contactNode)));
                }
            }
        }
    }

    private Map<Compartments, Long> numberOfNodesByCompartment(MutablePair<Double, SortedMap<Integer, Compartments>> dynamicState){
      return dynamicState.getValue()
              .values()
              .stream()
              .collect(Collectors.groupingBy(
                      Function.identity(),
                      Collectors.counting()));
    }

    private void reactionStep(MutablePair<Double, SortedMap<Integer, Compartments>> dynamicState,
                              NavigableMap<Double, ReactionSpecification> reactionHistory,
                              Map<Integer, List<Integer>> networkNeighbors, Parameters parameters){

        // next reaction time
        Double reactionTime = reactionsToCome.firstKey();

        // Modifying oldState into newState because of reaction
        SortedMap<Integer, Compartments> newState = new TreeMap<>();
        newState.putAll(dynamicState.getValue());

        if ( reactionsToCome.firstEntry().getValue().reactionType == ReactionType.Infection ){

            Integer sourceNode = reactionsToCome.firstEntry().getValue().reactionNodes.get(0);
            Integer targetNode = reactionsToCome.firstEntry().getValue().reactionNodes.get(1);
            reactionsToCome.remove(reactionTime);

            // Ensuring that target node is still susceptible
            if ( newState.get(targetNode) == Compartments.S ){
                newState.put(targetNode, Compartments.I);

                // Updating state
                dynamicState.setLeft(reactionTime);
                dynamicState.setRight(newState);

                // Updating reaction history
                reactionHistory.put(reactionTime,
                        new ReactionSpecification(ReactionType.Infection, Arrays.asList(sourceNode, targetNode)));

                // Assigning recovery and transmission times for newly infected node
                assignRecoveryTimeToInfectedNode(targetNode, dynamicState);
                assignTransmissionTimesToSourceNode(targetNode, dynamicState, networkNeighbors, parameters);
            }

        } else if ( reactionsToCome.firstEntry().getValue().reactionType == ReactionType.Recovery ) {
            Integer recoveringNode = reactionsToCome.firstEntry().getValue().reactionNodes.get(0);
            reactionsToCome.remove(reactionTime);

            newState.put(recoveringNode, Compartments.R);

            // Updating state
            dynamicState.setLeft(reactionTime);
            dynamicState.setRight(newState);

            // Updating reaction history
            reactionHistory.put(reactionTime,
                    new ReactionSpecification(ReactionType.Recovery, Arrays.asList(recoveringNode)));
        }
    }

    void reactionStepping(Integer experiment, MutablePair<Double, SortedMap<Integer, Compartments>> dynamicState,
                          NavigableMap<Double, ReactionSpecification> reactionHistory,
                          Map<Integer, List<Integer>> networkNeighbors, Parameters parameters) {

        // Printing IC
        printDynamicState(experiment, dynamicState);
        printSummarizedDynamicState(experiment, dynamicState);

        while ( reactionsToCome.size() > 0 ) {
//        for (int i = 1; i <= 10; i++){
            // Single step
            reactionStep(dynamicState, reactionHistory, networkNeighbors, parameters);

            // Printing new state to files
            if ( parameters.N <= 200 & experiment <= 10){
                printDynamicState(experiment, dynamicState);
            }
            printSummarizedDynamicState(experiment, dynamicState);
        }

    }


    private FileWriter writerSummarizedDynamicState;
    private FileWriter writerDynamicState;
    private FileWriter writerDynamicStateMinimal;
    private FileWriter writerReactionHistory;
    void openOutputFiles(String outputPath, Integer networkSize) {
        try {
            //------------- Summarized Dynamic State
            writerSummarizedDynamicState =
                    new FileWriter(outputPath + "summarizedDynamicState.csv");

            // File headers
            writerSummarizedDynamicState.append("iter, t, ");
            for (Compartments compartments : Compartments.values()) {
                writerSummarizedDynamicState.append(String.valueOf(compartments));

                if ( compartments == Compartments.values()[Compartments.values().length-1] ) {
                    writerSummarizedDynamicState.append("\n");
                } else {
                    writerSummarizedDynamicState.append(", ");
                }
            }

            //------------- Detailed Dynamic State
            writerDynamicState =
                    new FileWriter(outputPath + "dynamicState.csv");

            // File headers
            writerDynamicState.append("iter, t, ");
            for (Integer node = 1; node <= networkSize; node++) {
                writerDynamicState.append(node.toString());

                if ( node.equals(networkSize) ) {
                    writerDynamicState.append("\n");
                } else {
                    writerDynamicState.append(", ");
                }
            }

            //------------- Detailed Dynamic State Minimal
            writerDynamicStateMinimal =
                    new FileWriter(outputPath + "dynamicStateMinimal.csv");

            // File headers
            writerDynamicStateMinimal.append("iter, t, node, compartment\n");

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
            writerSummarizedDynamicState.close();
            writerDynamicState.close();
            writerDynamicStateMinimal.close();
            writerReactionHistory.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void printDynamicState(Integer experiment, MutablePair<Double, SortedMap<Integer, Compartments>> dynamicState){
        // Output file
        try  {
            // Writing to file
            writerDynamicState.append(experiment.toString());
            writerDynamicState.append(", ");
            writerDynamicState.append(String.valueOf(dynamicState.getKey()));
            writerDynamicState.append(", ");

            for (Iterator<Integer> it = dynamicState.getValue().keySet().iterator(); it.hasNext(); ) {
                writerDynamicState.append(String.valueOf(dynamicState.getValue().get(it.next())));

                if ( it.hasNext() ) {
                    writerDynamicState.append(", ");
                } else {
                    writerDynamicState.append("\n");
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    private void printSummarizedDynamicState(Integer experiment,
                                             MutablePair<Double, SortedMap<Integer, Compartments>> dynamicState) {

        // Counting number of nodes in each compartment
        Map<Compartments, Long> numbersByCompartment = numberOfNodesByCompartment(dynamicState);
        Arrays.stream(Compartments.values())
                .forEach(c -> numbersByCompartment.putIfAbsent(c, (long) 0));

        // Output file
        try {
            // Writing to file
            writerSummarizedDynamicState.append(experiment.toString());
            writerSummarizedDynamicState.append(", ");
            writerSummarizedDynamicState.append(String.valueOf(dynamicState.getKey()));
            writerSummarizedDynamicState.append(", ");

            for ( Compartments compartments : Compartments.values() ){
                    writerSummarizedDynamicState.append(String.valueOf(numbersByCompartment.get(compartments)));

                    if ( compartments == Compartments.values()[Compartments.values().length-1] ) {
                        writerSummarizedDynamicState.append("\n");
                    } else {
                        writerSummarizedDynamicState.append(", ");
                    }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    void printDynamicStateMinimal(Integer experiment,
                                  SortedMap<Integer, Compartments> initialState,
                                  NavigableMap<Double, ReactionSpecification> reactionHistory){
        // Output file
        try {
            // Initial state
            for (Map.Entry<Integer, Compartments> nodeState : initialState.entrySet()) {
                writerDynamicStateMinimal.append(experiment.toString());
                writerDynamicStateMinimal.append(", ");
                writerDynamicStateMinimal.append(simulationStartTime.toString());
                writerDynamicStateMinimal.append(", ");
                writerDynamicStateMinimal.append(nodeState.getKey().toString());
                writerDynamicStateMinimal.append(", ");
                writerDynamicStateMinimal.append(nodeState.getValue().toString());
                writerDynamicStateMinimal.append("\n");
            }

            // Writing to file
            for (Map.Entry<Double, ReactionSpecification> entry : reactionHistory.entrySet()) {
                writerDynamicStateMinimal.append(experiment.toString());
                writerDynamicStateMinimal.append(", ");
                writerDynamicStateMinimal.append(entry.getKey().toString());
                writerDynamicStateMinimal.append(", ");
                if (entry.getValue().reactionType == ReactionType.Infection) {
                    writerDynamicStateMinimal.append(entry.getValue().reactionNodes.get(1).toString());
                    writerDynamicStateMinimal.append(", ");
                    writerDynamicStateMinimal.append(Compartments.I.toString());
                } else if (entry.getValue().reactionType == ReactionType.Recovery) {
                    writerDynamicStateMinimal.append(entry.getValue().reactionNodes.get(0).toString());
                    writerDynamicStateMinimal.append(", ");
                    writerDynamicStateMinimal.append(Compartments.R.toString());
                }
                writerDynamicStateMinimal.append("\n");
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
