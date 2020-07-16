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
    private GammaDistribution gammaDistribution;

    private Double simulationStartTime;

    private NavigableMap<Double, ReactionSpecification> reactionsToCome = new TreeMap<>();   // list of potential reactions to come
    private Map<Integer, Double> recoveryTimesByNode = new HashMap<>();   // Future recovery Times by node


    Simulation(Double time, Parameters parameters){
        this.simulationStartTime = time;

        this.gammaDistribution  = new GammaDistribution(parameters.recov_scale, parameters.recov_shape);

    }

    void seedRNG(){
        this.gammaDistribution.reseedRandomGenerator(randomSeed);
        this.random = new Random(randomSeed);
    }

    private Double recoveryTime() {
        return gammaDistribution.sample();
    }

    void simulationSetUp(MutablePair<Double, SortedMap<Integer, Compartments>> initialState){
        Map<Integer, Compartments> ic = initialState.getValue();

        recoveryTimesByNode.clear();
        reactionsToCome.clear();

        // Adding recovery times to initially infected nodes
        for (Integer node : ic.keySet()){
            if (ic.get(node) == Compartments.I){
                Double recoveryTime = recoveryTime();
                recoveryTimesByNode.put(node, simulationStartTime + recoveryTime);

                reactionsToCome.put(recoveryTime,
                        new ReactionSpecification(ReactionType.Recovery, Arrays.asList(node)));
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

    private Double nextInfectionTime(Double lambda){
        if ( lambda == 0 ) {
            return Double.POSITIVE_INFINITY;
        } else {
            return Math.log(1 - random.nextDouble()) / ( - lambda );
        }
    }

    private Integer nextInfection(Map<Integer, Double> infectionRates, Double sumInfectionRates){

        Integer nextInfectionReaction;
        if ( sumInfectionRates != 0 ) {
            nextInfectionReaction = 0;
            Double runningSumW = 0.0;
            Double rnd = random.nextDouble();

            while (runningSumW < rnd) {
                runningSumW += infectionRates.get(nextInfectionReaction + 1) / sumInfectionRates;
                nextInfectionReaction++;
            }
        } else {
            nextInfectionReaction = null;
        }

        return nextInfectionReaction;
    }

    private void reactionStep(MutablePair<Double, SortedMap<Integer, Compartments>> dynamicState,
                              NavigableMap<Double, ReactionSpecification> reactionHistory,
                              Map<Integer, List<Integer>> networkNeighbors, Parameters parameters){

        // Reactions at current time
        Reactions reactions = new Reactions(dynamicState.getValue(), networkNeighbors, parameters);

        Double currentTime = dynamicState.getKey();
        Double nextReactionTime;
        Integer nextReaction;

        // Potential next infection
        nextReactionTime = currentTime + nextInfectionTime(reactions.sumReactionRates);
        nextReaction = nextInfection(reactions.reactionRates, reactions.sumReactionRates);

        if ( reactions.sumReactionRates != 0.0 | reactionsToCome.size() > 0) {

            if (nextReactionTime < reactionsToCome.firstEntry().getKey()) {  // New infection happens

                Integer targetNode = reactions.reactionNodes.get(nextReaction).get(1);

                // Modifying oldState into newState because of reaction
                SortedMap<Integer, Compartments> newState = new TreeMap<>();
                newState.putAll(dynamicState.getValue());
                newState.put(targetNode, Compartments.I);

                // Assigning Recovery of just-infected node
                recoveryTimesByNode.put(targetNode, nextReactionTime +
                        recoveryTime());
                reactionsToCome.put(recoveryTimesByNode.get(targetNode),
                        new ReactionSpecification(ReactionType.Recovery, Arrays.asList(targetNode)));


                // Updating reaction to History and state
                reactionHistory.put(nextReactionTime,
                        new ReactionSpecification(reactions.reactionType.get(nextReaction),
                                reactions.reactionNodes.get(nextReaction)));

                // Updating state
                dynamicState.setLeft(nextReactionTime);
                dynamicState.setRight(newState);

            } else {    // One Recovery happens before new infection

                Double recoveryTime = reactionsToCome.firstEntry().getKey();
                Integer recoveringNode = reactionsToCome.firstEntry().getValue().reactionNodes.get(0);
                reactionsToCome.remove(recoveryTime);

                // Modifying oldState into newState because of reaction
                SortedMap<Integer, Compartments> newState = new TreeMap<>();
                newState.putAll(dynamicState.getValue());
                newState.put(recoveringNode, Compartments.R);

                // Adding reaction to History and updating state
                reactionHistory.put(recoveryTime,
                        new ReactionSpecification(ReactionType.Recovery,
                                Arrays.asList(recoveringNode)));

                // Updating state
                dynamicState.setLeft(recoveryTime);
                dynamicState.setRight(newState);

            }
        }
    }

    void reactionStepping(Integer experiment, MutablePair<Double, SortedMap<Integer, Compartments>> dynamicState,
                          NavigableMap<Double, ReactionSpecification> reactionHistory,
                          Map<Integer, List<Integer>> networkNeighbors, Parameters parameters) {

        // Printing IC
        printDynamicState(experiment, dynamicState);
        printSummarizedDynamicState(experiment, dynamicState);

        Map<Compartments, Long> nodeComposition = numberOfNodesByCompartment(dynamicState);

        while ( nodeComposition.get(Compartments.I) != null ) {
//        for (int i = 1; i <= 10; i++){
            // Single step
            reactionStep(dynamicState, reactionHistory, networkNeighbors, parameters);

            // new node composition
            nodeComposition = numberOfNodesByCompartment(dynamicState);

            // Printing new state to files
            printDynamicState(experiment, dynamicState);
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
