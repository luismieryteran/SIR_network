import org.apache.commons.lang3.tuple.MutablePair;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;


class Simulation {
    private Integer randomSeed = 1234567890;
    private Random random = new Random(randomSeed);
    private Double simulationStartTime;

    Simulation(Double time){
        this.simulationStartTime = time;
    }

    private Map<Compartments, Long> numberOfNodesByCompartment(MutablePair<Double, SortedMap<Integer, Compartments>> dynamicState){
      return dynamicState.getValue()
              .values()
              .stream()
              .collect(Collectors.groupingBy(
                      Function.identity(),
                      Collectors.counting()));
    }

    private Double nextReactionTime(Double lambda){
        if ( lambda == 0 ) {
            return Double.POSITIVE_INFINITY;
        } else {
            return Math.log(1 - random.nextDouble()) / ( - lambda );
        }
    }

    private Integer nextReaction(Map<Integer, Double> reactionRates, Double sumReactionRates){

        Integer nextReaction;
        if ( sumReactionRates != 0 ) {
            nextReaction = 0;
            Double runningSumW = 0.0;
            Double rnd = random.nextDouble();

            while (runningSumW < rnd) {
                runningSumW += reactionRates.get(nextReaction + 1) / sumReactionRates;
                nextReaction++;
            }
        } else {
            nextReaction = null;
        }

        return nextReaction;
    }


    private void reactionStep(MutablePair<Double, SortedMap<Integer, Compartments>> dynamicState,
                              NavigableMap<Double, ReactionSpecification> reactionHistory,
                              Map<Integer, List<Integer>> networkNeighbors, Parameters parameters){

        // Reactions at current time
        Reactions reactions = new Reactions(dynamicState.getValue(), networkNeighbors, parameters);

        if ( reactions.sumReactionRates != 0.0 ) {
            // When and which reaction occurs
            Double nextReactionTime =
                    dynamicState.getKey() +
                            nextReactionTime(reactions.sumReactionRates);

            Integer nextReaction = nextReaction(reactions.reactionRates, reactions.sumReactionRates);

            // Modifying oldState into newState because of reaction
            SortedMap<Integer, Compartments> newState = new TreeMap<>();
            newState.putAll(dynamicState.getValue());

            if (reactions.reactionType.get(nextReaction) == ReactionType.Recovery) {
                newState.put(reactions.reactionNodes.get(nextReaction).get(0), Compartments.R);

            } else if (reactions.reactionType.get(nextReaction) == ReactionType.Infection) {
                newState.put(reactions.reactionNodes.get(nextReaction).get(1), Compartments.I);
            }

            reactionHistory.put(nextReactionTime,
                    new ReactionSpecification(reactions.reactionType.get(nextReaction),
                            reactions.reactionNodes.get(nextReaction)));
            dynamicState.setLeft(nextReactionTime);
            dynamicState.setRight(newState);
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


    private void printSummarizedDynamicState(Integer experiment, MutablePair<Double, SortedMap<Integer, Compartments>> dynamicState) {
        Map<Compartments, Long> numbersByCompartment =
                dynamicState.getValue()
                        .values()
                        .stream()
                        .collect(Collectors.groupingBy(
                                Function.identity(),
                                Collectors.counting()));

        for (Compartments state : Compartments.values()){
            numbersByCompartment.putIfAbsent(state, (long) 0);
        }

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
