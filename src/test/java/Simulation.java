import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

class Simulation {
    Path currentRelativePath = Paths.get("");

    private Random random = new Random(1234567890);

    private Double nextReactionTime(Double lambda){
        return  Math.log(1 - random.nextDouble()) / ( - lambda );
    }

    private Integer nextReaction(Map<Integer, Double> reactionRates, Double sumReactionRates){
        Integer nextReaction = 0;
        Double runningSumW = 0.0;
        Double rnd = random.nextDouble();

        while ( runningSumW < rnd ){
            runningSumW += reactionRates.get(nextReaction + 1) / sumReactionRates;
            nextReaction++;
        }

        return nextReaction;
    }


    private void reactionStep(NavigableMap<Double, Map<Integer, Object>> dynamicState,
                              NavigableMap<Double, ReactionSpecification> reactionHistory,
                              Map<Integer, List<Integer>> networkNeighbors, Parameters parameters){

        // Reactions at current time
        Reactions reactions = new Reactions(dynamicState.lastEntry().getValue(), networkNeighbors, parameters);

        if ( reactions.sumReactionRates != 0.0 ) {
            // When and which reaction occurs
            Double nextReactionTime =
                    dynamicState.lastEntry().getKey() +
                            nextReactionTime(reactions.sumReactionRates);

            Integer nextReaction = nextReaction(reactions.reactionRates, reactions.sumReactionRates);

            // Modifying oldState into newState because of reaction
            Map<Integer, Object> newState = new HashMap<>();
            newState.putAll(dynamicState.lastEntry().getValue());

            if (reactions.reactionType.get(nextReaction) == ReactionType.Recovery) {
                newState.put(reactions.reactionNodes.get(nextReaction).get(0), Compartments.R);

            } else if (reactions.reactionType.get(nextReaction) == ReactionType.Infection) {
                newState.put(reactions.reactionNodes.get(nextReaction).get(1), Compartments.I);
            }

            reactionHistory.put(nextReactionTime,
                    new ReactionSpecification(reactions.reactionType.get(nextReaction),
                            reactions.reactionNodes.get(nextReaction)));
            dynamicState.put(nextReactionTime, newState);
        }
    }


    void reactionStepping(NavigableMap<Double, Map<Integer, Object>> dynamicState,
                          NavigableMap<Double, ReactionSpecification> reactionHistory,
                          Map<Integer, List<Integer>> networkNeighbors, Parameters parameters) {

        Integer dynamicStateSize = 0;
        while ( dynamicStateSize != dynamicState.size() ) {
            dynamicStateSize = dynamicState.size();

            reactionStep(dynamicState, reactionHistory, networkNeighbors, parameters);
        }
    }


    void printDynamicState(NavigableMap<Double, Map<Integer, Object>> dynamicState){

        // Output file
        try (FileWriter writer =
                     new FileWriter(currentRelativePath.toAbsolutePath() +
                             "/src/test/output/dynamicState.csv")) {

            // File headers
            writer.append("t, ");
            Set<Integer> nodes = dynamicState.firstEntry().getValue().keySet();

            for (Iterator<Integer> it = nodes.iterator(); it.hasNext(); ) {
                writer.append(String.valueOf(it.next()));

                if ( it.hasNext() ) {
                    writer.append(", ");
                } else {
                    writer.append("\n");
                }
            }

            // Writing to file
            for (Map.Entry<Double, Map<Integer, Object>> entry : dynamicState.entrySet()) {

                writer.append(String.valueOf(entry.getKey()));
                writer.append(", ");
                for (Iterator<Integer> it = nodes.iterator(); it.hasNext(); ){
                    writer.append(String.valueOf(entry.getValue().get(it.next())));

                    if ( it.hasNext() ) {
                        writer.append(", ");
                    } else {
                        writer.append("\n");
                    }
                }
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    void printSummarizedDynamicState(NavigableMap<Double, Map<Integer, Object>> dynamicState) {
        NavigableMap<Double, Map<Object, Long>> summarizedDynamicState = new TreeMap<>();

        for (Double t : dynamicState.keySet()){
            Map<Object, Long> dynamicStateAtTimeToPrint = new LinkedHashMap<>();
            Map<Object, Long> dynamicStateAtTime =
                    dynamicState.get(t)
                            .values()
                            .stream()
                            .collect(Collectors.groupingBy(
                                    Function.identity(),
                                    Collectors.counting()));

            for (Object state : Compartments.values()){
                if ( dynamicStateAtTime.get(state) == null ){
                    dynamicStateAtTimeToPrint.
                            put(state, (long) 0);
                } else {
                    dynamicStateAtTimeToPrint.
                            put(state, dynamicStateAtTime.get(state));
                }
            }

            summarizedDynamicState.put(t, dynamicStateAtTimeToPrint);
        }


        // Output file
        try (FileWriter writer =
                     new FileWriter(currentRelativePath.toAbsolutePath() +
                             "/src/test/output/summarizedDynamicState.csv")) {

            // File headers
            writer.append("t, ");
            Set<Object> states = summarizedDynamicState.firstEntry().getValue().keySet();

            for (Iterator<Object> it = states.iterator(); it.hasNext(); ) {
                writer.append(String.valueOf(it.next()));

                if ( it.hasNext() ) {
                    writer.append(", ");
                } else {
                    writer.append("\n");
                }
            }

            // Writing to file
            for (Map.Entry<Double, Map<Object, Long>> entry : summarizedDynamicState.entrySet()) {

                writer.append(String.valueOf(entry.getKey()));
                writer.append(", ");
                for (Iterator<Object> it = states.iterator(); it.hasNext(); ){
                    writer.append(String.valueOf(entry.getValue().get(it.next())));

                    if ( it.hasNext() ) {
                        writer.append(", ");
                    } else {
                        writer.append("\n");
                    }
                }
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    void printReactionHistory(Map<Double, ReactionSpecification> reactionHistory){

        // Output file
        try (FileWriter writer =
                     new FileWriter(currentRelativePath.toAbsolutePath() +
                             "/src/test/output/reactionHistory.csv")) {

            // File header
            writer.append("t, ReactionType, ReactionNodes\n");

            // Writing to file
            for (Map.Entry<Double, ReactionSpecification> entry : reactionHistory.entrySet()) {

                writer.append(entry.getKey().toString());
                writer.append(", ");
                writer.append(entry.getValue().reactionType.toString());
                writer.append(", ");
                writer.append(entry.getValue().reactionNodes.toString().replace(",", ""));
                writer.append("\n");

            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
