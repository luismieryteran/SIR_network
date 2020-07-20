import javafx.util.Pair;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
import java.util.NavigableMap;

public class PrintOutput {

    static private FileWriter writerSummarizedNetworkState;
    static private FileWriter writerNetworkState;
    static private FileWriter writerNetworkStateMinimal;
    static private FileWriter writerReactionHistory;
    static void openOutputFiles(String outputPath, Integer networkSize) {
        try {
            //------------- Summarized Dynamic State
            writerSummarizedNetworkState =
                    new FileWriter(outputPath + "summarizedDynamicState.csv");

            // File headers
            writerSummarizedNetworkState.append("iter, t, ");
            for (Compartments compartments : Compartments.values()) {
                writerSummarizedNetworkState.append(String.valueOf(compartments));

                if ( compartments == Compartments.values()[Compartments.values().length-1] ) {
                    writerSummarizedNetworkState.append("\n");
                } else {
                    writerSummarizedNetworkState.append(", ");
                }
            }

            //------------- Detailed Dynamic State
            writerNetworkState =
                    new FileWriter(outputPath + "dynamicState.csv");

            // File headers
            writerNetworkState.append("iter, t, ");
            for (Integer node = 1; node <= networkSize; node++) {
                writerNetworkState.append(node.toString());

                if ( node.equals(networkSize) ) {
                    writerNetworkState.append("\n");
                } else {
                    writerNetworkState.append(", ");
                }
            }

            //------------- Detailed Dynamic State Minimal
            writerNetworkStateMinimal =
                    new FileWriter(outputPath + "dynamicStateMinimal.csv");

            // File headers
            writerNetworkStateMinimal.append("iter, t, node, compartment\n");

            //------------- Reaction History
            writerReactionHistory =
                    new FileWriter(outputPath + "reactionHistory.csv");
            // File header
            writerReactionHistory.append("iter, t, ReactionType, ReactionNodes\n");

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    static void closeOutputFiles(){
        try {
            writerSummarizedNetworkState.close();
            writerNetworkState.close();
            writerNetworkStateMinimal.close();
            writerReactionHistory.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    void printNetworkState(Integer experiment, NetworkState networkState){
        // Output file
        try  {
            // Writing to file
            writerNetworkState.append(experiment.toString());
            writerNetworkState.append(", ");
            writerNetworkState.append(String.valueOf(networkState.getTime()));
            writerNetworkState.append(", ");

            for (Iterator<Integer> it = networkState.getState().keySet().iterator(); it.hasNext(); ) {
                writerNetworkState.append(String.valueOf(networkState.getState().get(it.next())));

                if ( it.hasNext() ) {
                    writerNetworkState.append(", ");
                } else {
                    writerNetworkState.append("\n");
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    void printSummarizedNetworkState(Integer experiment,
                                             Pair<Double, Map<Compartments, Long>> summarizedNetworkState) {

        // Counting number of nodes in each compartment
//        Map<Compartments, Long> numbersByCompartment = numberOfNodesByCompartment(networkState);
//        Arrays.stream(Compartments.values())
//                .forEach(c -> numbersByCompartment.putIfAbsent(c, (long) 0));

        // Output file
        try {
            // Writing to file
            writerSummarizedNetworkState.append(experiment.toString());
            writerSummarizedNetworkState.append(", ");
            writerSummarizedNetworkState.append(String.valueOf(summarizedNetworkState.getKey()));
            writerSummarizedNetworkState.append(", ");

            for ( Compartments compartments : Compartments.values() ){
                writerSummarizedNetworkState.append(String.valueOf(summarizedNetworkState.getValue().get(compartments)));

                if ( compartments == Compartments.values()[Compartments.values().length-1] ) {
                    writerSummarizedNetworkState.append("\n");
                } else {
                    writerSummarizedNetworkState.append(", ");
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    void printNetworkStateMinimal(Integer experiment,
                                  NetworkState initialState,
                                  NavigableMap<Double, ReactionSpecification> reactionHistory){
        // Output file
        try {
            // Initial state
            for (Map.Entry<Integer, Compartments> nodeState : initialState.getState().entrySet()) {
                writerNetworkStateMinimal.append(experiment.toString());
                writerNetworkStateMinimal.append(", ");
                writerNetworkStateMinimal.append(initialState.getTime().toString());
                writerNetworkStateMinimal.append(", ");
                writerNetworkStateMinimal.append(nodeState.getKey().toString());
                writerNetworkStateMinimal.append(", ");
                writerNetworkStateMinimal.append(nodeState.getValue().toString());
                writerNetworkStateMinimal.append("\n");
            }

            // Writing to file
            for (Map.Entry<Double, ReactionSpecification> entry : reactionHistory.entrySet()) {
                writerNetworkStateMinimal.append(experiment.toString());
                writerNetworkStateMinimal.append(", ");
                writerNetworkStateMinimal.append(entry.getKey().toString());
                writerNetworkStateMinimal.append(", ");
                if (entry.getValue().reactionType == ReactionType.Infection) {
                    writerNetworkStateMinimal.append(entry.getValue().reactionNodes.get(1).toString());
                    writerNetworkStateMinimal.append(", ");
                    writerNetworkStateMinimal.append(Compartments.I.toString());
                } else if (entry.getValue().reactionType == ReactionType.Recovery) {
                    writerNetworkStateMinimal.append(entry.getValue().reactionNodes.get(0).toString());
                    writerNetworkStateMinimal.append(", ");
                    writerNetworkStateMinimal.append(Compartments.R.toString());
                }
                writerNetworkStateMinimal.append("\n");
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
