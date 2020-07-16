import java.util.*;

class Reactions {
    Map<Integer, Double> reactionRates = new HashMap<>();
    Map<Integer, ReactionType> reactionType = new HashMap<>();
    Map<Integer, List<Integer>> reactionNodes = new HashMap<>();
    Double sumReactionRates = 0.0;

    Reactions(Map<Integer, Compartments> state, Map<Integer, List<Integer>> networkNeighbors, Parameters parameters){

        for ( Map.Entry<Integer, Compartments> sourceNode : state.entrySet()){
            if ( sourceNode.getValue() == Compartments.I){

                // Cycling over contacts of node to consider transmission
                for (Integer contactNode : networkNeighbors.get(sourceNode.getKey())){
                    if ( state.get(contactNode) == Compartments.S) {
                        reactionRates.put(reactionRates.size() + 1, parameters.beta);
                        sumReactionRates += parameters.beta;

                        reactionType.put(reactionType.size() + 1, ReactionType.Infection);
                        reactionNodes.put(reactionNodes.size() + 1, Arrays.asList(sourceNode.getKey(), contactNode));
                    }
                }
            }
        }
    }
}
