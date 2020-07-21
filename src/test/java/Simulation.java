import javafx.util.Pair;
import org.apache.commons.math3.distribution.GammaDistribution;

import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;


class Simulation {
    private Integer randomSeed = 1234567890;
    private Random random;
    private GammaDistribution infectiousPeriodDistribution;

    private NetworkState initialState;

    private NavigableMap<Double, ReactionSpecification> reactionsToCome = new TreeMap<>();   // list of potential reactions to come
    private Map<Integer, Double> recoveryTimesByNode = new HashMap<>();   // Future recovery Times by node
    private Map<Integer, Double> infectionTimesByNode = new HashMap<>();   // Future infection Times by node

    Simulation(NetworkState initialState, Parameters parameters){
        this.initialState = initialState;

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
        infectionTimesByNode.clear();
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
        for (Integer targetNode : networkNeighbors.get(sourceNode)){
            if ( networkState.getState().get(targetNode) == Compartments.S ) {
                Double potentialInfectionTime = currentTime + transmissionTime(parameters.beta);

                // If transmission happens before recovery of source, proceed
                if ( potentialInfectionTime <= recoveryTimesByNode.get(sourceNode) ){

                    // Check if target node had already been marked for transmission
                    if ( !infectionTimesByNode.containsKey(targetNode) ) {
                        infectionTimesByNode.put(targetNode, potentialInfectionTime);
                        reactionsToCome.put(potentialInfectionTime,
                                new ReactionSpecification(ReactionType.Infection, Arrays.asList(sourceNode, targetNode)));
                    } else {
                        if ( potentialInfectionTime < infectionTimesByNode.get(targetNode) ){
                            // Replacing infection of target node with new infection time
                            reactionsToCome.remove(infectionTimesByNode.get(targetNode));

                            infectionTimesByNode.put(targetNode, potentialInfectionTime);
                            reactionsToCome.put(potentialInfectionTime,
                                    new ReactionSpecification(ReactionType.Infection, Arrays.asList(sourceNode, targetNode)));
                        }
                    }
                }
            }
        }
    }

    private Pair<Double, Map<Compartments, Long>> summarizedNetworkState(NetworkState networkState){
        Map<Compartments, Long> numberOfNodesByCompartment = networkState.getState()
                .values()
                .stream()
                .collect(Collectors.groupingBy(
                        Function.identity(),
                        Collectors.counting()));

        Arrays.stream(Compartments.values())
                .forEach(c -> numberOfNodesByCompartment.putIfAbsent(c, (long) 0));

        Pair<Double, Map<Compartments, Long>> summarizedNetworkState =
                new Pair<>(networkState.getTime(), numberOfNodesByCompartment);

        return summarizedNetworkState;
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


            // Target node is still susceptible by construction of reactionsToCome table
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

        PrintOutput printOutput = new PrintOutput();

        // Printing IC
        printOutput.printNetworkState(experiment, networkState);
        printOutput.printSummarizedNetworkState(experiment, summarizedNetworkState(networkState));

        while ( reactionsToCome.size() > 0 ) {
//        for (int i = 1; i <= 5; i++){
            // Single step
            reactionStep(networkState, reactionHistory, networkNeighbors, parameters);
            
            // Printing new state to files
            if ( parameters.N <= 200 & experiment <= 10){
                printOutput.printNetworkState(experiment, networkState);
            }

            // Counting number of nodes in each compartment
            printOutput.printSummarizedNetworkState(experiment, summarizedNetworkState(networkState));
        }

        // Printing to File
        printOutput.printNetworkStateMinimal(experiment, initialState, reactionHistory);
        printOutput.printReactionHistory(experiment, reactionHistory);

    }
}
