import javafx.util.Pair;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

class NetworkBuilder {
    List<Pair<Integer, Integer>> network = new ArrayList<>();
    Map<Integer, List<Integer>> networkNeighbors = new HashMap<>();


    void buildNetwork (int N, double p) {
        Random random = new Random(1234567890);

        // Initializing network with no contacts
        for (int i = 1; i <= N; i++){
            List<Integer> connectedNodes = new ArrayList<>();
            networkNeighbors.put(i, connectedNodes);
        }

        // Filling in contacts (Erdös-Renyi network)
        for (int i = 1; i <= N; i++){
            for (int j = i + 1; j <= N; j++){
                if ( random.nextFloat() <= p ) {
                    network.add(new Pair(i, j));

                    List<Integer> connectedNodes1 = networkNeighbors.get(i);
                    connectedNodes1.add(j);

                    // Symmetrizing connections
                    List<Integer> connectedNodes2 = networkNeighbors.get(j);
                    connectedNodes2.add(i);
                }
            }
        }

        // Capturing unconnected nodes and assigning them 1 connection
        for (int i = 1; i <= N; i++) {
            if (networkNeighbors.get(i).size() == 0) {
                Integer j = 1 + random.nextInt(N);

                List<Integer> connectedNodes1 = networkNeighbors.get(i);
                connectedNodes1.add(j);

                // Symmetrizing connections
                List<Integer> connectedNodes2 = networkNeighbors.get(j);
                connectedNodes2.add(i);

                network.add(new Pair(i, j));
            }
        }
        // Printing network to file
//        printNetwork(network);
    }

    void printNetwork(String outputPath, List<Pair<Integer, Integer>> network) throws IOException {

        // Output file
        try (FileWriter writer =
                     new FileWriter(outputPath + "network.csv")) {

            writer.append("Node1, Node2 \n");
            for (Pair linkedNodes : network) {

                writer.append(String.valueOf(linkedNodes.getKey()));
                writer.append(", ");
                writer.append(String.valueOf(linkedNodes.getValue()));
                writer.append("\n");
            }
            writer.close();
        }
    }
}
