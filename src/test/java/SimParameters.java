import java.nio.file.Path;
import java.nio.file.Paths;

class SimParameters {

    Path currentPath = Paths.get("");
    Integer numberExperiments = 1; // stochastic realizations
    String outputPath = currentPath.toAbsolutePath() + "/src/test/output/";

}
