import java.util.SortedMap;

public class NetworkState {
    public Double time;
    public SortedMap<Integer, Compartments> state;

    public NetworkState(Double time, SortedMap<Integer, Compartments> state){
        this.time = time;
        this.state = state;
    }




}
