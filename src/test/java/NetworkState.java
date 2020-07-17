import java.util.SortedMap;

public class NetworkState {
    private Double time;
    private SortedMap<Integer, Compartments> state;

    public NetworkState(Double time, SortedMap<Integer, Compartments> state){
        this.time = time;
        this.state = state;
    }

    public void setState(SortedMap<Integer, Compartments> state) {
        this.state = state;
    }

    public void setTime(Double time) {
        this.time = time;
    }

    public Double getTime() {
        return time;
    }

    public SortedMap<Integer, Compartments> getState() {
        return state;
    }
}
