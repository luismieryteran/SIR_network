import java.util.SortedMap;

class NetworkState {
    private Double time;
    private SortedMap<Integer, Compartments> state;

    NetworkState(Double time, SortedMap<Integer, Compartments> state){
        this.time = time;
        this.state = state;
    }

    void setState(SortedMap<Integer, Compartments> state) {
        this.state = state;
    }

    void setTime(Double time) {
        this.time = time;
    }

    Double getTime() {
        return time;
    }

    SortedMap<Integer, Compartments> getState() {
        return state;
    }
}
