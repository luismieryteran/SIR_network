import java.util.List;

class ReactionSpecification {
    ReactionType reactionType;
    List<Integer> reactionNodes;

    ReactionSpecification(ReactionType oneReactionType, List<Integer> oneReactionNodes){
        this.reactionType = oneReactionType;
        this.reactionNodes = oneReactionNodes;

    }

    public String toString(){
        return "[" + reactionType + " = " + reactionNodes + "]";
    }

}
