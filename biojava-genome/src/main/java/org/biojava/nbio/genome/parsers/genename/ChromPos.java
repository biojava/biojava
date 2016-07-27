package org.biojava.nbio.genome.parsers.genename;

/**
 * Created by ap3 on 27/10/2014.
 */
public class ChromPos {

    int pos;
    int phase;

    public int getPhase() {
        return phase;
    }

    public void setPhase(int phase) {
        this.phase = phase;
    }

    public int getPos() {
        return pos;
    }

    public void setPos(int pos) {
        this.pos = pos;
    }

    public ChromPos(int pos, int phase){
        this.pos = pos;
        this.phase = phase;
    }
}