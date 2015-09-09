package org.biojava.nbio.core.search.io.blast;


import org.biojava.nbio.core.search.io.Hsp;
import java.util.List;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * Designed by Paolo Pavan.
 * You may want to find my contacts on Github and LinkedIn for code info 
 * or discuss major changes.
 * https://github.com/paolopavan
 * 
 * @author Paolo Pavan
 */
public class BlastHitBuilder {
    private int hitNum;
    private String hitId;
    private String hitDef;
    private String hitAccession;
    private int hitLen;
    private Sequence hitSequence;
    private List<Hsp> hsps;

    public BlastHitBuilder() {
    }

    public BlastHitBuilder setHitNum(int hitNum) {
        this.hitNum = hitNum;
        return this;
    }

    public BlastHitBuilder setHitId(String hitId) {
        this.hitId = hitId;
        return this;
    }

    public BlastHitBuilder setHitDef(String hitDef) {
        this.hitDef = hitDef;
        return this;
    }

    public BlastHitBuilder setHitAccession(String hitAccession) {
        this.hitAccession = hitAccession;
        return this;
    }

    public BlastHitBuilder setHitLen(int hitLen) {
        this.hitLen = hitLen;
        return this;
    }
    
    public BlastHitBuilder setHitSequence(Sequence s) {
        this.hitSequence = s;
        return this;
    }

    public BlastHitBuilder setHsps(List<Hsp> hsps) {
        this.hsps = hsps;
        return this;
    }

    public BlastHit createBlastHit() {
        return new BlastHit(hitNum, hitId, hitDef, hitAccession, hitLen, hsps, hitSequence);
    }
    
}
