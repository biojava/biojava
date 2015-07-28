/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava.nbio.core.search.io.blast;


import org.biojava.nbio.core.search.io.Hsp;
import java.util.ArrayList;
import java.util.List;
import org.biojava.nbio.core.sequence.template.Sequence;


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
