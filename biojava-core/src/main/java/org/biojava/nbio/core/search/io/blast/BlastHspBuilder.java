/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.core.search.io.blast;

import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * Designed by Paolo Pavan.
 * You may want to find my contacts on Github and LinkedIn for code info 
 * or discuss major changes.
 * https://github.com/paolopavan
 * 
 * @author Paolo Pavan
 */

public class BlastHspBuilder<S extends Sequence<C>, C extends Compound> {
    private int hspNum;
    private double hspBitScore;
    private int hspScore;
    private double hspEvalue;
    private int hspQueryFrom;
    private int hspQueryTo;
    private int hspHitFrom;
    private int hspHitTo;
    private int hspQueryFrame;
    private int hspHitFrame;
    private int hspIdentity;
    private int hspPositive;
    private int hspGaps;
    private int hspAlignLen;
    private String hspQseq;
    private String hspHseq;
    private String hspIdentityString;
    private Double percentageIdentity;
    private Integer mismatchCount;

    public BlastHspBuilder() {
    }

    public BlastHspBuilder<S,C> setHspNum(int hspNum) {
        this.hspNum = hspNum;
        return this;
    }

    public BlastHspBuilder<S,C> setHspBitScore(double hspBitScore) {
        this.hspBitScore = hspBitScore;
        return this;
    }

    public BlastHspBuilder<S,C> setHspScore(int hspScore) {
        this.hspScore = hspScore;
        return this;
    }

    public BlastHspBuilder<S,C> setHspEvalue(double hspEvalue) {
        this.hspEvalue = hspEvalue;
        return this;
    }

    public BlastHspBuilder<S,C> setHspQueryFrom(int hspQueryFrom) {
        this.hspQueryFrom = hspQueryFrom;
        return this;
    }

    public BlastHspBuilder<S,C> setHspQueryTo(int hspQueryTo) {
        this.hspQueryTo = hspQueryTo;
        return this;
    }

    public BlastHspBuilder<S,C> setHspHitFrom(int hspHitFrom) {
        this.hspHitFrom = hspHitFrom;
        return this;
    }

    public BlastHspBuilder<S,C> setHspHitTo(int hspHitTo) {
        this.hspHitTo = hspHitTo;
        return this;
    }

    public BlastHspBuilder<S,C> setHspQueryFrame(int hspQueryFrame) {
        this.hspQueryFrame = hspQueryFrame;
        return this;
    }

    public BlastHspBuilder<S,C> setHspHitFrame(int hspHitFrame) {
        this.hspHitFrame = hspHitFrame;
        return this;
    }

    public BlastHspBuilder<S,C> setHspIdentity(int hspIdentity) {
        this.hspIdentity = hspIdentity;
        return this;
    }

    public BlastHspBuilder<S,C> setHspPositive(int hspPositive) {
        this.hspPositive = hspPositive;
        return this;
    }

    public BlastHspBuilder<S,C> setHspGaps(int hspGaps) {
        this.hspGaps = hspGaps;
        return this;
    }

    public BlastHspBuilder<S,C> setHspAlignLen(int hspAlignLen) {
        this.hspAlignLen = hspAlignLen;
        return this;
    }

    public BlastHspBuilder<S,C> setHspQseq(String hspQseq) {
        this.hspQseq = hspQseq;
        return this;
    }

    public BlastHspBuilder<S,C> setHspHseq(String hspHseq) {
        this.hspHseq = hspHseq;
        return this;
    }

    public BlastHspBuilder<S,C> setHspIdentityString(String hspIdentityString) {
        this.hspIdentityString = hspIdentityString;
        return this;
    }

    public BlastHspBuilder<S,C> setPercentageIdentity(Double percentageIdentity) {
        this.percentageIdentity = percentageIdentity;
        return this;
    }

    public BlastHspBuilder<S,C> setMismatchCount(Integer mismatchCount) {
        this.mismatchCount = mismatchCount;
        return this;
    }
    
    public BlastHsp<S,C> createBlastHsp() {
        return new BlastHsp<S,C>(hspNum, hspBitScore, hspScore, hspEvalue, hspQueryFrom, hspQueryTo, hspHitFrom, hspHitTo, hspQueryFrame, hspHitFrame, hspIdentity, hspPositive, hspGaps, hspAlignLen, hspQseq, hspHseq, hspIdentityString, percentageIdentity, mismatchCount);
    }
    
}
