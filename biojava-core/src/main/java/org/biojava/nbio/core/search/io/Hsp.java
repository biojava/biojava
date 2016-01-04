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
package org.biojava.nbio.core.search.io;

import java.util.ArrayList;
import java.util.List;
import org.biojava.nbio.core.alignment.SimpleAlignedSequence;
import org.biojava.nbio.core.alignment.SimpleSequencePair;
import org.biojava.nbio.core.alignment.template.AlignedSequence.Step;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This class models a search Hsp.
 * You will retrieve a list of this using iterator of a Hit
 * 
 * Designed by Paolo Pavan.
 * You may want to find my contacts on Github and LinkedIn for code info 
 * or discuss major changes.
 * https://github.com/paolopavan
 * 
 * @author Paolo Pavan
 */

public abstract class Hsp <S extends Sequence<C>, C extends Compound> {
    private static final Logger logger = LoggerFactory.getLogger(Hsp.class);
    private Integer hspNum;
    private Double hspBitScore;
    private Integer hspScore;
    private Double hspEvalue;
    private Integer hspQueryFrom;
    private Integer hspQueryTo;
    private Integer hspHitFrom;
    private Integer hspHitTo;
    private Integer hspQueryFrame;
    private Integer hspHitFrame;
    private Integer hspIdentity;
    private Integer hspPositive;
    private Integer hspGaps;
    private Integer hspAlignLen;
    private String hspQseq;
    private String hspHseq;
    private String hspIdentityString;
    private Double percentageIdentity = null;
    private Integer mismatchCount = null;
    private SimpleSequencePair<S, C> returnAln;

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 67 * hash + (this.hspQseq != null ? this.hspQseq.hashCode() : 0);
        hash = 67 * hash + (this.hspHseq != null ? this.hspHseq.hashCode() : 0);
        hash = 67 * hash + (this.hspIdentityString != null ? this.hspIdentityString.hashCode() : 0);
        return hash;
    }
    /**
     * Experimental.
     * Wants to implement conceptual comparisons of search results.
     * Fields unrelated to search are deliberately not considered.
     * 
     * In HSP case, alignment representation strings are considered.
     * @return true if HSP alignments are the same, 
     * false otherwise or if alignment strings are undetermined
     */
    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Hsp<?, ?> other = (Hsp<?, ?>) obj;
        if ((this.hspQseq == null) ? (other.hspQseq != null) : !this.hspQseq.equals(other.hspQseq)) {
            return false;
        }
        if ((this.hspHseq == null) ? (other.hspHseq != null) : !this.hspHseq.equals(other.hspHseq)) {
            return false;
        }
        if ((this.hspIdentityString == null) ? (other.hspIdentityString != null) : !this.hspIdentityString.equals(other.hspIdentityString)) {
            return false;
        }
        return true;
    }
    
    public SequencePair<S,C> getAlignment(){
        if (returnAln != null) return returnAln;
        
        SimpleAlignedSequence<S,C> alignedQuery, alignedHit;
        // queryFrom e hitTo?
        int numBefore, numAfter;
        
        alignedQuery = new SimpleAlignedSequence(getSequence(hspQseq), getAlignmentsSteps(hspQseq));
        alignedHit = new SimpleAlignedSequence(getSequence(hspHseq), getAlignmentsSteps(hspHseq));
        
        returnAln = new SimpleSequencePair<S, C>(alignedQuery, alignedHit);
        
        return returnAln;
    }
    
    private Sequence getSequence(String gappedSequenceString){
        if (gappedSequenceString == null) return null;
        
        Sequence returnSeq = null;
        String sequenceString = gappedSequenceString.replace("-", "");
        
        try {
            if (sequenceString.matches("^[ACTG]+$")) 
                returnSeq = new DNASequence(sequenceString, DNACompoundSet.getDNACompoundSet());
            else if (sequenceString.matches("^[ACUG]+$"))
                returnSeq = new RNASequence(sequenceString, DNACompoundSet.getDNACompoundSet());
            else
                returnSeq = new ProteinSequence(sequenceString, AminoAcidCompoundSet.getAminoAcidCompoundSet());
        } catch (CompoundNotFoundException ex) {
            logger.error("Unexpected error, could not find compound when creating Sequence object from Hsp", ex);
        }
        return returnSeq;
    }
    
    private List<Step> getAlignmentsSteps(String gappedSequenceString){
        List<Step> returnList = new ArrayList<Step>();
        
        for (char c: gappedSequenceString.toCharArray()){
            if (c=='-') returnList.add(Step.GAP); else returnList.add(Step.COMPOUND);
        }
        return returnList;
    }

    public int getHspNum() {
        return hspNum;
    }

    public double getHspBitScore() {
        return hspBitScore;
    }

    public int getHspScore() {
        return hspScore;
    }

    public double getHspEvalue() {
        return hspEvalue;
    }

    public int getHspQueryFrom() {
        return hspQueryFrom;
    }

    public int getHspQueryTo() {
        return hspQueryTo;
    }

    public int getHspHitFrom() {
        return hspHitFrom;
    }

    public int getHspHitTo() {
        return hspHitTo;
    }

    public int getHspQueryFrame() {
        return hspQueryFrame;
    }

    public int getHspHitFrame() {
        return hspHitFrame;
    }

    public int getHspIdentity() {
        return hspIdentity;
    }

    public int getHspPositive() {
        return hspPositive;
    }

    public int getHspGaps() {
        return hspGaps;
    }

    public int getHspAlignLen() {
        return hspAlignLen;
    }
    /**
     * HSP aligned query sequence string
     * @return 
     */
    public String getHspQseq() {
        return hspQseq;
    }
    /**
     * HSP aligned hit sequence string
     * @return 
     */
    public String getHspHseq() {
        return hspHseq;
    }
    /**
     * Identity string representing correspondence between aligned residues
     * @return 
     */
    public String getHspIdentityString() {
        return hspIdentityString;
    }

    public Double getPercentageIdentity() {
        if (percentageIdentity != null) return percentageIdentity;
        if (hspIdentity!= null && hspAlignLen != null) return (double)hspIdentity/hspAlignLen;
        return null;
    }

    public Integer getMismatchCount() {
        if (mismatchCount != null) return mismatchCount;
        if (hspIdentity!= null && hspAlignLen != null) return hspIdentity-hspAlignLen;
        return null;
    }

    public Hsp(int hspNum, double hspBitScore, int hspScore, double hspEvalue, int hspQueryFrom, int hspQueryTo, int hspHitFrom, int hspHitTo, int hspQueryFrame, int hspHitFrame, int hspIdentity, int hspPositive, int hspGaps, int hspAlignLen, String hspQseq, String hspHseq, String hspIdentityString, Double percentageIdentity, Integer mismatchCount) {
        this.hspNum = hspNum;
        this.hspBitScore = hspBitScore;
        this.hspScore = hspScore;
        this.hspEvalue = hspEvalue;
        this.hspQueryFrom = hspQueryFrom;
        this.hspQueryTo = hspQueryTo;
        this.hspHitFrom = hspHitFrom;
        this.hspHitTo = hspHitTo;
        this.hspQueryFrame = hspQueryFrame;
        this.hspHitFrame = hspHitFrame;
        this.hspIdentity = hspIdentity;
        this.hspPositive = hspPositive;
        this.hspGaps = hspGaps;
        this.hspIdentity = hspAlignLen;
        this.hspQseq = hspQseq;
        this.hspHseq = hspHseq;
        this.hspIdentityString = hspIdentityString;
        this.percentageIdentity = percentageIdentity; 
        this.mismatchCount = mismatchCount;
        
        // sanity check
        if (percentageIdentity != null && (percentageIdentity < 0 || percentageIdentity >1))
            throw new IllegalArgumentException("Percentage identity must be between 0 and 1");
        
    }
      
}
