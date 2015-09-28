package org.biojava.nbio.core.search.io;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
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
    
    /**
     * Experimental.
     * Wants to return an hashcode designed to allow conceptual comparisons of search results.
     * Wants to implement conceptual comparisons of search results.
     * Fields unrelated to search are deliberately not considered.
     * @return 
     */
    public int hashCode(){
        String allInOne = hspQseq+"\n"+hspIdentityString+"\n"+hspHseq;
        return allInOne.hashCode();
    }
    
    @Override
    public boolean equals(Object o){
        if (!(o instanceof Hsp)) return false;
        Hsp other = (Hsp)o;
        //if (this.getRepresentationString()==null || other.getRepresentationString()==null) return false;
        
        return o.hashCode() == this.hashCode();
    }
    
    public SequencePair<S,C> getAlignment(){
        if (returnAln != null) return returnAln;
        
        SimpleAlignedSequence alignedQuery, alignedHit;
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
            Logger.getLogger(Hsp.class.getName()).log(Level.SEVERE, null, ex);
        }
        return returnSeq;
    }
    
    private List<Step> getAlignmentsSteps(String gappedSequenceString){
        List<Step> returnList = new ArrayList();
        
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
