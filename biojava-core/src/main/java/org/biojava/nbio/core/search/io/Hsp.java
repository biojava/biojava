/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava.nbio.core.search.io;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.nbio.alignment.SimpleAlignedSequence;
import org.biojava.nbio.alignment.SimpleSequencePair;
import org.biojava.nbio.alignment.template.AlignedSequence.Step;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 *
 * @author pavanpa
 */
public abstract class Hsp <S extends Sequence<C>, C extends Compound> {
    private final int hspNum;
    private final double hspBitScore;
    private final int hspScore;
    private final double hspEvalue;
    private final int hspQueryFrom;
    private final int hspQueryTo;
    private final int hspHitFrom;
    private final int hspHitTo;
    private final int hspQueryFrame;
    private final int hspHitFrame;
    private final int hspIdentity;
    private final int hspPositive;
    private final int hspGaps;
    private final int hspAlignLen;
    private final String hspQseq;
    private final String hspHseq;
    private final String hspIdentityString;
    
    /**
     * Experimental.
     * Wants to return an hashcode designed to allow conceptual comparisons of search results.
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
        
        return o.hashCode() == this.hashCode();
    }
    
    public SequencePair<S,C> getAlignment(){
        SimpleSequencePair<S, C> returnAln;
        SimpleAlignedSequence alignedQuery, alignedHit;
        // queryFrom e hitTo?
        int numBefore, numAfter;
        
        alignedQuery = new SimpleAlignedSequence(getSequence(hspQseq), getAlignmentsSteps(hspQseq));
        alignedHit = new SimpleAlignedSequence(getSequence(hspHseq), getAlignmentsSteps(hspHseq));
        
        returnAln = new SimpleSequencePair<S, C>(alignedQuery, alignedHit);
        
        return returnAln;
    }
    
    private Sequence getSequence(String gappedSequenceString){
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
    
    

    public Hsp(int hspNum, double hspBitScore, int hspScore, double hspEvalue, int hspQueryFrom, int hspQueryTo, int hspHitFrom, int hspHitTo, int hspQueryFrame, int hspHitFrame, int hspIdentity, int hspPositive, int hspGaps, int hspAlignLen, String hspQseq, String hspHseq, String hspIdentityString) {
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
        this.hspAlignLen = hspAlignLen;
        this.hspQseq = hspQseq;
        this.hspHseq = hspHseq;
        this.hspIdentityString = hspIdentityString;
    }
      
}
