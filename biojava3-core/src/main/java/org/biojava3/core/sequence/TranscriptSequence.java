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
 * Created on DATE
 *
 */

package org.biojava3.core.sequence;

import java.util.LinkedHashMap;

import org.biojava3.core.sequence.transcription.TranscriptionEngine;

/**
 *
 * @author Scooter Willis
 */
public class TranscriptSequence extends DNASequence {

public enum Sense{
    POSITIVE, NEGATIVE, UNDEFINED
}

private final LinkedHashMap<String ,IntronSequence> intronSequenceHashMap = new LinkedHashMap<String,IntronSequence>();
private final LinkedHashMap<String ,ExonSequence> exonSequenceHashMap = new LinkedHashMap<String,ExonSequence>();
private Sense sense = Sense.UNDEFINED;



/**
 *
 * @param parentDNASequence
 * @param begin
 * @param end inclusive of end
 */
    public TranscriptSequence(DNASequence parentDNASequence, int begin, int end, Sense sense){
        setParentDNASequence(parentDNASequence);
        setBegin(begin);
        setEnd(end);
        this.sense = sense;
    }



    public Sense getSense(){
        return sense;
    }

    /**
     *
     * @param accession
     * @return
     */

    public IntronSequence removeIntron(String accession){
        return intronSequenceHashMap.remove(accession);
    }

    /**
     *
     * @param accession
     * @param begin
     * @param end
     * @return
     */

    public IntronSequence addIntron(AccessionID accession,int begin,int end){
        IntronSequence intronSequence = new IntronSequence(this,begin,end);
        intronSequence.setAccession(accession);
        intronSequenceHashMap.put(accession.toString(), intronSequence);
        return intronSequence;
    }

        /**
     *
     * @param accession
     * @return
     */

    public ExonSequence removeExon(String accession){
        return exonSequenceHashMap.remove(accession);
    }

    /**
     *
     * @param accession
     * @param begin
     * @param end
     * @return
     */

    public ExonSequence addExon(AccessionID accession,int begin,int end){
        ExonSequence exonSequence = new ExonSequence(this,begin,end);
        exonSequence.setAccession(accession);
        exonSequenceHashMap.put(accession.toString(), exonSequence);
        return exonSequence;
    }

    public RNASequence getRNACodingSequence(){
        StringBuilder sb = new StringBuilder();
        System.err.println("Need to sort exon list");
        for(ExonSequence exonSequence: exonSequenceHashMap.values()){
            sb.append(exonSequence.toString());
        }
        return new RNASequence(sb.toString());
    }

    public ProteinSequence getProteinSequence() {
      return getProteinSequence(TranscriptionEngine.getDefault());
    }

    public ProteinSequence getProteinSequence(TranscriptionEngine engine) {
        RNASequence rnaCodingSequence = getRNASequence();
        return rnaCodingSequence.getProteinSequence(engine);
    }

}
