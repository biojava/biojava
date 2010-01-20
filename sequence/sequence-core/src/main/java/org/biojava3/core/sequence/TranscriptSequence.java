/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence;

import org.biojava3.core.sequence.transcription.RNAProteinTranscription;
import org.biojava3.core.sequence.transcription.DefaultRNAProteinTranscription;
import java.util.LinkedHashMap;

/**
 *
 * @author Scooter
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
 * @param _parentDNASequence
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
     * @param _begin
     * @param _end
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
     * @param _begin
     * @param _end
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

    public ProteinSequence getProteinSequence() throws Exception{
         return getProteinSequence(new DefaultRNAProteinTranscription());
    }

    public ProteinSequence getProteinSequence(RNAProteinTranscription rnaProteinTranslator) throws Exception{
        RNASequence rnaCodingSequence = getRNASequence();
        String aminoAcidSequence = rnaProteinTranslator.translate(rnaCodingSequence);
        return new ProteinSequence(aminoAcidSequence);

    }

}
