/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence;

/**
 *
 * @author Scooter
 */
public class IntronSequence extends DNASequence{
private AccessionID accession;
public DNASequence parentGeneSequence = null;
int begin = -1;
int end = -1;

    public IntronSequence(TranscriptSequence _parentGeneSequence, int _begin, int _end){
        parentGeneSequence = _parentGeneSequence;
        begin = _begin;
        end = _end;
    }


        /**
     * @return the accession
     */
    public AccessionID getAccession() {
        return accession;
    }

    /**
     * @param accession the accession to set
     */
    public void setAccession(AccessionID accession) {
        this.accession = accession;
    }
}
