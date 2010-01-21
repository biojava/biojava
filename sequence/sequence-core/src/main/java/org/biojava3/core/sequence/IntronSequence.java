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

/**
 *
 * @author Scooter Willis
 */
public class IntronSequence extends DNASequence{
private AccessionID accession;
public DNASequence parentGeneSequence = null;
int begin = -1;
int end = -1;

    public IntronSequence(TranscriptSequence parentGeneSequence, int begin, int end){
        this.parentGeneSequence = parentGeneSequence;
        this.begin = begin;
        this.end = end;
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
