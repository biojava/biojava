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
 * Used to map the stop codon sequence on a gene
 * @author Scooter Willis
 */
public class StopCodonSequence extends DNASequence {

public DNASequence parentGeneSequence = null;


    public StopCodonSequence(TranscriptSequence parentGeneSequence, int begin, int end){
        this.parentGeneSequence = parentGeneSequence;
        setBioBegin(begin);
        setBioEnd(end);
    }


        @Override
    public int getLength() {
        return Math.abs(this.getBioEnd() - this.getBioBegin()) + 1;
    }
}
