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

import java.util.logging.Logger;

/**
 * A gene contains a collection of Exon sequences
 * @author Scooter Willis
 */
public class ExonSequence extends DNASequence {

    private static final Logger log = Logger.getLogger(ExonSequence.class.getName());

    /**
     * Need a parent gene sequence and the bioBegin and bioEnd. An Exon sequence doesn't actually imply what the
     * protein coding sequence will be. This is a little difficult to model and have it make sense.
     * A gene has a collection of Exon and Intron sequences where the Exon sequences will join up. A gene
     * sequences has a collection of different possible isoform proteins based on the transcription rules.
     * A TranscriptionSequence will contain CDSSequence where the CDSSequence will be contained in the ExonSequence.
     * Thus a ExonSequence is the union of overlapping CDSSequences.
     * @param parentGeneSequence
     * @param bioBegin
     * @param bioEnd
     */
    public ExonSequence(GeneSequence parentGeneSequence, int bioBegin, int bioEnd) {
        this.setParentSequence(parentGeneSequence);
        setBioBegin(bioBegin);
        setBioEnd(bioEnd);

    }

    @Override
    public int getLength() {
        return Math.abs(this.getBioEnd() - this.getBioBegin()) + 1;
    }


}
