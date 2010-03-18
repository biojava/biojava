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
import org.biojava3.core.sequence.TranscriptSequence.Sense;

/**
 *
 * @author Scooter Willis
 */
public class CDSSequence extends DNASequence {
    private static final Logger log = Logger.getLogger(CDSSequence.class.getName());
    public DNASequence parentGeneSequence = null;
    int phase = 0; // 0, 1, 2 http://www.sequenceontology.org/gff3.shtml
    Sense sense = Sense.UNDEFINED;

    public CDSSequence(ExonSequence parentGeneSequence, int begin, int end, int phase, Sense sense) {
        this.parentGeneSequence = parentGeneSequence;
        setBegin(begin);
        setEnd(end);
        this.phase = phase;
        this.sense = sense;
    }


    @Override
    public String toString(){
        String sequence = super.toString();
        if(sense == Sense.NEGATIVE){
            StringBuilder sb = new StringBuilder(sequence);
            sequence = sb.reverse().toString();
        }
        if(phase != 0){
            log.severe("Phase does not = 0 in exon sequence code needs to fix this");
        }
        return sequence;
    }

   



}
