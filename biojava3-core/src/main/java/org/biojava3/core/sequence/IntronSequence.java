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

public DNASequence parentGeneSequence = null;
int frame = 0; // not sure if this makes sense for an intron
 Strand sense = Strand.UNDEFINED;
    public IntronSequence(TranscriptSequence parentGeneSequence, int begin, int end, int frame,  Strand sense){
        this.parentGeneSequence = parentGeneSequence;
        setBegin(begin);
        setEnd(end);
        this.sense = sense;
    }



}
