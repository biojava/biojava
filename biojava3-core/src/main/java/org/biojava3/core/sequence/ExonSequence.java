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

import java.util.ArrayList;
import java.util.logging.Logger;
import org.biojava3.core.sequence.TranscriptSequence.Sense;

/**
 *
 * @author Scooter Willis
 */
public class ExonSequence extends DNASequence {
    private static final Logger log = Logger.getLogger(ExonSequence.class.getName());
    public DNASequence parentGeneSequence = null;    
    Sense sense = Sense.UNDEFINED;
    private final ArrayList<CDSSequence> cdsSequenceList = new ArrayList<CDSSequence>();

    public ExonSequence(TranscriptSequence parentGeneSequence, int begin, int end, Sense sense) {
        this.parentGeneSequence = parentGeneSequence;
        setBegin(begin);
        setEnd(end);
        this.sense = sense;
    }


        /**
     *
     * @param accession
     * @return
     */
    public CDSSequence removeCDS(String accession) {
        for (CDSSequence cdsSequence : cdsSequenceList) {
            if (cdsSequence.getAccession().getID().equals(accession)) {
                cdsSequenceList.remove(cdsSequence);
                return cdsSequence;
            }
        }
        return null;
    }

    /**
     *
     * @param accession
     * @param begin
     * @param end
     * @param phase 0,1,2
     * @return
     */
    public CDSSequence addCDS(AccessionID accession, int begin, int end, int phase) {
        CDSSequence cdsSequence = new CDSSequence(this, begin, end,phase, sense); //sense should be the same as parent
        cdsSequence.setAccession(accession);
        cdsSequenceList.add(cdsSequence);
        return cdsSequence;
    }


    @Override
    public String toString(){
        String sequence = super.toString();
        if(sense == Sense.NEGATIVE){
            StringBuilder sb = new StringBuilder(sequence);
            sequence = sb.reverse().toString();
        }

        return sequence;
    }

   



}
