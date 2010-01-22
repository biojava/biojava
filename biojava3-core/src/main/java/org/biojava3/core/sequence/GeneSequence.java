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

import java.util.Collection;
import java.util.LinkedHashMap;

import org.biojava3.core.sequence.TranscriptSequence.Sense;

/**
 *
 * @author Scooter Willis
 */
public class GeneSequence extends DNASequence{

private final LinkedHashMap<String ,TranscriptSequence> transcriptSequenceHashMap = new LinkedHashMap<String,TranscriptSequence>();

/**
 * 
 * @param parentDNASequence
 * @param begin
 * @param end inclusive of end
 */
    public GeneSequence(DNASequence parentDNASequence, int begin, int end){
        setParentDNASequence(parentDNASequence);
        setBegin(begin);
        setEnd(end);
    }

    public TranscriptSequence getTranscript(String accession){
        return transcriptSequenceHashMap.get(accession);
    }

    public Collection<TranscriptSequence> getTranscripts(){
        return transcriptSequenceHashMap.values();
    }

    public TranscriptSequence removeTranscript(String accession){


        return transcriptSequenceHashMap.remove(accession);
    }

    public TranscriptSequence addTranscript(AccessionID accession,int begin,int end, Sense sense){
        TranscriptSequence transcriptSequence = new TranscriptSequence(this,begin,end, sense);
        transcriptSequence.setAccession(accession);
        transcriptSequenceHashMap.put(accession.toString(), transcriptSequence);
        return transcriptSequence;
    }



}
