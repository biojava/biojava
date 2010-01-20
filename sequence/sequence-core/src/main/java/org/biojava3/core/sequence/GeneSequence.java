/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence;

import java.util.Collection;
import java.util.LinkedHashMap;
import org.biojava3.core.sequence.TranscriptSequence.Sense;

/**
 *
 * @author Scooter
 */
public class GeneSequence extends DNASequence{

private final LinkedHashMap<String ,TranscriptSequence> transcriptSequenceHashMap = new LinkedHashMap<String,TranscriptSequence>();

/**
 * 
 * @param _parentDNASequence
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
