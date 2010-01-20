/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.transcription;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.RNATools;
import org.biojava.bio.symbol.SymbolList;
import org.biojava3.core.sequence.RNASequence;

/**
 *
 * @author Scooter
 */
public class DefaultRNAProteinTranscription implements RNAProteinTranscription {

    public String translate(RNASequence rnaCodingSequence) throws Exception {
        String codingSequence = rnaCodingSequence.getString();
        SymbolList rnaSymbolList = RNATools.createRNA(codingSequence);

        //truncate to a length divisible by three.
        rnaSymbolList = rnaSymbolList.subList(1, rnaSymbolList.length() - (rnaSymbolList.length() % 3));
        SymbolList aminoAcidSymbolList = RNATools.translate(rnaSymbolList);

        return aminoAcidSymbolList.seqString();
    }
}
