package org.biojava3.core.util;

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.template.Sequence;

public class SequenceTools {

	protected static final String NUCLEOTIDE_LETTERS = "GCTAUXN";

    public static int percentNucleotideSequence(String sequence)
    {
            if (sequence == null || sequence.length() == 0) return 0;

            int l = sequence.length();
            int n =0;

            for (int i = 0; i < l; i++)
            {
                    if (NUCLEOTIDE_LETTERS.indexOf(sequence.charAt(i)) < 0)
                    {
                            continue;
                    }
                    n++;
            }
            return (100 * n) / l;
    }

    public static boolean isNucleotideSequence(String sequence)
    {
            if (sequence == null || sequence.length() == 0) return false;

            int l = sequence.length();
            for (int i = 0; i < l; i++)
            {
                    if (NUCLEOTIDE_LETTERS.indexOf(sequence.charAt(i)) < 0)
                    {
                            return false;
                    }
            }
            return true;
    }
    
    public Sequence<?> getSeqeunceFromString(String sequence){
    	
  
    	if( isNucleotideSequence(sequence)) {
    		return  new DNASequence(sequence);
    	} else {
    		return new ProteinSequence(sequence);
    	}
    	
    }
	
}
