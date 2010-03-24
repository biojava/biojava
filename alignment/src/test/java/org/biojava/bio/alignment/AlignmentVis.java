/**
 * 
 */
package org.biojava.bio.alignment;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SymbolList;

/**
 * @author draeger
 * 
 */
public class AlignmentVis {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		FiniteAlphabet alpha = DNATools.getDNA();
		short match = 0, replace = 4, insert = 2, delete = 2, ext = 1;
		SubstitutionMatrix mat = new SubstitutionMatrix(alpha, (short) -match, (short) -replace);
		AlignmentAlgorithm aa = new NeedlemanWunsch(match, replace, insert,
				delete, ext, mat);
		try {
			SymbolList s1 = new SimpleSymbolList(
					alpha.getTokenization("token"), args[0]);
			SymbolList s2 = new SimpleSymbolList(
					alpha.getTokenization("token"), args[1]);
			AlignmentPair ap = aa.pairwiseAlignment(s1, s2);
			
			System.out.println(ap.formatOutput(20));
			
		} catch (IllegalSymbolException e) {
			e.printStackTrace();
		} catch (BioException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
