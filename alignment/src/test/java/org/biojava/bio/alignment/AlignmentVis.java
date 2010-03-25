/**
 * 
 */
package org.biojava.bio.alignment;

import java.io.IOException;

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
		short match = -2, replace = 5, insert = 3, delete = 3, ext = 0;
		SubstitutionMatrix mat;
		try {
			mat = new SubstitutionMatrix(alpha, 
					"#"+
"# This matrix was created by Todd Lowe   12/10/92\n" +
"#\n" +
"# Uses ambiguous nucleotide codes, probabilities rounded to\n" +
"#  nearest integer\n" +
"#\n" +
"# Lowest score = -4, Highest score = 5\n" +
"#\n" +
"    A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N\n" +
"A   5  -4  -4  -4  -4   1   1  -4  -4   1  -4  -1  -1  -1  -2\n" +
"T  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2\n" +
"G  -4  -4   5  -4   1  -4   1  -4   1  -4  -1  -1  -4  -1  -2\n" +
"C  -4  -4  -4   5   1  -4  -4   1  -4   1  -1  -1  -1  -4  -2\n" +
"S  -4  -4   1   1  -1  -4  -2  -2  -2  -2  -1  -1  -3  -3  -1\n" +
"W   1   1  -4  -4  -4  -1  -2  -2  -2  -2  -3  -3  -1  -1  -1\n" +
"R   1  -4   1  -4  -2  -2  -1  -4  -2  -2  -3  -1  -3  -1  -1\n" +
"Y  -4   1  -4   1  -2  -2  -4  -1  -2  -2  -1  -3  -1  -3  -1\n" +
"K  -4   1   1  -4  -2  -2  -2  -2  -1  -4  -1  -3  -3  -1  -1\n" +
"M   1  -4  -4   1  -2  -2  -2  -2  -4  -1  -3  -1  -1  -3  -1\n" +
"B  -4  -1  -1  -1  -1  -3  -3  -1  -1  -3  -1  -2  -2  -2  -1\n" +
"V  -1  -4  -1  -1  -1  -3  -1  -3  -3  -1  -2  -1  -2  -2  -1\n" +
"H  -1  -1  -4  -1  -3  -1  -3  -1  -3  -1  -2  -2  -1  -2  -1\n" +  
"D  -1  -1  -1  -4  -3  -1  -1  -3  -1  -3  -2  -2  -2  -1  -1\n" +
"N  -2  -2  -2  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1\n", "Nuc.4.4");
		AlignmentAlgorithm aa = new NeedlemanWunsch(match, replace, insert,
				delete, ext, mat);
			SymbolList s1 = new SimpleSymbolList(
					alpha.getTokenization("token"), args[0]);
			SymbolList s2 = new SimpleSymbolList(
					alpha.getTokenization("token"), args[1]);
			AlignmentPair ap = aa.pairwiseAlignment(s1, s2);

			System.out.println(ap.formatOutput(60));
			
		} catch (IllegalSymbolException e) {
			e.printStackTrace();
		} catch (BioException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
