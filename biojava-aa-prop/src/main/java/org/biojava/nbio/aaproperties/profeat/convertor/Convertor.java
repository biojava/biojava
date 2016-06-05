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
 */
package org.biojava.nbio.aaproperties.profeat.convertor;

import org.biojava.nbio.core.sequence.ProteinSequence;

public abstract class Convertor {
	/**
	 * Based on Table 2 of http://nar.oxfordjournals.org/content/34/suppl_2/W32.full.pdf<br/>
	 * An abstract class to convert a protein sequence into representation of different attribute with each attribute having 3 groups.<br/>
	 * The seven different attributes are<p/>
	 * Hydrophobicity (Polar, Neutral, Hydrophobicity)<br/>
	 * Normalized van der Waals volume (Range 0 - 2.78, 2.95 - 4.0, 4.03 - 8.08)<br/>
	 * Polarity (Value 4.9 - 6.2, 8.0 - 9.2, 10.4 - 13.0)<br/>
	 * Polarizability (Value 0 - 1.08, 0.128 - 0.186, 0.219 - 0.409)<br/>
	 * Charge (Positive, Neutral, Negative)<br/>
	 * Secondary structure (Helix, Strand, Coil)<br/>
	 * Solvent accessibility (Buried, Exposed, Intermediate)<br/>
	 *
	 * @author kohchuanhock
	 * @version 2011.06.09
	 */
	public final static char group1 = '1';
	public final static char group2 = '2';
	public final static char group3 = '3';
	public final static char unknownGroup = '0';

	/**
	 * Returns the grouping of the amino acid character.
	 * The aminoAcid argument is preferably of non-ambiguous characters.
	 * Standard amino acids will be converted to '1', '2' or '3' depending on its grouping
	 * Non-standard amino acids are simply converted to '0'.
	 *
	 * @param aminoAcid
	 * 		an amino acid character preferably of non-ambiguous characters
	 * @return its grouping
	 */
	public abstract char convert(char aminoAcid);

	/**
	 * Returns the groupings of the attribute
	 * @return the groupings of the attribute
	 */
	public abstract String[] getGrouping();

	/**
	 * Return the attribute of the grouping
	 * @return the attribute of the grouping
	 */
	public abstract String getAttribute();

	/**
	 * Returns the converted sequence.
	 * The sequence argument must be a protein sequence consisting of preferably non-ambiguous characters only.
	 * Standard amino acids will be converted to '1', '2' or '3' depending on its grouping
	 * Non-standard amino acids are simply converted to '0'.
	 *
	 * @param sequence
	 * 		a protein sequence consisting of preferably non-ambiguous characters only
	 * @return the converted sequence
	 */
	public String convert(ProteinSequence sequence){
		String convertedSequence = "";
		String uppercaseSequence = sequence.getSequenceAsString().toUpperCase();
		for(int x = 0; x < uppercaseSequence.length(); x++){
			convertedSequence += String.valueOf(convert(uppercaseSequence.charAt(x)));
		}
		return convertedSequence;
	}

}
