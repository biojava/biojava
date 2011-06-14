package org.biojava3.aaproperties.convertor;

public abstract class Convertor {
	/**
	 * Based on Table 2 of http://nar.oxfordjournals.org/content/34/suppl_2/W32.full.pdf<br/>
	 * An abstract class to convert a protein sequence into representation of different attribute with each attribute having 3 groups.<br/>
	 * The seven different attributes are<p/>
	 * Hydrophobicity (Polar, Neutral, Hydrophobicity)<br/>
	 * Normalized van der Waals volume (Range 0 Ð 2.78, 2.95 Ð 4.0, 4.03 Ð 8.08)<br/>
	 * Polarity (Value 4.9 Ð 6.2, 8.0 Ð 9.2, 10.4 Ð 13.0)<br/>
	 * Polarizability (Value 0 Ð 1.08, 0.128 Ð 0.186, 0.219 Ð 0.409)<br/>
	 * Charge (Positive, Neutral, Negative)<br/>
	 * Secondary structure (Helix, Strand, Coil)<br/>
	 * Solvent accessibility (Buried, Exposed, Intermediate)<br/>
	 * 
	 * @author kohchuanhock
	 * @version 2011.05.09
	 */
	
	protected final char group1 = '1';
	protected final char group2 = '2';
	protected final char group3 = '3';
	protected final char unknownGroup = '0';
	
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
	public String convert(String sequence){
		String convertedSequence = "";
		String uppercaseSequence = sequence.toUpperCase();
		for(int x = 0; x < sequence.length(); x++){
			convertedSequence += convert(uppercaseSequence.charAt(x));
		}
		return convertedSequence;
	}
	
}
