package org.biojava.bio.structure;

import java.util.Arrays;
import java.util.List;

/** contains only the static declaration of which types of Groups are available
 *
 * @author Andreas Prlic
 * @since 1.7
 *
 */
public class GroupType {

	/** the type for amino acids
	 *
	 */
	public static final String AMINOACID 	= "amino";

	/** the type for hetero groups
	 *
	 */
	public static final String HETATM    	= "hetatm";

	/** the type for nucelotide groups
	 *
	 */
	public static final String NUCLEOTIDE   = "nucleotide";
	
	
	public static final List<String> WATERNAMES = Arrays.asList(new String[]{"HOH", "DOD",  "WAT"});
	   
}
