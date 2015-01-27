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
 * Created on December 19, 2013
 * Author: Douglas Myers-Turnbull
 */

package org.biojava.bio.structure;

import java.util.List;

/**
 * An identifier that <em>uniquely</em> identifies a whole {@link Structure} or arbitrary substructure,
 * including whole chains, {@link org.biojava.bio.structure.scop.ScopDomain ScopDomains}, and {@link org.biojava.bio.structure.cath.CathDomain CathDomains}.
 * @author dmyersturnbull
 */
public interface StructureIdentifier {

	/**
	 * The unique identifier, using the following formal specification:
	 * <pre>
	 * 		name     := pdbID
	 * 		               | pdbID '.' chainID
	 * 		               | pdbID '.' range
	 * 		               | scopID
	 * 		range         := '('? range (',' range)? ')'?
	 * 		               | chainID
	 * 		               | chainID '_' resNum '-' resNum
	 * 		pdbID         := [0-9][a-zA-Z0-9]{3}
	 * 		chainID       := [a-zA-Z0-9]
	 * 		scopID        := 'd' pdbID [a-z_][0-9_]
	 * 		cathID        := pdbID [A-Z][0-9]{2}
	 * 		resNum        := [-+]?[0-9]+[A-Za-z]?
	 * </pre>
	 * For example:
	 * <pre>
	 * 		1TIM                            #whole structure
	 * 		1tim                            #same as above
	 * 		4HHB.C                          #single chain
	 * 		3AA0.A,B                        #two chains
	 * 		d2bq6a1                         #SCOP domain
	 *      1cukA01                         #CATH domain
	 * 		4GCR.A_1-40                     #substructure
	 *      3iek.A_17-28,A_56-294,A_320-377 #substructure of 3 disjoint parts
	 * </pre>
	 * More options may be added to the specification at a future time.
	 */
	String getIdentifier();
	
	/**
	 * Returns the PDB identifier associated with this StructureIdentifier.
	 */
	String getPdbId();
	
	/**
	 * Returns the list of {@link ResidueRange ResidueRanges} that this StructureIdentifier defines.
	 * This is a unique representation.
	 */
	List<? extends ResidueRange> getResidueRanges();
	
	/**
	 * Returns a list of ranges of the form described in {@link #getIdentifier()}. For example:
	 * <pre>
	 * getRanges().get(0): 'A'
	 * getRanges().get(1): 'B_5-100'
	 * </pre>
	 * This is a unique representation.
	 */
	List<String> getRanges();
	
}
