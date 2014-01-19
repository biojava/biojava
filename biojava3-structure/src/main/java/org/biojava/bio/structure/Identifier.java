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

import java.net.URISyntaxException;

/**
 * A collection of utilities to create {@link StructureIdentifier StructureIdentifiers}.
 * @author dmyersturnbull
 */
public class Identifier {

	/**
	 * Loads a {@link StructureIdentifier} from the specified string.
	 * The type returned for any particular string can be considered relatively stable
	 * but should not be relied on.
	 */
	public static StructureIdentifier loadIdentifier(String id) {
		try {
			return new UriIdentifier(id);
		} catch (URISyntaxException e) {
			return new SubstructureIdentifier(id);
		}
	}
	
}
