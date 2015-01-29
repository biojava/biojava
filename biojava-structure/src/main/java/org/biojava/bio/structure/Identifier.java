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

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.cath.CathFactory;
import org.biojava.bio.structure.scop.ScopFactory;

/**
 * A collection of utilities to create {@link StructureIdentifier StructureIdentifiers}.
 * @author dmyersturnbull
 */
public class Identifier {

	private static final String CATH_PATTERN = "[0-9][a-z0-9]{3}.[0-9]{2}";
	private static final String SCOP_PATTERN = "d[0-9][a-zA-Z0-9]{3,4}([a-zA-Z][0-9_]|\\.[0-9]+)";

	/**
	 * Loads a {@link StructureIdentifier} from the specified string.
	 * The type returned for any particular string can be considered relatively stable
	 * but should not be relied on.
	 * 
	 */
	public static StructureIdentifier loadIdentifier(String id, AtomCache cache) {
		if (id.matches(CATH_PATTERN)) {
			return CathFactory.getCathDatabase().getDescriptionByCathId(id);
		} else if (id.matches(SCOP_PATTERN)) {
			return ScopFactory.getSCOP().getDomainByScopID(id);
		}
		try {
			return new SubstructureIdentifier(id, cache);
		} catch (Exception e) {
			throw new IllegalArgumentException("Couldn't understand id " + id, e);
		}
	}

}
