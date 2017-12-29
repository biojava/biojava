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
package org.biojava.nbio.structure;

import java.io.IOException;
import java.util.ArrayList;

import org.biojava.nbio.structure.align.util.AtomCache;

/**
 * A stub StructureIdentifier, representing the full structure in all cases.
 * @author Spencer Bliven
 *
 */
public class PassthroughIdentifier implements StructureIdentifier {

	private static final long serialVersionUID = -2773111624414448950L;
	
	private String identifier;
	public PassthroughIdentifier(String identifier) {
		this.identifier = identifier;
	}
	@Override
	public String getIdentifier() {
		return identifier;
	}

	/**
	 * @return A SubstructureIdentifier without ranges (e.g. including all residues)
	 */
	@Override
	public SubstructureIdentifier toCanonical() {
		return new SubstructureIdentifier(null, new ArrayList<ResidueRange>());
	}

	@Override
	public Structure reduce(Structure input) throws StructureException {
		return input;
	}
	/**
	 * Passthrough identifiers don't know how to load a structure
	 * @return null
	 */
	@Override
	public Structure loadStructure(AtomCache cache) throws StructureException,
			IOException {
		return null;
	}

}
