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
