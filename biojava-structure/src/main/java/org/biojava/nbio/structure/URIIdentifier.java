package org.biojava.nbio.structure;

import java.io.IOException;
import java.net.URI;
import java.net.URL;
import java.util.ArrayList;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A stub StructureIdentifier, representing the full structure in all cases.
 * @author Spencer Bliven
 *
 */
public class URIIdentifier implements StructureIdentifier {
	private static final Logger logger = LoggerFactory.getLogger(URIIdentifier.class);

	private URI uri;
	public URIIdentifier(URI uri) {
		this.uri = uri;
	}
	@Override
	public String getIdentifier() {
		return uri.toString();
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
		//TODO
		return null;
	}

}
