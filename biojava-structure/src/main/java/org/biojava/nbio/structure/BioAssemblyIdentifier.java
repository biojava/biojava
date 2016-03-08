package org.biojava.nbio.structure;

import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.nbio.structure.align.util.AtomCache;

public class BioAssemblyIdentifier implements StructureIdentifier {
	private String pdbCode;
	private int biolNr;

	public static final Pattern BIO_NAME_PATTERN = Pattern.compile("^(?:BIO:)([0-9][a-z0-9]{3})(?::([0-9]+))?$", Pattern.CASE_INSENSITIVE);

	public BioAssemblyIdentifier(String name) {
		Matcher match = BIO_NAME_PATTERN.matcher(name);
		if(! match.matches() ) {
			throw new IllegalArgumentException("Invalid BIO identifier");
		}
		pdbCode = match.group(1);
		if(match.group(2) != null) {
			biolNr = Integer.parseInt(match.group(2));
		} else {
			biolNr = 1;
		}
	}

	public BioAssemblyIdentifier(String pdbCode, int biolNr) {
		this.pdbCode = pdbCode;
		this.biolNr = biolNr;
	}

	@Override
	public String getIdentifier() {
		if( biolNr < 0) {
			return "BIO:"+pdbCode;
		} else {
			return String.format("BIO:%s:%d",pdbCode,biolNr);
		}
	}
	@Override
	public String toString() {
		return getIdentifier();
	}

	@Override
	public Structure loadStructure(AtomCache cache) throws StructureException,
			IOException {
		return cache.getBiologicalAssembly(pdbCode, biolNr);
	}

	@Override
	public SubstructureIdentifier toCanonical() throws StructureException {
		// null pdbCode indicates that the structure can't be loaded by AtomCache
		return new SubstructureIdentifier(null, new ArrayList<ResidueRange>());
	}

	@Override
	public Structure reduce(Structure input) throws StructureException {
		// Should be the full structure
		return input;
	}

}
