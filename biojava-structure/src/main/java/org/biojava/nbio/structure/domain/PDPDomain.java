package org.biojava.nbio.structure.domain;

import java.io.IOException;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.nbio.structure.ResidueRange;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.SubstructureIdentifier;
import org.biojava.nbio.structure.align.util.AtomCache;

public class PDPDomain implements StructureIdentifier {
	String identifier;
	SubstructureIdentifier canonical;
	public static final Pattern PDP_NAME_PATTERN = Pattern.compile("^(?:PDP:)([0-9][a-z0-9]{3})(\\w)(\\w)$",Pattern.CASE_INSENSITIVE);

	public PDPDomain(String pdpDomainName, List<ResidueRange> ranges) {
		this.identifier = pdpDomainName;
		Matcher matcher = PDP_NAME_PATTERN.matcher(identifier);
		if(!matcher.matches()) {
			throw new IllegalArgumentException("Malformed PDP domain name");
		}
		String pdbId = matcher.group(1);
		this.canonical = new SubstructureIdentifier(pdbId,ranges);
	}

	@Override
	public String getIdentifier() {
		return identifier;
	}

	public String getPdbId() {
		return canonical.getPdbId();
	}

	@Override
	public SubstructureIdentifier toCanonical() {
		return canonical;
	}

	@Override
	public Structure reduce(Structure input) throws StructureException {
		return canonical.reduce(input);
	}

	@Override
	public String toString() {
		return getIdentifier();
	}

	@Override
	public Structure loadStructure(AtomCache cache) throws StructureException,
			IOException {
		return canonical.loadStructure(cache);
	}
}
