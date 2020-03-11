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
package org.biojava.nbio.structure.domain;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;

import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PDPDomain implements StructureIdentifier {
	private static final long serialVersionUID = 6894463080739943026L;

	private final String identifier;
	private final SubstructureIdentifier canonical;

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
        return canonical.pdbId;
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
	public Structure loadStructure(AtomCache cache) throws StructureException {
		return canonical.loadStructure(cache);
	}
}
