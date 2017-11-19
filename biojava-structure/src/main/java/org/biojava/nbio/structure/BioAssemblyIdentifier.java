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
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.nbio.structure.align.util.AtomCache;

public class BioAssemblyIdentifier implements StructureIdentifier {

	private static final long serialVersionUID = -356206725119993449L;
	
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
		return cache.getBiologicalAssembly(pdbCode, biolNr, AtomCache.DEFAULT_BIOASSEMBLY_STYLE);
	}

	@Override
	public SubstructureIdentifier toCanonical() throws StructureException {
		return new SubstructureIdentifier(pdbCode, new ArrayList<ResidueRange>());
	}

	@Override
	public Structure reduce(Structure input) throws StructureException {
		// Should be the full structure
		return input;
	}

}
