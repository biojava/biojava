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

import java.io.IOException;
import java.util.List;

import org.biojava.bio.structure.align.util.AtomCache;

/**
 * An arbitrary collection of residues in a {@link Structure}.
 * @author dmyersturnbull
 */
public class SubstructureIdentifier implements StructureIdentifier {

	private final String pdbId;
	private final List<ResidueRange> ranges;
	
	public SubstructureIdentifier(String pdbId, List<ResidueRange> ranges) {
		this.pdbId = pdbId;
		this.ranges = ranges;
	}

	public SubstructureIdentifier(String id, AtomCache cache) throws IOException, StructureException {
		if (id.contains(".")) {
			this.pdbId = id.split("\\.")[0];
		} else {
			this.pdbId = id;
		}
		AtomPositionMap map = new AtomPositionMap(cache.getAtoms(pdbId));
		if (id.contains("_")) {
			String[] s = id.split("\\.");
			this.ranges = ResidueRange.parseMultiple(s[1], map);
		} else {
//			this.ranges = new ArrayList<ResidueRange>();
			this.ranges = map.getRanges();
		}
	}

	@Override
	public String getIdentifier() {
		if (ranges.isEmpty()) return pdbId;
		return pdbId + "." + ResidueRange.toString(ranges);
	}

	@Override
	public String getPdbId() {
		return pdbId;
	}

	@Override
	public List<ResidueRange> getResidueRanges() {
		return ranges;
	}

	@Override
	public List<String> getRanges() {
		return ResidueRange.toStrings(ranges);
	}

	@Override
	public String toString() {
		return getIdentifier();
	}
	
}
