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

import java.util.ArrayList;
import java.util.List;

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

	public SubstructureIdentifier(String id) {
		if (id.contains(".")) {
			this.pdbId = id;
			this.ranges = new ArrayList<ResidueRange>();
		} else {
			String[] s = id.split("\\.");
			this.pdbId = s[0];
			this.ranges = ResidueRange.parseMultiple(s[1]);
		}
	}

	@Override
	public String getIdentifier() {
		return ResidueRange.toString(ranges);
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

}
