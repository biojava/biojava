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

import java.util.List;

/**
 * An arbitrary collection of residues in a {@link Structure}.
 * @author dmyersturnbull
 */
public class Substructure implements StructureIdentifier {

	private String pdbId;
	private List<ResidueRange> ranges;
	
	public Substructure(String pdbId, List<ResidueRange> ranges) {
		this.pdbId = pdbId;
		this.ranges = ranges;
	}

	public Substructure(String id) {
		// TODO Auto-generated constructor stub
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
