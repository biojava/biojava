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
package org.biojavax.bio.phylo.io.nexus;

import org.biojava.bio.seq.io.ParseException;

/**
 * Builds Nexus distances blocks.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public class DistancesBlockBuilder extends NexusBlockBuilder.Abstract implements
		DistancesBlockListener {

	private DistancesBlock block;

	protected void addComment(final NexusComment comment) {
		this.block.addComment(comment);
	}

	protected NexusBlock startBlockObject() {
		this.block = new DistancesBlock();
		this.resetStatus();
		return this.block;
	}

	private void resetStatus() {
		// Nothing to do.
	}

	public void endBlock() {
		// Don't care.
	}

	public void endTokenGroup() {
		// Nothing to do.
	}

	public void addMatrixEntry(String taxa) {
		this.block.addMatrixEntry(taxa);
	}

	public void addTaxLabel(String taxLabel) throws ParseException {
		this.block.addTaxLabel(taxLabel);
	}

	public void appendMatrixData(String taxa, Object data) {
		this.block.appendMatrixData(taxa, data);
	}

	public void setDimensionsNChar(int dimensionsNChar) {
		this.block.setDimensionsNChar(dimensionsNChar);
	}

	public void setDimensionsNTax(int dimensionsNTax) {
		this.block.setDimensionsNTax(dimensionsNTax);
	}

	public void setInterleaved(boolean interleaved) {
		this.block.setInterleaved(interleaved);
	}

	public void setLabels(boolean labels) {
		this.block.setLabels(labels);
	}

	public void setMissing(String missing) {
		this.block.setMissing(missing);
	}

	public void setDiagonal(boolean diagonal) {
		this.block.setDiagonal(diagonal);
	}

	public void setTriangle(String triangle) {
		this.block.setTriangle(triangle);
	}

}
