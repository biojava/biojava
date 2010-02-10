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
 * Listens to events that represent Nexus distances blocks.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public interface DistancesBlockListener extends NexusBlockListener {

	/**
	 * Set the NTAX value.
	 * 
	 * @param dimensionsNTax
	 *            the NTAX value.
	 */
	public void setDimensionsNTax(int dimensionsNTax);

	/**
	 * Set the NCHAR value.
	 * 
	 * @param dimensionsNChar
	 *            the NCHAR value.
	 */
	public void setDimensionsNChar(int dimensionsNChar);

	public void setTriangle(final String triangle);

	public void setDiagonal(final boolean diagonal);

	public void setLabels(final boolean labels);

	public void setMissing(final String missing);

	public void setInterleaved(final boolean interleaved);

	/**
	 * Add a TAXLABEL. If it already exists, or is a number that refers to an
	 * index position that already exists, an exception is thrown.
	 * 
	 * @param taxLabel
	 *            the label to add.
	 * @throws ParseException
	 *             if the label cannot be added.
	 */
	public void addTaxLabel(final String taxLabel) throws ParseException;

	public void addMatrixEntry(final String taxa);

	public void appendMatrixData(final String taxa, final Object data);
}
