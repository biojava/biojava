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

import java.util.List;

import org.biojava.bio.seq.io.ParseException;

/**
 * Listens to events that represent Nexus characters blocks.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public interface CharactersBlockListener extends NexusBlockListener {

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

	public void setDataType(final String dataType);

	public void setRespectCase(final boolean respectCase);

	public void setMissing(final String missing);

	public void setGap(final String gap);

	public void addSymbol(final String symbol);

	public void addEquate(final String symbol, final List symbols);

	public void setMatchChar(final String matchChar);

	public void setLabels(final boolean labels);

	public void setTransposed(final boolean transposed);

	public void setInterleaved(final boolean interleaved);

	public void addItem(final String item);

	public void setStatesFormat(final String statesFormat);

	public void setTokens(final boolean tokens);

	public void setEliminateStart(final int eliminateStart);

	public void setEliminateEnd(final int eliminateEnd);

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

	public void addCharState(final String charState);

	public void setCharStateLabel(final String charState, final String label);

	public void addCharStateKeyword(final String charState, final String keyword);

	public void addCharLabel(final String charLabel);

	public void addState(final String state);

	public void addStateLabel(final String state, final String label);

	public void addMatrixEntry(final String taxa);

	public void appendMatrixData(final String taxa, final Object data);
}
