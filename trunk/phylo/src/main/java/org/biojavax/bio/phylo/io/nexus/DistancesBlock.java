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

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.seq.io.ParseException;

/**
 * Represents Nexus distances blocks.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public class DistancesBlock extends NexusBlock.Abstract {

	/**
	 * A constant representing the name of Distances blocks.
	 */
	public static final String DISTANCES_BLOCK = "DISTANCES";

	private int dimensionsNTax = 0;

	private int dimensionsNChar = 0;

	private String triangle = "LOWER";

	private boolean diagonal = true;

	private boolean labels = true;

	private String missing = "?";

	private boolean interleaved = false;

	private List taxLabels = new ArrayList();

	// values are lists, containing strings and nulls which are gaps
	private Map matrix = new LinkedHashMap();

	private List comments = new ArrayList();

	/**
	 * Delegates to NexusBlock.Abstract constructor using
	 * DistancesBlock.DISTANCES_BLOCK as the name.
	 */
	public DistancesBlock() {
		super(DistancesBlock.DISTANCES_BLOCK);
	}

	/**
	 * Set the NTAX value.
	 * 
	 * @param dimensionsNTax
	 *            the NTAX value.
	 */
	public void setDimensionsNTax(int dimensionsNTax) {
		this.dimensionsNTax = dimensionsNTax;
	}

	/**
	 * Get the NTAX value.
	 * 
	 * @return the NTAX value.
	 */
	public int getDimensionsNTax() {
		return this.dimensionsNTax;
	}

	/**
	 * Set the NCHAR value.
	 * 
	 * @param dimensionsNChar
	 *            the NCHAR value.
	 */
	public void setDimensionsNChar(int dimensionsNChar) {
		this.dimensionsNChar = dimensionsNChar;
	}

	/**
	 * Get the NCHAR value.
	 * 
	 * @return the NCHAR value.
	 */
	public int getDimensionsNChar() {
		return this.dimensionsNChar;
	}

	public void setTriangle(final String triangle) {
		this.triangle = triangle;
	}

	public void setDiagonal(final boolean diagonal) {
		this.diagonal = diagonal;
	}

	public boolean isDiagonal() {
		return this.diagonal;
	}

	public void setLabels(final boolean labels) {
		this.labels = labels;
	}

	public boolean isLabels() {
		return this.labels;
	}

	public void setMissing(final String missing) {
		this.missing = missing;
	}

	public String getMissing() {
		return this.missing;
	}

	public void setInterleaved(final boolean interleaved) {
		this.interleaved = interleaved;
	}

	public boolean isInterleaved() {
		return this.interleaved;
	}

	/**
	 * Add a TAXLABEL. If it already exists, or is a number that refers to an
	 * index position that already exists, an exception is thrown.
	 * 
	 * @param taxLabel
	 *            the label to add.
	 * @throws ParseException
	 *             if the label cannot be added.
	 */
	public void addTaxLabel(final String taxLabel) throws ParseException {
		if (this.taxLabels.contains(taxLabel))
			throw new ParseException("Duplicate taxa label: " + taxLabel);
		else
			try {
				// Try it as a number to see if it refers to
				// position we already have.
				final int i = Integer.parseInt(taxLabel);
				if (i <= this.taxLabels.size() + 1)
					throw new ParseException("Taxa label " + i
							+ " refers to already extant taxa position");
			} catch (NumberFormatException e) {
				// It is not a number, so ignore.
			} catch (ParseException e) {
				// Throw it.
				throw e;
			}
		this.taxLabels.add(taxLabel);
	}

	/**
	 * Removes the given TAXLABEL.
	 * 
	 * @param taxLabel
	 *            the label to remove.
	 */
	public void removeTaxLabel(final String taxLabel) {
		this.taxLabels.remove(taxLabel);
	}

	/**
	 * Checks to see if we contain the given TAXLABEL.
	 * 
	 * @param taxLabel
	 *            the label to check for.
	 * @return <tt>true</tt> if we already contain it.
	 */
	public boolean containsTaxLabel(final String taxLabel) {
		if (this.taxLabels.contains(taxLabel))
			return true;
		else
			try {
				// Try it as a number to see if it refers to
				// position we already have.
				final int i = Integer.parseInt(taxLabel);
				if (i <= this.taxLabels.size() + 1)
					return true;
			} catch (NumberFormatException e) {
				// It is not a number, so ignore.
			}
		return false;
	}

	/**
	 * Get the TAXLABEL values added so far.
	 * 
	 * @return this labels so far.
	 */
	public List getTaxLabels() {
		return this.taxLabels;
	}

	public void addMatrixEntry(final String taxa) {
		if (!this.matrix.containsKey(taxa))
			this.matrix.put(taxa, new ArrayList());
	}

	public void appendMatrixData(final String taxa, final Object data) {
		((List) this.matrix.get(taxa)).add(data);
	}

	public List getMatrixData(final String taxa) {
		return (List) this.matrix.get(taxa);
	}
	
	public Collection getMatrixLabels() {
		return Collections.unmodifiableSet(this.matrix.keySet());
	}

	/**
	 * Adds a comment.
	 * 
	 * @param comment
	 *            the comment to add.
	 */
	public void addComment(final NexusComment comment) {
		this.comments.add(comment);
	}

	/**
	 * Removes a comment.
	 * 
	 * @param comment
	 *            the comment to remove.
	 */
	public void removeComment(final NexusComment comment) {
		this.comments.remove(comment);
	}

	/**
	 * Returns all comments.
	 * 
	 * @return all the selected comments.
	 */
	public List getComments() {
		return this.comments;
	}

	protected void writeBlockContents(Writer writer) throws IOException {
		for (final Iterator i = this.comments.iterator(); i.hasNext();) {
			((NexusComment) i.next()).writeObject(writer);
			writer.write(NexusFileFormat.NEW_LINE);
		}
		writer.write(" DIMENSIONS ");
		if (!this.taxLabels.isEmpty())
			writer.write("NEWTAXA ");
		if (this.dimensionsNTax > 0)
			writer.write("NTAX=" + this.dimensionsNTax + " ");
		writer.write("NCHAR=" + this.dimensionsNChar + ";"
				+ NexusFileFormat.NEW_LINE);

		writer.write(" FORMAT TRIANGLE=");
		this.writeToken(writer, this.triangle);
		writer.write(this.diagonal ? " DIAGONAL" : " NODIAGONAL");
		writer.write(this.labels ? " LABELS" : " NOLABELS");
		writer.write(" MISSING=");
		this.writeToken(writer, this.missing);
		if (this.interleaved)
			writer.write(" INTERLEAVED");
		writer.write(";" + NexusFileFormat.NEW_LINE);

		if (this.taxLabels.size() > 0) {
			writer.write(" TAXLABELS");
			for (final Iterator i = this.taxLabels.iterator(); i.hasNext();) {
				writer.write(' ');
				this.writeToken(writer, (String) i.next());
			}
			writer.write(";" + NexusFileFormat.NEW_LINE);
		}

		// if statesformat=statespresent and items=1, bracket only multi values,
		// otherwise bracket all values
		// only space tokens if reallyUseTokens=true
		writer.write(" MATRIX" + NexusFileFormat.NEW_LINE);
		for (final Iterator i = this.matrix.entrySet().iterator(); i.hasNext();) {
			final Map.Entry entry = (Map.Entry) i.next();
			writer.write('\t');
			this.writeToken(writer, "" + entry.getKey());
			writer.write('\t');
			for (final Iterator j = ((List) entry.getValue()).iterator(); j
					.hasNext();) {
				final Object obj = j.next();
				if (obj instanceof String)
					this.writeToken(writer, (String) obj);
				if (j.hasNext())
					writer.write('\t');
			}
			writer.write(NexusFileFormat.NEW_LINE);
		}
		writer.write(";" + NexusFileFormat.NEW_LINE);

	}
}
