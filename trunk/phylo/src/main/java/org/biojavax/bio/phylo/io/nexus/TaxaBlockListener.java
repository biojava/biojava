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
 * Listens to events that represent Nexus taxa blocks.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public interface TaxaBlockListener extends NexusBlockListener {

	/**
	 * Set the DIMENSIONS NTAX value.
	 * 
	 * @param dimensionsNTax
	 *            the new value.
	 */
	public void setDimensionsNTax(int dimensionsNTax);

	/**
	 * Add another value after the TAXLABEL tag.
	 * 
	 * @param taxLabel
	 *            the new taxa to add.
	 * @throws ParseException
	 *             if the label is invalid.
	 */
	public void addTaxLabel(final String taxLabel) throws ParseException;
}
