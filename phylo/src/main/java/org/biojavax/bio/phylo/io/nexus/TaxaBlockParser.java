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
 * Parses Nexus taxa blocks.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public class TaxaBlockParser extends NexusBlockParser.Abstract {

	private boolean expectingDimension;

	private boolean expectingNTax;

	private boolean expectingNTaxEquals;

	private boolean expectingNTaxValue;

	private boolean expectingTaxLabel;

	private boolean expectingTaxLabelValue;

	/**
	 * Delegates to NexusBlockParser.Abstract.
	 * 
	 * @param blockListener
	 *            the listener to send parse events to.
	 */
	public TaxaBlockParser(TaxaBlockListener blockListener) {
		super(blockListener);
	}

	public void resetStatus() {
		this.expectingDimension = true;
		this.expectingNTax = false;
		this.expectingNTaxEquals = false;
		this.expectingNTaxValue = false;
		this.expectingTaxLabel = false;
		this.expectingTaxLabelValue = false;
	}

	public void parseToken(String token) throws ParseException {
		if (token.trim().length() == 0)
			return;
		else if (this.expectingDimension
				&& "DIMENSIONS".equalsIgnoreCase(token)) {
			this.expectingDimension = false;
			this.expectingNTax = true;
		} else if (this.expectingNTax && token.toUpperCase().startsWith("NTAX")) {
			this.expectingNTax = false;
			if (token.indexOf('=') >= 0) {
				final String[] parts = token.split("=");
				if (parts.length > 1) {
					this.expectingTaxLabel = true;
					try {
						((TaxaBlockListener) this.getBlockListener())
								.setDimensionsNTax(Integer.parseInt(parts[1]));
					} catch (NumberFormatException e) {
						throw new ParseException("Invalid NTAX value: "
								+ parts[1]);
					}
				} else
					this.expectingNTaxValue = true;
			} else
				this.expectingNTaxEquals = true;
		} else if (this.expectingNTaxEquals && token.startsWith("=")) {
			this.expectingNTaxEquals = false;
			final String[] parts = token.split("=");
			if (parts.length > 1) {
				this.expectingTaxLabel = true;
				try {
					((TaxaBlockListener) this.getBlockListener())
							.setDimensionsNTax(Integer.parseInt(parts[1]));
				} catch (NumberFormatException e) {
					throw new ParseException("Invalid NTAX value: " + parts[1]);
				}
			} else
				this.expectingNTaxValue = true;
		} else if (this.expectingNTaxValue) {
			this.expectingNTaxValue = false;
			try {
				((TaxaBlockListener) this.getBlockListener())
						.setDimensionsNTax(Integer.parseInt(token));
			} catch (NumberFormatException e) {
				throw new ParseException("Invalid NTAX value: " + token);
			}
			this.expectingTaxLabel = true;
		} else if (this.expectingTaxLabel
				&& "TAXLABELS".equalsIgnoreCase(token)) {
			this.expectingTaxLabel = false;
			this.expectingTaxLabelValue = true;
		} else if (this.expectingTaxLabelValue)
			// Use untoken version to preserve spaces.
			((TaxaBlockListener) this.getBlockListener()).addTaxLabel(token);
		else
			throw new ParseException("Found unexpected token " + token
					+ " in TAXA block");
	}

}
