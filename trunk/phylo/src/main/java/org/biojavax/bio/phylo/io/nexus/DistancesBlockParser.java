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

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.seq.io.ParseException;

/**
 * Parses Nexus distances blocks.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public class DistancesBlockParser extends NexusBlockParser.Abstract {

	private boolean expectingDimension;

	private boolean expectingNewTaxa;

	private boolean expectingNTax;

	private boolean expectingNTaxEquals;

	private boolean expectingNTaxValue;

	private boolean expectingNChar;

	private boolean expectingNCharEquals;

	private boolean expectingNCharValue;

	private boolean expectingFormat;

	private boolean expectingTaxLabel;

	private boolean expectingTaxLabelValue;

	private boolean expectingMatrix;

	private boolean expectingTriangle;

	private boolean expectingTriangleEquals;

	private boolean expectingTriangleContent;

	private boolean expectingDiagonal;

	private boolean expectingMissing;

	private boolean expectingMissingEquals;

	private boolean expectingMissingContent;

	private boolean expectingLabels;

	private boolean expectingInterleave;

	private boolean expectingMatrixKey;

	private boolean expectingMatrixContent;

	private String currentMatrixKey;

	private String matrixFirstLineKey;

	private List matrixSeenKeys = new ArrayList();

	private String triangleType;

	/**
	 * Delegates to NexusBlockParser.Abstract.
	 * 
	 * @param blockListener
	 *            the listener to send parse events to.
	 */
	public DistancesBlockParser(DistancesBlockListener blockListener) {
		super(blockListener);
	}

	public void resetStatus() {
		this.expectingDimension = true;
		this.expectingNewTaxa = false;
		this.expectingNTax = false;
		this.expectingNTaxEquals = false;
		this.expectingNTaxValue = false;
		this.expectingNChar = false;
		this.expectingNCharEquals = false;
		this.expectingNCharValue = false;
		this.expectingFormat = true;
		this.expectingDiagonal = false;
		this.expectingTaxLabel = true;
		this.expectingTaxLabelValue = false;
		this.expectingMatrix = true;
		this.expectingTriangle = false;
		this.expectingTriangleEquals = false;
		this.expectingTriangleContent = false;
		this.expectingMissing = false;
		this.expectingMissingEquals = false;
		this.expectingMissingContent = false;
		this.expectingLabels = false;
		this.expectingInterleave = false;
		this.expectingMatrixKey = false;
		this.expectingMatrixContent = false;
		this.currentMatrixKey = null;
		this.matrixFirstLineKey = null;
		this.triangleType = "LOWER";
		this.matrixSeenKeys.clear();
	}

	public boolean wantsBracketsAndBraces() {
		return this.expectingMatrixContent;
	}

	public void parseToken(String token) throws ParseException {
		if (this.expectingMatrixContent
				&& "\n".equals(token)) {
			// Special handling for new lines inside matrix data.
			this.expectingMatrixContent = false;
			this.expectingMatrixKey = true;
		} else if (token.trim().length() == 0)
			return;
		else if (this.expectingDimension
				&& "DIMENSIONS".equalsIgnoreCase(token)) {
			this.expectingDimension = false;
			this.expectingNewTaxa = true;
			this.expectingNChar = true;
		} else if (this.expectingNewTaxa && "NEWTAXA".equalsIgnoreCase(token)) {
			this.expectingNewTaxa = false;
			this.expectingNTax = true;
			this.expectingNChar = false;
		} else if (this.expectingNTax && token.toUpperCase().startsWith("NTAX")) {
			this.expectingNTax = false;
			if (token.indexOf('=') >= 0) {
				final String[] parts = token.split("=");
				if (parts.length > 1) {
					this.expectingNChar = true;
					try {
						((DistancesBlockListener) this.getBlockListener())
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
				this.expectingNChar = true;
				try {
					((DistancesBlockListener) this.getBlockListener())
							.setDimensionsNTax(Integer.parseInt(parts[1]));
				} catch (NumberFormatException e) {
					throw new ParseException("Invalid NTAX value: " + parts[1]);
				}
			} else
				this.expectingNTaxValue = true;
		} else if (this.expectingNTaxValue) {
			this.expectingNTaxValue = false;
			try {
				((DistancesBlockListener) this.getBlockListener())
						.setDimensionsNTax(Integer.parseInt(token));
			} catch (NumberFormatException e) {
				throw new ParseException("Invalid NTAX value: " + token);
			}
			this.expectingNChar = true;
		} else if (this.expectingNChar
				&& token.toUpperCase().startsWith("NCHAR")) {
			this.expectingNChar = false;
			if (token.indexOf('=') >= 0) {
				final String[] parts = token.split("=");
				if (parts.length > 1) {
					try {
						((DistancesBlockListener) this.getBlockListener())
								.setDimensionsNChar(Integer.parseInt(parts[1]));
					} catch (NumberFormatException e) {
						throw new ParseException("Invalid NCHAR value: "
								+ parts[1]);
					}
				} else
					this.expectingNCharValue = true;
			} else
				this.expectingNCharEquals = true;
		} else if (this.expectingNCharEquals && token.startsWith("=")) {
			this.expectingNCharEquals = false;
			final String[] parts = token.split("=");
			if (parts.length > 1) {
				try {
					((DistancesBlockListener) this.getBlockListener())
							.setDimensionsNChar(Integer.parseInt(parts[1]));
				} catch (NumberFormatException e) {
					throw new ParseException("Invalid NCHAR value: " + parts[1]);
				}
			} else
				this.expectingNCharValue = true;
		} else if (this.expectingNCharValue) {
			this.expectingNCharValue = false;
			try {
				((DistancesBlockListener) this.getBlockListener())
						.setDimensionsNChar(Integer.parseInt(token));
			} catch (NumberFormatException e) {
				throw new ParseException("Invalid NCHAR value: " + token);
			}
		}

		else if (this.expectingFormat && "FORMAT".equalsIgnoreCase(token)) {
			this.expectingFormat = false;
			this.expectingTriangle = true;
			this.expectingDiagonal = true;
			this.expectingMissing = true;
			this.expectingLabels = true;
			this.expectingInterleave = true;
		}

		else if (this.expectingTriangle
				&& token.toUpperCase().startsWith("TRIANGLE")) {
			this.expectingTriangle = false;

			if (token.indexOf("=") >= 0) {
				final String[] parts = token.split("=");
				if (parts.length > 1) {
					this.triangleType = parts[1];
					((DistancesBlockListener) this.getBlockListener())
							.setTriangle(parts[1]);
				} else
					this.expectingTriangleContent = true;
			} else
				this.expectingTriangleEquals = true;
		}

		else if (this.expectingTriangleEquals && token.startsWith("=")) {
			this.expectingTriangleEquals = false;
			if (token.length() > 1) {
				token = token.substring(1);
				this.triangleType = token;
				((DistancesBlockListener) this.getBlockListener())
						.setTriangle(token);
			} else
				this.expectingTriangleContent = true;
		}

		else if (this.expectingTriangleContent) {
			this.triangleType = token;
			((DistancesBlockListener) this.getBlockListener())
					.setTriangle(token);
			this.expectingTriangleContent = false;
		}

		else if (this.expectingDiagonal && "DIAGONAL".equalsIgnoreCase(token)) {
			((DistancesBlockListener) this.getBlockListener())
					.setDiagonal(true);
			this.expectingDiagonal = false;
		}

		else if (this.expectingDiagonal && "NODIAGONAL".equalsIgnoreCase(token)) {
			((DistancesBlockListener) this.getBlockListener())
					.setDiagonal(false);
			this.expectingDiagonal = false;
		}

		else if (this.expectingLabels && "LABELS".equalsIgnoreCase(token)) {
			((DistancesBlockListener) this.getBlockListener()).setLabels(true);
			this.expectingLabels = false;
		}

		else if (this.expectingLabels && "NOLABELS".equalsIgnoreCase(token)) {
			((DistancesBlockListener) this.getBlockListener()).setLabels(false);
			this.expectingLabels = false;
		}

		else if (this.expectingMissing
				&& token.toUpperCase().startsWith("MISSING")) {
			this.expectingMissing = false;

			if (token.indexOf("=") >= 0) {
				final String[] parts = token.split("=");
				if (parts.length > 1)
					((DistancesBlockListener) this.getBlockListener())
							.setMissing(parts[1]);
				else
					this.expectingMissingContent = true;
			} else
				this.expectingMissingEquals = true;
		}

		else if (this.expectingMissingEquals && token.startsWith("=")) {
			this.expectingMissingEquals = false;
			if (token.length() > 1)
				((DistancesBlockListener) this.getBlockListener())
						.setMissing(token.substring(1));
			else
				this.expectingMissingContent = true;
		}

		else if (this.expectingMissingContent) {
			((DistancesBlockListener) this.getBlockListener())
					.setMissing(token);
			this.expectingMissingContent = false;
		}

		else if (this.expectingInterleave
				&& "INTERLEAVE".equalsIgnoreCase(token)) {
			((DistancesBlockListener) this.getBlockListener())
					.setInterleaved(true);
			this.expectingInterleave = false;
		}

		else if (this.expectingTaxLabel && "TAXLABELS".equalsIgnoreCase(token)) {
			this.expectingFormat = false;
			this.expectingTriangle = false;
			this.expectingLabels = false;
			this.expectingDiagonal = false;
			this.expectingMissing = false;
			this.expectingInterleave = false;
			this.expectingTaxLabel = false;
			this.expectingTaxLabelValue = true;
		}

		else if (this.expectingMatrix && "MATRIX".equalsIgnoreCase(token)) {
			this.expectingFormat = false;
			this.expectingTriangle = false;
			this.expectingLabels = false;
			this.expectingDiagonal = false;
			this.expectingMissing = false;
			this.expectingInterleave = false;
			this.expectingTaxLabel = false;
			this.expectingTaxLabelValue = false;
			this.expectingMatrix = false;
			this.expectingMatrixKey = true;
		}

		else if (this.expectingTaxLabelValue)
			// Use untoken version to preserve spaces.
			((DistancesBlockListener) this.getBlockListener())
					.addTaxLabel(token);

		else if (this.expectingMatrixKey) {
			this.currentMatrixKey = token;
			// Use untoken version to preserve spaces.
			((DistancesBlockListener) this.getBlockListener())
					.addMatrixEntry(token);
			this.expectingMatrixKey = false;
			this.expectingMatrixContent = true;
			// Update first line info and set up stack for entry.
			if (!this.matrixSeenKeys.contains(token)) {
				if (this.triangleType.equalsIgnoreCase("UPPER"))
					for (int i = 0; i < this.matrixSeenKeys.size(); i++)
						((DistancesBlockListener) this.getBlockListener())
								.appendMatrixData(this.currentMatrixKey, null);
				this.matrixSeenKeys.add(token);
			}
			if (this.matrixFirstLineKey == null)
				this.matrixFirstLineKey = this.currentMatrixKey;
		}

		else if (this.expectingMatrixContent)
			((DistancesBlockListener) this.getBlockListener())
					.appendMatrixData(this.currentMatrixKey, token);

		else
			throw new ParseException("Found unexpected token " + token
					+ " in DISTANCES block");
	}
}
