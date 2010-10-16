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
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

import org.biojava.bio.seq.io.ParseException;

/**
 * Parses Nexus characters blocks.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public class CharactersBlockParser extends NexusBlockParser.Abstract {

	private boolean expectingDimension;

	private boolean expectingNewTaxa;

	private boolean expectingNTax;

	private boolean expectingNTaxEquals;

	private boolean expectingNTaxValue;

	private boolean expectingNChar;

	private boolean expectingNCharEquals;

	private boolean expectingNCharValue;

	private boolean expectingFormat;

	private boolean expectingEliminate;

	private boolean expectingTaxLabel;

	private boolean expectingTaxLabelValue;

	private boolean expectingCharStateLabel;

	private boolean expectingCharLabel;

	private boolean expectingStateLabel;

	private boolean expectingMatrix;

	private boolean expectingDataType;

	private boolean expectingDataTypeEquals;

	private boolean expectingDataTypeContent;

	private boolean expectingRespectCase;

	private boolean expectingMissing;

	private boolean expectingMissingEquals;

	private boolean expectingMissingContent;

	private boolean expectingGap;

	private boolean expectingGapEquals;

	private boolean expectingGapContent;

	private boolean expectingSymbols;

	private boolean expectingSymbolsEquals;

	private boolean expectingSymbolsContent;

	private boolean expectingEquate;

	private boolean expectingEquateEquals;

	private boolean expectingEquateContent;

	private boolean expectingMatchChar;

	private boolean expectingMatchCharEquals;

	private boolean expectingMatchCharContent;

	private boolean expectingLabels;

	private boolean expectingTranspose;

	private boolean expectingInterleave;

	private boolean expectingItems;

	private boolean expectingItemsEquals;

	private boolean expectingItemsContent;

	private boolean itemsInBrackets;

	private boolean expectingStatesFormat;

	private boolean expectingStatesFormatEquals;

	private boolean expectingStatesFormatContent;

	private boolean expectingTokens;

	private String specifiedDataType;

	private boolean tokenizedMatrix;

	private boolean expectingEliminateRange;

	private boolean expectingCharStateLabelKey;

	private boolean expectingCharStateLabelName;

	private boolean expectingCharStateLabelSynonym;

	private boolean expectingCharLabelValue;

	private boolean expectingStateLabelKey;

	private boolean expectingStateLabelContent;

	private boolean expectingMatrixKey;

	private boolean expectingMatrixContent;

	private String currentCharStateLabelKey;

	private String currentStateLabelKey;

	private String currentMatrixKey;

	private List currentMatrixBracket;

	private Map matrixStack = new HashMap();

	private int matrixFirstLineLength;

	private String matrixFirstLineKey;

	private int matrixPrependNulls;

	private boolean seenSymbol;

	/**
	 * Delegates to NexusBlockParser.Abstract.
	 * 
	 * @param blockListener
	 *            the listener to send parse events to.
	 */
	public CharactersBlockParser(CharactersBlockListener blockListener) {
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
		this.expectingFormat = false;
		this.expectingEliminate = false;
		this.expectingTaxLabel = false;
		this.expectingTaxLabelValue = false;
		this.expectingCharStateLabel = false;
		this.expectingCharLabel = false;
		this.expectingStateLabel = false;
		this.expectingMatrix = false;
		this.tokenizedMatrix = false;
		this.specifiedDataType = null;
		this.expectingDataType = false;
		this.expectingDataTypeEquals = false;
		this.expectingDataTypeContent = false;
		this.expectingRespectCase = false;
		this.expectingMissing = false;
		this.expectingMissingEquals = false;
		this.expectingMissingContent = false;
		this.expectingGap = false;
		this.expectingGapEquals = false;
		this.expectingGapContent = false;
		this.expectingSymbols = false;
		this.expectingSymbolsEquals = false;
		this.expectingSymbolsContent = false;
		this.expectingEquate = false;
		this.expectingEquateEquals = false;
		this.expectingEquateContent = false;
		this.expectingMatchChar = false;
		this.expectingMatchCharEquals = false;
		this.expectingMatchCharContent = false;
		this.expectingLabels = false;
		this.expectingTranspose = false;
		this.expectingInterleave = false;
		this.expectingItems = false;
		this.expectingItemsEquals = false;
		this.expectingItemsContent = false;
		this.itemsInBrackets = false;
		this.expectingStatesFormat = false;
		this.expectingStatesFormatEquals = false;
		this.expectingStatesFormatContent = false;
		this.expectingTokens = false;
		this.expectingEliminateRange = false;
		this.expectingCharStateLabelKey = false;
		this.expectingCharStateLabelName = false;
		this.expectingCharStateLabelSynonym = false;
		this.expectingCharLabelValue = false;
		this.expectingStateLabelKey = false;
		this.expectingStateLabelContent = false;
		this.expectingMatrixKey = false;
		this.expectingMatrixContent = false;
		this.currentCharStateLabelKey = null;
		this.currentStateLabelKey = null;
		this.currentMatrixKey = null;
		this.currentMatrixBracket = null;
		this.matrixStack.clear();
		this.matrixFirstLineKey = null;
		this.matrixFirstLineLength = 0;
		this.matrixPrependNulls = 0;
		this.seenSymbol = false;
	}

	public boolean wantsBracketsAndBraces() {
		return this.expectingMatrixContent;
	}

	public void parseToken(String token) throws ParseException {
		if (this.expectingMatrixContent && "\n".equals(token)) {
			// Special handling for new lines inside matrix data.
			if (this.currentMatrixBracket != null) {
				((CharactersBlockListener) this.getBlockListener())
						.appendMatrixData(this.currentMatrixKey,
								this.currentMatrixBracket);
				this.currentMatrixBracket = null;
			}
			this.expectingMatrixContent = false;
			this.expectingMatrixKey = true;
		} else if (this.expectingMatrixKey && "\n".equals(token)) {
			if (this.matrixFirstLineKey != null)
				this.matrixPrependNulls = this.matrixFirstLineLength;
		} else if (token.trim().length() == 0)
			return;
		else if (this.expectingDimension
				&& "DIMENSIONS".equalsIgnoreCase(token)) {
			this.expectingDimension = false;
			this.expectingNewTaxa = true;
			this.expectingNTax = true;
			this.expectingNChar = true;
		} else if (this.expectingNewTaxa && "NEWTAXA".equalsIgnoreCase(token)) {
			this.expectingNewTaxa = false;
			this.expectingNTax = true;
			this.expectingNChar = false;
		} else if (this.expectingNTax && token.toUpperCase().startsWith("NTAX")) {
			this.expectingNewTaxa = false;
			this.expectingNTax = false;
			if (token.indexOf('=') >= 0) {
				final String[] parts = token.split("=");
				if (parts.length > 1) {
					this.expectingNChar = true;
					try {
						((CharactersBlockListener) this.getBlockListener())
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
					((CharactersBlockListener) this.getBlockListener())
							.setDimensionsNTax(Integer.parseInt(parts[1]));
				} catch (NumberFormatException e) {
					throw new ParseException("Invalid NTAX value: " + parts[1]);
				}
			} else
				this.expectingNTaxValue = true;
		} else if (this.expectingNTaxValue) {
			this.expectingNTaxValue = false;
			try {
				((CharactersBlockListener) this.getBlockListener())
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
					this.expectingFormat = true;
					this.expectingEliminate = true;
					this.expectingTaxLabel = true;
					this.expectingCharStateLabel = true;
					this.expectingCharLabel = true;
					this.expectingStateLabel = true;
					this.expectingMatrix = true;
					try {
						((CharactersBlockListener) this.getBlockListener())
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
				this.expectingFormat = true;
				this.expectingEliminate = true;
				this.expectingTaxLabel = true;
				this.expectingCharStateLabel = true;
				this.expectingCharLabel = true;
				this.expectingStateLabel = true;
				this.expectingMatrix = true;
				try {
					((CharactersBlockListener) this.getBlockListener())
							.setDimensionsNChar(Integer.parseInt(parts[1]));
				} catch (NumberFormatException e) {
					throw new ParseException("Invalid NCHAR value: " + parts[1]);
				}
			} else
				this.expectingNCharValue = true;
		} else if (this.expectingNCharValue) {
			this.expectingNCharValue = false;
			try {
				((CharactersBlockListener) this.getBlockListener())
						.setDimensionsNChar(Integer.parseInt(token));
			} catch (NumberFormatException e) {
				throw new ParseException("Invalid NCHAR value: " + token);
			}
			this.expectingFormat = true;
			this.expectingEliminate = true;
			this.expectingTaxLabel = true;
			this.expectingCharStateLabel = true;
			this.expectingCharLabel = true;
			this.expectingStateLabel = true;
			this.expectingMatrix = true;
		}

		else if (this.expectingFormat && "FORMAT".equalsIgnoreCase(token)) {
			this.expectingFormat = false;
			this.expectingDataType = true;
			this.expectingRespectCase = true;
			this.expectingMissing = true;
			this.expectingGap = true;
			this.expectingSymbols = true;
			this.expectingEquate = true;
			this.expectingMatchChar = true;
			this.expectingLabels = true;
			this.expectingTranspose = true;
			this.expectingInterleave = true;
			this.expectingItems = true;
			this.expectingStatesFormat = true;
			this.expectingTokens = true;
		}

		else if (this.expectingDataType
				&& token.toUpperCase().startsWith("DATATYPE")) {
			this.expectingDataType = false;

			if (token.indexOf("=") >= 0) {
				final String[] parts = token.split("=");
				if (parts.length > 1) {
					this.specifiedDataType = parts[1];
					((CharactersBlockListener) this.getBlockListener())
							.setDataType(parts[1]);
				} else
					this.expectingDataTypeContent = true;
			} else
				this.expectingDataTypeEquals = true;
		}

		else if (this.expectingDataTypeEquals && token.startsWith("=")) {
			this.expectingDataTypeEquals = false;
			if (token.length() > 1) {
				token = token.substring(1);
				this.specifiedDataType = token;
				((CharactersBlockListener) this.getBlockListener())
						.setDataType(token);
			} else
				this.expectingDataTypeContent = true;
		}

		else if (this.expectingDataTypeContent) {
			this.specifiedDataType = token;
			((CharactersBlockListener) this.getBlockListener())
					.setDataType(token);
			this.expectingDataTypeContent = false;
		}

		else if (this.expectingRespectCase
				&& "RESPECTCASE".equalsIgnoreCase(token)) {
			((CharactersBlockListener) this.getBlockListener())
					.setRespectCase(true);
			this.expectingRespectCase = false;
		}

		else if (this.expectingMissing
				&& token.toUpperCase().startsWith("MISSING")) {
			this.expectingMissing = false;

			if (token.indexOf("=") >= 0) {
				final String[] parts = token.split("=");
				if (parts.length > 1)
					((CharactersBlockListener) this.getBlockListener())
							.setMissing(parts[1]);
				else
					this.expectingMissingContent = true;
			} else
				this.expectingMissingEquals = true;
		}

		else if (this.expectingMissingEquals && token.startsWith("=")) {
			this.expectingMissingEquals = false;
			if (token.length() > 1)
				((CharactersBlockListener) this.getBlockListener())
						.setMissing(token.substring(1));
			else
				this.expectingMissingContent = true;
		}

		else if (this.expectingMissingContent) {
			((CharactersBlockListener) this.getBlockListener())
					.setMissing(token);
			this.expectingMissingContent = false;
		}

		else if (this.expectingGap && token.toUpperCase().startsWith("GAP")) {
			this.expectingGap = false;

			if (token.indexOf("=") >= 0) {
				final String[] parts = token.split("=");
				if (parts.length > 1)
					((CharactersBlockListener) this.getBlockListener())
							.setGap(parts[1]);
				else
					this.expectingGapContent = true;
			} else
				this.expectingGapEquals = true;
		}

		else if (this.expectingGapEquals && token.startsWith("=")) {
			this.expectingGapEquals = false;
			if (token.length() > 1)
				((CharactersBlockListener) this.getBlockListener())
						.setGap(token.substring(1));
			else
				this.expectingGapContent = true;
		}

		else if (this.expectingGapContent) {
			((CharactersBlockListener) this.getBlockListener()).setGap(token);
			this.expectingGapContent = false;
		}

		else if (this.expectingSymbols
				&& token.toUpperCase().startsWith("SYMBOLS")) {
			this.expectingSymbols = false;

			if (token.indexOf("=") >= 0) {
				final String[] parts = token.split("=");
				if (parts.length > 1) {
					if (!parts[1].startsWith("\""))
						throw new ParseException(
								"Symbols string must start with '\"'");
					parts[1] = parts[1].substring(1);
					this.expectingSymbolsContent = true;
					if (parts[1].endsWith("\"")) {
						parts[1] = parts[1].substring(0, parts[1].length() - 1);
						this.expectingSymbolsContent = false;
					}
					for (int i = 0; i < parts[1].length(); i++)
						((CharactersBlockListener) this.getBlockListener())
								.addSymbol("" + parts[1].charAt(i));
				} else
					this.expectingSymbolsContent = true;
			} else
				this.expectingSymbolsEquals = true;
		}

		else if (this.expectingSymbolsEquals && token.startsWith("=")) {
			this.expectingSymbolsEquals = false;
			if (token.length() > 1) {
				token = token.substring(1);
				if (!token.startsWith("\""))
					throw new ParseException(
							"Symbols string must start with '\"'");
				token = token.substring(1);
				this.expectingSymbolsContent = true;

				if (token.endsWith("\"")) {
					token = token.substring(0, token.length() - 1);
					this.expectingSymbolsContent = false;
				}
				for (int i = 0; i < token.length(); i++)
					((CharactersBlockListener) this.getBlockListener())
							.addSymbol("" + token.charAt(i));
			} else
				this.expectingSymbolsContent = true;
		}

		else if (this.expectingSymbolsContent) {
			if (token.startsWith("\""))
				token = token.substring(1);
			if (token.endsWith("\"")) {
				token = token.substring(0, token.length() - 1);
				this.expectingSymbolsContent = false;
			}
			if (token.equals(""))
				this.expectingSymbolsContent = !this.seenSymbol;
			else {
				for (int i = 0; i < token.length(); i++)
					((CharactersBlockListener) this.getBlockListener())
							.addSymbol("" + token.charAt(i));
				this.seenSymbol = true;
			}
		}

		else if (this.expectingEquate
				&& token.toUpperCase().startsWith("EQUATE")) {
			this.expectingEquate = false;

			if (token.indexOf("=") >= 0) {
				final String[] parts = token.split("=");
				if (parts.length > 1) {
					if (!parts[1].startsWith("\""))
						throw new ParseException(
								"Symbols string must start with '\"'");
					parts[1] = parts[1].substring(1);
					this.expectingEquateContent = true;
					if (parts[1].endsWith("\"")) {
						parts[1] = parts[1].substring(0, parts[1].length() - 1);
						this.expectingEquateContent = false;
					}
					final String[] subParts = parts[1].split("=");
					final String symbol = subParts[0];
					final StringBuffer text = new StringBuffer();
					for (int i = 1; i < subParts.length; i++) {
						if (i >= 2)
							text.append('=');
						text.append(subParts[i]);
					}
					final List symbols = new ArrayList();
					if (text.charAt(0) == '(')
						symbols.addAll(Arrays.asList(text.substring(1,
								text.length() - 2).split("")));
					else
						symbols
								.addAll(Arrays
										.asList(text.toString().split("")));
					((CharactersBlockListener) this.getBlockListener())
							.addEquate(symbol, symbols);
				} else
					this.expectingEquateContent = true;
			} else
				this.expectingEquateEquals = true;
		}

		else if (this.expectingEquateEquals && token.startsWith("=")) {
			this.expectingEquateEquals = false;
			if (token.length() > 1) {
				token = token.substring(1);
				if (!token.startsWith("\""))
					throw new ParseException(
							"Symbols string must start with '\"'");
				token = token.substring(1);
				this.expectingEquateContent = true;

				if (token.endsWith("\"")) {
					token = token.substring(0, token.length() - 1);
					this.expectingEquateContent = false;
				}
				final String[] subParts = token.split("=");
				final String symbol = subParts[0];
				final StringBuffer text = new StringBuffer();
				for (int i = 1; i < subParts.length; i++) {
					if (i >= 2)
						text.append('=');
					text.append(subParts[i]);
				}
				final List symbols = new ArrayList();
				if (text.charAt(0) == '(')
					symbols.addAll(Arrays.asList(text.substring(1,
							text.length() - 2).split("")));
				else
					symbols.addAll(Arrays.asList(text.toString().split("")));
				((CharactersBlockListener) this.getBlockListener()).addEquate(
						symbol, symbols);
			} else
				this.expectingEquateContent = true;
		}

		else if (this.expectingEquateContent) {
			if (token.startsWith("\""))
				token = token.substring(1);
			if (token.endsWith("\"")) {
				token = token.substring(0, token.length() - 1);
				this.expectingEquateContent = false;
			}
			final String[] subParts = token.split("=");
			final String symbol = subParts[0];
			final StringBuffer text = new StringBuffer();
			for (int i = 1; i < subParts.length; i++) {
				if (i >= 2)
					text.append('=');
				text.append(subParts[i]);
			}
			final List symbols = new ArrayList();
			if (text.charAt(0) == '(')
				symbols.addAll(Arrays.asList(text.substring(1,
						text.length() - 2).split("")));
			else
				symbols.addAll(Arrays.asList(text.toString().split("")));
			((CharactersBlockListener) this.getBlockListener()).addEquate(
					symbol, symbols);
		}

		else if (this.expectingMatchChar
				&& token.toUpperCase().startsWith("MATCHCHAR")) {
			this.expectingMatchChar = false;

			if (token.indexOf("=") >= 0) {
				final String[] parts = token.split("=");
				if (parts.length > 1)
					((CharactersBlockListener) this.getBlockListener())
							.setMatchChar(parts[1]);
				else
					this.expectingMatchCharContent = true;
			} else
				this.expectingMatchCharEquals = true;
		}

		else if (this.expectingMatchCharEquals && token.startsWith("=")) {
			this.expectingMatchCharEquals = false;
			if (token.length() > 1)
				((CharactersBlockListener) this.getBlockListener())
						.setMatchChar(token.substring(1));
			else
				this.expectingMatchCharContent = true;
		}

		else if (this.expectingMatchCharContent) {
			((CharactersBlockListener) this.getBlockListener())
					.setMatchChar(token);
			this.expectingMatchCharContent = false;
		}

		else if (this.expectingLabels && "LABELS".equalsIgnoreCase(token)) {
			((CharactersBlockListener) this.getBlockListener()).setLabels(true);
			this.expectingLabels = false;
		}

		else if (this.expectingLabels && "NOLABELS".equalsIgnoreCase(token)) {
			((CharactersBlockListener) this.getBlockListener())
					.setLabels(false);
			this.expectingLabels = false;
		}

		else if (this.expectingTranspose && "TRANSPOSE".equalsIgnoreCase(token)) {
			((CharactersBlockListener) this.getBlockListener())
					.setTransposed(true);
			this.expectingTranspose = false;
		}

		else if (this.expectingInterleave
				&& token.toUpperCase().startsWith("INTERLEAVE")) {
			boolean interleaved = true;
			if (token.indexOf("=") >= 0) {
				final String[] parts = token.split("=");
				if (parts.length > 1) {
					if (!("YES".equalsIgnoreCase(parts[1]) || "TRUE".equalsIgnoreCase(parts[1]))) {
						interleaved = false;
					}
				}
			}
			((CharactersBlockListener) this.getBlockListener())
					.setInterleaved(interleaved);
			this.expectingInterleave = false;
		}

		else if (this.expectingItems && token.toUpperCase().startsWith("ITEMS")) {
			this.expectingItems = false;

			if (token.indexOf("=") >= 0) {
				final String[] parts = token.split("=");
				if (parts.length > 1) {
					if (parts[1].startsWith("(")) {
						parts[1] = parts[1].substring(1);
						this.itemsInBrackets = true;
						this.expectingItemsContent = true;
					}
					if (parts[1].endsWith(")")) {
						parts[1] = parts[1].substring(0, parts[1].length() - 1);
						this.itemsInBrackets = false;
						this.expectingItemsContent = false;
					}
					((CharactersBlockListener) this.getBlockListener())
							.setStatesFormat(parts[1]);
				} else
					this.expectingItemsContent = true;
			} else
				this.expectingItemsEquals = true;
		}

		else if (this.expectingItemsEquals && token.startsWith("=")) {
			this.expectingItemsEquals = false;
			if (token.length() > 1) {
				token = token.substring(1);
				if (token.startsWith("(")) {
					token = token.substring(1);
					this.itemsInBrackets = true;
					this.expectingItemsContent = true;
				}
				if (token.endsWith(")")) {
					token = token.substring(0, token.length() - 1);
					this.itemsInBrackets = false;
					this.expectingItemsContent = false;
				}
				((CharactersBlockListener) this.getBlockListener())
						.setStatesFormat(token);
			} else
				this.expectingItemsContent = true;
		}

		else if (this.expectingItemsContent) {
			if (token.startsWith("(")) {
				token = token.substring(1);
				this.itemsInBrackets = true;
				this.expectingItemsContent = true;
			}
			if (token.endsWith(")")) {
				token = token.substring(0, token.length() - 1);
				this.itemsInBrackets = false;
				this.expectingItemsContent = false;
			}
			((CharactersBlockListener) this.getBlockListener())
					.setStatesFormat(token);
			this.expectingItemsContent = this.itemsInBrackets;
		}

		else if (this.expectingStatesFormat
				&& token.toUpperCase().startsWith("STATESFORMAT")) {
			this.expectingStatesFormat = false;

			if (token.indexOf("=") >= 0) {
				final String[] parts = token.split("=");
				if (parts.length > 1)
					((CharactersBlockListener) this.getBlockListener())
							.setStatesFormat(parts[1]);
				else
					this.expectingStatesFormatContent = true;
			} else
				this.expectingStatesFormatEquals = true;
		}

		else if (this.expectingStatesFormatEquals && token.startsWith("=")) {
			this.expectingStatesFormatEquals = false;
			if (token.length() > 1)
				((CharactersBlockListener) this.getBlockListener())
						.setStatesFormat(token.substring(1));
			else
				this.expectingStatesFormatContent = true;
		}

		else if (this.expectingStatesFormatContent) {
			((CharactersBlockListener) this.getBlockListener())
					.setStatesFormat(token);
			this.expectingStatesFormatContent = false;
		}

		else if (this.expectingTokens && "TOKENS".equalsIgnoreCase(token)) {
			((CharactersBlockListener) this.getBlockListener()).setTokens(true);
			this.expectingTokens = false;
			this.tokenizedMatrix = true;
		}

		else if (this.expectingTokens && "NOTOKENS".equalsIgnoreCase(token)) {
			((CharactersBlockListener) this.getBlockListener())
					.setTokens(false);
			this.expectingTokens = false;
			this.tokenizedMatrix = false;
		}

		else if (this.expectingEliminate && "ELIMINATE".equalsIgnoreCase(token)) {
			this.expectingFormat = false;
			this.expectingDataType = false;
			this.expectingRespectCase = false;
			this.expectingMissing = false;
			this.expectingGap = false;
			this.expectingSymbols = false;
			this.expectingEquate = false;
			this.expectingMatchChar = false;
			this.expectingLabels = false;
			this.expectingTranspose = false;
			this.expectingInterleave = false;
			this.expectingItems = false;
			this.expectingStatesFormat = false;
			this.expectingTokens = false;
			this.expectingEliminate = false;
			this.expectingEliminateRange = true;
		}

		else if (this.expectingEliminateRange) {
			final String parts[] = token.split("-");
			if (parts.length != 2)
				throw new ParseException("Eliminate range " + token
						+ " not in form X-Y");
			try {
				final int eliminateStart = Integer.parseInt(parts[0]);
				final int eliminateEnd = Integer.parseInt(parts[1]);
				((CharactersBlockListener) this.getBlockListener())
						.setEliminateStart(eliminateStart);
				((CharactersBlockListener) this.getBlockListener())
						.setEliminateEnd(eliminateEnd);
			} catch (NumberFormatException e) {
				throw new ParseException("Values in eliminate range " + token
						+ " not parseable integers");
			}
			this.expectingEliminateRange = false;
		}

		else if (this.expectingTaxLabel && "TAXLABELS".equalsIgnoreCase(token)) {
			this.expectingFormat = false;
			this.expectingDataType = false;
			this.expectingRespectCase = false;
			this.expectingMissing = false;
			this.expectingGap = false;
			this.expectingSymbols = false;
			this.expectingEquate = false;
			this.expectingMatchChar = false;
			this.expectingLabels = false;
			this.expectingTranspose = false;
			this.expectingInterleave = false;
			this.expectingItems = false;
			this.expectingStatesFormat = false;
			this.expectingTokens = false;
			this.expectingEliminate = false;
			this.expectingEliminateRange = false;
			this.expectingTaxLabel = false;
			this.expectingTaxLabelValue = true;
		}

		else if (this.expectingCharStateLabel
				&& "CHARSTATELABELS".equalsIgnoreCase(token)) {
			this.expectingFormat = false;
			this.expectingDataType = false;
			this.expectingRespectCase = false;
			this.expectingMissing = false;
			this.expectingGap = false;
			this.expectingSymbols = false;
			this.expectingEquate = false;
			this.expectingMatchChar = false;
			this.expectingLabels = false;
			this.expectingTranspose = false;
			this.expectingInterleave = false;
			this.expectingItems = false;
			this.expectingStatesFormat = false;
			this.expectingTokens = false;
			this.expectingEliminate = false;
			this.expectingEliminateRange = false;
			this.expectingTaxLabel = false;
			this.expectingTaxLabelValue = false;
			this.expectingCharStateLabel = false;
			this.expectingCharStateLabelKey = true;
		}

		else if (this.expectingCharLabel
				&& "CHARLABELS".equalsIgnoreCase(token)) {
			this.expectingFormat = false;
			this.expectingDataType = false;
			this.expectingRespectCase = false;
			this.expectingMissing = false;
			this.expectingGap = false;
			this.expectingSymbols = false;
			this.expectingEquate = false;
			this.expectingMatchChar = false;
			this.expectingLabels = false;
			this.expectingTranspose = false;
			this.expectingInterleave = false;
			this.expectingItems = false;
			this.expectingStatesFormat = false;
			this.expectingTokens = false;
			this.expectingEliminate = false;
			this.expectingEliminateRange = false;
			this.expectingTaxLabel = false;
			this.expectingTaxLabelValue = false;
			this.expectingCharStateLabel = false;
			this.expectingCharStateLabelKey = false;
			this.expectingCharStateLabelName = false;
			this.expectingCharStateLabelSynonym = false;
			this.expectingCharLabel = false;
			this.expectingCharLabelValue = true;
		}

		else if (this.expectingStateLabel
				&& "STATELABELS".equalsIgnoreCase(token)) {
			this.expectingFormat = false;
			this.expectingDataType = false;
			this.expectingRespectCase = false;
			this.expectingMissing = false;
			this.expectingGap = false;
			this.expectingSymbols = false;
			this.expectingEquate = false;
			this.expectingMatchChar = false;
			this.expectingLabels = false;
			this.expectingTranspose = false;
			this.expectingInterleave = false;
			this.expectingItems = false;
			this.expectingStatesFormat = false;
			this.expectingTokens = false;
			this.expectingEliminate = false;
			this.expectingEliminateRange = false;
			this.expectingTaxLabel = false;
			this.expectingTaxLabelValue = false;
			this.expectingCharStateLabel = false;
			this.expectingCharStateLabelKey = false;
			this.expectingCharStateLabelName = false;
			this.expectingCharStateLabelSynonym = false;
			this.expectingCharLabel = false;
			this.expectingCharLabelValue = false;
			this.expectingStateLabel = false;
			this.expectingStateLabelKey = true;
		}

		else if (this.expectingMatrix && "MATRIX".equalsIgnoreCase(token)) {
			this.expectingFormat = false;
			this.expectingDataType = false;
			this.expectingRespectCase = false;
			this.expectingMissing = false;
			this.expectingGap = false;
			this.expectingSymbols = false;
			this.expectingEquate = false;
			this.expectingMatchChar = false;
			this.expectingLabels = false;
			this.expectingTranspose = false;
			this.expectingInterleave = false;
			this.expectingItems = false;
			this.expectingStatesFormat = false;
			this.expectingTokens = false;
			this.expectingEliminate = false;
			this.expectingEliminateRange = false;
			this.expectingTaxLabel = false;
			this.expectingTaxLabelValue = false;
			this.expectingCharStateLabel = false;
			this.expectingCharStateLabelKey = false;
			this.expectingCharStateLabelName = false;
			this.expectingCharStateLabelSynonym = false;
			this.expectingCharLabel = false;
			this.expectingCharLabelValue = false;
			this.expectingStateLabel = false;
			this.expectingStateLabelKey = false;
			this.expectingStateLabelContent = false;
			this.expectingMatrix = false;
			this.expectingMatrixKey = true;
		}

		else if (this.expectingTaxLabelValue)
			// Use untoken version to preserve spaces.
			((CharactersBlockListener) this.getBlockListener())
					.addTaxLabel(token);

		else if (this.expectingCharStateLabelKey) {
			this.currentCharStateLabelKey = token;
			// Use untoken version to preserve spaces.
			((CharactersBlockListener) this.getBlockListener())
					.addCharState(token);
			this.expectingCharStateLabelKey = false;
			this.expectingCharStateLabelName = true;
		}

		else if (this.expectingCharStateLabelName) {
			String actualName = token;
			String firstSynonym = null;
			if (token.indexOf("/") >= 0) {
				actualName = token.substring(0, token.indexOf("/"));
				if (token.indexOf("/") < token.length() - 2)
					firstSynonym = token.substring(token.indexOf("/") + 1);
			}
			final boolean skipSynonyms = actualName.endsWith(",")
					|| (firstSynonym != null && firstSynonym.endsWith(","));
			if (skipSynonyms) {
				if (firstSynonym != null)
					firstSynonym = firstSynonym.substring(0, firstSynonym
							.length() - 1);
				else
					actualName = actualName.substring(0,
							actualName.length() - 1);
			}
			// Use untoken version to preserve spaces.
			((CharactersBlockListener) this.getBlockListener())
					.setCharStateLabel(this.currentCharStateLabelKey,
							actualName);
			if (firstSynonym != null)
				((CharactersBlockListener) this.getBlockListener())
						.addCharStateKeyword(this.currentCharStateLabelKey,
								token);
			this.expectingCharStateLabelName = false;
			if (!skipSynonyms)
				this.expectingCharStateLabelSynonym = true;
			else
				this.expectingCharStateLabelKey = true;
		}

		else if (this.expectingCharStateLabelSynonym) {
			if (token.startsWith("/") && token.length() > 1)
				token = token.substring(1);
			final boolean skipSynonyms = token.endsWith(",");
			if (skipSynonyms)
				token = token.substring(0, token.length() - 1);
			if (!"/".equals(token))
				// Use untoken version to preserve spaces.
				((CharactersBlockListener) this.getBlockListener())
						.addCharStateKeyword(this.currentCharStateLabelKey,
								token);
			if (skipSynonyms) {
				this.expectingCharStateLabelSynonym = false;
				this.expectingCharStateLabelKey = true;
			}
		}

		else if (this.expectingCharLabelValue)
			// Use untoken version to preserve spaces.
			((CharactersBlockListener) this.getBlockListener())
					.addCharLabel(token);

		else if (this.expectingStateLabelKey) {
			final boolean skipContent = token.endsWith(",");
			if (skipContent)
				token = token.substring(0, token.length() - 1);
			this.currentStateLabelKey = token;
			// Use untoken version to preserve spaces.
			((CharactersBlockListener) this.getBlockListener()).addState(token);
			if (!skipContent) {
				this.expectingStateLabelKey = false;
				this.expectingStateLabelContent = true;
			}
		}

		else if (this.expectingStateLabelContent) {
			final boolean skipContent = token.endsWith(",");
			if (skipContent)
				token = token.substring(0, token.length() - 1);
			// Use untoken version to preserve spaces.
			((CharactersBlockListener) this.getBlockListener()).addStateLabel(
					this.currentStateLabelKey, token);
			if (skipContent) {
				this.expectingStateLabelKey = true;
				this.expectingStateLabelContent = false;
			}
		}

		else if (this.expectingMatrixKey) {
			this.currentMatrixKey = token;
			// Use untoken version to preserve spaces.
			((CharactersBlockListener) this.getBlockListener())
					.addMatrixEntry(token);
			this.expectingMatrixKey = false;
			this.expectingMatrixContent = true;
			// Update first line info and set up stack for entry.
			if (!this.matrixStack.containsKey(token)) {
				this.matrixStack.put(token, new Stack());
				if (this.matrixPrependNulls > 0)
					for (int i = 0; i < this.matrixPrependNulls; i++)
						((CharactersBlockListener) this.getBlockListener())
								.appendMatrixData(this.currentMatrixKey, null);
			}
			if (this.matrixFirstLineKey == null)
				this.matrixFirstLineKey = this.currentMatrixKey;
		}

		else if (this.expectingMatrixContent) {
			final Stack stack = (Stack) this.matrixStack
					.get(this.currentMatrixKey);
			if ("(".equals(token)) {
				final List newList = new ArrayList();
				if (!stack.isEmpty())
					((Collection) stack.peek()).add(newList);
				else
					((CharactersBlockListener) this.getBlockListener())
							.appendMatrixData(this.currentMatrixKey, newList);
				stack.push(newList);
			} else if ("{".equals(token)) {
				final Set newSet = new LinkedHashSet();
				if (!stack.isEmpty())
					((Collection) stack.peek()).add(newSet);
				else
					((CharactersBlockListener) this.getBlockListener())
							.appendMatrixData(this.currentMatrixKey, newSet);
				stack.push(newSet);
			} else if (")".equals(token) && !stack.isEmpty()
					&& (stack.peek() instanceof List)) {
				stack.pop();
				if (stack.isEmpty()
						&& this.currentMatrixKey
								.equals(this.matrixFirstLineKey))
					this.matrixFirstLineLength++;
			} else if ("}".equals(token) && !stack.isEmpty()
					&& (stack.peek() instanceof Set)) {
				stack.pop();
				if (stack.isEmpty()
						&& this.currentMatrixKey
								.equals(this.matrixFirstLineKey))
					this.matrixFirstLineLength++;
			} else {
				final boolean reallyUseTokens = (this.tokenizedMatrix || "CONTINUOUS"
						.equals(this.specifiedDataType))
						&& !("DNA".equals(this.specifiedDataType)
								|| "RNA".equals(this.specifiedDataType) || "NUCLEOTIDE"
								.equals(this.specifiedDataType));
				if (reallyUseTokens) {
					if (!stack.isEmpty())
						((Collection) stack.peek()).add(token);
					else {
						((CharactersBlockListener) this.getBlockListener())
								.appendMatrixData(this.currentMatrixKey, token);
						if (this.currentMatrixKey
								.equals(this.matrixFirstLineKey))
							this.matrixFirstLineLength++;
					}
				} else {
					final String[] toks = token.split("");
					for (int i = 0; i < toks.length; i++) {
						final String tok = toks[i];
						if (!stack.isEmpty())
							((Collection) stack.peek()).add(tok);
						else {
							((CharactersBlockListener) this.getBlockListener())
									.appendMatrixData(this.currentMatrixKey,
											tok);
							if (this.currentMatrixKey
									.equals(this.matrixFirstLineKey))
								this.matrixFirstLineLength++;
						}
					}
				}
			}
		}

		else
			throw new ParseException("Found unexpected token " + token
					+ " in CHARACTERS block");
	}
}
