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
import java.util.Set;

import org.biojava.bio.seq.io.ParseException;

/**
 * Represents Nexus characters blocks.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public class CharactersBlock extends NexusBlock.Abstract {

	/**
	 * A constant representing the name of Characters blocks.
	 */
	public static final String CHARACTERS_BLOCK = "CHARACTERS";

	private int dimensionsNTax = 0;

	private int dimensionsNChar = 0;

	private String dataType = "STANDARD";

	private boolean respectCase = false;

	private String missing = "?";

	private String gap;

	private List symbols = new ArrayList();

	private Map equate = new LinkedHashMap(); // values are lists

	private String matchChar;

	private boolean labels = true;

	private boolean transposed = false;

	private boolean interleaved = false;

	private List items = new ArrayList();

	private String statesFormat = "STATESPRESENT";

	private boolean tokens = false;

	private int eliminateStart = 0;

	private int eliminateEnd = 0;

	private List taxLabels = new ArrayList();

	// map containing two-value arrays where second value is list
	private Map charStateLabels = new LinkedHashMap();

	private List charLabels = new ArrayList();

	private Map stateLabels = new LinkedHashMap(); // values are lists

	// values are lists, containing nested mix of strings and lists and sets
	private Map matrix = new LinkedHashMap();

	private List comments = new ArrayList();

	/**
	 * Delegates to NexusBlock.Abstract constructor using
	 * CharactersBlock.CHARACTERS_BLOCK as the name.
	 */
	public CharactersBlock() {
		this(CharactersBlock.CHARACTERS_BLOCK);
	}

	/**
	 * For the DATA block subclass.
	 * 
	 * @param replacementLabel
	 *            the different label to use.
	 */
	protected CharactersBlock(final String replacementLabel) {
		super(replacementLabel);
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

	public void setDataType(final String dataType) {
		this.dataType = dataType;
	}

	public String getDataType() {
		return this.dataType;
	}

	public void setRespectCase(final boolean respectCase) {
		this.respectCase = respectCase;
	}

	public boolean isRespectCase() {
		return this.respectCase;
	}

	public void setMissing(final String missing) {
		this.missing = missing;
	}

	public String getMissing() {
		return this.missing;
	}

	public void setGap(final String gap) {
		this.gap = gap;
	}

	public String getGap() {
		return this.gap;
	}

	public void addSymbol(final String symbol) {
		if (!this.symbols.contains(symbol))
			this.symbols.add(symbol);
	}

	public void removeSymbol(final String symbol) {
		this.symbols.remove(symbol);
	}

	public void removeAllSymbols() {
		this.symbols.clear();
	}

	public List getSymbols() {
		return this.symbols;
	}

	public void addEquate(final String symbol, final List symbols) {
		this.equate.put(symbol, symbols);
	}

	public void removeEquate(final String symbol) {
		this.equate.remove(symbol);
	}

	public Map getEquates() {
		return this.equate;
	}

	public void setMatchChar(final String matchChar) {
		this.matchChar = matchChar;
	}

	public String getMatchChar() {
		return this.matchChar;
	}

	public void setLabels(final boolean labels) {
		this.labels = labels;
	}

	public boolean isLabels() {
		return this.labels;
	}

	public void setTransposed(final boolean transposed) {
		this.transposed = transposed;
	}

	public boolean isTransposed() {
		return this.transposed;
	}

	public void setInterleaved(final boolean interleaved) {
		this.interleaved = interleaved;
	}

	public boolean isInterleaved() {
		return this.interleaved;
	}

	public void addItem(final String item) {
		if (!this.items.contains(item))
			this.items.add(item);
	}

	public void removeItem(final String item) {
		this.items.remove(item);
	}

	public void removeAllItems() {
		this.items.clear();
	}

	public List getItems() {
		return this.items;
	}

	public void setStatesFormat(final String statesFormat) {
		this.statesFormat = statesFormat;
	}

	public String getStatesFormat() {
		return this.statesFormat;
	}

	public void setTokens(final boolean tokens) {
		this.tokens = tokens;
	}

	public boolean isTokens() {
		return this.tokens;
	}

	public void setEliminateStart(final int eliminateStart) {
		this.eliminateStart = eliminateStart;
	}

	public int getEliminateStart() {
		return this.eliminateStart;
	}

	public void setEliminateEnd(final int eliminateEnd) {
		this.eliminateEnd = eliminateEnd;
	}

	public int getEliminateEnd() {
		return this.eliminateEnd;
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

	public void addCharState(final String charState) {
		this.charStateLabels.put(charState, new Object[] { null,
				new ArrayList() });
	}

	public void setCharStateLabel(final String charState, final String label) {
		if (!this.charStateLabels.containsKey(charState))
			this.addCharState(charState);
		((Object[]) this.charStateLabels.get(charState))[0] = label;
	}

	public void addCharStateKeyword(final String charState, final String keyword) {
		if (!this.charStateLabels.containsKey(charState))
			this.addCharState(charState);
		((List) ((Object[]) this.charStateLabels.get(charState))[1])
				.add(keyword);
	}

	public String getCharStateLabel(final String charState) {
		return (String) (((Object[]) this.charStateLabels.get(charState))[0]);
	}

	public List getCharStateLabelKeywords(final String charState) {
		return (List) (((Object[]) this.charStateLabels.get(charState))[1]);
	}

	public void removeCharState(final String charState) {
		this.charStateLabels.remove(charState);
	}

	public Set getAllCharStates() {
		return this.charStateLabels.keySet();
	}

	public void addCharLabel(final String charLabel) {
		this.charLabels.add(charLabel);
	}

	public void removeCharLabel(final String charLabel) {
		this.charLabels.remove(charLabel);
	}

	public boolean containsCharLabel(final String charLabel) {
		return this.charLabels.contains(charLabel);
	}

	public List getCharLabels() {
		return this.charLabels;
	}

	public void addState(final String state) {
		this.stateLabels.put(state, new ArrayList());
	}

	public void addStateLabel(final String state, final String label) {
		if (!this.stateLabels.containsKey(state))
			this.addState(state);
		((List) this.stateLabels.get(state)).add(label);
	}

	public List getStateLabels(final String state) {
		return (List) this.stateLabels.get(state);
	}

	public void removeState(final String state) {
		this.stateLabels.remove(state);
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

		writer.write(" FORMAT DATATYPE=");
		this.writeToken(writer, this.dataType);
		if (this.respectCase && "STANDARD".equals(this.dataType))
			writer.write(" RESPECTCASE");
		writer.write(" MISSING=");
		this.writeToken(writer, this.missing);
		writer.write(" GAP=");
		this.writeToken(writer, this.gap);
		writer.write(" SYMBOLS=\"");
		if (this.symbols.isEmpty()) {
			this.symbols.add("0");
			this.symbols.add("1");
		}
		for (final Iterator i = this.symbols.iterator(); i.hasNext();)
			this.writeToken(writer, (String) i.next());
		writer.write('"');
		if (!this.equate.isEmpty()) {
			writer.write(" EQUATE=\"");
			for (final Iterator i = this.equate.entrySet().iterator(); i
					.hasNext();) {
				final Map.Entry entry = (Map.Entry) i.next();
				this.writeToken(writer, "" + entry.getKey());
				writer.write("=(");
				for (final Iterator j = ((List) entry.getValue()).iterator(); j
						.hasNext();)
					this.writeToken(writer, "" + j.next());
				writer.write(')');
				if (i.hasNext())
					writer.write(' ');
			}
			writer.write('"');
		}
		if (this.matchChar != null) {
			writer.write(" MATCHCHAR=");
			this.writeToken(writer, this.matchChar);
		}
		writer.write(this.labels ? " LABELS" : " NOLABELS");
		if (this.transposed)
			writer.write(" TRANSPOSED");
		// FIXME Output files, for now, are never interleaved.
		// if (this.interleaved)
		// writer.write(" INTERLEAVED");
		writer.write(" ITEMS=");
		if (this.items.isEmpty())
			this.items.add("STATES");
		if (this.items.size() > 1)
			writer.write('(');
		for (final Iterator i = this.items.iterator(); i.hasNext();) {
			this.writeToken(writer, "" + i.next());
			if (i.hasNext())
				writer.write(' ');
		}
		if (this.items.size() > 1)
			writer.write(')');
		writer.write(" STATESFORMAT=");
		this.writeToken(writer, this.statesFormat);
		final boolean reallyUseTokens = (this.tokens || "CONTINUOUS"
				.equals(this.dataType))
				&& !("DNA".equals(this.dataType) || "RNA".equals(this.dataType) || "NUCLEOTIDE"
						.equals(this.dataType));
		writer.write(reallyUseTokens ? " TOKENS" : " NOTOKENS");
		writer.write(";" + NexusFileFormat.NEW_LINE);

		if (this.eliminateStart > 0 && this.eliminateEnd > 0) {
			writer.write(" ELIMINATE " + this.eliminateStart + "-"
					+ this.eliminateEnd);
			writer.write(";" + NexusFileFormat.NEW_LINE);
		}

		if (this.taxLabels.size() > 0) {
			writer.write(" TAXLABELS");
			for (final Iterator i = this.taxLabels.iterator(); i.hasNext();) {
				writer.write(' ');
				this.writeToken(writer, (String) i.next());
			}
			writer.write(";" + NexusFileFormat.NEW_LINE);
		}

		if (!this.charStateLabels.isEmpty()
				&& !"CONTINUOUS".equals(this.dataType)) {
			writer.write(" CHARSTATELABELS" + NexusFileFormat.NEW_LINE);
			for (final Iterator i = this.charStateLabels.entrySet().iterator(); i
					.hasNext();) {
				final Map.Entry topEntry = (Map.Entry) i.next();
				writer.write('\t');
				this.writeToken(writer, "" + topEntry.getKey());
				writer.write('\t');
				final Object[] parts = (Object[]) topEntry.getValue();
				this.writeToken(writer, "" + parts[0]);
				final List names = (List) parts[1];
				if (!names.isEmpty()) {
					writer.write('/');
					for (final Iterator k = names.iterator(); k.hasNext();) {
						this.writeToken(writer, "" + k.next());
						if (k.hasNext())
							writer.write(' ');
					}
				}
				if (!i.hasNext())
					writer.write(';');
				else
					writer.write(',');
				writer.write(NexusFileFormat.NEW_LINE);
			}
		}

		if (!this.charLabels.isEmpty() && !this.transposed) {
			writer.write(" CHARLABELS" + NexusFileFormat.NEW_LINE);
			writer.write('\t');
			for (final Iterator i = this.charLabels.iterator(); i.hasNext();) {
				this.writeToken(writer, "" + i.next());
				if (i.hasNext())
					writer.write(' ');
			}
			writer.write(";" + NexusFileFormat.NEW_LINE);
		}

		if (!this.stateLabels.isEmpty() && !"CONTINUOUS".equals(this.dataType)) {
			writer.write(" STATELABELS" + NexusFileFormat.NEW_LINE);
			for (final Iterator i = this.stateLabels.entrySet().iterator(); i
					.hasNext();) {
				final Map.Entry topEntry = (Map.Entry) i.next();
				writer.write('\t');
				this.writeToken(writer, "" + topEntry.getKey());
				writer.write('\t');
				final List names = (List) topEntry.getValue();
				for (final Iterator k = names.iterator(); k.hasNext();) {
					this.writeToken(writer, "" + k.next());
					if (k.hasNext())
						writer.write(' ');
				}
				if (!i.hasNext())
					writer.write(';');
				else
					writer.write(',');
				writer.write(NexusFileFormat.NEW_LINE);
			}
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
				this.writeMatrixEntry(writer, j.next(), reallyUseTokens);
				if (reallyUseTokens && j.hasNext())
					writer.write(' ');
			}
			writer.write(NexusFileFormat.NEW_LINE);
		}
		writer.write(";" + NexusFileFormat.NEW_LINE);

	}

	private void writeMatrixEntry(final Writer writer, final Object obj,
			final boolean reallyUseTokens) throws IOException {
		if (obj == null)
			this.writeToken(writer, this.missing);
		else if (obj instanceof String)
			this.writeToken(writer, (String) obj);
		else if (obj instanceof List) {
			writer.write('(');
			for (final Iterator k = ((List) obj).iterator(); k.hasNext();) {
				this.writeMatrixEntry(writer, k.next(), reallyUseTokens);
				if (k.hasNext() && reallyUseTokens)
					writer.write(' ');
			}
			writer.write(')');
		} else if (obj instanceof Set) {
			writer.write('{');
			for (final Iterator k = ((Set) obj).iterator(); k.hasNext();) {
				this.writeMatrixEntry(writer, k.next(), reallyUseTokens);
				if (k.hasNext() && reallyUseTokens)
					writer.write(' ');
			}
			writer.write('}');
		}
	}
}
