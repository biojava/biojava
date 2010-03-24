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

package org.biojava.bio.alignment;

import java.io.Serializable;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.BioError;
import org.biojava.bio.symbol.AbstractSymbolList;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * A simple implementation of an Alignment.
 * <p>
 * This is a simple-stupid implementation that is made from a set of
 * same-lengthed SymbolList objects each with an associated label. It does not
 * handle differently lengthed sequences and doesn't contain any gap-editing
 * concepts.
 * 
 * @author Matthew Pocock
 * @author Greg Cox
 * @author Nimesh Singh
 */
public class SimpleAlignment extends AbstractSymbolList implements Alignment,
		Serializable {
	private static final long serialVersionUID = -1760075176220928440L;

	/**
	 * 
	 */
	private LinkedHashMap<String, SymbolList> labelToSymbolList;
	/**
	 * 
	 */
	private List<String> labels;
	/**
	 * 
	 */
	private Alphabet alphabet;
	/**
	 * 
	 */
	private int length;
	/**
	 * 
	 */
	private int score;

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#finalize()
	 */
	protected void finalize() throws Throwable {
		super.finalize();
		// System.err.println("Finalizing a SimpleAlignement");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.biojava.bio.symbol.SymbolList#length()
	 */
	public int length() {
		return length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.biojava.bio.symbol.SymbolList#getAlphabet()
	 */
	public Alphabet getAlphabet() {
		return alphabet;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.biojava.bio.symbol.SymbolList#symbolAt(int)
	 */
	public Symbol symbolAt(int index) {
		try {
			if (labels.size() == 1) {
				return symbolAt(labels.get(0), index);
			} else {
				return alphabet.getSymbol(new ColAsList<Symbol>(index));
			}
		} catch (IllegalSymbolException ire) {
			throw new BioError(

			"Somehow my crossproduct alphabet is incompatible with column "
					+ index, ire);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.biojava.bio.alignment.Alignment#getLabels()
	 */
	public List<String> getLabels() {
		return labels;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.biojava.bio.alignment.Alignment#symbolAt(java.lang.Object, int)
	 */
	public Symbol symbolAt(String label, int column) {
		return symbolListForLabel(label).symbolAt(column);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.biojava.bio.alignment.Alignment#subAlignment(java.util.Set,
	 * org.biojava.bio.symbol.Location)
	 */
	public Alignment subAlignment(Set<String> labels, Location loc)
			throws NoSuchElementException {
		Map<String, SymbolList> labelsToResList = new LinkedHashMap<String, SymbolList>();
		Iterator<String> i;
		if (labels != null) {
			i = labels.iterator();
		} else {
			i = getLabels().iterator();
		}
		while (i.hasNext()) {
			String label = i.next();
			SymbolList sym = symbolListForLabel(label);
			if (loc != null) {
				sym = loc.symbols(sym);
			}
			labelsToResList.put(label, sym);
		}
		return new SimpleAlignment(labelsToResList);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.biojava.bio.alignment.Alignment#symbolListForLabel(java.lang.String)
	 */
	public SymbolList symbolListForLabel(String label)
			throws NoSuchElementException {
		SymbolList rl = labelToSymbolList.get(label);
		if (rl == null) {
			throw new NoSuchElementException(
					"No symbol list associated with label " + label);
		}
		return rl;
	}

	/**
	 * Generate an alignment from a list of SymbolLists.
	 * <p>
	 * The SymbolLists must all be of the same length.
	 * 
	 * @param labelToResList
	 *            the label-to-symbol list mapping
	 * @throws IllegalArgumentException
	 *             if the SymbolLists are not the same length
	 */
	public SimpleAlignment(Map<String, SymbolList> labelToResList)
			throws IllegalArgumentException {
		if (labelToResList.isEmpty()) {
			throw new IllegalArgumentException(
					"Can't create an alignment with no sequences");
		}

		this.labels = Collections.unmodifiableList(new ArrayList<String>(
				labelToResList.keySet()));
		this.labelToSymbolList = new LinkedHashMap<String, SymbolList>(
				labelToResList);

		int length = -1;
		List<Alphabet> alphaList = new ArrayList<Alphabet>();
		for (Iterator<String> li = labels.iterator(); li.hasNext();) {
			String label = li.next();
			try {
				SymbolList rl = symbolListForLabel(label);
				alphaList.add(rl.getAlphabet());
				if (length == -1) {
					length = rl.length();
				} else {
					if (rl.length() != length) {
						StringBuffer sb = new StringBuffer();
						for (Iterator<String> labI = labels.iterator(); labI
								.hasNext();) {
							String lab = labI.next();
							sb.append("\n\t" + lab + " ("
									+ symbolListForLabel(lab).length() + ")");
						}
						throw new IllegalArgumentException(
								"All SymbolLists must be the same length: "
										+ sb.substring(0));
					}
				}
			} catch (NoSuchElementException nsee) {
				if (labelToSymbolList.containsKey(label)) {
					throw new IllegalArgumentException(
							"The symbol list associated with " + label
									+ " is null");
				} else {
					throw new BioError(
							"Something is screwey - map is lying about key/values",
							nsee);
				}
			}
		}

		this.alphabet = AlphabetManager.getCrossProductAlphabet(alphaList);
		this.length = length;
	}

	/**
	 * 
	 */
	public Iterator<SymbolList> symbolListIterator() {
		return new Alignment.SymbolListIterator(this);
	}

	/**
	 * Makes a column of the alignment behave like a list.
	 * 
	 * @author Matthew Pocock
	 */
	private final class ColAsList<T extends Symbol> extends AbstractList<T>
			implements Serializable {
		/**
		 * Generated serial version identifier.
		 */
		private static final long serialVersionUID = 8254702569039851040L;
		private final int col;

		public ColAsList(int col) {
			this.col = col;
		}

		// protected ColAsList() {
		// this.col = 0;
		// }

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.AbstractList#get(int)
		 */
		@SuppressWarnings("unchecked")
		public T get(int indx) {
			return (T) symbolAt(labels.get(indx), col);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.AbstractCollection#size()
		 */
		public int size() {
			return labels.size();
		}
	}

	/**
	 * 
	 * @return
	 */
	public int getScore() {
		return score;
	}

	/**
	 * 
	 * @param score
	 */
	public void setScore(int score) {
		this.score = score;
	}

}
