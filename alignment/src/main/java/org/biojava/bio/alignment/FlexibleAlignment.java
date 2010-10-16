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

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.TreeSet;
import java.util.Vector;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.GappedSymbolList;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.SimpleGappedSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;

/**
 * <p>
 * FlexibleAlignment is a class which implements UnequalLengthAlignment,
 * ARAlignment and EditableAlignment <b>It places no restriction on where any
 * sequence can be in the alignment so there could be gaps in the alignment. You
 * tell it where to put the sequence, it will do it. I think I will be adding an
 * Exception NonContinuousAlignmentException. STILL UNDER CONSTRUCTION.
 * seqString does not work because there it does not seem to support
 * tokenization 'token' this is true for SimpleAlignment too.</b>
 * 
 * @author David Waring
 * @author Matthew Pocock
 */
public class FlexibleAlignment extends AbstractULAlignment implements
		ARAlignment, EditableAlignment {

	protected Map<Object, AlignmentElement> data;
	protected List<String> labelOrder;
	protected Location alignmentRange;
	List<Alphabet> alphaList = new ArrayList<Alphabet>();

	/**
	 * construct this object with the reference sequence which can either be a
	 * gappedSymbolList or not label in all cases refers to an object that holds
	 * the display name (generally just a String). since more than one sequence
	 * in an alignment could have the same name this works as long as the labels
	 * are different objects even though they may hold the same name.
	 */

	public FlexibleAlignment(List<AlignmentElement> seqList)
			throws BioException {
		data = new Hashtable<Object, AlignmentElement>();
		labelOrder = new Vector<String>();
		alignmentRange = new RangeLocation(1, 1);

		int k = 0;
		// go through the list make sure that all seqs are GappedSymbolLists
		for (Iterator<AlignmentElement> i = seqList.iterator(); i.hasNext();) {
			AlignmentElement ae = i.next();
			String label = ae.getLabel();
			Location loc = ae.getLoc();
			SymbolList seq = ae.getSymbolList();
			alphaList.add(seq.getAlphabet());
			if (!(seq instanceof GappedSymbolList)) {
				seq = new SimpleGappedSymbolList(seq);
				ae = new SimpleAlignmentElement(label, seq, loc);
			}
			data.put(label, ae);
			labelOrder.add(label);
			int min = lesser(alignmentRange.getMin(), loc.getMin());
			int max = greater(alignmentRange.getMax(), loc.getMax());
			alignmentRange = new RangeLocation(min, max);
			k++;
		}
		this.alphabet = AlphabetManager.getCrossProductAlphabet(alphaList);
		try {
			resetRange();
		} catch (ChangeVetoException e) {
			throw new BioError("Should not have a problem here");
		}

	}

	private int getOrder(Object label) throws Exception {
		for (int i = 0; i < labelOrder.size(); i++) {
			if (labelOrder.get(i).equals(label))
				return i;
		}
		throw new Exception("did not find label");
	}

	/**
	 * add a new a alignment usings a location to the reference sequence. This
	 * should either contain no gaps or it should be relative to a reference
	 * sequence that already has the gaps added
	 */

	public synchronized void addSequence(AlignmentElement ae)
			throws ChangeVetoException, BioException {
		ChangeSupport cs;
		ChangeEvent cevt;

		String label = ae.getLabel();
		SymbolList seq = ae.getSymbolList();
		Location loc = ae.getLoc();

		// give the listeners a change to veto this
		// create a new change event ->the EDIT is a static final variable of
		// type ChangeType in SymbolList interface
		cevt = new ChangeEvent(this, ARAlignment.ADD_LABEL, label);
		cs = getChangeSupport(ARAlignment.ADD_LABEL);

		// let the listeners know what we want to do
		cs.firePreChangeEvent(cevt);
		if (!(seq instanceof GappedSymbolList)) {
			seq = new SimpleGappedSymbolList(seq);
			ae = new SimpleAlignmentElement(label, seq, loc);
		}
		data.put(label, ae);
		labelOrder.add(label);
		alphaList.add(seq.getAlphabet());
		this.alphabet = AlphabetManager.getCrossProductAlphabet(alphaList);

		int min = lesser(alignmentRange.getMin(), loc.getMin());
		int max = greater(alignmentRange.getMax(), loc.getMax());
		alignmentRange = new RangeLocation(min, max);
		resetRange();

		cs.firePostChangeEvent(cevt);

	}

	public synchronized void removeSequence(Object label)
			throws ChangeVetoException {

		ChangeSupport cs;
		ChangeEvent cevt;

		// give the listeners a change to veto this
		// create a new change event ->the EDIT is a static final variable of
		// type ChangeType in SymbolList interface
		cevt = new ChangeEvent(this, ARAlignment.REMOVE_LABEL, label);
		cs = getChangeSupport(ARAlignment.REMOVE_LABEL);

		// let the listeners know what we want to do
		cs.firePreChangeEvent(cevt);
		try {
			alphaList.remove(getOrder(label));
			this.alphabet = AlphabetManager.getCrossProductAlphabet(alphaList);
		} catch (Throwable e) {
			e.printStackTrace();
		}
		data.remove(label);
		labelOrder.remove(label);
		resetRange();
		cs.firePostChangeEvent(cevt);

	}

	// ///////////////////////
	// methods from Interface UnequalLengthAlignment
	// ///////////////////////

	/**
	 * The location of an individual SymbolList relative to overall Alignment
	 */
	public Location locInAlignment(Object label) throws NoSuchElementException {
		return getAE(label).getLoc();
	}

	public List<Object> getLabelsAt(int column)
			throws IndexOutOfBoundsException {
		if (column < 1 || column > this.length())
			throw new IndexOutOfBoundsException();
		List<Object> labelList = new ArrayList<Object>();
		Location loc;
		Object label;
		for (Iterator<Object> labelIterator = data.keySet().iterator(); labelIterator
				.hasNext();) {
			label = labelIterator.next();
			loc = getAE(label).getLoc();
			if (loc.contains(column))
				labelList.add(label);
		}
		return labelList;
	}

	// ///////////////////////
	// methods from Interface Alignment
	// //////////////////////

	public synchronized int length() {
		return alignmentRange.getMax() - alignmentRange.getMin() + 1;
	}

	public Alphabet getAlphabet() {
		return alphabet;
	}

	/**
	 * getLabels will return a list of labels in left to right order
	 */

	public synchronized List<String> getLabels() {
		TreeSet<String> sorted = new TreeSet<String>(
				new LeftRightLocationComparator<String>());
		sorted.addAll(labelOrder);
		return new Vector<String>(sorted);
	}

	/**
	 * This gets the symbol for an individual sequence at position in the
	 * overall alignment If the sequence is not aligned at that location it
	 * returns null
	 */

	public synchronized Symbol symbolAt(String label, int column)
			throws NoSuchElementException, IndexOutOfBoundsException {

		SymbolList seq = symbolListForLabel(label);
		int cloc = posInSeq(label, column);
		Symbol symbol = null;
		// debug (label.toString() + " " + column + ":" + cloc);
		if (seq == null) {
			// debug("seq is null");
		}

		try {
			symbol = seq.symbolAt(cloc);
		} catch (IndexOutOfBoundsException e) {
			// leave symbol == null
		}
		return symbol;
	}

	/**
	 * 
	 * @param label
	 * @return
	 * @throws NoSuchElementException
	 */
	public synchronized SymbolList symbolListForLabel(String label)
			throws NoSuchElementException {
		return getAE(label).getSymbolList();
	}

	// methods from interface EditableAlignment

	public synchronized void edit(Object label, Edit edit)
			throws ChangeVetoException {
		throw new BioError("Not implemented yet");
	}

	/**
	 * loc in this case is the Alignment Location
	 */
	public synchronized void shiftAtAlignmentLoc(Object label, Location loc,
			int offset) throws ChangeVetoException,
			IllegalAlignmentEditException, IndexOutOfBoundsException {

		Location sourceLoc = locInSeq(label, loc);
		shiftAtSequenceLoc(label, sourceLoc, offset);

	}

	/**
	 * loc in this case is the SymbolList Location
	 */
	public synchronized void shiftAtSequenceLoc(Object label, Location loc,
			int offset) throws ChangeVetoException,
			IllegalAlignmentEditException, IndexOutOfBoundsException {

		ChangeSupport csgap;
		ChangeEvent cegap;
		ChangeSupport csloc;
		ChangeEvent celoc;
		celoc = new ChangeEvent(this, EditableAlignment.LOCATION, label);
		csloc = getChangeSupport(EditableAlignment.LOCATION);
		cegap = new ChangeEvent(this, EditableAlignment.GAPS, label);
		csgap = getChangeSupport(EditableAlignment.GAPS);

		int caseValue = 0;
		int absOffset = Math.abs(offset);
		Location seqLoc = locInAlignment(label);
		AlignmentElement ae = getAE(label);
		SymbolList seq = ae.getSymbolList();
		Location newLoc;
		int min = loc.getMin();
		int max = loc.getMax();
		if (min < 1 || max > seq.length()) {
			throw new IndexOutOfBoundsException();
		}
		if (offset == 0) {
			return;
		}
		if (offset > 1)
			caseValue += 1;
		if (min == 1)
			caseValue += 2;
		if (max == seq.length())
			caseValue += 4;

		switch (caseValue) {

		case 0: // internal shift to left
			if (!allGaps(seq, min + offset, min - 1)) {
				throw new IllegalAlignmentEditException();
			}
			csgap.firePreChangeEvent(cegap);

			((GappedSymbolList) seq).addGapsInView(max + 1, absOffset);
			removeGaps((GappedSymbolList) seq, min - absOffset, absOffset);
			csgap.firePostChangeEvent(cegap);
			break;

		case 1: // internal shift to right
			if (!allGaps(seq, max + 1, max + offset)) {
				throw new IllegalAlignmentEditException();
			}
			csgap.firePreChangeEvent(cegap);
			removeGaps((GappedSymbolList) seq, max + 1, offset);
			((GappedSymbolList) seq).addGapsInView(min, offset);
			csgap.firePostChangeEvent(cegap);
			break;

		case 2: // left end shift to left
			csgap.firePreChangeEvent(cegap);
			csloc.firePreChangeEvent(celoc);
			((GappedSymbolList) seq).addGapsInView(max + 1, absOffset);
			newLoc = new RangeLocation(seqLoc.getMin() - absOffset, seqLoc
					.getMax());
			ae.setLoc(newLoc);
			resetRange();
			csloc.firePostChangeEvent(celoc);
			csgap.firePostChangeEvent(cegap);
			break;

		case 3: // left end shift to right
			if (!allGaps(seq, max + 1, max + offset)) {
				throw new IllegalAlignmentEditException();
			}
			csgap.firePreChangeEvent(cegap);
			csloc.firePreChangeEvent(celoc);
			removeGaps((GappedSymbolList) seq, max + 1, offset);
			newLoc = new RangeLocation(seqLoc.getMin() + offset, seqLoc
					.getMax());
			ae.setLoc(newLoc);
			resetRange();
			csloc.firePostChangeEvent(celoc);
			csgap.firePostChangeEvent(cegap);
			break;

		case 4: // right end shift to left
			if (!allGaps(seq, min - absOffset, min - 1)) {
				throw new IllegalAlignmentEditException();
			}
			csgap.firePreChangeEvent(cegap);
			csloc.firePreChangeEvent(celoc);
			removeGaps((GappedSymbolList) seq, min - absOffset, absOffset);
			newLoc = new RangeLocation(seqLoc.getMin(), seqLoc.getMax()
					+ offset);
			ae.setLoc(newLoc);
			resetRange();
			csloc.firePostChangeEvent(celoc);
			csgap.firePostChangeEvent(cegap);
			break;

		case 5: // right end shift to right
			csgap.firePreChangeEvent(cegap);
			csloc.firePreChangeEvent(celoc);
			((GappedSymbolList) seq).addGapsInView(min, offset);
			newLoc = new RangeLocation(seqLoc.getMin(), seqLoc.getMax()
					+ offset);
			ae.setLoc(newLoc);
			resetRange();
			csloc.firePostChangeEvent(celoc);
			csgap.firePostChangeEvent(cegap);
			break;

		case 6: // whole seq shift to left
			debug("Shifting all to left " + absOffset);
			shift(label, offset);
			break;

		case 7: // whole seq shift to right
			debug("Shifting all to right " + absOffset);
			shift(label, offset);
			break;

		default:
			debug("OOOPS something is wrong " + loc.toString() + " "
					+ absOffset);
			return;
		}
	}

	/**
	 * because there is a bug in GappedSymbolList
	 */

	public synchronized void removeGaps(GappedSymbolList seq, int start,
			int length) {
		try {
			// seq.removeGaps (start , length);
			// because there is a bug in GappedSymbolList we do it one at a time
			for (int i = 1; i <= length; i++) {
				seq.removeGap(start);
			}
		} catch (IllegalSymbolException e) {
			throw new BioError("We should have tested for this already");
		}
	}

	/**
	 * make sure that all Symbols in this range are gaps
	 */

	protected synchronized boolean allGaps(SymbolList seq, int start, int end) {

		Symbol gs = seq.getAlphabet().getGapSymbol();
		for (int i = start; i <= end; i++) {
			if (!(seq.symbolAt(i).equals(gs))) {
				return false;
			}
		}
		return true;
	}

	/**
	 * check that begining is at 1 otherwise shift everything over
	 */
	protected synchronized void resetRange() throws ChangeVetoException {

		int min = 0;// just for the compiler
		int max = 0;// just for the compiler
		int lMin;
		int lMax;
		int count = 1;
		// get the current range from all labels
		for (Iterator<String> i = getLabels().iterator(); i.hasNext();) {
			Object label = i.next();
			lMin = locInAlignment(label).getMin();
			lMax = locInAlignment(label).getMax();
			if (count == 1) {
				min = lMin;
			} else {
				min = lesser(min, lMin);
			}
			if (count == 1) {
				max = lMax;
			} else {
				max = greater(max, lMax);
			}
			count++;
		}
		alignmentRange = new RangeLocation(min, max);
		if (min != 1) {
			int offset = 1 - alignmentRange.getMin();
			shiftAll(offset);
			alignmentRange = new RangeLocation(
					alignmentRange.getMin() + offset, alignmentRange.getMax()
							+ offset);
		}
	}

	protected synchronized void shiftAll(int offset) throws ChangeVetoException {
		List<String> lList = getLabels();
		for (Iterator<String> i = lList.iterator(); i.hasNext();) {
			Object label = i.next();
			shift(label, offset);
		}
	}

	/**
	 * moves the whole sequence
	 */

	protected synchronized void shift(Object label, int offset)
			throws ChangeVetoException {
		ChangeSupport csloc;
		ChangeEvent celoc;
		celoc = new ChangeEvent(this, EditableAlignment.LOCATION, label);
		csloc = getChangeSupport(EditableAlignment.LOCATION);
		Location oLoc = locInAlignment(label);
		Location nLoc = new RangeLocation(oLoc.getMin() + offset, oLoc.getMax()
				+ offset);
		csloc.firePreChangeEvent(celoc);
		debug("shifting " + label.toString());
		getAE(label).setLoc(nLoc);
		resetRange();
		debug("shifted " + label);
		csloc.firePostChangeEvent(celoc);
	}

	// utility methods
	protected int greater(int x, int y) {
		int greatest = (x > y) ? x : y;
		return greatest;
	}

	protected int lesser(int x, int y) {
		int least = (x < y) ? x : y;
		return least;
	}

	protected AlignmentElement getAE(Object label)
			throws NoSuchElementException {
		if (!(data.containsKey(label)))
			;
		return data.get(label);
	}

	/**
	 * get the position in the sequence corresponding to the postion within the
	 * alignment
	 */

	protected synchronized int posInSeq(Object label, int column)
			throws NoSuchElementException, IndexOutOfBoundsException {
		if (column < 1 || column > this.length()) {
			throw new IndexOutOfBoundsException();
		}
		Location loc = locInAlignment(label);
		return (column - loc.getMin() + 1);
	}

	protected synchronized Location locInSeq(Object label, Location viewLoc)
			throws NoSuchElementException, IndexOutOfBoundsException {
		int min = posInSeq(label, viewLoc.getMin());
		int max = posInSeq(label, viewLoc.getMax());
		return new RangeLocation(min, max);
	}
}
