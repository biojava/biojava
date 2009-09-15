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

package org.biojava.bio.symbol;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.BioError;
import org.biojava.utils.AssertionFailure;

/**
 * This implementation of GappedSymbolList wraps a SymbolList, allowing you to
 * insert gaps. Please note that this is a view onto another SymbolList. Gaps
 * created and removed are only in the view not the underlying original. This
 * means that any gaps present in the original cannot be manipulated in this
 * view. To manipulate the original you would need to use Edit objects.
 * 
 * @author Matthew Pocock
 * @since 1.3
 */
public class SimpleGappedSymbolList extends AbstractSymbolList implements
		GappedSymbolList, Serializable {
	/**
	 * An aligned block.
	 * <p>
	 * The alignment is actualy stoored as a list of these objects. Each block
	 * is contiguous with the next in the source fields, but may be gapped in
	 * the view fields.
	 * 
	 * @author Matthew Pocock
	 */
	protected static final class Block implements Serializable {
		private static final long serialVersionUID = 4090888450309921439L;

		public int sourceStart;
		public int sourceEnd;
		public int viewStart;
		public int viewEnd;

		public Block(Block block) {
			this.sourceStart = block.sourceStart;
			this.sourceEnd = block.sourceEnd;
			this.viewStart = block.viewStart;
			this.viewEnd = block.viewEnd;

			// fixme: should be using 1.4 assertion syntax
			// assert isSane() : "Block is a sane shape: " + this.toString();
			assert isSane() : "Block is a sane shape: " + this.toString();
		}

		public Block(int sourceStart, int sourceEnd, int viewStart, int viewEnd) {
			this.sourceStart = sourceStart;
			this.sourceEnd = sourceEnd;
			this.viewStart = viewStart;
			this.viewEnd = viewEnd;
		}

		public boolean isSane() {
			return (viewEnd - viewStart) == (sourceEnd - sourceStart);
		}

		public String toString() {
			return "Block: source=" + sourceStart + "," + sourceEnd + " view="
					+ viewStart + "," + viewEnd;
		}
	}

	private static final long serialVersionUID = 4258621048300499709L;

	/**
	 * The Alphabet - the same as source but guaranteed to include the gap
	 * character.
	 */
	private final Alphabet alpha;

	/**
	 * The SymbolList to view.
	 */
	private final SymbolList source;

	/**
	 * The list of ungapped blocks that align between source and this view.
	 */
	private final List<Block> blocks;

	/**
	 * The total length of the alignment - necessary to allow leading & trailing
	 * gaps.
	 */
	private int length;

	/**
	 * Create a new SimpleGappedSymbolList that will view source, inheriting all
	 * existing gaps.
	 * 
	 * @param gappedSource
	 *            the underlying sequence
	 */
	public SimpleGappedSymbolList(GappedSymbolList gappedSource) {
		this.source = gappedSource.getSourceSymbolList();
		this.alpha = this.source.getAlphabet();
		this.blocks = new ArrayList<Block>();
		this.length = this.source.length();
		Block b = new Block(1, length, 1, length);
		blocks.add(b);
		int n = 1;
		for (int i = 1; i <= this.length(); i++)
			if (this.alpha.getGapSymbol().equals(gappedSource.symbolAt(i)))
				this.addGapInSource(n);
			else
				n++;
	}

	/**
	 * Create a new SimpleGappedSymbolList that will view source.
	 * 
	 * @param source
	 *            the underlying sequence
	 */
	public SimpleGappedSymbolList(SymbolList source) {
		this.source = source;
		this.alpha = source.getAlphabet();
		this.blocks = new ArrayList<Block>();
		this.length = source.length();
		Block b = new Block(1, length, 1, length);
		blocks.add(b);
	}
	
	public SimpleGappedSymbolList(Alphabet alpha) {
		this.alpha = alpha;
		this.source = new SimpleSymbolList(this.alpha);
		this.blocks = new ArrayList<Block>();
		this.length = source.length();
		// Block b = new Block(1, length, 1, length);
		// blocks.add(b);
	}

	public void addGapInSource(int pos) throws IndexOutOfBoundsException {
		addGapsInSource(pos, 1);
	}

	public void addGapInView(int pos) throws IndexOutOfBoundsException {
		addGapsInView(pos, 1);
	}

	public void addGapsInSource(int pos, int length) {
		if (pos < 1 || pos > (length() + 1)) {
			throw new IndexOutOfBoundsException(
					"Attempted to add a gap outside of this sequence (1.."
							+ length() + ") at " + pos);
		}

		int i = blocks.size() / 2;
		int imin = 0;
		int imax = blocks.size() - 1;

		do {
			Block b = blocks.get(i);
			if (b.sourceStart < pos && b.sourceEnd >= pos) { // found block that
				// need
				// splitting
				// with gaps
				Block c = new Block(pos, b.sourceEnd, sourceToView(b, pos),
						b.viewEnd);
				b.viewEnd = sourceToView(b, pos - 1);
				b.sourceEnd = pos - 1;
				blocks.add(i + 1, c);
				renumber(i + 1, length);
				return;
			} else {
				if (b.sourceStart < pos) {
					imin = i + 1;
					i = imin + (imax - imin) / 2;
				} else {
					imax = i - 1;
					i = imin + (imax - imin) / 2;
				}
			}
		} while (imin <= imax);

		// extending an already existing run of gaps;
		if (i < blocks.size()) {
			Block b = blocks.get(i);
			if (b.sourceStart <= pos) {
				i--;
			}
		}
		renumber(i + 1, length);

		assert isSane() : "Data corrupted: " + blocks;
	}

	public void addGapsInView(int pos, int length)
			throws IndexOutOfBoundsException {
		if (pos < 1 || pos > (length() + 1)) {
			throw new IndexOutOfBoundsException(
					"Attempted to add a gap outside of this sequence (1.."
							+ length() + ") at " + pos);
		}

		int i = blocks.size() / 2;
		int imin = 0;
		int imax = blocks.size() - 1;

		do {
			Block b = blocks.get(i);
			if (b.viewStart < pos && b.viewEnd >= pos) { // found block that
				// need splitting
				// with gaps
				Block c = new Block(viewToSource(b, pos), b.sourceEnd, pos,
						b.viewEnd);
				b.viewEnd = pos - 1;
				b.sourceEnd = viewToSource(b, pos - 1);
				blocks.add(i + 1, c);
				renumber(i + 1, length);
				return;
			} else {
				if (b.viewStart < pos) {
					imin = i + 1;
					i = imin + (imax - imin) / 2;
				} else {
					imax = i - 1;
					i = imin + (imax - imin) / 2;
				}
			}
		} while (imin <= imax);

		// extending an already existing run of gaps;
		if (i < blocks.size()) {
			Block b = blocks.get(i);
			if (pos <= b.viewStart) {
				i--;
			}
		} else {
			i--;
		}
		renumber(i + 1, length);
	}

	/**
	 * Get list of the un-gapped region of the SymbolList.
	 * <p>
	 * The gapped symbol list can be represented as a list of ungapped regions.
	 * Given a list of start-stop pairs in the ungapped coordinate system each
	 * with a corresponding pair of start-stop pairs in the gapped view, the
	 * entire gapped list can be reconstructed.
	 * 
	 * @return a List of Block instances
	 */
	public List<Block> BlockIterator() {
		return Collections.unmodifiableList(blocks);
	}

	/**
	 * Debugging method
	 */
	public void dumpBlocks() {
		for (Iterator<Block> i = blocks.iterator(); i.hasNext();) {
			Block b = i.next();
			System.out.println(b.sourceStart + ".." + b.sourceEnd + "\t"
					+ b.viewStart + ".." + b.viewEnd);
		}
	}

	public int firstNonGap() {
		int first = (blocks.get(0)).viewStart;
		return first;
	}

	/**
	 * Translates a Location from the gapped view into the underlying sequence.
	 * End points that are in gaps are moved 'inwards' to shorten the location.
	 * 
	 * @since 1.3
	 */

	public Location gappedToLocation(Location l) {
		if (l.isContiguous()) {
			return gappedToBlock(l);
		} else {
			List<Location> lblocks = new ArrayList<Location>();
			for (Iterator<Location> i = l.blockIterator(); i.hasNext();) {
				lblocks.add(gappedToBlock((Location) i.next()));
			}
			return LocationTools.union(lblocks);
		}
	}

	public Alphabet getAlphabet() {
		return alpha;
	}

	public SymbolList getSourceSymbolList() {
		return source;
	}

	public Location getUngappedLocation() {
		List<Location> locList = new ArrayList<Location>(blocks.size());
		for (Iterator<Block> i = blocks.iterator(); i.hasNext();) {
			Block block = i.next();
			locList.add(new RangeLocation(block.viewStart, block.viewEnd));
		}

		return LocationTools.union(locList);
	}

	public int lastNonGap() {
		int last = (blocks.get(blocks.size() - 1)).viewEnd;
		return last;
	}

	public int length() {
		return length;
	}

	/**
	 * Translate a Location onto the gapped view, splitting blocks if necessary
	 * 
	 * @since 1.3
	 */

	public Location locationToGapped(Location l) {
		if (l.isContiguous()) {
			return blockToGapped(l);
		} else {
			List<Location> lblocks = new ArrayList<Location>();
			for (Iterator<Location> i = l.blockIterator(); i.hasNext();)
				lblocks.add(blockToGapped(i.next()));
			return LocationTools.union(lblocks);
		}
	}

	public void removeGap(int pos) throws IndexOutOfBoundsException,
			IllegalSymbolException {
		if (pos < 1 || pos > length()) {
			throw new IndexOutOfBoundsException(
					"Attempted to remove gap outside of this sequence (1.."
							+ length() + ") at " + pos);
		}
		int i = findViewGap(pos);
		if (i == -2) {
			throw new IllegalSymbolException(
					"Attempted to remove a gap at a non-gap index: " + pos
							+ " -> " + symbolAt(pos).getName());
		}

		if (i == -1 || i == (blocks.size() - 1)) { // at the beginning or the
			// end
			renumber(i + 1, -1);
		} else { // internal
			Block l = blocks.get(i);
			Block r = blocks.get(i + 1);
			renumber(i + 1, -1);
			int gap = r.viewStart - l.viewEnd;
			if (gap == 1) { // removing the last gap - need to join blocks l & r
				// together
				l.sourceEnd = r.sourceEnd;
				l.viewEnd = r.viewEnd;
				blocks.remove(i + 1);
			}
		}

		assert isSane() : "Data corrupted: " + blocks;
	}

	public void removeGaps(int pos, int length)
			throws IndexOutOfBoundsException, IllegalSymbolException {
		int end = pos + length - 1;
		if (pos < 1 || length < 1 || end > length()) {
			throw new IndexOutOfBoundsException(
					"Attempted to remove gap outside of this sequence (1.."
							+ length() + ") at (" + pos + ".." + end + ")");
		}
		int i = findViewGap(pos);
		if (i == -2) {
			throw new IllegalSymbolException(
					"Attempted to remove a gap at a non-gap index: " + pos
							+ " -> " + symbolAt(pos).getName());
		}

		if (i == -1) { // removing track at the beginning
			Block b = blocks.get(0);
			if (b.viewStart <= end) {
				throw new IllegalSymbolException(
						"Attempted to remove some non-gap characters at ("
								+ pos + ".." + end + ")");
			}
			renumber(i + 1, -length);
		} else if (i == (blocks.size() - 1)) { // removing trailing gaps
			this.length -= length;
		} else { // removing internal gaps
			Block l = blocks.get(i);
			Block r = blocks.get(i + 1);
			int gap = r.viewStart - l.viewEnd;
			if (gap < length) {
				throw new IllegalSymbolException("Removing " + length
						+ " gaps from + " + i + " but there are only " + gap
						+ " gaps there: " + blocks);
			}

			renumber(i + 1, -length);
			if (gap == length) { // deleted an entire gapped region
				l.sourceEnd = r.sourceEnd;
				l.viewEnd = r.viewEnd;
				blocks.remove(i + 1);
			}
		}

		assert isSane() : "Data corrupted: removeGaps(" + pos + "," + length
				+ ") " + blocks;
	}

	public final int sourceToView(int indx) throws IndexOutOfBoundsException {
		if (indx < 1 || indx > source.length()) {
			throw new IndexOutOfBoundsException(
					"Attempted to address source sequence (1.." + length()
							+ ") at " + indx);
		}
		int j = findSourceBlock(indx);
		return sourceToView(blocks.get(j), indx);
	}

	public Symbol symbolAt(int indx) throws IndexOutOfBoundsException {
		if (indx > length() || indx < 1) {
			throw new IndexOutOfBoundsException(
					"Attempted to read outside of this sequence (1.."
							+ length() + ") at " + indx);
		}
		int i = findViewBlock(indx);
		if (i < 0) {
			if ((indx < firstNonGap()) || (indx > lastNonGap())) {
				return Alphabet.EMPTY_ALPHABET.getGapSymbol();
			} else {
				return getAlphabet().getGapSymbol();
			}
		} else {
			try {
				Block b = blocks.get(i);
				return source.symbolAt(b.sourceStart - b.viewStart + indx);
			} catch (IndexOutOfBoundsException e) {
				throw new AssertionFailure(
						"Internal book-keeping error fetching index: " + indx
								+ " of " + length() + " blocks: " + blocks, e);
			}
		}
	}

	public final int viewToSource(int indx) throws IndexOutOfBoundsException {
		if (indx < 1 || indx > length()) {
			throw new IndexOutOfBoundsException(
					"Attempted to address sequence (1.." + length() + ") at "
							+ indx);
		}
		int j = findViewBlock(indx);
		if (j < 0) {
			int nj = -j - 1;
			// System.out.println("Converted: " + indx + ":" + j + " to " + nj);
			if (nj < blocks.size()) {
				Block b = blocks.get(nj);
				// System.out.println("Has a following block: " + b);
				return -b.sourceStart;
			} else {
				// System.out.println("Has no following block");
				return -(source.length() + 1);
			}
		} else {
			return viewToSource(blocks.get(j), indx);
		}
	}

	private Location blockToGapped(Location l) {
		int start = l.getMin();
		int end = l.getMax();

		List<Location> lblocks = new ArrayList<Location>();

		while (start <= end) {
			int containingBlockIdx = findSourceBlock(start);
			Block containingBlock = blocks.get(containingBlockIdx);
			int gbstart = start - containingBlock.sourceStart
					+ containingBlock.viewStart;
			int gbend = Math
					.min(end - start + gbstart, containingBlock.viewEnd);
			lblocks.add(new RangeLocation(gbstart, gbend));
			start = start + gbend - gbstart + 1;
		}

		return LocationTools.union(lblocks);
	}

	private Location gappedToBlock(Location l) {
		int start = l.getMin();
		int end = l.getMax();

		int startBlockI = findViewBlock(start);
		int endBlockI = findViewBlock(end);

		if (startBlockI < 0) { // in a gap
			int sb = -startBlockI - 1;
			if (sb == blocks.size()) {
				start = Integer.MAX_VALUE;
			} else {
				Block startBlock = blocks.get(sb);

				start = startBlock.sourceStart;
			}
		} else {
			Block startBlock = blocks.get(startBlockI);
			start = start - startBlock.viewStart + startBlock.sourceStart;
		}

		if (endBlockI < 0) { // in a gap
			int eb = -endBlockI - 1;
			if (eb == 0) {
				end = Integer.MIN_VALUE;
			} else {
				Block endBlock = blocks.get(eb - 1);

				end = endBlock.sourceEnd;
			}
		} else {
			Block endBlock = blocks.get(endBlockI);
			end = end - endBlock.viewEnd + endBlock.sourceEnd;
		}

		if (start > end) {
			return Location.empty;
		} else {
			return new RangeLocation(start, end);
		}
	}

	/**
	 * Finds the index of the block containing the source coordinate indx.
	 * 
	 * @param indx
	 *            the index to find
	 * @return the index of the Block containing indx
	 */
	protected final int findSourceBlock(int indx) {
		int i = blocks.size() / 2;
		int imin = 0;
		int imax = blocks.size() - 1;

		do {
			Block b = blocks.get(i);
			if (b.sourceStart <= indx && b.sourceEnd >= indx) {
				return i;
			} else {
				if (b.sourceStart < indx) {
					imin = i + 1;
					i = imin + (imax - imin) / 2;
				} else {
					imax = i - 1;
					i = imin + (imax - imin) / 2;
				}
			}
		} while (imin <= imax);

		throw new BioError(
				"Something is screwed. Could not find source block containing index "
						+ indx + " in sequence of length " + source.length());
	}

	/**
	 * Finds the index of the Block before the gap at indx within the following
	 * gap.
	 * 
	 * @param indx
	 *            the index to find within a gap
	 * @return the index of the block with indx in the gap
	 */
	protected final int findSourceGap(int indx) {
		int i = blocks.size() / 2;
		int imin = 0;
		int imax = blocks.size() - 1;

		do {
			Block b = blocks.get(i);
			if (b.sourceStart <= indx && b.sourceEnd >= indx) {
				return -1;
			} else {
				if (b.sourceStart < indx) {
					imin = i + 1;
					i = imin + (imax - imin) / 2;
				} else {
					imax = i - 1;
					i = imin + (imax - imin) / 2;
				}
			}
		} while (imin <= imax);

		Block b = blocks.get(i);
		if (b.viewEnd < indx) {
			return i;
		} else {
			return i - 1;
		}
	}

	/**
	 * Finds the index of the Block containing indx within the view ranges.
	 * <p>
	 * If indx is not within a view block, then it is the index of a gap. The
	 * method will return -(indx+1) where indx is the block emediately following
	 * the gap.
	 * 
	 * @param indx
	 *            the index to find within a view range.
	 * @return the index of the block containing index or one less than the
	 *         negative of the index of the block following the gap
	 */
	protected final int findViewBlock(int indx) {
		int i = blocks.size() / 2;
		int imin = 0;
		int imax = blocks.size() - 1;
		// System.out.println("Searching for " + indx);

		Block b;
		do {
			// System.out.println(imin + " < " + i + " < " + imax);
			b = blocks.get(i);
			// System.out.println("Checking " + b.viewStart + ".." + b.viewEnd);
			if (b.viewStart <= indx && b.viewEnd >= indx) {
				// System.out.println("hit");
				return i;
			} else {
				if (b.viewStart < indx) {
					// System.out.println("Too low");
					imin = i + 1;
					i = imin + (imax - imin) / 2;
				} else {
					// System.out.println("Too high");
					imax = i - 1;
					i = imin + (imax - imin) / 2;
				}
			}
		} while (imin <= imax);

		if (i >= blocks.size()) {
			return -blocks.size() - 1;
		}

		if (blocks.get(i) != b) {
			if (blocks.get(i - 1) == b) {
				--i;
			} else if (blocks.get(i + 1) == b) {
				++i;
			}
		}

		// System.out.println("Finding block for: " + indx + " in " + blocks);
		// System.out.println("\ti=" + i);
		// System.out.println("\t" + blocks.get(i));
		if (indx > b.viewStart) { // block before gap - move to block after gap
			// System.out.println("\tAdvancing i");
			i++;
			if (i < blocks.size()) {
				// System.out.println("\t" + blocks.get(i));
			}
		}
		return -i - 1;
	}

	/**
	 * Finds the index of the Block before the gap at indx within the view
	 * range.
	 * <p>
	 * If indx is in-fact a real symbol, then there will be no Block before it.
	 * In this case, the method returns -2. It returns -1 if indx is within the
	 * leading gaps and blocks.size()-1 if it is within the trailing gaps.
	 * 
	 * @param indx
	 *            the index to find within a view range
	 * @return the index of the block with indx in the following gap
	 */
	protected final int findViewGap(int indx) {
		int i = blocks.size() / 2;
		int imin = 0;
		int imax = blocks.size() - 1;

		do {
			Block b = blocks.get(i);
			if (b.viewStart <= indx && b.viewEnd >= indx) {
				return -2;
			} else {
				if (b.viewStart < indx) {
					imin = i + 1;
					i = imin + (imax - imin) / 2;
				} else {
					imax = i - 1;
					i = imin + (imax - imin) / 2;
				}
			}
		} while (imin <= imax);

		if (i < blocks.size()) {
			Block b = blocks.get(i);
			if (b.viewEnd < indx) {
				return i;
			} else {
				return i - 1;
			}
		} else {
			return i - 1;
		}
	}

	protected boolean isSane() {
		for (Iterator<Block> i = blocks.iterator(); i.hasNext();)
			if (!i.next().isSane())
				return false;
		return true;
	}

	/**
	 * Renumber the view indexes from block, adding delta to each offset.
	 * <p>
	 * This adjusts viewStart and viewEnd to be += delta for each block
	 * i->blocks.size(), and sets the total length to += delta.
	 * 
	 * @param i
	 *            the first
	 */
	protected final void renumber(int i, int delta) {
		for (int j = i; j < blocks.size(); j++) {
			Block b = blocks.get(j);
			b.viewStart += delta;
			b.viewEnd += delta;
		}
		length += delta;
	}

	/**
	 * Coordinate conversion from source to view.
	 * 
	 * @param b
	 *            the block containing indx
	 * @param indx
	 *            the index to project
	 * @return the position of indx projected from source to view
	 */
	protected final int sourceToView(Block b, int indx) {
		return indx - b.sourceStart + b.viewStart;
	}

	/**
	 * Coordinate conversion from view to source.
	 * 
	 * @param b
	 *            the block containing indx
	 * @param indx
	 *            the index to project
	 * @return the position of indx projected from view to source
	 */
	protected final int viewToSource(Block b, int indx) {
		return indx - b.viewStart + b.sourceStart;
	}
}
