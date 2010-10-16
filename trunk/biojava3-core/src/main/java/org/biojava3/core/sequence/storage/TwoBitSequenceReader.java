package org.biojava3.core.sequence.storage;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Implementation of the 2bit encoding. This will default to the following
 * encodings:
 *
 * <ul>
 * <li>0 - T</li>
 * <li>1 - C</li>
 * <li>2 - A</li>
 * <li>3 - G</li>
 * </ul>
 *
 * We also do not support case sensitive encodings therefore if you pass a
 * lowercased a this will be treated as if it is an uppercase A and we will
 * erase that information.
 *
 * @author ayates
 */
public class TwoBitSequenceReader<C extends NucleotideCompound> extends BitSequenceReader<C> {

	public TwoBitSequenceReader(Sequence<C> sequence) {
		super(new TwoBitArrayWorker<C>(sequence), sequence.getAccession());
	}

	public TwoBitSequenceReader(String sequence, CompoundSet<C> compoundSet) {
		this(sequence, compoundSet, new AccessionID("Unknown"));
	}

	public TwoBitSequenceReader(String sequence, CompoundSet<C> compoundSet, AccessionID accession) {
		super(new TwoBitArrayWorker<C>(sequence, compoundSet), accession);
	}

	public TwoBitSequenceReader(TwoBitArrayWorker<C> worker) {
		super(worker, new AccessionID("unknown"));
	}

	public TwoBitSequenceReader(TwoBitArrayWorker<C> worker, AccessionID accession) {
		super(worker, accession);
	}

	/**
	 * Extension of the BitArrayWorker which provides the 2bit implementation
	 * code. This is intended to work with the 4 basic nucelotide types. If you
	 * require a different version of the encoding used here then extend
	 * and override as required.
	 *
	 * @param <C> Must extend NucleotideCompound
	 */
	public static class TwoBitArrayWorker<C extends NucleotideCompound> extends BitArrayWorker<C> {

		public TwoBitArrayWorker(CompoundSet<C> compoundSet, int length) {
			super(compoundSet, length);
		}

		public TwoBitArrayWorker(CompoundSet<C> compoundSet, int[] sequence) {
			super(compoundSet, sequence);
		}

		public TwoBitArrayWorker(Sequence<C> sequence) {
			super(sequence);
		}

		public TwoBitArrayWorker(String sequence, CompoundSet<C> compoundSet) {
			super(sequence, compoundSet);
		}

		/**
		 * Masking value used for extracting the right most 2 bits from a byte
		 */
		private final static byte MASK = (byte) ((int) Math.pow(2, 0) | (int) Math.pow(2, 1));

		@Override
		protected byte bitMask() {
			return MASK;
		}

		@Override
		protected int compoundsPerDatatype() {
			return 16;
		}

		/**
		 * Returns a Map which encodes TCAG into positions 0,1,2,3.
		 */
		@Override
		@SuppressWarnings("serial")
		protected Map<C, Integer> generateCompoundsToIndex() {
			final CompoundSet<C> cs = getCompoundSet();
			return new HashMap<C, Integer>() {

				{
					put(cs.getCompoundForString("T"), 0);
					put(cs.getCompoundForString("C"), 1);
					put(cs.getCompoundForString("A"), 2);
					put(cs.getCompoundForString("G"), 3);
					put(cs.getCompoundForString("t"), 0);
					put(cs.getCompoundForString("c"), 1);
					put(cs.getCompoundForString("a"), 2);
					put(cs.getCompoundForString("g"), 3);
				}
			};
		}

		/**
		 * Returns a List which encodes TCAG into positions 0,1,2,3.
		 */
		@Override
		protected List<C> generateIndexToCompounds() {
			CompoundSet<C> cs = getCompoundSet();
			List<C> result = new ArrayList<C>();
			result.add( cs.getCompoundForString("T"));


			result.add( cs.getCompoundForString("C"));
			result.add( cs.getCompoundForString("A"));
			result.add( cs.getCompoundForString("G"));
			return result;
		}
	}

}
