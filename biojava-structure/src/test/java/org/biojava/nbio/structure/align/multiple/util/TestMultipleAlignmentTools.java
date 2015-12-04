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
package org.biojava.nbio.structure.align.multiple.util;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.junit.Test;

/**
 * 
 * @author Spencer Bliven
 *
 */
public class TestMultipleAlignmentTools {

	/**
	 * Override the toString method to give shorter output for errors
	 * @author blivens
	 *
	 */
	public static class NamedBlock extends BlockImpl {
		private static final long serialVersionUID = 5060618718423340848L;
		private String name;
		public NamedBlock(String name, BlockSet bs) {
			super(bs);
			this.name = name;
		}
		@Override
		public String toString() {
			return String.format("Block[%s]", name);
		}
	}
	
	@SuppressWarnings("unchecked")
	@Test
	public void testSortBlocks() {
		
		// Sample alignment with blocks out of order
		
		// Row0: Already sorted
		// Row1: Unsorted
		// Row2: Test nulls at start
		// Row3: Test all nulls
		// Row4: Test overlapping ranges
		// Row5: Test only fist element used
		MultipleAlignment align = new MultipleAlignmentImpl();
		BlockSet bs = new BlockSetImpl(align);
		
		Block one = new NamedBlock("1",bs);
		one.setAlignRes(Arrays.asList(
				Arrays.asList( 10, 11, 12),
				Arrays.asList( 20, 21, 22),
				Arrays.asList( null, 21, 22),
				Arrays.asList( 20, 21, 22),
				Arrays.asList( 20, 21, 22),
				Arrays.asList( 20, 21, 22)
				));
		Block two = new NamedBlock("2",bs);
		two.setAlignRes(Arrays.asList(
				Arrays.asList( 20, 21, 22),
				Arrays.asList( 10, 11, 12),
				Arrays.asList( 10, null, 12),
				Arrays.asList( (Integer)null,null,null),
				Arrays.asList( 10, 11, 12),
				Arrays.asList( 10, 11, 12)
				));
		Block three = new NamedBlock("3",bs);
		three.setAlignRes(Arrays.asList(
				Arrays.asList( 30, 31, 32),
				Arrays.asList( 40, 41, 42),
				Arrays.asList( 40, 41, null),
				Arrays.asList( 40, 41, 42),
				Arrays.asList( 40, 41, 42),
				Arrays.asList( 20, 41, 42)
				));
		Block four = new NamedBlock("4",bs);
		four.setAlignRes(Arrays.asList(
				Arrays.asList( 40, 41, 42),
				Arrays.asList( 30, 31, 32),
				Arrays.asList( null, 31, 32),
				Arrays.asList( (Integer)null,null,null),
				Arrays.asList( 30, 51, 52),
				Arrays.asList( 30, 31, 32)
				));
		
		List<Block> blocks;
		int index;
		List<Block> expected;
		
		index = 0;
		blocks = align.getBlocks();
		MultipleAlignmentTools.sortBlocks(blocks, index);
		expected = Arrays.asList(one,two,three,four);
		assertEquals("Bad comparison of row "+index, expected,blocks);

		index = 1;
		blocks = align.getBlocks();
		MultipleAlignmentTools.sortBlocks(blocks, index);
		expected = Arrays.asList(two,one,four,three);
		assertEquals("Bad comparison of row "+index, expected,blocks);

		index = 2;
		blocks = align.getBlocks();
		MultipleAlignmentTools.sortBlocks(blocks, index);
		expected = Arrays.asList(two,one,four,three);
		assertEquals("Bad comparison of row "+index, expected,blocks);

		index = 3;
		blocks = align.getBlocks();
		try {
			MultipleAlignmentTools.sortBlocks(blocks, index);
			fail("Row "+index+" should throw NPE");
		} catch(NullPointerException e) {}

		index = 4;
		blocks = align.getBlocks();
		MultipleAlignmentTools.sortBlocks(blocks, index);
		expected = Arrays.asList(two,one,four,three);
		assertEquals("Bad comparison of row "+index, expected,blocks);

		index = 5;
		blocks = align.getBlocks();
		MultipleAlignmentTools.sortBlocks(blocks, index);
		expected = Arrays.asList(two,one,three,four);
		assertEquals("Bad comparison of row "+index, expected,blocks);

	}

}
