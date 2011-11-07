package org.biojava.bio.structure.align.util;

import org.biojava.bio.structure.align.model.AFPChain;

/**
 * Some utility methods for analyzing and manipulating AFPChains.
 * 
 * @author Spencer Bliven
 *
 */
public class AlignmentTools {
	/**
	 * Checks that the alignment given by afpChain is sequential. This means
	 * that the residue indices of both proteins increase monotonically as
	 * a function of the alignment position (ie both proteins are sorted).
	 * 
	 * This will return false for circularly permuted alignments or other
	 * non-topological alignments. It will also return false for cases where
	 * the alignment itself is sequential but it is not stored in the afpChain
	 * in a sorted manner.
	 * 
	 * Since algorithms which create non-sequential alignments split the
	 * alignment into multiple blocks, some computational time can be saved
	 * by only checking block boundaries for sequentiality. Setting
	 * <tt>checkWithinBlocks</tt> to <tt>true</tt> makes this function slower,
	 * but detects AFPChains with non-sequential blocks.
	 * 
	 * Note that this method should give the same results as
	 * {@link AFPChain#isSequentialAlignment()}. However, the AFPChain version
	 * relies on the StructureAlignment algorithm correctly setting this
	 * parameter, which is sadly not always the case.
	 * 
	 * @param afpChain An alignment
	 * @param checkWithinBlocks Indicates whether individual blocks should be
	 * 	checked for sequentiality
	 * @return True if the alignment is sequential.
	 */
	public static boolean isSequentialAlignment(AFPChain afpChain, boolean checkWithinBlocks) {
		int[][][] optAln = afpChain.getOptAln();
		int[] alnLen = afpChain.getOptLen();
		int blocks = afpChain.getBlockNum();
		
		if(blocks < 1) return true; //trivial case
		if ( alnLen[0] < 1) return true;
		
		// Check that blocks are sequential
		if(checkWithinBlocks) {
			for(int block = 0; block<blocks; block++) {
				if(alnLen[block] < 1 ) continue; //skip empty blocks
				
				int prevRes1 = optAln[block][0][0];
				int prevRes2 = optAln[block][1][0];
				
				for(int pos = 1; pos<alnLen[block]; pos++) {
					int currRes1 = optAln[block][0][pos];
					int currRes2 = optAln[block][1][pos];
					
					if(currRes1 < prevRes1) {
						return false;
					}
					if(currRes2 < prevRes2) {
						return false;
					}
					
					prevRes1 = currRes1;
					prevRes2 = currRes2;
				}
			}
		}

		// Check that blocks are sequential
		int prevRes1 = optAln[0][0][alnLen[0]-1];
		int prevRes2 = optAln[0][1][alnLen[0]-1];
		
		for(int block = 1; block<blocks;block++) {
			if(alnLen[block] < 1 ) continue; //skip empty blocks

			if(optAln[block][0][0]<prevRes1) {
				return false;
			}
			if(optAln[block][1][0]<prevRes2) {
				return false;
			}
			
			prevRes1 = optAln[block][0][alnLen[block]-1];
			prevRes2 = optAln[block][1][alnLen[block]-1];
		}
		
		return true;
	}
}
