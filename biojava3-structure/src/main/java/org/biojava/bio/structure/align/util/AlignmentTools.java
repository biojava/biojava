package org.biojava.bio.structure.align.util;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.Map.Entry;

import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;

/**
 * Some utility methods for analyzing and manipulating AFPChains.
 * 
 * @author Spencer Bliven
 *
 */
public class AlignmentTools {
	static final boolean debug = false;

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
	
	/**
	 * Creates a Map specifying the alignment as a mapping between residue indices
	 * of protein 1 and residue indices of protein 2.
	 * 
	 * <p>For example,<pre>
	 * 1234
	 * 5678</pre>
	 * becomes<pre>
	 * 1->5
	 * 2->6
	 * 3->7
	 * 4->8</pre>
	 * 
	 * @param afpChain An alignment
	 * @return A mapping from aligned residues of protein 1 to their partners in protein 2.
	 * @throws StructureException If afpChain is not one-to-one
	 */
	public static Map<Integer, Integer> alignmentAsMap(AFPChain afpChain) throws StructureException {
		Map<Integer,Integer> map = new HashMap<Integer,Integer>();
		
		int[][][] optAln = afpChain.getOptAln();
		int[] optLen = afpChain.getOptLen();
		for(int block = 0; block < afpChain.getBlockNum(); block++) {
			for(int pos = 0; pos < optLen[block]; pos++) {
				int res1 = optAln[block][0][pos];
				int res2 = optAln[block][1][pos];
				if(map.containsKey(res1)) {
					throw new StructureException(String.format("Residue %d aligned to both %d and %d.", res1,map.get(res1),res2));
				}
				map.put(res1,res2);
			}
		}
		return map;
	}
	

	/**
	 * Helper for {@link #getSymmetryOrder(Map, Map, int, float)} with a true
	 * identity function (X->X).
	 * 
	 * <p>This method should only be used in cases where the two proteins
	 * aligned have identical numbering, as for self-alignments. See
	 * {@link #getSymmetryOrder(AFPChain, int, float)} for a way to guess
	 * the sequential correspondence between two proteins.
	 * 
	 * @param alignment
	 * @param maxSymmetry
	 * @param minimumMetricChange
	 * @return
	 */
	public static int getSymmetryOrder(Map<Integer, Integer> alignment,
			final int maxSymmetry, final float minimumMetricChange) {
		return getSymmetryOrder(alignment, new IdentityMap<Integer>(), maxSymmetry, minimumMetricChange);
	}
	/**
	 * Tries to detect symmetry in an alignment.
	 * 
     * <p>Conceptually, an alignment is a function f:A->B between two sets of
     * integers. The function may have simple topology (meaning that if two
     * elements of A are close, then their images in B will also be close), or
     * may have more complex topology (such as a circular permutation). This
     * function checks <i>alignment</i> against a reference function
     * <i>identity</i>, which should have simple topology. It then tries to 
     * determine the symmetry order of <i>alignment</i> relative to
     * <i>identity</i>, up to a maximum order of <i>maxSymmetry</i>.
     * 
	 * 
	 * <p><strong>Details</strong><br/>
	 * Considers the offset (in number of residues) which a residue moves
	 * after undergoing <i>n</i> alternating transforms by alignment and
	 * identity. If <i>n</i> corresponds to the intrinsic order of the alignment,
	 * this will be small. This algorithm tries increasing values of <i>n</i>
	 * and looks for abrupt decreases in the root mean squared offset.
	 * If none are found at <i>n</i><=maxSymmetry, the alignment is reported as
	 * non-symmetric.
	 * 
	 * @param alignment The alignment to test for symmetry
	 * @param identity An alignment with simple topology which approximates
	 *  the sequential relationship between the two proteins. Should map in the
	 *  reverse direction from alignment.
	 * @param maxSymmetry Maximum symmetry to consider. High values increase
	 *  the calculation time and can lead to overfitting. 
	 * @param minimumMetricChange Percent decrease in root mean squared offsets
	 *  in order to declare symmetry. 0.4f seems to work well for CeSymm.
	 * @return The order of symmetry of alignment, or -1 if no order <= 
	 *  maxSymmetry is found.
	 *  
	 * @see IdentityMap For a simple identity function
	 */
	public static int getSymmetryOrder(Map<Integer, Integer> alignment, Map<Integer,Integer> identity,
			final int maxSymmetry, final float minimumMetricChange) {
		List<Integer> preimage = new ArrayList<Integer>(alignment.keySet()); // currently unmodified
		List<Integer> image = new ArrayList<Integer>(preimage);
				
		int bestSymmetry = -1;
		double bestMetric = Double.POSITIVE_INFINITY; //lower is better
		boolean foundSymmetry = false;
		
		if(debug) {
			System.out.println("Symm\tPos\tDelta");
		}
		
		for(int n=1;n<=maxSymmetry;n++) {
			int deltasSq = 0;
			int numDeltas = 0;
			// apply alignment
			for(int i=0;i<image.size();i++) {
				Integer pre = image.get(i);
				Integer intermediate = (pre==null?null: alignment.get(pre));
				Integer post = (intermediate==null?null: identity.get(intermediate));
				image.set(i, post);
				
				if(post != null) {
					int delta = post-preimage.get(i);
					
					deltasSq += delta*delta;
					numDeltas++;

					if(debug) {
						System.out.format("%d\t%d\t%d\n",n,preimage.get(i),delta);
					}
				}
				
			}
			
			// Metrics: RMS compensates for the trend of smaller numDeltas with higher order
			// Not normalizing by numDeltas favors smaller orders
			
			double metric = Math.sqrt((double)deltasSq/numDeltas); // root mean squared distance
			//double metric = Math.sqrt((double)deltasSq); // root squared distance
			
			//System.out.format("%d\t%f\n",n,metric);
			
			if(!foundSymmetry && metric < bestMetric * minimumMetricChange) {
				// n = 1 is never the best symmetry
				if(bestMetric < Double.POSITIVE_INFINITY) {
					foundSymmetry = true;
				}
				bestSymmetry = n;
				bestMetric = metric;
			}
			
			// When debugging need to loop over everything. Unneeded in production
			if(!debug && foundSymmetry) {
				break;
			}
			
		}
		if(foundSymmetry) {
			return bestSymmetry;
		} else {
			return -1;
		}
	}
	

	/**
	 * Guesses the order of symmetry in an alignment
	 * 
	 * <p>Uses {@link #getSymmetryOrder(Map alignment, Map identity, int, float)}
	 * to determine the the symmetry order. For the identity alignment, sorts
	 * the aligned residues of each protein sequentially, then defines the ith
	 * residues of each protein to be equivalent.
	 */
	public static int getSymmetryOrder(AFPChain afpChain, int maxSymmetry, float minimumMetricChange) throws StructureException {
		// alignment comes from the afpChain alignment
		Map<Integer,Integer> alignment = AlignmentTools.alignmentAsMap(afpChain);

		// Now construct identity to map aligned residues in sequential order
		Map<Integer, Integer> identity = guessSequentialAlignment(alignment, true);
		
		
		return AlignmentTools.getSymmetryOrder(alignment,
				identity,
				maxSymmetry, minimumMetricChange);
	}

	/**
	 * Takes a potentially non-sequential alignment and guesses a sequential
	 * version of it. Residues from each structure are sorted sequentially and
	 * then compared directly.
	 * 
	 * <p>The results of this method are consistent with what one might expect
	 * from an identity function, and are therefore useful with
	 * {@link #getSymmetryOrder(Map, Map identity, int, float)}.
	 * <ul>
	 *  <li>Perfect self-alignments will have the same pre-image and image,
	 *      so will map X->X</li>
	 *  <li>Gaps and alignment errors will cause errors in the resulting map,
	 *      but only locally. Errors do not propagate through the whole 
	 *      alignment.</li>
	 * </ul>
	 * 
	 * <h4>Example:</h4>
	 * A non sequential alignment, represented schematically as
	 * <pre>
	 * 12456789
	 * 78912345</pre>
	 * would result in a map
	 * <pre>
	 * 12456789
	 * 12345789</pre>
	 * @param afpChain The non-sequential input alignment
	 * @param inverseAlignment If false, map from structure1 to structure2. If
	 *  true, generate the inverse of that map.
	 * @return A mapping from sequential residues of one protein to those of the other
	 * @throws IllegalArgumentException if the input alignment is not one-to-one.
	 */
	public static Map<Integer, Integer> guessSequentialAlignment(
			Map<Integer,Integer> alignment, boolean inverseAlignment) {
		Map<Integer,Integer> identity = new HashMap<Integer,Integer>();
		
		SortedSet<Integer> aligned1 = new TreeSet<Integer>();
		SortedSet<Integer> aligned2 = new TreeSet<Integer>();
		
		for(Entry<Integer,Integer> pair : alignment.entrySet()) {
			aligned1.add(pair.getKey());
			if( !aligned2.add(pair.getValue()) )
				throw new IllegalArgumentException("Alignment is not one-to-one for residue "+pair.getValue()+" of the second structure.");
		}
		
		Iterator<Integer> it1 = aligned1.iterator();
		Iterator<Integer> it2 = aligned2.iterator();
		while(it1.hasNext()) {
			if(inverseAlignment) { // 2->1
				identity.put(it2.next(),it1.next());
			} else { // 1->2
				identity.put(it1.next(),it2.next());
			}
		}
		return identity;
	}
	
	/**
	 * Retrieves the optimum alignment from an AFPChain and returns it as a
	 * java collection. The result is indexed in the same way as 
	 * {@link AFPChain#getOptAln()}, but has the correct size().
	 * <pre>
	 * List<List<List<Integer>>> aln = getOptAlnAsList(AFPChain afpChain);
	 * aln.get(blockNum).get(structureNum={0,1}).get(pos)</pre>
	 * 
	 * @param afpChain
	 * @return
	 */
	public static List<List<List<Integer>>> getOptAlnAsList(AFPChain afpChain) {
		int[][][] optAln = afpChain.getOptAln();
		int[] optLen = afpChain.getOptLen();
		List<List<List<Integer>>> blocks = new ArrayList<List<List<Integer>>>(afpChain.getBlockNum());
		for(int blockNum=0;blockNum<afpChain.getBlockNum();blockNum++) {
			//TODO could improve speed an memory by wrapping the arrays with
			// an unmodifiable list, similar to Arrays.asList(...) but with the
			// correct size parameter.
			List<Integer> align1 = new ArrayList<Integer>(optLen[blockNum]);
			List<Integer> align2 = new ArrayList<Integer>(optLen[blockNum]);
			for(int pos=0;pos<optLen[blockNum];pos++) {
				align1.add(optAln[blockNum][0][pos]);
				align2.add(optAln[blockNum][1][pos]);
			}
			List<List<Integer>> block = new ArrayList<List<Integer>>(2);
			block.add(align1);
			block.add(align2);
			blocks.add(block);
		}
		
		return blocks;
	}



	/**
	 * A Map<K,V> can be viewed as a function from K to V. This class represents
	 * the identity function. Getting a value results in the value itself.
	 * 
	 * <p>The class is a bit inconsistent when representing its contents. On
	 * the one hand, containsKey(key) is true for all objects. However, 
	 * attempting to iterate through the values returns an empty set.
	 * 
	 * @author Spencer Bliven
	 *
	 * @param <K>
	 * @param <V>
	 */
	public static class IdentityMap<K> extends AbstractMap<K,K> {
		public IdentityMap() {}

		/**
		 * @param key
		 * @return the key
		 * @throws ClassCastException if key is not of type K
		 */
		@SuppressWarnings("unchecked")
		@Override
		public K get(Object key) {
			return (K)key;
		}
		
		/**
		 * Always returns the empty set
		 */
		@Override
		public Set<java.util.Map.Entry<K, K>> entrySet() {
			return Collections.emptySet();
		}
		
		@Override
		public boolean containsKey(Object key) {
			return true;
		}
	}
}
