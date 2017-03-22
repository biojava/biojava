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
package org.biojava.nbio.structure.align.util;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.AFPTwister;
import org.biojava.nbio.structure.align.ce.CECalculator;
import org.biojava.nbio.structure.align.fatcat.FatCatFlexible;
import org.biojava.nbio.structure.align.fatcat.FatCatRigid;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.xml.AFPChainXMLParser;
import org.biojava.nbio.structure.geometry.Matrices;
import org.biojava.nbio.structure.geometry.SuperPositions;
import org.biojava.nbio.structure.jama.Matrix;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.Writer;
import java.util.*;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.Matrix4d;

/**
 * Methods for analyzing and manipulating AFPChains and for
 * other pairwise alignment utilities. <p>
 * Current methods: replace optimal alignment, create new AFPChain,
 * format conversion, update superposition, etc.
 *
 * @author Spencer Bliven
 * @author Aleix Lafita
 *
 */
public class AlignmentTools {

	private static final Logger logger = LoggerFactory.getLogger(AlignmentTools.class);


	public static boolean debug = false;

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

		if( afpChain.getAlnLength() < 1 ) {
			return map;
		}
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
	 * Applies an alignment k times. Eg if alignmentMap defines function f(x),
	 * this returns a function f^k(x)=f(f(...f(x)...)).
	 *
	 * @param <T>
	 * @param alignmentMap The input function, as a map (see {@link AlignmentTools#alignmentAsMap(AFPChain)})
	 * @param k The number of times to apply the alignment
	 * @return A new alignment. If the input function is not automorphic
	 *  (one-to-one), then some inputs may map to null, indicating that the
	 *  function is undefined for that input.
	 */
	public static <T> Map<T,T> applyAlignment(Map<T, T> alignmentMap, int k) {
		return applyAlignment(alignmentMap, new IdentityMap<T>(), k);
	}

	/**
	 * Applies an alignment k times. Eg if alignmentMap defines function f(x),
	 * this returns a function f^k(x)=f(f(...f(x)...)).
	 *
	 * To allow for functions with different domains and codomains, the identity
	 * function allows converting back in a reasonable way. For instance, if
	 * alignmentMap represented an alignment between two proteins with different
	 * numbering schemes, the identity function could calculate the offset
	 * between residue numbers, eg I(x) = x-offset.
	 *
	 * When an identity function is provided, the returned function calculates
	 * f^k(x) = f(I( f(I( ... f(x) ... )) )).
	 *
	 * @param <S>
	 * @param <T>
	 * @param alignmentMap The input function, as a map (see {@link AlignmentTools#alignmentAsMap(AFPChain)})
	 * @param identity An identity-like function providing the isomorphism between
	 *  the codomain of alignmentMap (of type <T>) and the domain (type <S>).
	 * @param k The number of times to apply the alignment
	 * @return A new alignment. If the input function is not automorphic
	 *  (one-to-one), then some inputs may map to null, indicating that the
	 *  function is undefined for that input.
	 */
	public static <S,T> Map<S,T> applyAlignment(Map<S, T> alignmentMap, Map<T,S> identity, int k) {

		// This implementation simply applies the map k times.
		// If k were large, it would be more efficient to do this recursively,
		// (eg f^4 = (f^2)^2) but k will usually be small.

		if(k<0) throw new IllegalArgumentException("k must be positive");
		if(k==1) {
			return new HashMap<S,T>(alignmentMap);
		}
		// Convert to lists to establish a fixed order
		List<S> preimage = new ArrayList<S>(alignmentMap.keySet()); // currently unmodified
		List<S> image = new ArrayList<S>(preimage);

		for(int n=1;n<k;n++) {
			// apply alignment
			for(int i=0;i<image.size();i++) {
				S pre = image.get(i);
				T intermediate = (pre==null?null: alignmentMap.get(pre));
				S post = (intermediate==null?null: identity.get(intermediate));
				image.set(i, post);
			}
		}



		Map<S, T> imageMap = new HashMap<S,T>(alignmentMap.size());

		//TODO handle nulls consistently.
		// assure that all the residues in the domain are valid keys
		/*
		for(int i=0;i<preimage.size();i++) {
			S pre = preimage.get(i);
			T intermediate = (pre==null?null: alignmentMap.get(pre));
			S post = (intermediate==null?null: identity.get(intermediate));
			imageMap.put(post, null);
		}
		 */
		// now populate with actual values
		for(int i=0;i<preimage.size();i++) {
			S pre = preimage.get(i);

			// image is currently f^k-1(x), so take the final step
			S preK1 = image.get(i);
			T postK = (preK1==null?null: alignmentMap.get(preK1));
			imageMap.put(pre,postK);

		}
		return imageMap;
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
	 * @return The order of symmetry of alignment, or 1 if no order <=
	 *  maxSymmetry is found.
	 *
	 * @see IdentityMap For a simple identity function
	 */
	public static int getSymmetryOrder(Map<Integer, Integer> alignment, Map<Integer,Integer> identity,
									   final int maxSymmetry, final float minimumMetricChange) {
		List<Integer> preimage = new ArrayList<Integer>(alignment.keySet()); // currently unmodified
		List<Integer> image = new ArrayList<Integer>(preimage);

		int bestSymmetry = 1;
		double bestMetric = Double.POSITIVE_INFINITY; //lower is better
		boolean foundSymmetry = false;

		if(debug) {
			logger.trace("Symm\tPos\tDelta");
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
						logger.debug("%d\t%d\t%d\n",n,preimage.get(i),delta);
					}
				}

			}

			// Metrics: RMS compensates for the trend of smaller numDeltas with higher order
			// Not normalizing by numDeltas favors smaller orders

			double metric = Math.sqrt((double)deltasSq/numDeltas); // root mean squared distance

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
			return 1;
		}
	}


	/**
	 * Guesses the order of symmetry in an alignment
	 *
	 * <p>Uses {@link #getSymmetryOrder(Map alignment, Map identity, int, float)}
	 * to determine the the symmetry order. For the identity alignment, sorts
	 * the aligned residues of each protein sequentially, then defines the ith
	 * residues of each protein to be equivalent.
	 *
	 * <p>Note that the selection of the identity alignment here is <i>very</i>
	 * naive, and only works for proteins with very good coverage. Wherever
	 * possible, it is better to construct an identity function explicitly
	 * from a sequence alignment (or use an {@link IdentityMap} for internally
	 * symmetric proteins) and use {@link #getSymmetryOrder(Map, Map, int, float)}.
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
	 * @param alignment The non-sequential input alignment
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

	/**
	 * Fundamentally, an alignment is just a list of aligned residues in each
	 * protein. This method converts two lists of ResidueNumbers into an
	 * AFPChain.
	 *
	 * <p>Parameters are filled with defaults (often null) or sometimes
	 * calculated.
	 *
	 * <p>For a way to modify the alignment of an existing AFPChain, see
	 * {@link AlignmentTools#replaceOptAln(AFPChain, Atom[], Atom[], Map)}
	 * @param ca1 CA atoms of the first protein
	 * @param ca2 CA atoms of the second protein
	 * @param aligned1 A list of aligned residues from the first protein
	 * @param aligned2 A list of aligned residues from the second protein.
	 *  Must be the same length as aligned1.
	 * @return An AFPChain representing the alignment. Many properties may be
	 *  null or another default.
	 * @throws StructureException if an error occured during superposition
	 * @throws IllegalArgumentException if aligned1 and aligned2 have different
	 *  lengths
	 * @see AlignmentTools#replaceOptAln(AFPChain, Atom[], Atom[], Map)
	 */
	public static AFPChain createAFPChain(Atom[] ca1, Atom[] ca2,
										  ResidueNumber[] aligned1, ResidueNumber[] aligned2 ) throws StructureException {
		//input validation
		int alnLen = aligned1.length;
		if(alnLen != aligned2.length) {
			throw new IllegalArgumentException("Alignment lengths are not equal");
		}

		AFPChain a = new AFPChain(AFPChain.UNKNOWN_ALGORITHM);
		try {
			a.setName1(ca1[0].getGroup().getChain().getStructure().getName());
			if(ca2[0].getGroup().getChain().getStructure() != null) {
				// common case for cloned ca2
				a.setName2(ca2[0].getGroup().getChain().getStructure().getName());
			}
		} catch(Exception e) {
			// One of the structures wasn't fully created. Ignore
		}
		a.setBlockNum(1);
		a.setCa1Length(ca1.length);
		a.setCa2Length(ca2.length);

		a.setOptLength(alnLen);
		a.setOptLen(new int[] {alnLen});


		Matrix[] ms = new Matrix[a.getBlockNum()];
		a.setBlockRotationMatrix(ms);
		Atom[] blockShiftVector = new Atom[a.getBlockNum()];
		a.setBlockShiftVector(blockShiftVector);

		String[][][] pdbAln = new String[1][2][alnLen];
		for(int i=0;i<alnLen;i++) {
			pdbAln[0][0][i] = aligned1[i].getChainName()+":"+aligned1[i];
			pdbAln[0][1][i] = aligned2[i].getChainName()+":"+aligned2[i];
		}

		a.setPdbAln(pdbAln);

		// convert pdbAln to optAln, and fill in some other basic parameters
		AFPChainXMLParser.rebuildAFPChain(a, ca1, ca2);

		return a;

		// Currently a single block. Split into several blocks by sequence if needed
		//		return AlignmentTools.splitBlocksByTopology(a,ca1,ca2);
	}

	/**
	 *
	 * @param a
	 * @param ca1
	 * @param ca2
	 * @return
	 * @throws StructureException if an error occurred during superposition
	 */
	public static AFPChain splitBlocksByTopology(AFPChain a, Atom[] ca1, Atom[] ca2) throws StructureException {
		int[][][] optAln = a.getOptAln();
		int blockNum = a.getBlockNum();
		int[] optLen = a.getOptLen();

		// Determine block lengths
		// Split blocks if residue indices don't increase monotonically
		List<Integer> newBlkLen = new ArrayList<Integer>();
		boolean blockChanged = false;
		for(int blk=0;blk<blockNum;blk++) {
			int currLen=1;
			for(int pos=1;pos<optLen[blk];pos++) {
				if( optAln[blk][0][pos] <= optAln[blk][0][pos-1]
						|| optAln[blk][1][pos] <= optAln[blk][1][pos-1] )
				{
					//start a new block
					newBlkLen.add(currLen);
					currLen = 0;
					blockChanged = true;
				}
				currLen++;
			}
			if(optLen[blk] < 2 ) {
				newBlkLen.add(optLen[blk]);
			} else {
				newBlkLen.add(currLen);
			}
		}

		// Check if anything needs to be split
		if( !blockChanged ) {
			return a;
		}

		// Split blocks
		List<int[][]> blocks = new ArrayList<int[][]>( newBlkLen.size() );

		int oldBlk = 0;
		int pos = 0;
		for(int blkLen : newBlkLen) {
			if( blkLen == optLen[oldBlk] ) {
				assert(pos == 0); //should be the whole block
				// Use the old block
				blocks.add(optAln[oldBlk]);
			} else {
				int[][] newBlock = new int[2][blkLen];
				assert( pos+blkLen <= optLen[oldBlk] ); // don't overrun block
				for(int i=0; i<blkLen;i++) {
					newBlock[0][i] = optAln[oldBlk][0][pos + i];
					newBlock[1][i] = optAln[oldBlk][1][pos + i];
				}
				pos += blkLen;
				blocks.add(newBlock);

				if( pos == optLen[oldBlk] ) {
					// Finished this oldBlk, start the next
					oldBlk++;
					pos = 0;
				}
			}
		}

		// Store new blocks
		int[][][] newOptAln = blocks.toArray(new int[0][][]);
		int[] newBlockLens = new int[newBlkLen.size()];
		for(int i=0;i<newBlkLen.size();i++) {
			newBlockLens[i] = newBlkLen.get(i);
		}

		return replaceOptAln(a, ca1, ca2, blocks.size(), newBlockLens, newOptAln);
	}

	/**
	 * It replaces an optimal alignment of an AFPChain and calculates all the new alignment scores and variables.
	 */
	public static AFPChain replaceOptAln(int[][][] newAlgn, AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws StructureException {

		//The order is the number of groups in the newAlgn
		int order = newAlgn.length;

		//Calculate the alignment length from all the subunits lengths
		int[] optLens = new int[order];
		for(int s=0;s<order;s++) {
			optLens[s] = newAlgn[s][0].length;
		}
		int optLength = 0;
		for(int s=0;s<order;s++) {
			optLength += optLens[s];
		}

		//Create a copy of the original AFPChain and set everything needed for the structure update
		AFPChain copyAFP = (AFPChain) afpChain.clone();

		//Set the new parameters of the optimal alignment
		copyAFP.setOptLength(optLength);
		copyAFP.setOptLen(optLens);
		copyAFP.setOptAln(newAlgn);

		//Set the block information of the new alignment
		copyAFP.setBlockNum(order);
		copyAFP.setBlockSize(optLens);
		copyAFP.setBlockResList(newAlgn);
		copyAFP.setBlockResSize(optLens);
		copyAFP.setBlockGap(calculateBlockGap(newAlgn));

		//Recalculate properties: superposition, tm-score, etc
		Atom[] ca2clone = StructureTools.cloneAtomArray(ca2); // don't modify ca2 positions
		AlignmentTools.updateSuperposition(copyAFP, ca1, ca2clone);

		//It re-does the sequence alignment strings from the OptAlgn information only
		copyAFP.setAlnsymb(null);
		AFPAlignmentDisplay.getAlign(copyAFP, ca1, ca2clone);

		return copyAFP;
	}

	/**
	 * Takes an AFPChain and replaces the optimal alignment based on an alignment map
	 *
	 * <p>Parameters are filled with defaults (often null) or sometimes
	 * calculated.
	 *
	 * <p>For a way to create a new AFPChain, see
	 * {@link AlignmentTools#createAFPChain(Atom[], Atom[], ResidueNumber[], ResidueNumber[])}
	 *
	 * @param afpChain The alignment to be modified
	 * @param alignment The new alignment, as a Map
	 * @throws StructureException if an error occurred during superposition
	 * @see AlignmentTools#createAFPChain(Atom[], Atom[], ResidueNumber[], ResidueNumber[])
	 */
	public static AFPChain replaceOptAln(AFPChain afpChain, Atom[] ca1, Atom[] ca2,
										 Map<Integer, Integer> alignment) throws StructureException {

		// Determine block lengths
		// Sort ca1 indices, then start a new block whenever ca2 indices aren't
		// increasing monotonically.
		Integer[] res1 = alignment.keySet().toArray(new Integer[0]);
		Arrays.sort(res1);
		List<Integer> blockLens = new ArrayList<Integer>(2);
		int optLength = 0;
		Integer lastRes = alignment.get(res1[0]);
		int blkLen = lastRes==null?0:1;
		for(int i=1;i<res1.length;i++) {
			Integer currRes = alignment.get(res1[i]); //res2 index
			assert(currRes != null);// could be converted to if statement if assertion doesn't hold; just modify below as well.
			if(lastRes<currRes) {
				blkLen++;
			} else {
				// CP!
				blockLens.add(blkLen);
				optLength+=blkLen;
				blkLen = 1;
			}
			lastRes = currRes;
		}
		blockLens.add(blkLen);
		optLength+=blkLen;

		// Create array structure for alignment
		int[][][] optAln = new int[blockLens.size()][][];
		int pos1 = 0; //index into res1
		for(int blk=0;blk<blockLens.size();blk++) {
			optAln[blk] = new int[2][];
			blkLen = blockLens.get(blk);
			optAln[blk][0] = new int[blkLen];
			optAln[blk][1] = new int[blkLen];
			int pos = 0; //index into optAln
			while(pos<blkLen) {
				optAln[blk][0][pos]=res1[pos1];
				Integer currRes = alignment.get(res1[pos1]);
				optAln[blk][1][pos]=currRes;
				pos++;
				pos1++;
			}
		}
		assert(pos1 == optLength);

		// Create length array
		int[] optLens = new int[blockLens.size()];
		for(int i=0;i<blockLens.size();i++) {
			optLens[i] = blockLens.get(i);
		}

		return replaceOptAln(afpChain, ca1, ca2, blockLens.size(), optLens, optAln);
	}

	/**
	 * @param afpChain Input afpchain. UNMODIFIED
	 * @param ca1
	 * @param ca2
	 * @param optLens
	 * @param optAln
	 * @return A NEW AfpChain based off the input but with the optAln modified
	 * @throws StructureException if an error occured during superposition
	 */
	public static AFPChain replaceOptAln(AFPChain afpChain, Atom[] ca1, Atom[] ca2,
										 int blockNum, int[] optLens, int[][][] optAln) throws StructureException {
		int optLength = 0;
		for( int blk=0;blk<blockNum;blk++) {
			optLength += optLens[blk];
		}

		//set everything
		AFPChain refinedAFP = (AFPChain) afpChain.clone();
		refinedAFP.setOptLength(optLength);
		refinedAFP.setBlockSize(optLens);
		refinedAFP.setOptLen(optLens);
		refinedAFP.setOptAln(optAln);
		refinedAFP.setBlockNum(blockNum);

		//TODO recalculate properties: superposition, tm-score, etc
		Atom[] ca2clone = StructureTools.cloneAtomArray(ca2); // don't modify ca2 positions
		AlignmentTools.updateSuperposition(refinedAFP, ca1, ca2clone);

		AFPAlignmentDisplay.getAlign(refinedAFP, ca1, ca2clone);
		return refinedAFP;
	}


	/**
	 * After the alignment changes (optAln, optLen, blockNum, at a minimum),
	 * many other properties which depend on the superposition will be invalid.
	 *
	 * This method re-runs a rigid superposition over the whole alignment
	 * and repopulates the required properties, including RMSD (TotalRMSD) and
	 * TM-Score.
	 * @param afpChain
	 * @param ca1
	 * @param ca2 Second set of ca atoms. Will be modified based on the superposition
	 * @throws StructureException
	 * @see {@link CECalculator#calc_rmsd(Atom[], Atom[], int, boolean)}
	 *  contains much of the same code, but stores results in a CECalculator
	 *  instance rather than an AFPChain
	 */
	public static void updateSuperposition(AFPChain afpChain, Atom[] ca1, 
			Atom[] ca2) throws StructureException {

		//Update ca information, because the atom array might also be changed
		afpChain.setCa1Length(ca1.length);
		afpChain.setCa2Length(ca2.length);

		//We need this to get the correct superposition
		int[] focusRes1 = afpChain.getFocusRes1();
		int[] focusRes2 = afpChain.getFocusRes2();
		if (focusRes1 == null) {
			focusRes1 = new int[afpChain.getCa1Length()];
			afpChain.setFocusRes1(focusRes1);
		}
		if (focusRes2 == null) {
			focusRes2 = new int[afpChain.getCa2Length()];
			afpChain.setFocusRes2(focusRes2);
		}

		if (afpChain.getNrEQR() == 0) return;

		// create new arrays for the subset of atoms in the alignment.
		Atom[] ca1aligned = new Atom[afpChain.getOptLength()];
		Atom[] ca2aligned = new Atom[afpChain.getOptLength()];
		
		fillAlignedAtomArrays(afpChain, ca1, ca2, ca1aligned, ca2aligned);

		//Superimpose the two structures in correspondance to the new alignment
		Matrix4d trans = SuperPositions.superpose(Calc.atomsToPoints(ca1aligned),
				Calc.atomsToPoints(ca2aligned));

		Matrix matrix = Matrices.getRotationJAMA(trans);
		Atom shift = Calc.getTranslationVector(trans);

		Matrix[] blockMxs = new Matrix[afpChain.getBlockNum()];
		Arrays.fill(blockMxs, matrix);
		afpChain.setBlockRotationMatrix(blockMxs);
		Atom[] blockShifts = new Atom[afpChain.getBlockNum()];
		Arrays.fill(blockShifts, shift);
		afpChain.setBlockShiftVector(blockShifts);

		for (Atom a : ca2aligned) {
			Calc.rotate(a, matrix);
			Calc.shift(a, shift);
		}

		//Calculate the RMSD and TM score for the new alignment
		double rmsd = Calc.rmsd(ca1aligned, ca2aligned);
		double tmScore = Calc.getTMScore(ca1aligned, ca2aligned, ca1.length, ca2.length);
		afpChain.setTotalRmsdOpt(rmsd);
		afpChain.setTMScore(tmScore);
		
		int[] blockLens = afpChain.getOptLen();
		int[][][] optAln = afpChain.getOptAln();

		//Calculate the RMSD and TM score for every block of the new alignment
		double[] blockRMSD = new double[afpChain.getBlockNum()];
		double[] blockScore = new double[afpChain.getBlockNum()];
		for (int k=0; k<afpChain.getBlockNum(); k++){
			//Create the atom arrays corresponding to the aligned residues in the block
			Atom[] ca1block = new Atom[afpChain.getOptLen()[k]];
			Atom[] ca2block = new Atom[afpChain.getOptLen()[k]];
			int position=0;
			for(int i=0;i<blockLens[k];i++) {
				int pos1 = optAln[k][0][i];
				int pos2 = optAln[k][1][i];
				Atom a1 = ca1[pos1];
				Atom a2 = (Atom) ca2[pos2].clone();
				ca1block[position] = a1;
				ca2block[position] = a2;
				position++;
			}
			if (position != afpChain.getOptLen()[k]){
				logger.warn("AFPChainScorer getTMScore: Problems reconstructing block alignment! nr of loaded atoms is " + position + " but should be " + afpChain.getOptLen()[k]);
				// we need to resize the array, because we allocated too many atoms earlier on.
				ca1block = (Atom[]) resizeArray(ca1block, position);
				ca2block = (Atom[]) resizeArray(ca2block, position);
			}
			//Superimpose the two block structures
			Matrix4d transb = SuperPositions.superpose(Calc.atomsToPoints(ca1block),
					Calc.atomsToPoints(ca2block));

			blockMxs[k] = Matrices.getRotationJAMA(trans);
			blockShifts[k] = Calc.getTranslationVector(trans);

			Calc.transform(ca2block, transb);

			//Calculate the RMSD and TM score for the block
			double rmsdb = Calc.rmsd(ca1block, ca2block);
			double tmScoreb = Calc.getTMScore(ca1block, ca2block, ca1.length, ca2.length);
			blockRMSD[k] = rmsdb;
			blockScore[k] = tmScoreb;
		}
		afpChain.setOptRmsd(blockRMSD);
		afpChain.setBlockRmsd(blockRMSD);
		afpChain.setBlockScore(blockScore);
	}

	/**
	 * Reallocates an array with a new size, and copies the contents
	 * of the old array to the new array.
	 * @param oldArray  the old array, to be reallocated.
	 * @param newSize   the new array size.
	 * @return          A new array with the same contents.
	 */
	public static Object resizeArray (Object oldArray, int newSize) {
		int oldSize = java.lang.reflect.Array.getLength(oldArray);
		@SuppressWarnings("rawtypes")
		Class elementType = oldArray.getClass().getComponentType();
		Object newArray = java.lang.reflect.Array.newInstance(
				elementType,newSize);
		int preserveLength = Math.min(oldSize,newSize);
		if (preserveLength > 0)
			System.arraycopy (oldArray,0,newArray,0,preserveLength);
		return newArray;
	}

	/**
	 * Print an alignment map in a concise representation. Edges are given
	 * as two numbers separated by '>'. They are chained together where possible,
	 * or separated by spaces where disjoint or branched.
	 *
	 * <p>Note that more concise representations may be possible.</p>
	 *
	 * Examples:
	 * <li>1>2>3>1</li>
	 * <li>1>2>3>2 4>3</li>
	 *
	 * @param alignment The input function, as a map (see {@link AlignmentTools#alignmentAsMap(AFPChain)})
	 * @param identity An identity-like function providing the isomorphism between
	 *  the codomain of alignment (of type <T>) and the domain (type <S>).
	 * @return
	 */
	public static <S,T> String toConciseAlignmentString(Map<S,T> alignment, Map<T,S> identity) {
		// Clone input to prevent changes
		Map<S,T> alig = new HashMap<S,T>(alignment);

		// Generate inverse alignment
		Map<S,List<S>> inverse = new HashMap<S,List<S>>();
		for(Entry<S,T> e: alig.entrySet()) {
			S val = identity.get(e.getValue());
			if( inverse.containsKey(val) ) {
				List<S> l = inverse.get(val);
				l.add(e.getKey());
			} else {
				List<S> l = new ArrayList<S>();
				l.add(e.getKey());
				inverse.put(val,l);
			}
		}

		StringBuilder str = new StringBuilder();

		while(!alig.isEmpty()){
			// Pick an edge and work upstream to a root or cycle
			S seedNode = alig.keySet().iterator().next();
			S node = seedNode;
			if( inverse.containsKey(seedNode)) {
				node = inverse.get(seedNode).iterator().next();
				while( node != seedNode && inverse.containsKey(node)) {
					node = inverse.get(node).iterator().next();
				}
			}

			// Now work downstream, deleting edges as we go
			seedNode = node;
			str.append(node);

			while(alig.containsKey(node)) {
				S lastNode = node;
				node = identity.get( alig.get(lastNode) );

				// Output
				str.append('>');
				str.append(node);

				// Remove edge
				alig.remove(lastNode);
				List<S> inv = inverse.get(node);
				if(inv.size() > 1) {
					inv.remove(node);
				} else {
					inverse.remove(node);
				}
			}
			if(!alig.isEmpty()) {
				str.append(' ');
			}
		}

		return str.toString();
	}

	/**
	 * @see #toConciseAlignmentString(Map, Map)
	 */
	public static <T> String toConciseAlignmentString(Map<T, T> alignment) {
		return toConciseAlignmentString(alignment, new IdentityMap<T>());
	}

	/**
	 * @see #toConciseAlignmentString(Map, Map)
	 */
	public static Map<Integer, Integer> fromConciseAlignmentString(String string) {
		Map<Integer, Integer> map = new HashMap<Integer, Integer>();
		boolean matches = true;
		while (matches) {
			Pattern pattern = Pattern.compile("(\\d+)>(\\d+)");
			Matcher matcher = pattern.matcher(string);
			matches = matcher.find();
			if (matches) {
				Integer from = Integer.parseInt(matcher.group(1));
				Integer to = Integer.parseInt(matcher.group(2));
				map.put(from, to);
				string = string.substring(matcher.end(1) + 1);
			}
		}
		return map;
	}

	/**
	 * Method that calculates the number of gaps in each subunit block of an optimal AFP alignment.
	 *
	 * INPUT: an optimal alignment in the format int[][][].
	 * OUTPUT: an int[] array of <order> length containing the gaps in each block as int[block].
	 */
	public static int[] calculateBlockGap(int[][][] optAln){

		//Initialize the array to be returned
		int [] blockGap = new int[optAln.length];

		//Loop for every block and look in both chains for non-contiguous residues.
		for (int i=0; i<optAln.length; i++){
			int gaps = 0; //the number of gaps in that block
			int last1 = 0; //the last residue position in chain 1
			int last2 = 0; //the last residue position in chain 2
			//Loop for every position in the block
			for (int j=0; j<optAln[i][0].length; j++){
				//If the first position is evaluated initialize the last positions
				if (j==0){
					last1 = optAln[i][0][j];
					last2 = optAln[i][1][j];
				}
				else{
					//If one of the positions or both are not contiguous increment the number of gaps
					if (optAln[i][0][j] > last1+1 || optAln[i][1][j] > last2+1){
						gaps++;
						last1 = optAln[i][0][j];
						last2 = optAln[i][1][j];
					}
					//Otherwise just set the last position to the current one
					else{
						last1 = optAln[i][0][j];
						last2 = optAln[i][1][j];
					}
				}
			}
			blockGap[i] = gaps;
		}
		return blockGap;
	}

	/**
	 * Creates a simple interaction format (SIF) file for an alignment.
	 *
	 * The SIF file can be read by network software (eg Cytoscape) to analyze
	 * alignments as graphs.
	 *
	 * This function creates a graph with residues as nodes and two types of edges:
	 *   1. backbone edges, which connect adjacent residues in the aligned protein
	 *   2. alignment edges, which connect aligned residues
	 *
	 * @param out Stream to write to
	 * @param afpChain alignment to write
	 * @param ca1 First protein, used to generate node names
	 * @param ca2 Second protein, used to generate node names
	 * @param backboneInteraction Two-letter string used to identify backbone edges
	 * @param alignmentInteraction Two-letter string used to identify alignment edges
	 * @throws IOException
	 */
	public static void alignmentToSIF(Writer out,AFPChain afpChain,
									  Atom[] ca1,Atom[] ca2, String backboneInteraction,
									  String alignmentInteraction) throws IOException {

		//out.write("Res1\tInteraction\tRes2\n");
		String name1 = afpChain.getName1();
		String name2 = afpChain.getName2();
		if(name1==null) name1=""; else name1+=":";
		if(name2==null) name2=""; else name2+=":";

		// Print alignment edges
		int nblocks = afpChain.getBlockNum();
		int[] blockLen = afpChain.getOptLen();
		int[][][] optAlign = afpChain.getOptAln();
		for(int b=0;b<nblocks;b++) {
			for(int r=0;r<blockLen[b];r++) {
				int res1 = optAlign[b][0][r];
				int res2 = optAlign[b][1][r];

				ResidueNumber rn1 = ca1[res1].getGroup().getResidueNumber();
				ResidueNumber rn2 = ca2[res2].getGroup().getResidueNumber();

				String node1 = name1+rn1.getChainName()+rn1.toString();
				String node2 = name2+rn2.getChainName()+rn2.toString();

				out.write(String.format("%s\t%s\t%s\n",node1, alignmentInteraction, node2));
			}
		}

		// Print first backbone edges
		ResidueNumber rn = ca1[0].getGroup().getResidueNumber();
		String last = name1+rn.getChainName()+rn.toString();
		for(int i=1;i<ca1.length;i++) {
			rn = ca1[i].getGroup().getResidueNumber();
			String curr = name1+rn.getChainName()+rn.toString();
			out.write(String.format("%s\t%s\t%s\n",last, backboneInteraction, curr));
			last = curr;
		}

		// Print second backbone edges, if the proteins differ
		// Do some quick checks for whether the proteins differ
		// (Not perfect, but should detect major differences and CPs.)
		if(!name1.equals(name2) ||
				ca1.length!=ca2.length ||
				(ca1.length>0 && ca1[0].getGroup()!=null && ca2[0].getGroup()!=null &&
						!ca1[0].getGroup().getResidueNumber().equals(ca2[0].getGroup().getResidueNumber()) ) ) {
			rn = ca2[0].getGroup().getResidueNumber();
			last = name2+rn.getChainName()+rn.toString();
			for(int i=1;i<ca2.length;i++) {
				rn = ca2[i].getGroup().getResidueNumber();
				String curr = name2+rn.getChainName()+rn.toString();
				out.write(String.format("%s\t%s\t%s\n",last, backboneInteraction, curr));
				last = curr;
			}
		}
	}



	/** get an artificial List of chains containing the Atoms and groups.
	 * Does NOT rotate anything.
	 * @param ca
	 * @return a list of Chains that is built up from the Atoms in the ca array
	 * @throws StructureException
	 */
	public static final List<Chain> getAlignedModel(Atom[] ca){

		List<Chain> model = new ArrayList<Chain>();
		for ( Atom a: ca){

			Group g = a.getGroup();
			Chain parentC = g.getChain();

			Chain newChain = null;
			for ( Chain c :  model) {
				if ( c.getId().equals(parentC.getId())){
					newChain = c;
					break;
				}
			}
			if ( newChain == null){

				newChain = new ChainImpl();

				newChain.setId(parentC.getId());

				model.add(newChain);
			}

			newChain.addGroup(g);

		}

		return model;
	}


	/** Get an artifical Structure containing both chains.
	 * Does NOT rotate anything
	 * @param ca1
	 * @param ca2
	 * @return a structure object containing two models, one for each set of Atoms.
	 * @throws StructureException
	 */
	public static final Structure getAlignedStructure(Atom[] ca1, Atom[] ca2) throws StructureException{

		/* Previous implementation commented

		Structure s = new StructureImpl();


		List<Chain>model1 = getAlignedModel(ca1);
		List<Chain>model2 = getAlignedModel(ca2);
		s.addModel(model1);
		s.addModel(model2);

		return s;*/

		Structure s = new StructureImpl();

		List<Chain>model1 = getAlignedModel(ca1);
		s.addModel(model1);
		List<Chain> model2 = getAlignedModel(ca2);
		s.addModel(model2);

		return s;
	}

	/** Rotate the Atoms/Groups so they are aligned for the 3D visualisation
	 *
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @return an array of Groups that are transformed for 3D display
	 * @throws StructureException
	 */
	public static Group[] prepareGroupsForDisplay(AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws StructureException{


		if ( afpChain.getBlockRotationMatrix().length == 0 ) {
			// probably the alignment is too short!
			System.err.println("No rotation matrix found to rotate 2nd structure!");
			afpChain.setBlockRotationMatrix(new Matrix[]{Matrix.identity(3, 3)});
			afpChain.setBlockShiftVector(new Atom[]{new AtomImpl()});
		}

		// List of groups to be rotated according to the alignment
		Group[] twistedGroups = new Group[ ca2.length];

		//int blockNum = afpChain.getBlockNum();

		int i = -1;

		// List of groups from the structure not included in ca2 (e.g. ligands)
		// Will be rotated according to first block
		List<Group> hetatms2 = StructureTools.getUnalignedGroups(ca2);

		if (  (afpChain.getAlgorithmName().equals(FatCatRigid.algorithmName) ) || (afpChain.getAlgorithmName().equals(FatCatFlexible.algorithmName) ) ){

			for (Atom a: ca2){
				i++;
				twistedGroups[i]=a.getGroup();

			}

			twistedGroups = AFPTwister.twistOptimized(afpChain, ca1, ca2);

			//} else  if  (( blockNum == 1 ) || (afpChain.getAlgorithmName().equals(CeCPMain.algorithmName))) {
		} else {

			Matrix m   =  afpChain.getBlockRotationMatrix()[ 0];
			Atom shift =  afpChain.getBlockShiftVector()   [ 0 ];

			shiftCA2(afpChain, ca2, m,shift, twistedGroups);

		}

		if ( afpChain.getBlockNum() > 0){

			// Superimpose ligands relative to the first block
			if( hetatms2.size() > 0 ) {

				if ( afpChain.getBlockRotationMatrix().length > 0 ) {

					Matrix m1      = afpChain.getBlockRotationMatrix()[0];
					//m1.print(3,3);
					Atom   vector1 = afpChain.getBlockShiftVector()[0];
					//System.out.println("shift vector:" + vector1);

					for ( Group g : hetatms2){
						Calc.rotate(g, m1);
						Calc.shift(g,vector1);
					}
				}
			}
		}

		return twistedGroups;
	}

	/** only shift CA positions.
	 *
	 */
	public static void shiftCA2(AFPChain afpChain, Atom[] ca2,  Matrix m, Atom shift, Group[] twistedGroups) {

		int i = -1;
		for (Atom a: ca2){
			i++;
			Group g = a.getGroup();

			Calc.rotate(g,m);
			Calc.shift(g, shift);

			if (g.hasAltLoc()){
				for (Group alt: g.getAltLocs()){
					for (Atom alta : alt.getAtoms()){
						if ( g.getAtoms().contains(alta))
							continue;
						Calc.rotate(alta,m);
						Calc.shift(alta,shift);
					}
				}
			}
			twistedGroups[i]=g;
		}
	}
	
	/**
	 * Fill the aligned Atom arrays with the equivalent residues in the afpChain.
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @param ca1aligned
	 * @param ca2aligned
	 */
	public static void fillAlignedAtomArrays(AFPChain afpChain, Atom[] ca1, 
			Atom[] ca2, Atom[] ca1aligned, Atom[] ca2aligned) {
		
		int pos=0;
		int[] blockLens = afpChain.getOptLen();
		int[][][] optAln = afpChain.getOptAln();
		assert(afpChain.getBlockNum() <= optAln.length);

		for (int block=0; block < afpChain.getBlockNum(); block++) {
			for(int i=0;i<blockLens[block];i++) {
				int pos1 = optAln[block][0][i];
				int pos2 = optAln[block][1][i];
				Atom a1 = ca1[pos1];
				Atom a2 = (Atom) ca2[pos2].clone();
				ca1aligned[pos] = a1;
				ca2aligned[pos] = a2;
				pos++;
			}
		}

		// this can happen when we load an old XML serialization which did not support modern ChemComp representation of modified residues.
		if (pos != afpChain.getOptLength()){
			logger.warn("AFPChainScorer getTMScore: Problems reconstructing alignment! nr of loaded atoms is " + pos + " but should be " + afpChain.getOptLength());
			// we need to resize the array, because we allocated too many atoms earlier on.
			ca1aligned = (Atom[]) resizeArray(ca1aligned, pos);
			ca2aligned = (Atom[]) resizeArray(ca2aligned, pos);
		}
		
	}
	
	/**
	 * Find the alignment position with the highest atomic distance between the
	 * equivalent atomic positions of the arrays and remove it from the
	 * alignment.
	 * 
	 * @param afpChain
	 *            original alignment, will be modified
	 * @param ca1
	 *            atom array, will not be modified
	 * @param ca2
	 *            atom array, will not be modified
	 * @return the original alignment, with the alignment position at the
	 *         highest distance removed
	 * @throws StructureException
	 */
	public static AFPChain deleteHighestDistanceColumn(AFPChain afpChain,
			Atom[] ca1, Atom[] ca2) throws StructureException {

		int[][][] optAln = afpChain.getOptAln();

		int maxBlock = 0;
		int maxPos = 0;
		double maxDistance = Double.MIN_VALUE;

		for (int b = 0; b < optAln.length; b++) {
			for (int p = 0; p < optAln[b][0].length; p++) {
				Atom ca2clone = ca2[optAln[b][1][p]];
				Calc.rotate(ca2clone, afpChain.getBlockRotationMatrix()[b]);
				Calc.shift(ca2clone, afpChain.getBlockShiftVector()[b]);
				
				double distance = Calc.getDistance(ca1[optAln[b][0][p]],
						ca2clone);
				if (distance > maxDistance) {
					maxBlock = b;
					maxPos = p;
					maxDistance = distance;
				}
			}
		}

		return deleteColumn(afpChain, ca1, ca2, maxBlock, maxPos);
	}

	/**
	 * Delete an alignment position from the original alignment object.
	 * 
	 * @param afpChain
	 *            original alignment, will be modified
	 * @param ca1
	 *            atom array, will not be modified
	 * @param ca2
	 *            atom array, will not be modified
	 * @param block
	 *            block of the alignment position
	 * @param pos
	 *            position index in the block
	 * @return the original alignment, with the alignment position removed
	 * @throws StructureException
	 */
	public static AFPChain deleteColumn(AFPChain afpChain, Atom[] ca1,
			Atom[] ca2, int block, int pos) throws StructureException {

		// Check validity of the inputs
		if (afpChain.getBlockNum() <= block) {
			throw new IndexOutOfBoundsException(String.format(
					"Block index requested (%d) is higher than the total number of AFPChain blocks (%d).",
					block, afpChain.getBlockNum()));
		}
		if (afpChain.getOptAln()[block][0].length <= pos) {
			throw new IndexOutOfBoundsException(String.format(
					"Position index requested (%d) is higher than the total number of aligned position in the AFPChain block (%d).",
					block, afpChain.getBlockSize()[block]));
		}
		
		int[][][] optAln = afpChain.getOptAln();

		int[] newPos0 = new int[optAln[block][0].length - 1];
		int[] newPos1 = new int[optAln[block][1].length - 1];

		int position = 0;
		for (int p = 0; p < optAln[block][0].length; p++) {

			if (p == pos)
				continue;

			newPos0[position] = optAln[block][0][p];
			newPos1[position] = optAln[block][1][p];

			position++;
		}

		optAln[block][0] = newPos0;
		optAln[block][1] = newPos1;

		return AlignmentTools.replaceOptAln(optAln, afpChain, ca1, ca2);
	}
}
