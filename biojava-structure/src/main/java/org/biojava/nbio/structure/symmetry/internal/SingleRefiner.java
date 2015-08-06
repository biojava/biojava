package org.biojava.nbio.structure.symmetry.internal;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.TreeSet;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.align.util.AtomCache;

/**
 * Creates a refined alignment with the CE alternative self-alignment.
 * Needs the order of symmetry and assumes that the last subunit aligns
 * with the first, being thus a CLOSE symmetry.
 * 
 * @author Spencer Bliven
 * @author Aleix Lafita
 * 
 */
public class SingleRefiner implements Refiner {

	private int order;
	
	public SingleRefiner(int order) {
		this.order = order;
	}
	
	@Override
	public AFPChain refine(List<AFPChain> afpAlignments, Atom[] atoms)
			throws RefinerFailedException, StructureException {
		
		if (order == 1)	throw new RefinerFailedException(
				"Cannot refine an order 1 symmetry.");
		
		return refineSymmetry(afpAlignments.get(0), atoms, atoms, order);
	}
	
	/**
	 * Refines a CE-Symm alignment so that it is perfectly symmetric.
	 *
	 * The resulting alignment will have a one-to-one correspondance between
	 * aligned residues of each symmetric part.
	 * @param afpChain Input alignment from CE-Symm
	 * @param k Symmetry order. This can be guessed by {@link CeSymm#getSymmetryOrder(AFPChain)}
	 * @return The refined alignment
	 * @throws StructureException
	 */
	public static AFPChain refineSymmetry(AFPChain afpChain, Atom[] ca1, Atom[] ca2, int k) throws StructureException {
		// The current alignment
		Map<Integer, Integer> alignment = AlignmentTools.alignmentAsMap(afpChain);

		// Do the alignment
		Map<Integer, Integer> refined = refineSymmetry(alignment, k);
		
		//Substitute and partition the alignment
		AFPChain refinedAFP = AlignmentTools.replaceOptAln(afpChain, ca1, ca2, refined);
		return partitionAFPchain(refinedAFP, ca1, ca2, k);
	}

	/**
	 * Refines a CE-Symm alignment so that it is perfectly symmetric.
	 *
	 * The resulting alignment will have a one-to-one correspondance between
	 * aligned residues of each symmetric part.
	 * @param alignment The input alignment, as a map. This will be modified.
	 * @param k Symmetry order. This can be guessed by {@link CeSymm#getSymmetryOrder(AFPChain)}
	 * @return A modified map with the refined alignment
	 * @throws StructureException
	 */
	public static Map<Integer, Integer> refineSymmetry(Map<Integer, Integer> alignment,int k) throws StructureException {

		// Store scores
		Map<Integer, Double> scores = null;
		scores = initializeScores(alignment,scores, k);

		// Store eligible residues
		// Eligible if:
		//  1. score(x)>0
		//  2. f^K-1(x) is defined
		//	3. score(f^K-1(x))>0

		TreeSet<Integer> forwardLoops = new TreeSet<Integer>();
		TreeSet<Integer> backwardLoops = new TreeSet<Integer>();


		List<Integer> eligible = null;
		eligible = initializeEligible(alignment,scores,eligible,k,forwardLoops,backwardLoops);

		/* For future heap implementation
		Comparator<Integer> scoreComparator = new Comparator<Integer>() {
			@Override public int compare(Integer o1, Integer o2) {
				if(scores.containsKey(o1)) {
					if(scores.containsKey(o2)) {
						// If both have defined scores, compare the scores
						return scores.get(o1).compareTo(scores.get(o2));
					} else {
						// o2 has infinite score, so o1 < o2
						return -1;
					}
				} else {
					//o1 has infinite score
					if(scores.containsKey(o2)) {
						// o1 > o2
						return 1;
					} else {
						//both undefined
						return 0;
					}
				}
			}
		};
		PriorityQueue<Integer> heap = new PriorityQueue<Integer>(alignment.size(), scoreComparator);
		 */
		//int step = 0;
		while (!eligible.isEmpty()) {
			//System.out.format("Step %d: %s%n", ++step, AlignmentTools.toConciseAlignmentString(alignment));

			// Find eligible residue with lowest scores
			Integer bestRes = null;
			double bestResScore = Double.POSITIVE_INFINITY;
			for(Integer res : eligible) {
				Double score = scores.get(res);
				if (score != null && score < bestResScore) {
					bestResScore = score;
					bestRes = res;
				}
			}

			// Find f^k-1(bestRes)
			Integer resK1 = bestRes;
			for (int i = 0; i < k - 1; i++) {
				assert (resK1 != null);
				resK1 = alignment.get(resK1);

				// Update scores
				scores.put(resK1, 0.0);
			}
			scores.put(bestRes, 0.0);

			// Modify alignment
			alignment.put(resK1, bestRes);

			scores = initializeScores(alignment, scores, k);

			Map<Integer, Double> virginScores = initializeScores(alignment, null, k);
			if (scores.size() != virginScores.size()) {
				System.out.println("Size missmatch");
			} else {
				for (Integer key : scores.keySet()) {
					if (!virginScores.containsKey(key) || !scores.get(key).equals(virginScores.get(key))) {
						System.out.format("Mismatch %d -> %f/%f%n", key, scores.get(key), virginScores.get(key));
					}
				}
			}

			// Update eligible
			// TODO only update residues which could become ineligible
			eligible = initializeEligible(alignment, scores, eligible, k, forwardLoops, backwardLoops);

			// System.out.format("Modifying %d -> %d. %d now eligible.%n", resK1,bestRes,eligible.size());
		}
		//System.out.format("Step %d: %s%n", ++step, AlignmentTools.toConciseAlignmentString(alignment));

		// Remove remaining edges
		Iterator<Integer> alignmentIt = alignment.keySet().iterator();
		while (alignmentIt.hasNext()) {
			Integer res = alignmentIt.next();
			Double score = scores.get(res);
			if (score == null || score > 0.0) {
				alignmentIt.remove();
			}
		}
		//System.out.format("Step %d: %s%n", ++step, AlignmentTools.toConciseAlignmentString(alignment));

		return alignment;
	}

	/**
	 * Helper method to initialize eligible residues.
	 *
	 * Eligible if:
	 *  1. score(x)>0
	 *  2. f^K-1(x) is defined
	 *  3. score(f^K-1(x))>0
	 *  4. For all y, score(y)==0 implies sign(f^K-1(x)-y) == sign(x-f(y) )
	 * @param alignment The alignment with respect to which to calculate eligibility
	 * @param scores An up-to-date map from residues to their scores
	 * @param eligible Starting list of eligible residues. If null will be generated.
	 * @param k
	 * @param backwardLoops
	 * @param forwardLoops
	 * @return eligible after modification
	 */
	private static List<Integer> initializeEligible(Map<Integer, Integer> alignment,
			Map<Integer, Double> scores, List<Integer> eligible, int k, NavigableSet<Integer> forwardLoops, NavigableSet<Integer> backwardLoops) {
		// Eligible if:
		// 1. score(x)>0
		// 2. f^K-1(x) is defined
		// 3. score(f^K-1(x))>0
		// 4. For all y, score(y)==0 implies sign(f^K-1(x)-y) == sign(x-f(y) )
		// 5. Not in a loop of length less than k

		// Assume all residues are eligible to start
		if(eligible == null) {
			eligible = new LinkedList<Integer>(alignment.keySet());
		}

		// Precalculate f^K-1(x)
		// Map<Integer, Integer> alignK1 = AlignmentTools.applyAlignment(alignment, k-1);
		Map<Integer, Integer> alignK1 = applyAlignmentAndCheckCycles(alignment, k - 1, eligible);

		// Remove ineligible residues
		Iterator<Integer> eligibleIt = eligible.iterator();
		while(eligibleIt.hasNext()) {
			Integer res = eligibleIt.next();

			//  2. f^K-1(x) is defined
			if(!alignK1.containsKey(res)) {
				eligibleIt.remove();
				continue;
			}
			Integer k1 = alignK1.get(res);
			if(k1 == null) {
				eligibleIt.remove();
				continue;
			}

			//  1. score(x)>0
			Double score = scores.get(res);
			if(score == null || score <= 0.0) {
				eligibleIt.remove();

				// res is in a loop. Add it to the proper set
				if(res <= alignment.get(res)) {
					//forward
					forwardLoops.add(res);
				} else {
					//backward
					backwardLoops.add(res);
				}

				continue;
			}
			//	3. score(f^K-1(x))>0
			Double scoreK1 = scores.get(k1);
			if(scoreK1 == null || scoreK1 <= 0.0) {
				eligibleIt.remove();
				continue;
			}
		}


		// Now that loops are up-to-date, check for loop crossings
		eligibleIt = eligible.iterator();
		while(eligibleIt.hasNext()) {
			Integer res = eligibleIt.next();

			//4. For all y, score(y)==0 implies sign(f^K-1(x)-y) == sign(x-f(y) )
			//Test equivalent: All loop edges should be properly ordered wrt edge f^k-1(x) -> x

			Integer src = alignK1.get(res);

			if( src < res  ) {
				//forward
				// get interval [a,b) containing res
				Integer a = forwardLoops.floor(src);
				Integer b = forwardLoops.higher(src);

				// Ineligible unless f(a) < res < f(b)
				if(a != null && alignment.get(a) > res ) {
					eligibleIt.remove();
					continue;
				}
				if(b != null && alignment.get(b) < res ) {
					eligibleIt.remove();
					continue;
				}
			}
		}

		return eligible;
	}


	/**
	 * Like {@link AlignmentTools#applyAlignment(Map, int)}, returns a map of k applications of alignmentMap. However,
	 * it also sets loops of size less than k as ineligible.
	 *
	 * @param alignmentMap
	 *            f(x)
	 * @param k
	 * @param eligible
	 *            Eligible residues. Residues from small cycles are removed.
	 * @return f^k(x)
	 */
	private static Map<Integer, Integer> applyAlignmentAndCheckCycles(Map<Integer, Integer> alignmentMap, int k, List<Integer> eligible) {

		// Convert to lists to establish a fixed order (avoid concurrent modification)
		List<Integer> preimage = new ArrayList<Integer>(alignmentMap.keySet()); // currently unmodified
		List<Integer> image = new ArrayList<Integer>(preimage);

		for (int n = 1; n <= k; n++) {
			// apply alignment
			for (int i = 0; i < image.size(); i++) {
				final Integer pre = image.get(i);
				final Integer post = (pre == null ? null : alignmentMap.get(pre));
				image.set(i, post);

				// Make cycles ineligible
				if (post != null && post.equals(preimage.get(i))) {
					eligible.remove(preimage.get(i)); // Could be O(n) with List impl
				}
			}
		}

		Map<Integer, Integer> imageMap = new HashMap<Integer, Integer>(alignmentMap.size());

		// now populate with actual values
		for (int i = 0; i < preimage.size(); i++) {
			Integer pre = preimage.get(i);
			Integer postK = image.get(i);
			imageMap.put(pre, postK);
		}
		return imageMap;
	}

	/**
	 * Calculates all scores for an alignment
	 * @param alignment
	 * @param scores A mapping from residues to scores, which will be updated or
	 * 	created if null
	 * @return scores
	 */
	private static Map<Integer, Double> initializeScores(Map<Integer, Integer> alignment,
			Map<Integer, Double> scores, int k) {
		if(scores == null) {
			scores = new HashMap<Integer, Double>(alignment.size());
		} else {
			scores.clear();
		}
		Map<Integer,Integer> alignK = AlignmentTools.applyAlignment(alignment, k);

		// calculate input range
		int maxPre = Integer.MIN_VALUE;
		int minPre = Integer.MAX_VALUE;
		for(Integer pre : alignment.keySet()) {
			if(pre>maxPre) maxPre = pre;
			if(pre<minPre) minPre = pre;
		}

		for(Integer pre : alignment.keySet()) {
			Integer image = alignK.get(pre);

			// Use the absolute error score, |x - f^k(x)|
			double score = scoreAbsError(pre,image,minPre,maxPre);
			scores.put(pre, score);
		}
		return scores;
	}



	/**
	 * Calculate the score for a residue, specifically the Absolute Error
	 * 	score(x) = |x-f^k(x)|
	 *
	 * Also includes a small bias based on residue number, for uniqueness..
	 * @param pre x
	 * @param image f^k(x)
	 * @param minPre lowest possible residue number
	 * @param maxPre highest possible residue number
	 * @return
	 */
	private static double scoreAbsError(Integer pre, Integer image,int minPre,int maxPre) {
		// Use the absolute error score, |x - f^k(x)|
		double error;
		if(image == null) {
			error = Double.POSITIVE_INFINITY;
		} else {
			error = Math.abs(pre - image);
		}

		//TODO favor lower degree-in

		// Add fractional portion relative to sequence position, for uniqueness
		if(error > 0)
			error += (double)(pre-minPre)/(1+maxPre-minPre);

		return error;
	}
	
	public static void main(String[] args) {
		try {
			String name;

			name = "1itb.A"; // b-trefoil, C3
			//name = "1tim.A"; // tim-barrel, C8
			//name = "d1p9ha_"; // not rotational symmetry
			name = "3HKE.A"; // very questionable alignment
			//name = "d1jlya1"; // C3 with minimum RSSE at C6
			//name = "1YOX(A:95-160)";
			name = "3jut.A"; // b-trefoil FGF-1, C3
			name = "2jaj.A"; //C5
			name = "d1poqa_"; //superantigen Ypm
			name = "1hiv"; // hiv protease

			// Write a SIF (graph) file
			boolean writeSIF = false;
			String path = "/Users/blivens/dev/bourne/symmetry/refinement/";
			if( !(new File(path)).exists() ) {
				writeSIF = false;
			}

			AtomCache cache = new AtomCache();
			Atom[] ca1 = cache.getAtoms(name);
			Atom[] ca2 = cache.getAtoms(name);

			StructureAlignmentFactory.addAlgorithm(new CeSymm());
			CeSymm ce = (CeSymm) StructureAlignmentFactory.getAlgorithm(CeSymm.algorithmName);

			// CE-Symm alignment
			long startTime = System.currentTimeMillis();
			AFPChain afpChain = ce.align(ca1, ca2);
			long alignTime = System.currentTimeMillis()-startTime;

			afpChain.setName1(name);
			afpChain.setName2(name);

			System.out.format("Alignment took %dms%n", alignTime);

			startTime = System.currentTimeMillis();
			int symm = new SequenceFunctionOrderDetector().calculateOrder(afpChain, ca1);
			long orderTime = System.currentTimeMillis()-startTime;
			System.out.println("Symmetry="+symm);

			System.out.format("Finding order took %dms%n", orderTime);

			//Output SIF file
			String filename = path+name+".sif";
			Writer out = null;
			if(writeSIF) {
				System.out.println("Writing alignment to "+filename);
				out = new FileWriter(filename);
				alignmentToSIF(out, afpChain, ca1, ca2, "bb","ur");
			}

			//Refine alignment
			startTime = System.currentTimeMillis();
			AFPChain refinedAFP = refineSymmetry(afpChain, ca1, ca2, symm);
			long refineTime = System.currentTimeMillis()-startTime;

			refinedAFP.setAlgorithmName(refinedAFP.getAlgorithmName()+"-refine");

			//Output refined to SIF
			System.out.format("Refinement took %dms%n", refineTime);
			if(writeSIF) {
				alignmentToSIF(out, refinedAFP, ca1, ca2,"bb","rr");
				out.close();
			}

			//display jmol of alignment
			System.out.println("Original rmsd:"+afpChain.getTotalRmsdOpt());
			System.out.println("New rmsd:"+refinedAFP.getTotalRmsdOpt());

		} catch (Exception e) {
			e.printStackTrace();
		}
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
	private static void alignmentToSIF(Writer out,
			AFPChain afpChain, Atom[] ca1,Atom[] ca2, String backboneInteraction, String alignmentInteraction) throws IOException {
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

				String node1 = name1+rn1.getChainId()+rn1.toString();
				String node2 = name2+rn2.getChainId()+rn2.toString();

				out.write(String.format("%s\t%s\t%s\n",node1, alignmentInteraction, node2));
			}
		}

		// Print first backbone edges
		ResidueNumber rn = ca1[0].getGroup().getResidueNumber();
		String last = name1+rn.getChainId()+rn.toString();
		for(int i=1;i<ca1.length;i++) {
			rn = ca1[i].getGroup().getResidueNumber();
			String curr = name1+rn.getChainId()+rn.toString();
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
			last = name2+rn.getChainId()+rn.toString();
			for(int i=1;i<ca2.length;i++) {
				rn = ca2[i].getGroup().getResidueNumber();
				String curr = name2+rn.getChainId()+rn.toString();
				out.write(String.format("%s\t%s\t%s\n",last, backboneInteraction, curr));
				last = curr;
			}
		}

	}
	
	/**
	 *  Partitions an afpChain alignment into order blocks of aligned residues.
	 */
	private static AFPChain partitionAFPchain(AFPChain afpChain, Atom[] ca1, Atom[] ca2, int order) throws StructureException{
		
		int[][][] newAlgn = new int[order][][];
		int subunitLen = afpChain.getOptLength()/order;
		
		//Extract all the residues considered in the first chain of the alignment
		List<Integer> alignedRes = new ArrayList<Integer>();
		for (int su=0; su<afpChain.getBlockNum(); su++){
			for (int i=0; i<afpChain.getOptLen()[su]; i++){
				alignedRes.add(afpChain.getOptAln()[su][0][i]);
			}
		}
		
		//Build the new alignment
		for (int su=0; su<order; su++){
			newAlgn[su] = new int[2][];
			newAlgn[su][0] = new  int[subunitLen];
			newAlgn[su][1] = new  int[subunitLen];
			for (int i=0; i<subunitLen; i++){
				newAlgn[su][0][i] = alignedRes.get(subunitLen*su+i);
				newAlgn[su][1][i] = alignedRes.get((subunitLen*(su+1)+i)%alignedRes.size());
			}
		}
		
		return AlignmentTools.replaceOptAln(newAlgn, afpChain, ca1, ca2);
	}
}
