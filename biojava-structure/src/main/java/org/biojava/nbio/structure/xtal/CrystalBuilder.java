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
package org.biojava.nbio.structure.xtal;


import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.contact.AtomContactSet;
import org.biojava.nbio.structure.contact.StructureInterface;
import org.biojava.nbio.structure.contact.StructureInterfaceList;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;
import java.util.*;


/**
 * A class containing methods to find interfaces in a given crystallographic Structure by
 * reconstructing the crystal lattice through application of symmetry operators
 *
 * @author Jose Duarte
 *
 */

public class CrystalBuilder {

	// Default number of cell neighbors to try in interface search (in 3 directions of space).
	// In the search, only bounding box overlaps are tried, thus there's not so much overhead in adding
	// more cells. We actually tested it and using numCells from 1 to 10 didn't change runtimes at all.
	// Examples with interfaces in distant neighbor cells:
	//   2nd neighbors: 3hz3, 1wqj, 2de3, 1jcd
	//   3rd neighbors: 3bd3, 1men, 2gkp, 1wui
	//   5th neighbors: 2ahf, 2h2z
	//   6th neighbors: 1was (in fact interfaces appear only at 5th neighbors for it)
	// Maybe this could be avoided by previously translating the given molecule to the first cell,
	// BUT! some bona fide cases exist, e.g. 2d3e: it is properly placed at the origin but the molecule
	// is enormously long in comparison with the dimensions of the unit cell, some interfaces come at the 7th neighbor.
	// After a scan of the whole PDB (Oct 2013) using numCells=50, the highest one was 4jgc with
	// interfaces up to the 11th neighbor. Other high ones (9th neighbors) are 4jbm and 4k3t.
	// We set the default value to 12 based on that (having not seen any difference in runtime)
	public static final int DEF_NUM_CELLS = 12;

	/**
	 * Default maximum distance between two chains to be considered an interface.
	 * @see #getUniqueInterfaces(double)
	 */
	public static final double DEFAULT_INTERFACE_DISTANCE_CUTOFF = 5.5;

	public static final Matrix4d IDENTITY = new Matrix4d(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);


	/**
	 * Whether to consider HETATOMs in contact calculations
	 */
	private static final boolean INCLUDE_HETATOMS = true;

	private Structure structure;
	private PDBCrystallographicInfo crystallographicInfo;
	private int numPolyChainsAu;
	private int numOperatorsSg;
	private Map<String,Matrix4d> chainNcsOps = null;
	private Map<String,String> chainOrigNames = null;

	private static final Logger logger = LoggerFactory.getLogger(CrystalBuilder.class);

	private int numCells;

	private ArrayList<CrystalTransform> visitedCrystalTransforms;
	private Map<String,Map<Matrix4d,StructureInterface>> visitedNcsChainPairs = null;

	private boolean searchBeyondAU;
	private Matrix4d[] ops;

	/**
	 * Special constructor for NCS-aware CrystalBuilder.
	 * The output list of interfaces will be pre-clustered by NCS-equivalence.
	 * Run {@link CrystalBuilder#expandNcsOps(Structure, Map, Map)} first to extend the AU
	 * and get the equivalence information.
	 * @param structure
	 *          NCS-extended structure
	 * @param chainOrigNames
	 *          chain names mapped to the original chain names (pre-NCS extension)
	 * @param chainNcsOps
	 *          chain names mapped to the ncs operators that was used to generate them
	 * @since 5.0.0
	 */
	public CrystalBuilder(Structure structure, Map<String,String> chainOrigNames, Map<String,Matrix4d> chainNcsOps) {
		this(structure);
		this.chainOrigNames = chainOrigNames;
		this.chainNcsOps = chainNcsOps;
	}

	public CrystalBuilder(Structure structure) {
		this.structure = structure;
		this.crystallographicInfo = structure.getCrystallographicInfo();
		this.numPolyChainsAu = structure.getPolyChains().size();

		this.searchBeyondAU = false;
		if (structure.isCrystallographic()) {

			this.searchBeyondAU = true;

			// we need to check space group not null for the cases where the entry is crystallographic but
			// the space group is not a standard one recognized by biojava, e.g. 1mnk (SG: 'I 21')
			if (this.crystallographicInfo.isNonStandardSg()) {
				logger.warn("Space group is non-standard, will only calculate asymmetric unit interfaces.");
				this.searchBeyondAU = false;
			}

			// just in case we still check for space group null (a user pdb file could potentially be crystallographic and have no space group)
			if (this.crystallographicInfo.getSpaceGroup() == null) {
				logger.warn("Space group is null, will only calculate asymmetric unit interfaces.");
				this.searchBeyondAU = false;
			}

			// we need to check crystal cell not null for the rare cases where the entry is crystallographic but
			// the crystal cell is not given, e.g. 2i68, 2xkm, 4bpq
			if (this.crystallographicInfo.getCrystalCell() == null) {
				logger.warn("Could not find a crystal cell definition, will only calculate asymmetric unit interfaces.");
				this.searchBeyondAU = false;
			}

			// check for cases like 4hhb that are in a non-standard coordinate frame convention, see https://github.com/eppic-team/owl/issues/4
			if (this.crystallographicInfo.isNonStandardCoordFrameConvention()) {
				logger.warn("Non-standard coordinate frame convention, will only calculate asymmetric unit interfaces.");
				this.searchBeyondAU = false;
			}
		}

		if (this.searchBeyondAU) {
			// explore the crystal
			this.numOperatorsSg = this.crystallographicInfo.getSpaceGroup().getMultiplicity();
			this.ops = this.crystallographicInfo.getTransformationsOrthonormal();
		} else {
			// look for contacts within structure as given
			this.numOperatorsSg = 1;
			this.ops = new Matrix4d[1];
			this.ops[0] = new Matrix4d(IDENTITY);
		}

		this.numCells = DEF_NUM_CELLS;

	}


	/**
	 * @return true if this CrystalBuilder is NCS-aware.
	 * @since 5.0.0
	 */
	public boolean hasNcsOps() {
		return chainNcsOps != null;
	}

	/**
	 * Set the number of neighboring crystal cells that will be used in the search for contacts
	 * @param numCells
	 */
	public void setNumCells(int numCells) {
		this.numCells = numCells;
	}

	private void initialiseVisited() {
		visitedCrystalTransforms = new ArrayList<>();
		if(this.hasNcsOps()) {
			visitedNcsChainPairs = new HashMap<>();
		}
	}

	/**
	 * Returns the list of unique interfaces that the given Structure has upon
	 * generation of all crystal symmetry mates. An interface is defined as any pair of chains
	 * that contact, i.e. for which there is at least a pair of atoms (one from each chain) within
	 * the default cutoff distance.
	 * @return
	 * @see #DEFAULT_INTERFACE_DISTANCE_CUTOFF
	 */
	public StructureInterfaceList getUniqueInterfaces() {
		return getUniqueInterfaces(DEFAULT_INTERFACE_DISTANCE_CUTOFF);
	}

	/**
	 * Returns the list of unique interfaces that the given Structure has upon
	 * generation of all crystal symmetry mates. An interface is defined as any pair of chains
	 * that contact, i.e. for which there is at least a pair of atoms (one from each chain) within
	 * the given cutoff distance.
	 * @param cutoff the distance cutoff for 2 chains to be considered in contact
	 * @return
	 */
	public StructureInterfaceList getUniqueInterfaces(double cutoff) {


		StructureInterfaceList set = new StructureInterfaceList();

		// certain structures in the PDB are not macromolecules (contain no polymeric chains at all), e.g. 1ao2
		// with the current mmCIF parsing, those will be empty since purely non-polymeric chains are removed
		// see commit e9562781f23da0ebf3547146a307d7edd5741090
		if (numPolyChainsAu==0) {
			logger.warn("No chains present in the structure! No interfaces will be calculated");
			return set;
		}



		// initialising the visited ArrayList for keeping track of symmetry redundancy
		initialiseVisited();



		// the isCrystallographic() condition covers 3 cases:
		// a) entries with expMethod X-RAY/other diffraction and defined crystalCell (most usual case)
		// b) entries with expMethod null but defined crystalCell (e.g. PDB file with CRYST1 record but no expMethod annotation)
		// c) entries with expMethod not X-RAY (e.g. NMR) and defined crystalCell (NMR entries do have a dummy CRYST1 record "1 1 1 90 90 90 P1")
		// d) isCrystallographic will be false if the structure is crystallographic but the space group was not recognized


		calcInterfacesCrystal(set, cutoff);

		return set;
	}

	/**
	 * Calculate interfaces between original asymmetric unit and neighboring
	 * whole unit cells, including the original full unit cell i.e. i=0,j=0,k=0
	 * @param set
	 * @param cutoff
	 */
	private void calcInterfacesCrystal(StructureInterfaceList set, double cutoff) {


		// initialising debugging vars
		long start = -1;
		long end = -1;
		int trialCount = 0;
		int skippedRedundant = 0;
		int skippedAUsNoOverlap = 0;
		int skippedChainsNoOverlap = 0;
		int skippedSelfEquivalent = 0;

		// The bounding boxes of all AUs of the unit cell
		UnitCellBoundingBox bbGrid = new UnitCellBoundingBox(numOperatorsSg, numPolyChainsAu);;
		// we calculate all the bounds of each of the asym units, those will then be reused and translated
		bbGrid.setBbs(structure, ops, INCLUDE_HETATOMS);


		// if not crystallographic there's no search to do in other cells, only chains within "AU" will be checked
		if (!searchBeyondAU) numCells = 0;

		boolean verbose = logger.isDebugEnabled();

		if (verbose) {
			trialCount = 0;
			start= System.currentTimeMillis();
			int neighbors = (2*numCells+1)*(2*numCells+1)*(2*numCells+1)-1;
			int auTrials = (numPolyChainsAu*(numPolyChainsAu-1))/2;
			int trials = numPolyChainsAu*numOperatorsSg*numPolyChainsAu*neighbors;
			logger.debug("Chain clash trials within original AU: "+auTrials);
			logger.debug(
					"Chain clash trials between the original AU and the neighbouring "+neighbors+
					" whole unit cells ("+numCells+" neighbours)" +
					"(2x"+numPolyChainsAu+"chains x "+numOperatorsSg+"AUs x "+neighbors+"cells) : "+trials);
			logger.debug("Total trials: "+(auTrials+trials));
		}

		List<Chain> polyChains = structure.getPolyChains();

		for (int a=-numCells;a<=numCells;a++) {
			for (int b=-numCells;b<=numCells;b++) {
				for (int c=-numCells;c<=numCells;c++) {

					Point3i trans = new Point3i(a,b,c);
					Vector3d transOrth = new Vector3d(a,b,c);
					if (a!=0 || b!=0 || c!=0) {
						// we avoid doing the transformation for 0,0,0 (in case it's not crystallographic)
						this.crystallographicInfo.getCrystalCell().transfToOrthonormal(transOrth);
					}

					UnitCellBoundingBox bbGridTrans = bbGrid.getTranslatedBbs(transOrth);

					for (int n=0;n<numOperatorsSg;n++) {

						// short-cut strategies
						// 1) we skip first of all if the bounding boxes of the AUs don't overlap
						if (!bbGrid.getAuBoundingBox(0).overlaps(bbGridTrans.getAuBoundingBox(n), cutoff)) {
							skippedAUsNoOverlap++;
							continue;
						}

						// 2) we check if we didn't already see its equivalent symmetry operator partner
						CrystalTransform tt = new CrystalTransform(this.crystallographicInfo.getSpaceGroup(), n);
						tt.translate(trans);
						if (isRedundantTransform(tt)) {
							skippedRedundant++;
							continue;
						}
						addVisitedTransform(tt);


						boolean selfEquivalent = false;

						// 3) an operator can be "self redundant" if it is the inverse of itself (involutory, e.g. all pure 2-folds with no translation)
						if (tt.isEquivalent(tt)) {
							logger.debug("Transform "+tt+" is equivalent to itself, will skip half of i-chains to j-chains comparisons");
							// in this case we can't skip the operator, but we can skip half of the matrix comparisons e.g. j>i
							// we set a flag and do that within the loop below
							selfEquivalent = true;
						}

						StringBuilder builder = null;
						if (verbose) builder = new StringBuilder(String.valueOf(tt)).append(" ");

						// Now that we know that boxes overlap and operator is not redundant, we have to go to the details
						int contactsFound = 0;

						for (int j=0;j<numPolyChainsAu;j++) {

							for (int i=0;i<numPolyChainsAu;i++) { // we only have to compare the original asymmetric unit to every full cell around

								if(selfEquivalent && (j>i)) {
									// in case of self equivalency of the operator we can safely skip half of the matrix
									skippedSelfEquivalent++;
									continue;
								}
								// special case of original AU, we don't compare a chain to itself
								if (n==0 && a==0 && b==0 && c==0 && i==j) continue;

								// before calculating the AtomContactSet we check for overlap, then we save putting atoms into the grid
								if (!bbGrid.getChainBoundingBox(0,i).overlaps(bbGridTrans.getChainBoundingBox(n,j),cutoff)) {
									skippedChainsNoOverlap++;
									if (verbose) {
										builder.append(".");
									}
									continue;
								}

								trialCount++;

								// finally we've gone through all short-cuts and the 2 chains seem to be close enough:
								// we do the calculation of contacts
								Chain chaini = polyChains.get(i);
								Chain chainj = polyChains.get(j);

								if (n!=0 || a!=0 || b!=0 || c!=0) {
									Matrix4d mJCryst = new Matrix4d(ops[n]);
									translate(mJCryst, transOrth);
									chainj = (Chain)chainj.clone();
									Calc.transform(chainj,mJCryst);
								}

								StructureInterface interf = calcContacts(chaini, chainj, cutoff, tt, builder);
								if (interf == null) {
									continue;
								}

								contactsFound++;
								if(this.hasNcsOps()) {
									StructureInterface interfNcsRef = findNcsRef(interf);
									set.addNcsEquivalent(interf,interfNcsRef);
								} else {
									set.add(interf);
								}
							}
						}

						if( verbose ) {
							if (a==0 && b==0 && c==0 && n==0)
								builder.append(" "+contactsFound+"("+(numPolyChainsAu*(numPolyChainsAu-1))/2+")");
							else if (selfEquivalent)
								builder.append(" "+contactsFound+"("+(numPolyChainsAu*(numPolyChainsAu+1))/2+")");
							else
								builder.append(" "+contactsFound+"("+numPolyChainsAu*numPolyChainsAu+")");

							logger.debug(builder.toString());
						}
					}
				}
			}
		}

		end = System.currentTimeMillis();
		logger.debug("\n"+trialCount+" chain-chain clash trials done. Time "+(end-start)/1000+"s");
		logger.debug("  skipped (not overlapping AUs)       : "+skippedAUsNoOverlap);
		logger.debug("  skipped (not overlapping chains)    : "+skippedChainsNoOverlap);
		logger.debug("  skipped (sym redundant op pairs)    : "+skippedRedundant);
		logger.debug("  skipped (sym redundant self op)     : "+skippedSelfEquivalent);
		logger.debug("Found "+set.size()+" interfaces.");
	}


	/**
	 * Checks whether given interface is NCS-redundant, i.e., an identical interface between NCS copies of
	 * these molecules has already been seen, and returns this (reference) interface.
	 *
	 * @param interf
	 *          StructureInterface
	 * @return  already seen interface that is NCS-equivalent to interf,
	 *          null if such interface is not found.
	 */
	private StructureInterface findNcsRef(StructureInterface interf) {
		if (!this.hasNcsOps()) {
			return null;
		}
		String chainIName = interf.getMoleculeIds().getFirst();
		String iOrigName = chainOrigNames.get(chainIName);

		String chainJName = interf.getMoleculeIds().getSecond();
		String jOrigName = chainOrigNames.get(chainJName);

		Matrix4d mJCryst;
		if(this.searchBeyondAU) {
			mJCryst = interf.getTransforms().getSecond().getMatTransform();
			mJCryst = crystallographicInfo.getCrystalCell().transfToOrthonormal(mJCryst);
		} else {
			mJCryst = IDENTITY;
		}

		// Let X1,...Xn be the original coords, before NCS transforms (M1...Mk)
		// current chain i: M_i * X_i
		// current chain j: Cn * M_j * X_j

		// transformation to bring chain j near X_i: M_i^(-1) * Cn * M_j
		// transformation to bring chain i near X_j: (Cn * M_j)^(-1) * M_i = (M_i^(-1) * Cn * M_j)^(-1)

		Matrix4d mChainIInv = new Matrix4d(chainNcsOps.get(chainIName));
		mChainIInv.invert();

		Matrix4d mJNcs = new Matrix4d(chainNcsOps.get(chainJName));

		Matrix4d j2iNcsOrigin = new Matrix4d(mChainIInv);
		j2iNcsOrigin.mul(mJCryst);
		//overall transformation to bring current chainj from its NCS origin to i's
		j2iNcsOrigin.mul(mJNcs);

		//overall transformation to bring current chaini from its NCS origin to j's
		Matrix4d i2jNcsOrigin = new Matrix4d(j2iNcsOrigin);
		i2jNcsOrigin.invert();

		String matchChainIdsIJ = iOrigName + jOrigName;
		String matchChainIdsJI = jOrigName + iOrigName;

		// same original chain names
		Optional<Matrix4d> matchDirect =
				visitedNcsChainPairs.computeIfAbsent(matchChainIdsIJ, k-> new HashMap<>()).entrySet().stream().
					map(r->r.getKey()).
					filter(r->r.epsilonEquals(j2iNcsOrigin,0.01)).
					findFirst();

		Matrix4d matchMatrix = matchDirect.orElse(null);
		String matchChainIds = matchChainIdsIJ;

		if(matchMatrix == null) {
			// reversed original chain names with inverted transform
			Optional<Matrix4d> matchInverse =
					visitedNcsChainPairs.computeIfAbsent(matchChainIdsJI, k-> new HashMap<>()).entrySet().stream().
					map(r->r.getKey()).
					filter(r->r.epsilonEquals(i2jNcsOrigin,0.01)).
					findFirst();
			matchMatrix = matchInverse.orElse(null);
			matchChainIds = matchChainIdsJI;
		}

		StructureInterface matchInterface = null;

		if (matchMatrix == null) {
			visitedNcsChainPairs.get(matchChainIdsIJ).put(j2iNcsOrigin,interf);
		} else {
			matchInterface = visitedNcsChainPairs.get(matchChainIds).get(matchMatrix);
		}

		return matchInterface;
	}

	private StructureInterface calcContacts(Chain chaini, Chain chainj, double cutoff, CrystalTransform tt, StringBuilder builder) {
		// note that we don't consider hydrogens when calculating contacts
		AtomContactSet graph = StructureTools.getAtomsInContact(chaini, chainj, cutoff, INCLUDE_HETATOMS);

		if (graph.size()>0) {
			if (builder != null) builder.append("x");

			CrystalTransform transf = new CrystalTransform(this.crystallographicInfo.getSpaceGroup());
			StructureInterface interf = new StructureInterface(
					StructureTools.getAllAtomArray(chaini), StructureTools.getAllAtomArray(chainj),
					chaini.getName(), chainj.getName(),
					graph,
					transf, tt);

			return interf;

		} else {
			if (builder != null) builder.append("o");
			return null;
		}
	}

	private void addVisitedTransform(CrystalTransform tt) {
		visitedCrystalTransforms.add(tt);
	}

	/**
	 * Checks whether given transformId/translation is symmetry redundant
	 * Two transformations are symmetry redundant if their matrices (4d) multiplication gives the identity, i.e.
	 * if one is the inverse of the other.
	 * @param tt
	 * @return
	 */
	private boolean isRedundantTransform(CrystalTransform tt) {

		Iterator<CrystalTransform> it = visitedCrystalTransforms.iterator();
		while (it.hasNext()) {
			CrystalTransform v = it.next();

			if (tt.isEquivalent(v)) {

				logger.debug("Skipping redundant transformation: "+tt+", equivalent to "+v);

				// there's only 1 possible equivalent partner for each visited matrix
				// (since the equivalent is its inverse matrix and the inverse matrix is unique)
				// thus once the partner has been seen, we don't need to check it ever again
				it.remove();

				return true;
			}
		}

		return false;
	}

	public void translate(Matrix4d m, Vector3d translation) {
		m.m03 = m.m03+translation.x;
		m.m13 = m.m13+translation.y;
		m.m23 = m.m23+translation.z;
	}

	/**
	 * Apply the NCS operators in the given Structure adding new chains as needed.
	 * All chains are (re)assigned ids of the form: original_chain_id+ncs_operator_index+"n".
	 * @param structure
	 *          the structure to expand
	 * @param chainOrigNames
	 *          new chain names mapped to the original chain names
	 * @param chainNcsOps
	 *          new chain names mapped to the ncs operators that was used to generate them
	 * @since 5.0.0
	 */
	public static void expandNcsOps(Structure structure, Map<String,String> chainOrigNames, Map<String,Matrix4d> chainNcsOps) {
		PDBCrystallographicInfo xtalInfo = structure.getCrystallographicInfo();
		if (xtalInfo ==null) return;

		if (xtalInfo.getNcsOperators()==null || xtalInfo.getNcsOperators().length==0)
			return;

		List<Chain> chainsToAdd = new ArrayList<>();

		Matrix4d identity = new Matrix4d();
		identity.setIdentity();

		Matrix4d[] ncsOps = xtalInfo.getNcsOperators();

		for (Chain c:structure.getChains()) {
			String cOrigId = c.getId();
			String cOrigName = c.getName();

			for (int iOperator = 0; iOperator < ncsOps.length; iOperator++) {
				Matrix4d m = ncsOps[iOperator];

				Chain clonedChain = (Chain)c.clone();
				String newChainId = cOrigId+(iOperator+1)+"n";
				String newChainName = cOrigName+(iOperator+1)+"n";
				clonedChain.setId(newChainId);
				clonedChain.setName(newChainName);

				setChainIdsInResidueNumbers(clonedChain, newChainName);
				Calc.transform(clonedChain, m);

				chainsToAdd.add(clonedChain);
				c.getEntityInfo().addChain(clonedChain);

				chainOrigNames.put(newChainName,cOrigName);
				chainNcsOps.put(newChainName,m);
			}

			chainNcsOps.put(cOrigName,identity);
			chainOrigNames.put(cOrigName,cOrigName);
		}

		chainsToAdd.forEach(structure::addChain);
	}

	/**
	 * Auxiliary method to reset chain ids of residue numbers in a chain.
	 * Used when cloning chains and resetting their ids: one needs to take care of
	 * resetting the ids within residue numbers too.
	 * @param c
	 * @param newChainName
	 */
	private static void setChainIdsInResidueNumbers(Chain c, String newChainName) {
		for (Group g:c.getAtomGroups()) {
			g.setResidueNumber(newChainName, g.getResidueNumber().getSeqNum(), g.getResidueNumber().getInsCode());
		}
		for (Group g:c.getSeqResGroups()) {
			if (g.getResidueNumber()==null) continue;
			g.setResidueNumber(newChainName, g.getResidueNumber().getSeqNum(), g.getResidueNumber().getInsCode());
		}
	}

}
