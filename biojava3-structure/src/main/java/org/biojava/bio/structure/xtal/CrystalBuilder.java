package org.biojava.bio.structure.xtal;


import java.util.ArrayList;
import java.util.Iterator;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;

import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.PDBCrystallographicInfo;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.contact.AtomContactSet;
import org.biojava.bio.structure.contact.StructureInterface;
import org.biojava.bio.structure.contact.StructureInterfaceList;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;



/**
 * A class containing methods to find interfaces in a given crystallographic Structure by
 * reconstructing the crystal lattice through application of symmetry operators
 *
 * @author duarte_j
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
	private static final int DEF_NUM_CELLS = 12;

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
	private int numChainsAu;
	private int numOperatorsSg;

	private static final Logger logger = LoggerFactory.getLogger(CrystalBuilder.class);

	private int numCells;

	private ArrayList<CrystalTransform> visited;
	
	private boolean isCrystallographic;



	public CrystalBuilder(Structure structure) {
		this.structure = structure;
		this.crystallographicInfo = structure.getCrystallographicInfo();

		this.numChainsAu = structure.getChains().size();
		this.numOperatorsSg = 1;
		this.isCrystallographic = false;
				
		
		if (this.crystallographicInfo.getSpaceGroup()==null) {
			logger.warn("Could not find a space group, will only calculate asymmetric unit interfaces.");
		}
		// we need to check space group not null for the cases where the entry is crystallographic but 
		// the space group is not a standard one recognized by biojava, e.g. 1mnk (SG: 'I 21')
		if (structure.isCrystallographic() && this.crystallographicInfo.getSpaceGroup()!=null) {
			this.numOperatorsSg = this.crystallographicInfo.getSpaceGroup().getMultiplicity();
			this.isCrystallographic = true;
		}

		this.numCells = DEF_NUM_CELLS;

	}

	/**
	 * Set the number of neighboring crystal cells that will be used in the search for contacts
	 * @param numCells
	 */
	public void setNumCells(int numCells) {
		this.numCells = numCells;
	}

	private void initialiseVisited() {
		visited = new ArrayList<CrystalTransform>();
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

		// initialising the visited ArrayList for keeping track of symmetry redundancy
		initialiseVisited();



		// the isCrystallographic() condition covers 3 cases:
		// a) entries with expMethod X-RAY/other diffraction and defined crystalCell (most usual case)
		// b) entries with expMethod null but defined crystalCell (e.g. PDB file with CRYST1 record but no expMethod annotation)
		// c) entries with expMethod not X-RAY (e.g. NMR) and defined crystalCell (NMR entries do have a dummy CRYST1 record "1 1 1 90 90 90 P1")
		// d) isCrystallographic will be false if the structure is crystallographic but the space group was not recognized


		calcInterfacesCrystal(set, cutoff, isCrystallographic);


		return set;
	}

	/**
	 * Calculate interfaces between original asymmetric unit and neighboring
	 * whole unit cells, including the original full unit cell i.e. i=0,j=0,k=0
	 * @param set
	 * @param cutoff
	 */
	private void calcInterfacesCrystal(StructureInterfaceList set, double cutoff, boolean isCrystallographic) {


		// initialising debugging vars
		long start = -1;
		long end = -1;
		int trialCount = 0;
		int skippedRedundant = 0;
		int skippedAUsNoOverlap = 0;
		int skippedChainsNoOverlap = 0;
		int skippedSelfEquivalent = 0;


		Matrix4d[] ops = null;
		if (isCrystallographic) {
			ops = structure.getCrystallographicInfo().getTransformationsOrthonormal();
		} else {
			ops = new Matrix4d[1];
			ops[0] = new Matrix4d(IDENTITY);
		}

		// The bounding boxes of all AUs of the unit cell
		UnitCellBoundingBox bbGrid = new UnitCellBoundingBox(numOperatorsSg, numChainsAu);;
		// we calculate all the bounds of each of the asym units, those will then be reused and translated
		bbGrid.setBbs(structure, ops, INCLUDE_HETATOMS);


		// if not crystallographic there's no search to do in other cells, only chains within "AU" will be checked
		if (!isCrystallographic) numCells = 0;

		boolean verbose = logger.isDebugEnabled();

		if (verbose) {
			trialCount = 0;
			start= System.currentTimeMillis();
			int neighbors = (2*numCells+1)*(2*numCells+1)*(2*numCells+1)-1;
			int auTrials = (numChainsAu*(numChainsAu-1))/2;
			int trials = numChainsAu*numOperatorsSg*numChainsAu*neighbors;
			logger.debug("Chain clash trials within original AU: "+auTrials);
			logger.debug(
					"Chain clash trials between the original AU and the neighbouring "+neighbors+
					" whole unit cells ("+numCells+" neighbours)" +
					"(2x"+numChainsAu+"chains x "+numOperatorsSg+"AUs x "+neighbors+"cells) : "+trials);
			logger.debug("Total trials: "+(auTrials+trials));
		}


		for (int a=-numCells;a<=numCells;a++) {
			for (int b=-numCells;b<=numCells;b++) {
				for (int c=-numCells;c<=numCells;c++) {

					Point3i trans = new Point3i(a,b,c);
					Vector3d transOrth = new Vector3d(a,b,c);
					if (a!=0 || b!=0 || c!=0)
						// we avoid doing the transformation for 0,0,0 (in case it's not crystallographic)
						this.crystallographicInfo.getCrystalCell().transfToOrthonormal(transOrth);
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
						if (isRedundant(tt)) {
							skippedRedundant++;
							continue;
						}
						addVisited(tt);


						boolean selfEquivalent = false;

						// 3) an operator can be "self redundant" if it is the inverse of itself (involutory, e.g. all pure 2-folds with no translation)
						if (tt.isEquivalent(tt)) {
							logger.debug("Transform "+tt+" is equivalent to itself, will skip half of i-chains to j-chains comparisons");
							// in this case we can't skip the operator, but we can skip half of the matrix comparisons e.g. j>i
							// we set a flag and do that within the loop below
							selfEquivalent = true;
						}

						StringBuilder builder = null;
						if (verbose) builder = new StringBuilder(tt+" ");

						// Now that we know that boxes overlap and operator is not redundant, we have to go to the details
						int contactsFound = 0;

						for (int j=0;j<numChainsAu;j++) {
							for (int i=0;i<numChainsAu;i++) { // we only have to compare the original asymmetric unit to every full cell around

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
								Chain chainj = null;
								Chain chaini = structure.getChain(i);

								if (n==0 && a==0 && b==0 && c==0) {
									chainj = structure.getChain(j);
								} else {
									chainj = (Chain)structure.getChain(j).clone();
									Matrix4d m = new Matrix4d(ops[n]);
									translate(m, transOrth);
									Calc.transform(chainj,m);
								}

								StructureInterface interf = calcContacts(chaini, chainj, cutoff, tt, builder);

								if (interf!=null) {
									contactsFound++;
									set.add(interf);
								}
							}
						}

						if( verbose ) {
							if (a==0 && b==0 && c==0 && n==0)
								builder.append(" "+contactsFound+"("+(numChainsAu*(numChainsAu-1))/2+")");
							else if (selfEquivalent)
								builder.append(" "+contactsFound+"("+(numChainsAu*(numChainsAu+1))/2+")");
							else
								builder.append(" "+contactsFound+"("+numChainsAu*numChainsAu+")");

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

	private StructureInterface calcContacts(Chain chaini, Chain chainj, double cutoff, CrystalTransform tt, StringBuilder builder) {
		// note that we don't consider hydrogens when calculating contacts
		AtomContactSet graph = StructureTools.getAtomsInContact(chaini, chainj, cutoff, INCLUDE_HETATOMS);

		if (graph.size()>0) {
			if (builder != null) builder.append("x");

			CrystalTransform transf = new CrystalTransform(this.crystallographicInfo.getSpaceGroup());
			StructureInterface interf = new StructureInterface(
					StructureTools.getAllAtomArray(chaini), StructureTools.getAllAtomArray(chainj),
					chaini.getChainID(), chainj.getChainID(),
					graph,
					transf, tt);

			return interf;

		} else {
			if (builder != null) builder.append("o");
			return null;
		}
	}

	private void addVisited(CrystalTransform tt) {
		visited.add(tt);
	}

	/**
	 * Checks whether given transformId/translation is symmetry redundant
	 * Two transformations are symmetry redundant if their matrices (4d) multiplication gives the identity, i.e.
	 * if one is the inverse of the other.
	 * @param tt
	 * @return
	 */
	private boolean isRedundant(CrystalTransform tt) {

		Iterator<CrystalTransform> it = visited.iterator();
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

//	/**
//	 * If NCS operators are given in MTRIX records, a bigger AU has to be constructed based on those.
//	 * Later they have to be removed with {@link #removeExtraChains()}
//	 */
//	private void constructFullStructure() {
//
//		if (this.crystallographicInfo.getNcsOperators()==null ||
//			this.crystallographicInfo.getNcsOperators().length==0) {
//			// normal case: nothing to do
//			return;
//		}
//
//		// first we store the original chains in a new list to be able to restore the structure to its original state afterwards
//		origChains = new ArrayList<Chain>();
//		for (Chain chain:structure.getChains()) {
//			origChains.add(chain);
//		}
//
//		// if we are here, it means that the NCS operators exist and we have to complete the given AU by applying them
//		Matrix4d[] ncsOps = this.crystallographicInfo.getNcsOperators();
//
//		if (verbose)
//			System.out.println(ncsOps.length+" NCS operators found, generating new AU...");
//
//
//		for (int i=0;i<ncsOps.length;i++) {
//			Structure transformedStruct = (Structure)structure.clone();
//			Calc.transform(transformedStruct, ncsOps[i]);
//
//			for (Chain chain: transformedStruct.getChains()) {
//				// we assign a new AU id (chain ID) consisting in original chain ID + an operator index from 1 to n
//				chain.setChainID(chain.getChainID()+(i+1));
//				structure.addChain(chain);
//			}
//		}
//
//		// now we have more chains in AU, we have to update the value
//		this.numChainsAu = structure.getChains().size();
//	}
//
//	/**
//	 * Removes the extra chains that were added to the original structure in {@link #constructFullStructure()}
//	 */
//	private void removeExtraChains() {
//		structure.setChains(origChains);
//	}
}
