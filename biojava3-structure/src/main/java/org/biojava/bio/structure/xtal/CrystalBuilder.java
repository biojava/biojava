package org.biojava.bio.structure.xtal;


import java.util.ArrayList;
import java.util.Iterator;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;

import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainInterface;
import org.biojava.bio.structure.ChainInterfaceList;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.contact.AtomContactSet;



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
	 * Whether to consider HETATOMs in contact calculations
	 */
	private static final boolean INCLUDE_HETATOMS = true;
	
	private Structure pdb;
	private SpaceGroup sg;
	private int numChainsAu;
	private int numOperatorsSg;
	
	private boolean debug;
	
	private int numCells;
	
	private ArrayList<CrystalTransform> visited;

	/**
	 * The bounding boxes of all AUs of the unit cell
	 */
	private UnitCellBoundingBox bbGrid;
	
	// debugging vars
	private long start; 
	private long end;
	private int trialCount;	
	private int skippedRedundant;
	private int skippedAUsNoOverlap;
	private int skippedChainsNoOverlap;
	private int skippedSelfEquivalent;
	

	
	public CrystalBuilder(Structure pdb) {
		this.pdb = pdb;
		this.numChainsAu = pdb.getChains().size();
		this.sg = (pdb.getCrystallographicInfo()==null)?null:pdb.getCrystallographicInfo().getSpaceGroup();
		this.numOperatorsSg = 1;
		if (sg!=null) {
			this.numOperatorsSg = sg.getMultiplicity();
		}
		this.debug = false;
		this.numCells = DEF_NUM_CELLS;		
	}
	
	/**
	 * Set the debug flag for verbose output to stdout
	 * @param debug
	 */
	public void setDebug(boolean debug) {
		//TODO use logger for debug output once we have slf4j in BJ
		this.debug = debug;
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
	 * Returns a list of all unique interfaces that the given Structure has upon 
	 * generation of all crystal symmetry mates. An interface is defined as any pair of chains 
	 * that contact, i.e. for which there is at least a pair of atoms (one from each chain) within 
	 * the given cutoff distance.
	 * @param cutoff the distance cutoff for 2 chains to be considered in contact
	 * @return
	 */
	public ChainInterfaceList getUniqueInterfaces(double cutoff) {	

		
		ChainInterfaceList set = new ChainInterfaceList();
		
		// initialising the visited ArrayList for keeping track of symmetry redundancy
		initialiseVisited();
		
		
		
		// initialising debugging vars
		start = -1; 
		end = -1;
		trialCount = 0;
		skippedRedundant = 0;
		skippedAUsNoOverlap = 0;
		skippedChainsNoOverlap = 0;
		skippedSelfEquivalent = 0;
		
		bbGrid = new UnitCellBoundingBox(numOperatorsSg, numChainsAu);
		
		bbGrid.setOriginalAuBbs(pdb, INCLUDE_HETATOMS);		
		
		// we can always calculate contacts within AU (be it crystallographic or not)
		calcInterfacesWithinAu(set, cutoff);
		
		// this condition covers 3 cases:
		// a) entries with expMethod X-RAY/other diffraction and defined crystalCell (most usual case)
		// b) entries with expMethod null but defined crystalCell (e.g. PDB file with CRYST1 record but no expMethod annotation) 
		// c) entries with expMethod not X-RAY (e.g. NMR) and defined crystalCell (NMR entries do have a dummy CRYST1 record "1 1 1 90 90 90 P1")
		if (    pdb.getCrystallographicInfo()!=null && 
				pdb.getCrystallographicInfo().getCrystalCell()!=null && 
				pdb.isCrystallographic()) {
						

			// we can only do this for crystallographic structures
			calcInterfacesCrystal(set, cutoff);
			
			
			if (debug) {
				end = System.currentTimeMillis();
				System.out.println("\n"+trialCount+" chain-chain clash trials done. Time "+(end-start)/1000+"s");
				System.out.println("  skipped (not overlapping AUs)       : "+skippedAUsNoOverlap);
				System.out.println("  skipped (not overlapping chains)    : "+skippedChainsNoOverlap);
				System.out.println("  skipped (sym redundant op pairs)    : "+skippedRedundant);
				System.out.println("  skipped (sym redundant self op)     : "+skippedSelfEquivalent);

				System.out.println("Found "+set.size()+" interfaces.");
			}
		}
		
		return set;
	}
	
	/**
	 * Calculate interfaces within asymmetric unit
	 * @param set
	 * @param cutoff
	 */
	private void calcInterfacesWithinAu(ChainInterfaceList set, double cutoff) {
		
		
		if (debug) {
			trialCount = 0;
			start= System.currentTimeMillis();			
			System.out.println("\nInterfaces within asymmetric unit (total possible trials "+(numChainsAu*(numChainsAu-1))/2+")");
			System.out.print("[ 0-( 0, 0, 0)] "); // printing header for dots line to have same format as in calcInterfacesCrystal
		}
		
		int contactsFound = 0;
		
		int i = -1;
		for (Chain chaini:pdb.getChains()) {
			i++;
			int j = -1;
			for (Chain chainj:pdb.getChains()) {
				j++;
				if (j<=i) continue;
				
				// before calculating the AtomContactSet we check for overlap, then we save putting atoms into the grid
				if (!bbGrid.getChainBoundingBox(0,i).overlaps(bbGrid.getChainBoundingBox(0,j), cutoff)) { 
					if (debug) {
						skippedChainsNoOverlap++;
						System.out.print(".");
					}
					continue;
				}
				
				if (debug) trialCount++;
				
				AtomContactSet graph = StructureTools.getAtomsInContact(chaini, chainj, cutoff, INCLUDE_HETATOMS); 
				if (graph.size()>0) {
					contactsFound++;
					if (debug) System.out.print("x");					
					
					CrystalTransform transf = new CrystalTransform(sg);
					ChainInterface interf = new ChainInterface(chaini,chainj,graph,transf,transf);
					
					set.add(interf);
					
				} else {
					if (debug) System.out.print("o");
				}
			}
		}
		if (debug) {
			end = System.currentTimeMillis();
			System.out.println(" "+contactsFound+"("+(numChainsAu*(numChainsAu-1))/2+")");
			System.out.println("\n"+trialCount+" chain-chain clash trials done. Time "+(end-start)/1000+"s");
		}

		
		
	}
	
	/**
	 * Calculate interfaces between original asymmetric unit and neighboring 
	 * whole unit cells, including the original full unit cell i.e. i=0,j=0,k=0 
	 * @param set
	 * @param cutoff
	 */
	private void calcInterfacesCrystal(ChainInterfaceList set, double cutoff) {

		// both arrays below are of size numOperatorsSg (multiplicity of space group)
		// generate complete unit cell, by applying all SG operators
		Structure[] cell = getUnitCell();
		// we calculate all the bounds of each of the asym units, those will then be reused and translated
		bbGrid.setAllNonAuBbs(cell, INCLUDE_HETATOMS);
		
		if (debug) {
			trialCount = 0;
			start= System.currentTimeMillis();
			int neighbors = (2*numCells+1)*(2*numCells+1)*(2*numCells+1)-1;
			int trials = numChainsAu*numOperatorsSg*numChainsAu*neighbors;
			System.out.println("\nInterfaces between the original asym unit and the neighbouring "+neighbors+" whole unit cells ("+numCells+" neighbours)" +
					"(2x"+numChainsAu+"chains x "+numOperatorsSg+"AUs x "+neighbors+"cells = "+trials+" total possible trials)");
		}


		for (int i=-numCells;i<=numCells;i++) {
			for (int j=-numCells;j<=numCells;j++) {
				for (int k=-numCells;k<=numCells;k++) {
					
					Point3i trans = new Point3i(i,j,k);
					Vector3d transOrth = new Vector3d(i,j,k);
					pdb.getCrystallographicInfo().getCrystalCell().transfToOrthonormal(transOrth);
					UnitCellBoundingBox bbGridTrans = bbGrid.getTranslatedBbs(transOrth);

					for (int au=0;au<numOperatorsSg;au++) { 
						if (au==0 && i==0 && j==0 && k==0) continue; // that would be the original au 

						// short-cut strategies
						// 1) we skip first of all if the bounding boxes of the AUs don't overlap
						if (!bbGrid.getAuBoundingBox(0).overlaps(bbGridTrans.getAuBoundingBox(au), cutoff)) {
							if (debug) skippedAUsNoOverlap++;
							continue;
						}

						// 2) we check if we didn't already see its equivalent symmetry operator partner 													
						CrystalTransform tt = new CrystalTransform(sg,au);
						tt.translate(trans);
						if (isRedundant(tt)) { 								
							if (debug) skippedRedundant++;								
							continue;
						}
						addVisited(tt);
						
						
						boolean selfEquivalent = false;
						
						// now we copy and actually translate the AU if we saw it does overlap and the sym op was not redundant
						Structure jAsym = cell[au].clone();
						Calc.translate(jAsym, transOrth);
						

						// 3) an operator can be "self redundant" if it is the inverse of itself (involutory, e.g. all pure 2-folds with no translation)						
						if (tt.isEquivalent(tt)) { 
							if (debug) 
								System.out.println("Transform "+tt+" is equivalent to itself, will skip half of i-chains to j-chains comparisons");
							// in this case we can't skip the operator, but we can skip half of the matrix comparisons e.g. j>i
							// we set a flag and do that within the loop below
							selfEquivalent = true;
						}
						
						if (debug) System.out.print(tt+" ");
						
						// Now that we know that boxes overlap and operator is not redundant, we have to go to the details 
						int contactsFound = 0;
												
						int jIdx = -1;
						for (Chain chainj:jAsym.getChains()) {
							jIdx++;
							int iIdx = -1;
							for (Chain chaini:pdb.getChains()) { // we only have to compare the original asymmetric unit to every full cell around
								iIdx++;
								if(selfEquivalent && (jIdx>iIdx)) {
									// in case of self equivalency of the operator we can safely skip half of the matrix
									skippedSelfEquivalent++;
									continue;
								}
								// before calculating the AtomContactSet we check for overlap, then we save putting atoms into the grid
								if (!bbGrid.getChainBoundingBox(0,iIdx).overlaps(bbGridTrans.getChainBoundingBox(au,jIdx),cutoff)) {
									if (debug) {
										skippedChainsNoOverlap++;
										System.out.print(".");
									}
									continue;
								}
								if (debug) trialCount++;

								AtomContactSet graph = StructureTools.getAtomsInContact(chaini, chainj, cutoff, INCLUDE_HETATOMS);
								if (graph.size()>0) {
									contactsFound++;										
									if (debug) System.out.print("x");
									
									CrystalTransform transf = new CrystalTransform(sg);
									ChainInterface interf = new ChainInterface(chaini,chainj,graph,transf,tt);

									set.add(interf);
									
								} else {
									if (debug) System.out.print("o");
								}
							}
						}
						if (debug) {
							if (selfEquivalent) 								
								System.out.println(" "+contactsFound+"("+(numChainsAu*(numChainsAu+1))/2+")");							
							else
								System.out.println(" "+contactsFound+"("+numChainsAu*numChainsAu+")");
						}
					}
				}
			}
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

				if (debug) System.out.println("Skipping redundant transformation: "+tt+", equivalent to "+v);
				
				// there's only 1 possible equivalent partner for each visited matrix 
				// (since the equivalent is its inverse matrix and the inverse matrix is unique)
				// thus once the partner has been seen, we don't need to check it ever again
				it.remove();
				
				return true;
			}
		}
		
		return false;
	}
	
	/**
	 * Generates all symmetry-related objects from this asym unit and returns the whole
	 * unit cell (this asymmetric unit plus the symmetry-related objects). 
	 * @return
	 */
	private Structure[] getUnitCell() {

		Structure[] aus = new Structure[numOperatorsSg];
		aus[0] = pdb;

		int i = 1;
		for (Matrix4d m:pdb.getCrystallographicInfo().getTransformationsOrthonormal()) {
			
			Structure sym = pdb.clone();
			
			Calc.transform(sym, m); 

			aus[i] = sym;
			
			i++;
			
		}
		
		return aus;
	}
	
	
}
