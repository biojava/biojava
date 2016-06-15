package demo;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.StructureIO;

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
 * created at Sep 19, 2013
 * Author: Andreas Prlic
 */

public class DemoShowLargeAssembly {

	public static void main(String[] args){

		// This loads the PBCV-1 virus capsid, one of, if not the biggest biological assembly in terms on nr. of atoms.
		// The 1m4x.pdb1.gz file has 313 MB (compressed)
		// This Structure requires a minimum of 9 GB of memory to be able to be loaded in memory.

		String pdbId = "1M4X";

		Structure bigStructure = readStructure(pdbId,1);

		// let's take a look how much memory this consumes currently

		Runtime r = Runtime.getRuntime();

		// let's try to trigger the Java Garbage collector
		r.gc();

		System.out.println("Memory consumption after " + pdbId +
				" structure has been loaded into memory:");

		String mem = String.format("Total %dMB, Used %dMB, Free %dMB, Max %dMB",
				r.totalMemory() / 1048576,
				(r.totalMemory() - r.freeMemory()) / 1048576,
				r.freeMemory() / 1048576,
				r.maxMemory() / 1048576);

		System.out.println(mem);

		// 9693 atoms in the asymmetric unit * 1680 copies per assembly = 16284240 atoms
		System.out.println("# atoms: " + StructureTools.getNrAtoms(bigStructure));

	}
	/** Load a specific biological assembly for a PDB entry
	 *
	 * @param pdbId .. the PDB ID
	 * @param bioAssemblyId .. the first assembly has the bioAssemblyId 1
	 * @return a Structure object or null if something went wrong.
	 */
	public static Structure  readStructure(String pdbId, int bioAssemblyId) {

		// pre-computed files use lower case PDB IDs
		pdbId = pdbId.toLowerCase();

		// we just need this to track where to store PDB files
		// this checks the PDB_DIR property (and uses a tmp location if not set)
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		FileParsingParameters p = cache.getFileParsingParams();

		// some bio assemblies are large, we want an all atom representation and avoid
		// switching to a Calpha-only representation for large molecules
		// note, this requires several GB of memory for some of the largest assemblies, such a 1MX4
		p.setAtomCaThreshold(Integer.MAX_VALUE);

		// parse remark 350
		p.setParseBioAssembly(true);

		// download missing files

		Structure structure = null;
		try {
			structure = StructureIO.getBiologicalAssembly(pdbId,bioAssemblyId);
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}
		return structure;
	}
}
