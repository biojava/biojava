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
 * Created on May 17, 2010
 * Author: Andreas Prlic
 *
 */

package demo;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.MMCIFFileReader;
import org.biojava.nbio.structure.io.StructureProvider;

import java.util.List;

/** An example of how to read MMcif files
 *
 * @author Andreas Prlic
 *
 */
public class DemoMMCIFReader
{

	public static void main(String[] args){

		DemoMMCIFReader demo = new DemoMMCIFReader();

		demo.loadSimple();

		demo.loadFromDirectAccess();

	}

	/** 
	 * A basic example how to load an mmCif file and get a Structure object
	 *
	 */
	public void loadSimple(){
		String pdbId = "4hhb";

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);

		StructureIO.setAtomCache(cache);

		try {
			Structure s = StructureIO.getStructure(pdbId);

			System.out.println(pdbId + " has nr atoms: " + StructureTools.getNrAtoms(s));

		} catch (Exception e){
			e.printStackTrace();
		}
	}


	/**
	 * An example demonstrating how to directly use the mmCif file parsing classes. This could potentially be used
	 * to use the parser to populate a data-structure that is different from the biojava-structure data model.
	 *
	 */
	public void loadFromDirectAccess(){
		String pdbId = "1A4W";

		StructureProvider pdbreader = new MMCIFFileReader();

		try {
			Structure s = pdbreader.getStructureById(pdbId);
			
			System.out.println("Getting chain H of 1A4W");

			List<Chain> hs = s.getNonPolyChainsByPDB("H");

			Chain h = hs.get(0);
			List<Group> ligands = h.getAtomGroups();

			System.out.println("These ligands have been found in chain " + h.getName());

			for (Group l:ligands){
				System.out.println(l);
			}

			System.out.println("Accessing QWE directly: ");
			Group qwe = s.getNonPolyChainsByPDB("H").get(2).getGroupByPDB(new ResidueNumber("H",373,null));

			System.out.println(qwe.getChemComp());

			System.out.println(h.getSeqResSequence());
			System.out.println(h.getAtomSequence());
			System.out.println(h.getAtomGroups(GroupType.HETATM));

			System.out.println("Entities: " + s.getEntityInfos());

		} catch (Exception e) {
			e.printStackTrace();
		}


	}
}
