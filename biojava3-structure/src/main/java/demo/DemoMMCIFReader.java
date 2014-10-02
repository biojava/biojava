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

import java.util.List;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupType;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.MMCIFFileReader;
import org.biojava.bio.structure.io.StructureIOFile;
import org.biojava3.structure.StructureIO;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/** An example of how to read MMcif files
 * 
 * @author Andreas Prlic
 * 
 */
public class DemoMMCIFReader {

	private static final Logger logger = LoggerFactory.getLogger(DemoMMCIFReader.class);

	public static void main(String[] args){

		DemoMMCIFReader demo = new DemoMMCIFReader();

		demo.loadSimple();

		//demo.loadFromDirectAccess();

	}

	/** A basic example how to load an mmCif file and get a Structure object
	 *  
	 */
	public void loadSimple(){
		String pdbId = "4hhb";

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);

		StructureIO.setAtomCache(cache);

		try {
			Structure s = StructureIO.getStructure(pdbId);

			logger.info(pdbId + " has nr atoms: " + StructureTools.getNrAtoms(s));

		} catch (Exception e){
			logger.error("Exception: ", e);
		}
	}


	/** An example demonstrating how to directly use the mmCif file parsing classes. This could potentially be used
	 * to use the parser to populate a data-structure that is different from the biojava-structure data model.
	 * 
	 */
	public void loadFromDirectAccess(){
		String pdbId = "1A4W";

		StructureIOFile pdbreader = new MMCIFFileReader();

		try {
			pdbreader.setAutoFetch(true);
			Structure s = pdbreader.getStructureById(pdbId);

			Chain h = s.getChainByPDB("H");

			List<Group> ligands = h.getAtomLigands();

			logger.info("These ligands have been found in chain {}", h.getChainID());

			for (Group l:ligands){
				logger.info("{}", l);
			}

			logger.info("Accessing QWE directly: ");
			Group qwe = h.getGroupByPDB(new ResidueNumber("H",373,null));

			logger.info("{}", qwe.getChemComp());

			logger.info("{}", h.getSeqResSequence());
			logger.info("{}", h.getAtomSequence());
			logger.info("{}", h.getAtomGroups(GroupType.HETATM));

			logger.info("Compounds: {}", s.getCompounds());
		} catch (Exception e) {
			logger.error("Exception: ", e);
		}
	}
}
