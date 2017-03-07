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
 * Created on Jan 27, 2010
 * Author: Andreas Prlic
 *
 */

package demo;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.core.util.InputStreamProvider;



/** Example for how to load protein structures (from PDB files).
 *
 * @author Andreas Prlic
 *
 */
public class DemoLoadStructure
{

	public static void main(String[] args){

		DemoLoadStructure demo  = new DemoLoadStructure();

		demo.loadStructureIO();

		//demo.basicLoad();

		//demo.loadStructureFromCache();
	}

	public void loadStructureIO(){
		try {
			Structure s1 = StructureIO.getStructure("1gav");
			System.out.println(s1.getPDBCode() + " asym unit has nr atoms:");
			System.out.println(StructureTools.getNrAtoms(s1));


			Chain chain1 = s1.getChainByIndex(0);

			System.out.println("First chain: " + chain1);

			System.out.println("Chain " + chain1.getName() + " has the following sequence mismatches:");
			for (SeqMisMatch mm : chain1.getSeqMisMatches()){
				System.out.println(mm);
			}

			Structure s2 = StructureIO.getBiologicalAssembly("1gav");
			System.out.println(s2.getPDBCode() + " biological assembly has nr atoms:");
			System.out.println(StructureTools.getNrAtoms(s2));

		} catch (Exception e){
			e.printStackTrace();
		}

	}


	public void basicLoad(){
		try {

			PDBFileReader reader = new PDBFileReader();

			// the path to the local PDB installation
			reader.setPath("/tmp");

			// configure the parameters of file parsing

			FileParsingParameters params = new FileParsingParameters();

			// should the ATOM and SEQRES residues be aligned when creating the internal data model?
			params.setAlignSeqRes(true);

			// should secondary structure get parsed from the file
			params.setParseSecStruc(false);

			reader.setFileParsingParameters(params);

			Structure structure = reader.getStructureById("4hhb");

			System.out.println(structure);

			Chain c = structure.getPolyChainByPDB("C");


			System.out.print(c);

			System.out.println(c.getEntityInfo());


		} catch (Exception e){
			e.printStackTrace();
		}

	}

	public void loadStructureFromCache(){
		String pdbId = "4hhb";
		String chainName = "4hhb.A";
		String entityName = "4hhb:0";

		// we can set a flag if the file should be cached in memory
		// this will enhance IO massively if the same files have to be accessed over and over again.
		// since this is a soft cache, no danger of memory leak
		// this is actually not necessary to provide, since the default is "true" if the AtomCache is being used.
		System.setProperty(InputStreamProvider.CACHE_PROPERTY, "true");

		AtomCache cache = new AtomCache();

		try {
			System.out.println("======================");
			Structure s = cache.getStructure(pdbId);

			System.out.println("Full Structure:" + s);

			Atom[] ca = cache.getAtoms(chainName);
			System.out.println("got " + ca.length + " CA atoms");

			Structure firstEntity = cache.getStructure(entityName);
			System.out.println("First entity: " + firstEntity);

		} catch (Exception e){
			e.printStackTrace();
		}

	}


}
