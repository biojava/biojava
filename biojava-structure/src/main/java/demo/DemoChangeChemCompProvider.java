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
package demo;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.io.mmcif.AllChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.ChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;

import java.util.List;


/**
 *  This demo shows how to use an alternative ChemCompProvider. The default mechanism in BioJava is to access chemical components
 *  by using the {@link DownloadChemCompProvider}. It fetches and locally caches chemical component definitions as they are encountered during file parsing.
 *  It can be enabled by using the {@link FileParsingParameters#setLoadChemCompInfo(boolean)} method.
 *
 * The {@link AllChemCompProvider} downloads and unpacks all chemcomps. It is slower and requires more memory than the default {@link DownloadChemCompProvider},
 * but it avoids network access to the FTP site, if a new chemcomp is detected, that has not been downloaded yet.
 *
 * Since all chemcomps will be kept in memory, the standard memory that is available to a JVM will not be sufficient
 * in order to run this demo. Please start with -Xmx200M
 *
 * @author Andreas Prlic
 *
 */
public class DemoChangeChemCompProvider {

	public static void main(String[] args){
		String pdbId = "1O1G";

		boolean loadChemComp = true;

		//////
		// no need to change anything below here...
		//////

		PDBFileReader reader = new PDBFileReader();

		// Set the system wide property where PDB and chemcomp files are being cached.
		// you can set the path to the local PDB installation either like this
//		reader.setPath(PDB_PATH);
		// or via
		// by setting the PDB_PATH environmental variable or system property
		// when running the demo (e.g. -DPDB_DIR=/path/to/pdb)

		if ( loadChemComp) {

			// The AllChemCompProvider loads all chem comps at startup.
			// This is slow (13 sec on my laptop) and requires more
			// memory than the default DownloadChemCompProvider.
			// In contrast to it it keeps all definitions in memory.
			ChemCompProvider all = new AllChemCompProvider();

			ChemCompGroupFactory.setChemCompProvider(all);
		}

		DemoChangeChemCompProvider demo = new DemoChangeChemCompProvider();

		// run the demo
		demo.basicLoad(reader,loadChemComp, pdbId);

	}

	public void basicLoad(PDBFileReader reader, boolean loadChemComp, String pdbId){

		try {
			// configure the parameters of file parsing

			FileParsingParameters params = new FileParsingParameters();

			// should the ATOM and SEQRES residues be aligned when creating the internal data model?
			// only do this if you need to work with SEQRES sequences. If all you need are ATOMs, then
			// set it to false to have quicker file loading.
			params.setAlignSeqRes(true);

			//
			// should secondary structure get parsed from the file
			params.setParseSecStruc(false);

			reader.setFileParsingParameters(params);

			Structure struc = reader.getStructureById(pdbId);

			printStructure(struc);


		} catch (Exception e){
			e.printStackTrace();
		}

	}

	private void printStructure(Structure struc) {

		System.out.println(struc);

		//Chain c = struc.getChainByPDB("C");
		String pdbid = struc.getPDBCode();
		for (int i = 0; i < struc.nrModels(); i++) {

			// loop chain
			for (Chain ch : struc.getModel(i)) {
				if (! ch.getName().equals("A") )
					continue;
				System.out.println(pdbid + ">>>" + ch.getName() + ">>>"
						+ ch.getAtomSequence());
				System.out.println(pdbid + ">>>" + ch.getName() + ">>>"
						+ ch.getSeqResSequence());
				// Test the getAtomGroups() and getSeqResGroups() method

				List<Group> group = ch.getSeqResGroups();
				int seqPos = 0;
				for (Group gp : group) {
					System.out.println(ch.getName() + ":"+seqPos + ":" + gp.getResidueNumber() + ":"
							+ gp.getPDBName() + " " + gp.getType());
					seqPos++;
				}
			}
		}


	}
}
