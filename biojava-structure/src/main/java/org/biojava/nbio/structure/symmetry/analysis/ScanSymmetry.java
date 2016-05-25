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
package org.biojava.nbio.structure.symmetry.analysis;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.List;
import java.util.Set;

import org.biojava.nbio.structure.PDBCrystallographicInfo;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.AllChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.rcsb.GetRepresentatives;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.core.Subunits;
import org.biojava.nbio.structure.symmetry.misc.ProteinComplexSignature;
import org.biojava.nbio.structure.symmetry.utils.BlastClustReader;
import org.biojava.nbio.structure.xtal.SpaceGroup;


public class ScanSymmetry implements Runnable {
	private AtomCache cache = null;
	private static String RESULT_DIR = "/Users/peter/Results/ScanSymmetry/";


	public ScanSymmetry () {
		initializeCache();
	}

	public static void main(String[] args) {
		new ScanSymmetry().run();
	}

	@Override
	public void run() {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());

		System.out.println("Reading blastclust files");

		BlastClustReader reader95 = new BlastClustReader(95);
		BlastClustReader reader30 = new BlastClustReader(30);


		PrintWriter out = null;
		PrintWriter error = null;

		try {
			out = new PrintWriter(new FileWriter(RESULT_DIR + timeStamp + "_symm.csv"));
			error = new PrintWriter(new FileWriter(RESULT_DIR + timeStamp + "_error.txt"));
		} catch (IOException e1) {
			e1.printStackTrace();
			System.exit(-1);
		}


		long t1 = System.nanoTime();

		int success = 0;
		int proteins = 0;
		int failure = 0;

		String header = "pdbId,bioassembly,local,pseudostoichiometric,stoichiometry,pseudosymmetric,pointgroup,order," +
				"lowSymmetry,minidentity,maxidentity,subunitrmsd,rmsd,tm,minrmsd,maxrmsd,mintm,maxtm,rmsdintra,tmintra,symdeviation,subunits,nucleiacids,cacount,time,signature95,stoich95,signature30,stoich30,spacegroup";
		out.println(header);

		QuatSymmetryParameters parameters = new QuatSymmetryParameters();

		Set<String> set = GetRepresentatives.getAll();
		// pr testing
		set.clear();
		set.add("4HHB");

		// set skip to true to restart calculation with a specified PDB ID
		boolean skip = false;
		String restartId = "10MH";

		for (String pdbId: set) {
			if (skip && pdbId.equals(restartId)) {
				skip = false;
			}
			if (skip) {
				continue;
			}

			System.out.println("------------- " + pdbId  + "-------------");

			StructureIO.setAtomCache(cache);

			
			
			List<Structure> structures = null;
			try {
				structures = StructureIO.getBiologicalAssemblies(pdbId);
			} catch (StructureException|IOException e) {
				e.printStackTrace();
				error.println(pdbId + ": " + e.getMessage());
				error.flush();
				continue;
			}
			
			int i = 0;
			for (Structure structure : structures) {

				// note: before biojava 5.0 refactoring, the structures without bioassemblies would
				// use i=0 as the identifier for the default bioassembly (the asymmetric unit). Now
				// identifier i=1 is used - JD 2016-05-17
				i++;				

				long ts1 = System.nanoTime();

				try {
					SpaceGroup spaceGroup =null;
					//float resolution = 0.0f;
					if (structure != null) {
						PDBCrystallographicInfo info = structure.getCrystallographicInfo();
						if (info != null) {
							spaceGroup = info.getSpaceGroup();
						}
						//PDBHeader pdbHeader = structure.getPDBHeader();
						//resolution = pdbHeader.getResolution();
					}
					QuatSymmetryDetector detector = new QuatSymmetryDetector(structure, parameters);

					if (detector.hasProteinSubunits()) {
						long ts2 = System.nanoTime();

						int time = Math.round((ts2-ts1)/1000000.0f);

						// save global symmetry results
						List<QuatSymmetryResults> globalResults = detector.getGlobalSymmetry();
						printToCsv(reader95, reader30, out, pdbId,
								i, time, globalResults, spaceGroup);

						// save local symmetry results
						for (List<QuatSymmetryResults> localResults: detector.getLocalSymmetries()) {
							printToCsv(reader95, reader30, out, pdbId,
									i, time, localResults, spaceGroup);
						}
						proteins++;
					}
					success++;
					out.flush();
				} catch (Exception e) {
					failure++;
					e.printStackTrace();
					error.println(pdbId + "[" + i + "]: " + e.getMessage());
					error.flush();
				}
			}
		}
		long t2 = System.nanoTime();

		System.out.println("PDBs succeeded: " + success);
		System.out.println("PDBs failed   : " + failure);
		System.out.println("Proteins      : " + proteins);
		System.out.println("Total structure: " + set.size());
		System.out.println("Cpu time: " + (t2-t1)/1000000 + " ms.");

		out.close();
		error.close();
	}

	private void printToCsv(BlastClustReader reader95,
			BlastClustReader reader30, PrintWriter out, String pdbId,
			int bioAssemblyId, int time, List<QuatSymmetryResults> resultsList, SpaceGroup spaceGroup) {

		for (QuatSymmetryResults results: resultsList) {
			ProteinComplexSignature s95 = new ProteinComplexSignature(pdbId, results.getSubunits().getChainIds(), reader95);
			String signature95 = s95.getComplexSignature();
			String stoich95 = s95.getComplexStoichiometry();
			ProteinComplexSignature s30 = new ProteinComplexSignature(pdbId, results.getSubunits().getChainIds(), reader30);
			String signature30 = s30.getComplexSignature();
			String stoich30 = s30.getComplexStoichiometry();
			int order = 1;
			if (!results.getSymmetry().equals("H")) {
				order = results.getRotationGroup().getOrder();
			}

			out.println("PDB" + pdbId +"," + bioAssemblyId + "," + results.isLocal() +
					"," + results.getSubunits().isPseudoStoichiometric() +
					"," + results.getSubunits().getStoichiometry() +
					"," + results.getSubunits().isPseudoSymmetric() +
					"," + results.getSymmetry() +
					"," + order +
					"," + isLowSymmetry(results) +
					"," + Math.round(results.getSubunits().getMinSequenceIdentity()*100.0) +
					"," + Math.round(results.getSubunits().getMaxSequenceIdentity()*100.0) +
					"," + (float) results.getScores().getRmsdCenters() +
					"," + (float) results.getScores().getRmsd() +
					"," + (float) results.getScores().getTm() +
					"," + (float) results.getScores().getMinRmsd() +
					"," + (float) results.getScores().getMaxRmsd() +
					"," + (float) results.getScores().getMinTm() +
					"," + (float) results.getScores().getMaxTm() +
					"," + (float) results.getScores().getRmsdIntra() +
					"," + (float) results.getScores().getTmIntra() +
					"," + (float) results.getScores().getSymDeviation() +
					"," + results.getSubunits().getSubunitCount() +
					"," + results.getNucleicAcidChainCount() +
					"," + results.getSubunits().getCalphaCount() +
					"," + time +
					"," + signature95 +
					"," + stoich95 +
					"," + signature30 +
					"," + stoich30 +
					"," + spaceGroup);
		}
	}

	private boolean isLowSymmetry(QuatSymmetryResults results) {
		return getMinFold(results.getSubunits()) > 1 && results.getRotationGroup() != null && results.getRotationGroup().getPointGroup().equals("C1");
	}

	private int getMinFold(Subunits subunits) {
		if (subunits.getFolds().size() > 1) {
			return subunits.getFolds().get(1);
		}
		return subunits.getFolds().get(0);
	}

	private void initializeCache() {
		cache = new AtomCache();
		FileParsingParameters params = cache.getFileParsingParams();
		cache.setUseMmCif(true);
		params.setParseCAOnly(true);
//		MmCifBiolAssemblyProvider mmcifProvider = new MmCifBiolAssemblyProvider();
//		BioUnitDataProviderFactory.setBioUnitDataProvider(mmcifProvider.getClass().getCanonicalName());
		ChemCompGroupFactory.setChemCompProvider(new AllChemCompProvider());
//		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());
	}
}
