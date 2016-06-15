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
package org.biojava.nbio.structure.align;

/**
 * Simple Callable Class that calculates a pairwise alignment in a different
 * thread, so that multiple pairwise alignments can be run in parallel
 * (examples: all-to-all alignments, DB search alignments).
 * Adapted to a more general implementation since 4.1.0, because before it
 * was thought for DB search only.
 *
 * @author Aleix Lafita
 *
 */
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.client.PdbPair;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.SynchronizedOutFile;
import org.biojava.nbio.structure.align.xml.AFPChainXMLConverter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.concurrent.Callable;
import java.util.zip.GZIPOutputStream;

public class CallableStructureAlignment implements  Callable<AFPChain> {

	private final static Logger logger = LoggerFactory.getLogger(
			CallableStructureAlignment.class);

	//Structure information
	private PdbPair pair;
	private AtomCache cache;
	private Atom[] ca1;
	private Atom[] ca2;

	//File output information - for DB searches
	private SynchronizedOutFile outFile;
	private File outFileDir;

	//Algorithm information
	private String algorithmName;
	private ConfigStrucAligParams params;

	/**
	 * Default constructor. Used in DB search.
	 * Instantiates an empty object, everything has to be set independently.
	 */
	public CallableStructureAlignment() {}

	/**
	 * Constructor for all-to-all alignment calculation.
	 * Used for MultipleMC seed alignment calculation, for example.
	 *
	 * @param ca1 Atoms to align of the first structure
	 * @param ca2 Atoms to align of the second structure
	 * @param algorithmName the pairwise aligner algorithm to use, a new
	 * 			instance will be created for each thread.
	 * @param params parameter bean for the alignment.
	 */
	public CallableStructureAlignment(Atom[] ca1, Atom[] ca2, String algorithmName, ConfigStrucAligParams params){
		this.ca1 = ca1;
		this.ca2 = ca2;
		this.algorithmName = algorithmName;
		this.params = params;
	}

	@Override
	public AFPChain call() throws Exception {

		//Prepare the alignment algorithm
		StructureAlignment algorithm = StructureAlignmentFactory.getAlgorithm(algorithmName);
		if (params!=null) algorithm.setParameters(params);

		AFPChain afpChain = null;
		try {
			//Download the Atoms if they are not provided from the outisde (DB searches usually)
			if (ca1 == null) {
				Structure structure1 = cache.getStructure(pair.getName1());
				ca1 =  StructureTools.getRepresentativeAtomArray(structure1);
			} else ca1 = StructureTools.cloneAtomArray(ca1);

			Structure structure2 = null;
			if (ca2 == null) {
				structure2 = cache.getStructure(pair.getName2());
				ca2 = StructureTools.getRepresentativeAtomArray(structure2);
			} else ca2 = StructureTools.cloneAtomArray(ca2);

			afpChain = algorithm.align(ca1, ca2);

			if (pair!=null){
				afpChain.setName1(pair.getName1());
				afpChain.setName2(pair.getName2());
			}

			//Do not output anything if there is no File information
			if (outFile != null && outFileDir != null){
				String desc2 = structure2.getPDBHeader().getDescription();
				if ( desc2 == null)
					desc2="";
				afpChain.setDescription2(desc2);
				String result = afpChain.toDBSearchResult();
				logger.info("{}", result);

				outFile.write(result);

				String xml = AFPChainXMLConverter.toXML(afpChain, ca1, ca2);
				writeXML(outFileDir,pair.getName1(), pair.getName2(), xml);
			}

		} catch ( Exception e){
			logger.error("Exception: ", e);
		}
		return afpChain;
	}

	public PdbPair getPair() {
		return pair;
	}

	public void setPair(PdbPair pair) {
		this.pair = pair;
	}

	public AtomCache getCache() {
		return cache;
	}

	public void setCache(AtomCache cache) {
		this.cache = cache;
	}

	public SynchronizedOutFile getOutFile() {
		return outFile;
	}

	public void setOutFile(SynchronizedOutFile outFile) {
		this.outFile = outFile;
	}

	public Atom[] getCa1() {
		return ca1;
	}

	public void setCa1(Atom[] ca1) {
		this.ca1 = ca1;
	}

	private void writeXML(File outFileF, String name1, String name2, String xml)
	{
		try{
			// Create file
			File newF = new File(outFileF, "dbsearch_" +name1+"_" + name2+".xml.gz");

			FileOutputStream fstream = new FileOutputStream(newF);

			GZIPOutputStream gz = new GZIPOutputStream(fstream);
			OutputStreamWriter writer = new OutputStreamWriter(gz);
			writer.write(xml);
			writer.close();
		} catch (Exception e){//Catch exception if any
			logger.error("Exception: ", e);
		}
	}

	public void setOutputDir(File outFileF) {
		this.outFileDir = outFileF;
	}

	public void setAlgorithmName(String algorithmName) {
		this.algorithmName = algorithmName;
	}

	public void setParameters(ConfigStrucAligParams parameters) {
		this.params = parameters;
	}
}
