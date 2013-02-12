package org.biojava.bio.structure.align;

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
 * Created on Feb 11, 2013
 * Author: Andreas Prlic 
 *
 */

import java.io.File;
import java.util.SortedSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.atomic.AtomicBoolean;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.ce.CeSideChainMain;
import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.bio.structure.align.client.FarmJobParameters;
import org.biojava.bio.structure.align.client.JFatCatClient;
import org.biojava.bio.structure.align.client.PdbPair;
import org.biojava.bio.structure.align.client.StructureName;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.SynchronizedOutFile;
import org.biojava.bio.structure.domain.DomainProvider;
import org.biojava.bio.structure.domain.DomainProviderFactory;
import org.biojava.bio.structure.domain.RemoteDomainProvider;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava3.core.util.ConcurrencyTools;


/** Performs a multi threaded database search for an input protein structure
 * 
 * @author Andreas Prlic
 *
 */

public class MultiThreadedDBSearch {

	AtomicBoolean interrupted  ;

	StructureAlignment algorithm;
	
	String outFile;
	
	String name1;
	
	int nrCPUs;
	
	AtomCache cache;
	File resultList;
	SortedSet<String> representatives;
	
	boolean domainSplit;
	
	Structure structure1;
	
	String customFile1;
	String customChain1;
	
	public MultiThreadedDBSearch(String name, Structure structure,
			String outFile, 
			StructureAlignment algorithm,
			int nrCPUs, boolean domainSplit){
		
		interrupted = new AtomicBoolean(false);
		this.name1= name;
		this.structure1 = structure;
		this.outFile = outFile;
		this.algorithm = algorithm;
		this.nrCPUs = nrCPUs;
		this.domainSplit = domainSplit;
		cache  = new AtomCache();
		
		String serverLocation = FarmJobParameters.DEFAULT_SERVER_URL;
		if ( representatives == null){
			SortedSet<String> repre = JFatCatClient.getRepresentatives(serverLocation,40);
			System.out.println("got  " + repre.size() + " representatives for comparison");
			representatives = repre;
		}
	}
	
	
	public String getCustomFile1() {
		return customFile1;
	}


	/** set the file path for a custom, user provided file, not a standard PDB file.
	 * 
	 * @param customFile1
	 */
	public void setCustomFile1(String customFile1) {
		this.customFile1 = customFile1;
	}


	
	public String getCustomChain1() {
		return customChain1;
	}

	/** sets a chain in a custom, user provided file
	 * 
	 * @param customChain1
	 */
	public void setCustomChain1(String customChain1) {
		this.customChain1 = customChain1;
	}


	public AtomCache getAtomCache() {
		return cache;
	}

	public void setAtomCache(AtomCache cache) {
		this.cache = cache;
	}



	public StructureAlignment getAlgorithm() {
		return algorithm;
	}

	public void setAlgorithm(StructureAlignment algo) {
		this.algorithm = algo;
	}
	
	
	public String getOutFile() {
		return outFile;
	}


	public void setOutFile(String outFile) {
		this.outFile = outFile;
	}


	public static String getLegend(String algorithmName){

		if ( algorithmName.equalsIgnoreCase(CeMain.algorithmName) || 
				algorithmName.equalsIgnoreCase(CeSideChainMain.algorithmName) ||
				algorithmName.equalsIgnoreCase(CeCPMain.algorithmName)) {
			return "# name1\tname2\tscore\tz-score\trmsd\tlen1\tlen2\tcov1\tcov2\t%ID\tDescription\t " ;
		}

		// looks like a FATCAT alignment

		return "# name1\tname2\tscore\tprobability\trmsd\tlen1\tlen2\tcov1\tcov2\t%ID\tDescription\t " ;

	}
	
	
	
	public File getResultFile() {
		return resultList;
	}


	public void setResultFile(File resultList) {
		this.resultList = resultList;
	}


	public void run(){

		checkLocalFiles();
		
		if ( interrupted.get())
			return;
		
		String header = "# algorithm:" + algorithm.getAlgorithmName();
		String legend = getLegend(algorithm.getAlgorithmName());
		


		File outFileF = new File(outFile);
		if ( ! outFileF.isDirectory()){
			System.err.println( outFileF.getAbsolutePath() + " is not a directory, can't create result files in there... ");
			interrupt();
			cleanup();
		}

		if ( name1 == null) 
			name1 = "CUSTOM";

		SynchronizedOutFile out ;

		resultList = new File(outFileF,"results_" + name1 + ".out");

		try {

			out = new SynchronizedOutFile(resultList);

			out.write(header);
			out.write(AFPChain.newline);
			out.write(legend);
			out.write(AFPChain.newline);

			if ( name1.equals("CUSTOM")) {

				String config1 = "#param:file1=" + customFile1;
				out.write(config1);
				out.write(AFPChain.newline);

				if ( customChain1 != null) {
				String config2 = "#param:chain1=" + customChain1;
				out.write(config2);
				out.write(AFPChain.newline);
				}

			}

			if ( algorithm.getAlgorithmName().startsWith("jCE")){
				ConfigStrucAligParams params = algorithm.getParameters();
				if ( params instanceof CeParameters){
					CeParameters ceParams = (CeParameters) params;
					if ( ceParams.getScoringStrategy() != CeParameters.DEFAULT_SCORING_STRATEGY) {
						String scoring = "#param:scoring=" + ceParams.getScoringStrategy();
						out.write(scoring);
						out.write(AFPChain.newline);
					}
				}
			}

		} catch (Exception e){
			System.err.println("Error while loading representative structure " + name1);
			e.printStackTrace();
			interrupt();
			cleanup();
			return;
		}



		DomainProvider domainProvider = DomainProviderFactory.getDomainProvider();

		ConcurrencyTools.setThreadPoolSize(nrCPUs);

		Atom[] ca1 = StructureTools.getAtomCAArray(structure1);

		int nrJobs = 0;
		for (String repre : representatives){

			if( domainSplit ) {
				SortedSet<String> domainNames = domainProvider.getDomainNames(repre);
				//System.out.println(repre +" got domains: " +domainNames);
				if( domainNames == null || domainNames.size()==0){
					// no domains found, use whole chain.
					submit(name1, repre, ca1, algorithm, outFileF, out, cache);
					nrJobs++;
					continue;
				}
				//System.out.println("got " + domainNames.size() + " for " + repre);
				for( String domain : domainNames){
					submit(name1, domain, ca1, algorithm, outFileF, out, cache);
					nrJobs++;
				}
			} else {
				submit(name1, repre, ca1, algorithm, outFileF, out, cache);
				nrJobs++;
			}			

		}


		ThreadPoolExecutor  pool = ConcurrencyTools.getThreadPool();
		System.out.println(pool.getPoolSize());

		long startTime = System.currentTimeMillis();

		try {
			while ( pool.getCompletedTaskCount() < nrJobs-1  ) {
				//long now = System.currentTimeMillis();
				//System.out.println( pool.getCompletedTaskCount() + " " + (now-startTime)/1000 + " " + pool.getPoolSize() + " " + pool.getActiveCount()  + " " + pool.getTaskCount()  );
				//				if ((now-startTime)/1000 > 60) {
				//					
				//					interrupt();
				//					System.out.println("completed: " + pool.getCompletedTaskCount());
				//				}

				if ( interrupted.get())
					break;

				Thread.sleep(2000);

			}
			out.close();
		}
		catch (Exception e){
			e.printStackTrace();
			interrupt();
			cleanup();
		}

		if (domainProvider instanceof RemoteDomainProvider){
			RemoteDomainProvider remote = (RemoteDomainProvider) domainProvider;
			remote.flushCache();
		}
		long now = System.currentTimeMillis();
		System.out.println("Calculation took : " + (now-startTime)/1000 + " sec.");
		System.out.println( pool.getCompletedTaskCount() + " "  + pool.getPoolSize() + " " + pool.getActiveCount()  + " " + pool.getTaskCount()  );
	}
	
	

	private void checkLocalFiles() {
		System.out.println("checking local PDB installation in directory: " + cache.getPath());
		DomainProvider domainProvider = DomainProviderFactory.getDomainProvider();
		
		for (String repre : representatives){

			if ( interrupted.get())
				return;
			
			if( domainSplit ) {
				SortedSet<String> domainNames = domainProvider.getDomainNames(repre);
				//System.out.println(repre +" got domains: " +domainNames);
				if( domainNames == null || domainNames.size()==0){
					// no domains found, use whole chain.
					//submit(name1, repre, ca1, algorithm, outFileF, out, cache);
					checkFile(repre);
					continue;
				}
				//System.out.println("got " + domainNames.size() + " for " + repre);
				for( String domain : domainNames){
					//submit(name1, domain, ca1, algorithm, outFileF, out, cache);
					checkFile(domain);
				}
			} else {
				//submit(name1, repre, ca1, algorithm, outFileF, out, cache);
				checkFile(repre);
			}			

		}

		System.out.println("done checking local files...");
		
	}


	private void checkFile(String repre) {
		
		StructureName name = new StructureName(repre);
		
		PDBFileReader reader = new PDBFileReader();
		reader.setAutoFetch(true);
		reader.setPath(cache.getPath());
		reader.setFileParsingParameters(cache.getFileParsingParams());
		reader.setPdbDirectorySplit(cache.isSplit());
		if ( ! reader.checkFileExists(name.getPdbId()))
			reader.downloadPDB(name.getPdbId());
		
	}


	private void submit(String name12, String repre, Atom[] ca1, StructureAlignment algorithm , File outFileF , SynchronizedOutFile out , AtomCache cache ) {
		CallableStructureAlignment ali = new CallableStructureAlignment();

		PdbPair pair = new PdbPair(name1, repre);
		try {
			ali.setCa1(ca1);
		} catch (Exception e){
			e.printStackTrace();
			ConcurrencyTools.shutdown();
			return;
		}
		ali.setAlgorithmName(algorithm.getAlgorithmName());
		ali.setParameters(algorithm.getParameters());
		ali.setPair(pair);
		ali.setOutFile(out);			
		ali.setOutputDir(outFileF);
		ali.setCache(cache);

		ConcurrencyTools.submit(ali);

	}
	
	
	/** stops what is currently happening and does not continue
	 * 
	 *
	 */
	public void interrupt() {
		interrupted.set(true);
		ExecutorService pool = ConcurrencyTools.getThreadPool();
		pool.shutdownNow();
		DomainProvider domainProvider = DomainProviderFactory.getDomainProvider();
		if (domainProvider instanceof RemoteDomainProvider){
			RemoteDomainProvider remote = (RemoteDomainProvider) domainProvider;
			remote.flushCache();
		}

	}
	
	public void cleanup()
	{
	
		structure1 = null;
		

	}

}
