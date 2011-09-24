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
 * Created on Nov 5, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.gui;


import java.io.File;

import java.util.SortedSet;
import java.util.concurrent.ExecutorService;

import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.logging.Logger;


import org.biojava.bio.structure.Atom;

import org.biojava.bio.structure.Structure;

import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.CallableStructureAlignment;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeSideChainMain;
import org.biojava.bio.structure.align.client.FarmJobParameters;
import org.biojava.bio.structure.align.client.JFatCatClient;
import org.biojava.bio.structure.align.client.PdbPair;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.SynchronizedOutFile;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.domain.DomainProvider;
import org.biojava.bio.structure.domain.DomainProviderFactory;
import org.biojava.bio.structure.domain.RemoteDomainProvider;

import org.biojava3.core.util.ConcurrencyTools;

public class AlignmentCalcDB implements AlignmentCalculationRunnable {
	public static Logger logger =  Logger.getLogger("org.biojava");

	AtomicBoolean interrupted ;


	String name1;

	Structure structure1;

	AlignmentGui parent;

	UserConfiguration config;
	SortedSet<String> representatives;

	String outFile;
	File resultList;
	int nrCPUs;
Boolean domainSplit ;
	public AlignmentCalcDB(AlignmentGui parent, Structure s1,  String name1, UserConfiguration config,String outFile, Boolean domainSplit) {

		this.parent= parent;

		structure1 = s1;

		this.name1 = name1;

		this.config = config;
		//this.representatives = representatives;
		interrupted = new AtomicBoolean(false);
		this.outFile = outFile;
		this.domainSplit = domainSplit;
	}

	public void run() {
		
		AtomCache cache = new AtomCache(config);
		StructureAlignment algorithm = parent.getStructureAlignment();	
		String serverLocation = FarmJobParameters.DEFAULT_SERVER_URL;
		if ( representatives == null){
			SortedSet<String> repre = JFatCatClient.getRepresentatives(serverLocation,40);
			System.out.println("got  " + repre.size() + " representatives for comparison");
			representatives = repre;
		}

		String header = "# algorithm:" + algorithm.getAlgorithmName(); 

		String legend = "# name1\tname2\tscore\tprobability\trmsd\tlen1\tlen2\tcov1\tcov2\t%ID\t " ;
		if (    algorithm.getAlgorithmName().equalsIgnoreCase(CeMain.algorithmName) || 
				algorithm.getAlgorithmName().equalsIgnoreCase(CeSideChainMain.algorithmName)){
			legend =  "# name1\tname2\tscore\tz-score\trmsd\tlen1\tlen2\tcov1\tcov2\t%ID\t " ;
		}


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

				String config1 = "#param:file1=" + parent.getDBSearch().getPDBUploadPanel().getFilePath1();
				out.write(config1);
				out.write(AFPChain.newline);

				String config2 = "#param:chain1=" + parent.getDBSearch().getPDBUploadPanel().getChain1();
				out.write(config2);
				out.write(AFPChain.newline);

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
		for (String repre : representatives){

			if( domainSplit ) {
				SortedSet<String> domainNames = domainProvider.getDomainNames(repre);
				//System.out.println(repre +" got domains: " +domainNames);
				if( domainNames == null || domainNames.size()==0){
					// no domains found, use whole chain.
					submit(name1, repre, ca1, algorithm, outFileF, out, cache);
					continue;
				}
				//System.out.println("got " + domainNames.size() + " for " + repre);
				for( String domain : domainNames){
					submit(name1, domain, ca1, algorithm, outFileF, out, cache);
				}
			} else {
				submit(name1, repre, ca1, algorithm, outFileF, out, cache);
			}			

		}


		ThreadPoolExecutor  pool = ConcurrencyTools.getThreadPool();
		System.out.println(pool.getPoolSize());

		long startTime = System.currentTimeMillis();

		try {
			while ( pool.getCompletedTaskCount() < representatives.size()-1  ) {
				//long now = System.currentTimeMillis();
				//System.out.println( pool.getCompletedTaskCount() + " " + (now-startTime)/1000 + " " + pool.getPoolSize() + " " + pool.getActiveCount()  + " " + pool.getTaskCount()  );
//				if ((now-startTime)/1000 > 60) {
//					
//					interrupt();
//					System.out.println("completed: " + pool.getCompletedTaskCount());
//				}

				if ( interrupted.get())
					break;

				Thread.sleep(1000);

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
//		if ((now-startTime)/1000 > 30) {
		

		//		try {
		//			out.flush();
		//			out.close();
		//		} catch (Exception e) {
		//			e.printStackTrace();
		//		}
		parent.notifyCalcFinished();
		DBResultTable table = new DBResultTable();
		table.show(resultList,config);
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

	}

	public void cleanup()
	{
		parent.notifyCalcFinished();

		parent=null;
		// cleanup...

		structure1 = null;
		config = null;

	}

	public void setNrCPUs(int useNrCPUs) {
		nrCPUs = useNrCPUs;

	}

	public synchronized boolean isInterrupted() {
		return interrupted.get();
	}





}
