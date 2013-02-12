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


import java.util.concurrent.atomic.AtomicBoolean;
import java.util.logging.Logger;


import org.biojava.bio.structure.Structure;


import org.biojava.bio.structure.align.MultiThreadedDBSearch;
import org.biojava.bio.structure.align.StructureAlignment;

import org.biojava.bio.structure.align.util.AtomCache;

import org.biojava.bio.structure.align.util.UserConfiguration;

public class AlignmentCalcDB implements AlignmentCalculationRunnable {
	public static Logger logger =  Logger.getLogger("org.biojava");

	AtomicBoolean interrupted ;


	String name1;

	Structure structure1;

	AlignmentGui parent;

	UserConfiguration config;


	String outFile;

	int nrCPUs;
	Boolean domainSplit ;

	StructureAlignment customAlgorithm;

	MultiThreadedDBSearch job = null;

	public StructureAlignment getAlgorithm() {
		return customAlgorithm;
	}

	public void setAlgorithm(StructureAlignment algo) {
		this.customAlgorithm = algo;
	}

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

		StructureAlignment algorithm = null;

		if ( parent != null )
			algorithm = parent.getStructureAlignment();
		else {
			algorithm = customAlgorithm;
		}


		if ( name1.startsWith("file:/"))
			name1= "CUSTOM";
		
		job = new MultiThreadedDBSearch(name1,structure1, outFile, algorithm, nrCPUs, domainSplit);

		AtomCache cache = new AtomCache(config);
		System.out.println("using cache: " + cache.getPath());
		System.out.println("name1: " + name1);
		System.out.println("structure:" + structure1.getName());
		job.setAtomCache(cache);

		if ( name1.equals("CUSTOM")) {
			job.setCustomFile1(parent.getDBSearch().getPDBUploadPanel().getFilePath1());
			job.setCustomChain1(parent.getDBSearch().getPDBUploadPanel().getChain1());
		}

		job.run();
		
		File resultList = job.getResultFile();
		//		if ((now-startTime)/1000 > 30) {


		//		try {
		//			out.flush();
		//			out.close();
		//		} catch (Exception e) {
		//			e.printStackTrace();
		//		}
		if ( parent != null ) {
			parent.notifyCalcFinished();
			if ( resultList != null) {
				DBResultTable table = new DBResultTable();
				table.show(resultList,config);
			}
		}

	}





	/** stops what is currently happening and does not continue
	 * 
	 *
	 */
	public void interrupt() {
		interrupted.set(true);
		if ( job != null)
			job.interrupt();



	}

	public void cleanup()
	{
		parent.notifyCalcFinished();

		parent=null;
		// cleanup...

		structure1 = null;
		config = null;

		if ( job != null)
			job.cleanup();

	}

	public void setNrCPUs(int useNrCPUs) {
		nrCPUs = useNrCPUs;

	}

	public synchronized boolean isInterrupted() {
		return interrupted.get();
	}





}
