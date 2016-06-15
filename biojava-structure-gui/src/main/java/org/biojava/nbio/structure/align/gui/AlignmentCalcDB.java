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

package org.biojava.nbio.structure.align.gui;


import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.MultiThreadedDBSearch;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.scop.ScopFactory;

import java.io.File;
import java.util.concurrent.atomic.AtomicBoolean;
//import org.slf4j.Logger;
//import org.slf4j.LoggerFactory;

public class AlignmentCalcDB implements AlignmentCalculationRunnable {


	public static String SCOP_VERSION =  "1.75";

	//private static final Logger logger = LoggerFactory.getLogger(AlignmentCalcDB.class);

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

		System.out.println("AlignmentCalcDB: Using SCOP version " + SCOP_VERSION);
		ScopFactory.setScopDatabase(SCOP_VERSION);

	}



	@Override
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
	@Override
	public void interrupt() {
		interrupted.set(true);
		if ( job != null)
			job.interrupt();



	}

	@Override
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

	@Override
	public void setNrCPUs(int useNrCPUs) {
		nrCPUs = useNrCPUs;

	}

	public synchronized boolean isInterrupted() {
		return interrupted.get();
	}





}
