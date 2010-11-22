/*
 *                  BioJava development code
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
 * Created on Jul 16, 2006
 *
 */
package org.biojava.bio.structure.gui.util;


import java.util.logging.Logger;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.ClusterAltAligs;
import org.biojava.bio.structure.align.StructurePairAligner;
import org.biojava.bio.structure.align.gui.AlignmentGui;
import org.biojava.bio.structure.align.pairwise.AlternativeAlignment;



/** A class that obtains two structures via DAS and aligns them
 *  This is done in a separate thread.
 *  It is possible to register Event listeners to get notification of when the download has finished.
 *  
 * @author Andreas Prlic
 * @since 1.7
 * @version %I% %G%
 */
public class AlignmentCalc implements Runnable {


	public static Logger logger =  Logger.getLogger("org.biojava");

	boolean interrupted = false;

	String pdb1;
	String pdb2;
	String chain1;
	String chain2;

	Structure structure1;
	Structure structure2;

	AlignmentGui parent;


	/** requests an alignment of pdb1 vs pdb 2.
	 * Chain 1 and chain2 are optional.
	 * If they are empty strings, they are ignored
	 * @param parent the alignment gui frame that interacts with this class          
	 * @param s1 structure 1
	 * @param s2 structure 2
	 */
	public AlignmentCalc(AlignmentGui parent, Structure s1, Structure s2 ) {

		this.parent= parent;

		structure1 = s1;
		structure2 = s2;

	}

	public void run() {

		// both structure have been downloaded, now calculate the alignment ...


		StructurePairAligner aligner = new StructurePairAligner();
		aligner.setDebug(true);
		try {
			aligner.align(structure1,structure2);
		} catch (StructureException e){
			logger.warning(e.getMessage());

		}



		AlternativeAlignment[] aligs = aligner.getAlignments();
		//cluster similar results together 
		ClusterAltAligs.cluster(aligs);
		showAlignment(aligner,aligs);

		//logger.info("done!");

		parent.notifyCalcFinished();

	}



	private void showAlignment(StructurePairAligner alignment, AlternativeAlignment[] aligs) {
		AlternativeAlignmentFrame frame = new AlternativeAlignmentFrame(structure1, structure2);
		frame.setStructurePairAligner(alignment);
		frame.setAlternativeAlignments(aligs);
		frame.pack();
		frame.setVisible(true);



	}

	/** stops what is currently happening and does not continue
	 * 
	 *
	 */
	public void interrupt() {
		interrupted = true;
	}

	public void cleanup() {

		parent.notifyCalcFinished();

		parent=null;
		// cleanup...

		structure1 = null;
		structure2 = null;

	}



}




