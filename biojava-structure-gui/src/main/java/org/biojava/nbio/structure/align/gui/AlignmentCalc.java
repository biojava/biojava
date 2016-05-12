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
package org.biojava.nbio.structure.align.gui;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/** A class that obtains two structures via DAS and aligns them
 *  This is done in a separate thread.
 *  It is possible to register Event listeners to get notification of when the download has finished.
 *
 * @author Andreas Prlic
 * @since 1.7
 * @version %I% %G%
 */
public class AlignmentCalc implements AlignmentCalculationRunnable {

	private static final Logger logger = LoggerFactory.getLogger(AlignmentCalc.class);

	boolean interrupted = false;

	String pdb1;
	String pdb2;

	String name1;
	String name2;

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
	public AlignmentCalc(AlignmentGui parent, Structure s1, Structure s2 , String name1, String name2) {

		this.parent= parent;

		structure1 = s1;
		structure2 = s2;

		this.name1 = name1;
		this.name2 = name2;

	}

	@Override
	public void run() {

		// both structure have been downloaded, now calculate the alignment ...

		StructureAlignment algorithm = parent.getStructureAlignment();
		//StructurePairAligner aligner = new StructurePairAligner();
		//aligner.setDebug(true);
		try {

			Atom[] ca1 = StructureTools.getRepresentativeAtomArray(structure1);
			Atom[] ca2 = StructureTools.getRepresentativeAtomArray(structure2);

			//System.out.println("ca1 size:" + ca1.length + " ca2 size: " + ca2.length);
			AFPChain afpChain = algorithm.align(ca1, ca2);

			afpChain.setName1(name1);
			afpChain.setName2(name2);

			StructureAlignmentJmol jmol =   StructureAlignmentDisplay.display(afpChain, ca1, ca2);

			String title = jmol.getTitle();
			ConfigStrucAligParams params = algorithm.getParameters();
			if ( params != null)
				title += " " + algorithm.getParameters().toString();
			jmol.setTitle(title);

			DisplayAFP.showAlignmentPanel(afpChain,ca1,ca2,jmol);

			System.out.println(afpChain.toCE(ca1,ca2));

		} catch (StructureException e){
			e.printStackTrace();
			logger.warn(e.getMessage());
		}



		//logger.info("done!");

		parent.notifyCalcFinished();

	}




	/** stops what is currently happening and does not continue
	 *
	 *
	 */
	@Override
	public void interrupt() {
		interrupted = true;
	}

	@Override
	public void cleanup() {

		parent.notifyCalcFinished();

		parent=null;
		// cleanup...

		structure1 = null;
		structure2 = null;

	}

	/** does not do anything here...
	 *
	 */
	@Override
	public void setNrCPUs(int useNrCPUs) {
		// TODO Auto-generated method stub
		//
	}



}




