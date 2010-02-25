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
package org.biojava.bio.structure.align.gui;


import java.util.List;
import java.util.logging.Logger;



import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.gui.AlignmentGui;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava.bio.structure.Calc;


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

	public void run() {

		// both structure have been downloaded, now calculate the alignment ...

		StructureAlignment algorithm = parent.getStructureAlignment();  
		//StructurePairAligner aligner = new StructurePairAligner();
		//aligner.setDebug(true);
		try {

			Atom[] ca1 = StructureTools.getAtomCAArray(structure1);
			Atom[] ca2 = StructureTools.getAtomCAArray(structure2);

			System.out.println("ca1 size:" + ca1.length + " ca2 size: " + ca2.length);
			AFPChain afpChain = algorithm.align(ca1, ca2);

			afpChain.setName1(name1);
			afpChain.setName2(name2);

			//System.out.println(afpChain);

			List<Group> hetatms = structure1.getChain(0).getAtomGroups("hetatm");
			List<Group> nucs    = structure1.getChain(0).getAtomGroups("nucleotide");

			List<Group> hetatms2 = structure2.getChain(0).getAtomGroups("hetatm");
			List<Group> nucs2    = structure2.getChain(0).getAtomGroups("nucleotide");
			
			StructureAlignmentJmol jmol;
			if(afpChain.getAlgorithmName().startsWith("jCE")) {
				//rotate all ca2 together
				Matrix m = afpChain.getBlockRotationMatrix()[0];
				Atom s = afpChain.getBlockShiftVector()[0];
				for(Atom a : ca2 ) {
					Calc.rotate(a.getParent(),m);
					Calc.shift(a.getParent(),s);
				}
				jmol = new StructureAlignmentJmol(afpChain, ca1, ca2);
			}
			else {
				jmol = StructureAlignmentDisplay.display(afpChain, ca1, ca2, hetatms, nucs,hetatms2, nucs2);
			}

			//String result = afpChain.toFatcat(ca1, ca2);

			//String rot = afpChain.toRotMat();

			DisplayAFP.showAlignmentImage(afpChain,ca1,ca2,jmol);           

			System.out.println(afpChain.toCE(ca1,ca2));
			/*
			JFrame frame = new JFrame();
			frame.setTitle(algorithm.getAlgorithmName()+ " : " + name1 + " vs. " + name2);
			frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

			String html = "<html><body><pre>"+result+"</pre></body></html>";

			//JTextPane tp = new JTextPane();
			JEditorPane tp = new JEditorPane("text/html", html);
			tp.setEditable(false);
			JScrollPane js = new JScrollPane();
			js.getViewport().add(tp);


			frame.getContentPane().add(js);
			frame.pack();      
			frame.setVisible(true);*/

			//aligner.align(structure1,structure2);
		} catch (StructureException e){
			logger.warning(e.getMessage());

		}



		//logger.info("done!");

		parent.notifyCalcFinished();

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




