/*
 *                    PDB web development code
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
 *
 * Created on Jul 28, 2009
 * Created by ap3
 *
 */

package org.biojava.bio.structure.align.gui;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.helper.JointFragments;
import org.biojava.bio.structure.align.model.AFP;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava.bio.structure.gui.ScaleableMatrixPanel;

/**
 * Displays the dot plot trace for an alignment.
 * 
 * This class adapts ScaleableMatrixPanel, which uses code from org.biojava.bio.structure.align.pairwise,
 * with more BioJava-friendly methods based off AFPChains.
 * 
 * @author Spencer Bliven
 *
 */
public class DotPlotPanel extends ScaleableMatrixPanel {

	private static final long serialVersionUID = -7641953255857483895L;

	/**
	 * 
	 * @param alignment The alignment to plot
	 * @param background [Optional]A matrix of 'background colors' over which to draw the alignment.
	 * 
	 *	Originally designed as a matrix of RMSD values between AFPs, so it is colorized 
	 *	accordingly from red (0) to black (>10). 
	 *
	 *  If this set to null, the background is set to black.
	 */
	public DotPlotPanel(AFPChain alignment ){
		super();
		
		final double defaultBackground = 100.;
				
		// Convert the AFPChain alignment into the MatrixPanel format
		AlternativeAlignment[] aligns = new AlternativeAlignment[alignment.getBlockNum()];
		int alignNumber = 0;
		
		//One alternative alignment for each block
		int[][][] optAln = alignment.getOptAln(); // [block #][{0,1} chain index][pos]

		for(;alignNumber < optAln.length;alignNumber++) {
			List<int[]> alignPairs = new ArrayList<int[]>();
			for(int pos = 0; pos<optAln[alignNumber][0].length; pos++ ) {
				alignPairs.add( new int[] {
						optAln[alignNumber][0][pos],
						optAln[alignNumber][1][pos] }
				);
			}
			JointFragments frag = new JointFragments();
			frag.setIdxlist(alignPairs);
			aligns[alignNumber] = new AlternativeAlignment();
			aligns[alignNumber].apairs_from_idxlst(frag);

		}
		
		/* TODO After the AFPSet is fixed in CeMain#filterDuplicateAFPs, maybe include this again
		//add alignment for the AFPs
		List<AFP> afps = alignment.getAfpSet();
		List<int[]> alignPairs = new ArrayList<int[]>();
		for(AFP afp : afps) {
			int start1 = afp.getP1();
			int start2 = afp.getP2();
			for(int i=0;i<afp.getFragLen();i++) {
				alignPairs.add( new int[] { start1+i, start2+i } );
			}
		}
		JointFragments frag = new JointFragments();
		frag.setIdxlist(alignPairs);
		aligns[alignNumber] = new AlternativeAlignment();
		aligns[alignNumber].apairs_from_idxlst(frag);
		*/

		
		/* AFP boxes are unnecessary.
		// Calculate FragmentPairs based on alignment.
		// These are displayed as a small box around the start of each alignment.
		FragmentPair[] pairs = new FragmentPair[afps.size()];
		for(int i=0;i<pairs.length;i++) {
			AFP afp = afps.get(i);
			pairs[i] = new FragmentPair(afp.getFragLen(), afp.getP1(), afp.getP2());
			pairs[i].setRms(afp.getRmsd());
		}
		
		this.setFragmentPairs(pairs);
		*/

		
		// Now the alignments have been build; add it
		this.setAlternativeAligs(aligns);
		this.setSelectedAlignmentPos(0); //color white, not red

		Matrix background = alignment.getDistanceMatrix();
		//Fill with default black background if none given
		if(background == null) {
			background = new Matrix(alignment.getCa1Length(),alignment.getCa2Length());
			for(int i=0;i<background.getRowDimension();i++)
				for(int j=0;j<background.getColumnDimension(); j++) {
					background.set(i, j, defaultBackground);
				}
		}
		
		// Set parameters
		this.setMatrix(background);
	}

	/**
	 * Helper function to create and display a JFrame with a single DotPlotPanel
	 * 
	 * @param afpChain
	 * @param background
	 */
	private static JFrame showDotPlotJFrame(AFPChain afpChain ) {
		
		DotPlotPanel dotplot = new DotPlotPanel(afpChain);			
				
		//Create JFrame
		
		String title = String.format("Dot plot of %s vs. %s", afpChain.getName1(),afpChain.getName2());

		// Create window
		JFrame frame = new JFrame(title);
		frame.addWindowListener(new WindowAdapter(){
			public void windowClosing(WindowEvent e){
				JFrame f = (JFrame) e.getSource();
				f.setVisible(false);
				f.dispose();
			}
		});

		
		frame.getContentPane().add(dotplot);

		frame.pack();
		frame.setVisible(true);
		
		return frame;
	}	
	
	public static void main(String[] args) {

//		String name2= "1k5j.A"; //16-68,73-119
//		String name1= "1lrh.A"; //80-127,37-79
		
		String name1= "1iu9.A";
		String name2= "1h0r.A";
		
		// Hard case
//		String name1= "1uiz.A";
//		String name2= "1xxa.C";

		AtomCache cache = new AtomCache();


		try {
			CeMain ceA = (CeMain) StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);

			CeParameters params = (CeParameters) ceA.getParameters();
			params.setMaxGapSize(0);
			
			Atom[] ca1 = cache.getAtoms(name1);
			Atom[] ca2 = cache.getAtoms(name2);
			
			// Create initial alignment
			AFPChain afpChain = ceA.align(ca1,ca2);
			afpChain.setName1(name1);
			afpChain.setName2(name2);
			for ( AFP afpI : afpChain.getAfpSet()){
				System.out.println(afpI);
			}
			
			/*
			// Get background distances
			CECalculator calculator = ceA.getCECalculator();
			int winSize = params.getWinSize();
			int winSizeComb1 = (winSize-1)*(winSize-2)/2;	
			double[][] m = calculator.initSumOfDistances(ca1.length, ca2.length, params.getWinSize(), winSizeComb1, ca1, ca2);
			//double[][] m = calculator.getMatMatrix();
			Matrix mat = new Matrix(m);
			
			//Find range
			double min = mat.get(0, 0);
			double max = min;
			for(int r=0;r<mat.getRowDimension();r++) {
				for(int c=0;c<mat.getColumnDimension();c++) {
					double y = mat.get(r,c);
					if(y<min)
						min = y;
					if(y>max)
						max = y;
				}
			}
			System.out.format("[%f, %f]\n", min, max);
			*/
			
			//afpChain.setDistanceMatrix(mat);
			showDotPlotJFrame(afpChain);
			
			//StructureAlignmentJmol jmol = new StructureAlignmentJmol(afpChain, ca1, ca2);
//			jmol.setStructure(cache.getStructure(name1));

			
			//////////////////////////
			// Now make it circular
			ceA = (CeMain) StructureAlignmentFactory.getAlgorithm(CeCPMain.algorithmName);
	
			System.out.format("Aligning %s[%d] with %s[%d] with CPs\n",name1,ca1.length,name2,ca2.length);
			afpChain = ceA.align(ca1,ca2);
			afpChain.setName1(name1);
			afpChain.setName2(name2+"-"+name2);
			for ( AFP afpI : afpChain.getAfpSet()){
				System.out.println(afpI);
			}
			
			/*/ Reuse mat from the non-cp case, for simplicity

			// Get background distances
			Atom[] ca2clone = new Atom[ca2.length*2];
			int pos = 0;
			for (Atom a : ca2){
				Group g = (Group)a.getParent().clone(); // works because each group has only a CA atom

				ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

				pos++;
			}
			for (Atom a : ca2){
				Group g = (Group)a.getParent().clone();

				ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

				pos++;
			}
			m = calculator.initSumOfDistances(ca1.length, ca2clone.length, params.getWinSize(), winSizeComb1, ca1, ca2clone);
			//m = calculator.getMatMatrix();
			mat = new Matrix(m);/*ca2.length,ca1.length);
			for(int i=0;i<ca2.length;i++)
				for(int j=0;j<ca1.length;j++) {
					mat.set(i, j, m[i][j]);
				}
			*/

			showDotPlotJFrame(afpChain);
			
			
		} catch (StructureException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}

