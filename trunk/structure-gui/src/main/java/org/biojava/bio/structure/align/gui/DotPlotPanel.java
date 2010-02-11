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
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CECalculator;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.helper.JointFragments;
import org.biojava.bio.structure.align.model.AFP;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.bio.structure.align.pairwise.FragmentPair;
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
	DotPlotPanel(AFPChain alignment, Matrix background ){
		
		final double defaultBackground = 100.;
		
		List<AFP> afps = alignment.getAfpSet();
		
		// Generate the lists of equivalent amino acids based on alignment
		List<int[]> alignPairs = new ArrayList<int[]>(afps.size()*afps.get(0).getFragLen() );
		for(AFP afp : afps) {
			int start1 = afp.getP1();
			int start2 = afp.getP2();
			for(int i=0;i<afp.getFragLen();i++) {
				alignPairs.add( new int[] { start1+i, start2+i } );
			}
		}
		JointFragments frag = new JointFragments();
		frag.setIdxlist(alignPairs);
		AlternativeAlignment[] aligns = new AlternativeAlignment[1];
		aligns[0] = new AlternativeAlignment();
		aligns[0].apairs_from_idxlst(frag);

		
		// Calculate FragmentPairs based on alignment.
		// These are displayed as a small box around the start of each alignment.
		FragmentPair[] pairs = new FragmentPair[afps.size()];
		for(int i=0;i<pairs.length;i++) {
			AFP afp = afps.get(i);
			pairs[i] = new FragmentPair(afp.getFragLen(), afp.getP1(), afp.getP2());
			pairs[i].setRms(afp.getRmsd());
		}
		
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
		this.setAlternativeAligs(aligns);
		this.setSelectedAlignmentPos(0);
		this.setFragmentPairs(pairs);
	}

	public static void main(String[] args) {
		try {
			StructureAlignment ce = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);

			CeMain ceA = (CeMain) ce;
			
		
			
			String name1= "1cdg.A";
			String name2= "1tim.A";

			AtomCache cache = new AtomCache("/tmp/", true);

			Atom[] ca1 = cache.getAtoms(name1, true);
			Atom[] ca2 = cache.getAtoms(name2, true);
			
			AFPChain afpChain = ce.align(ca1,ca2);
			for ( AFP afpI : afpChain.getAfpSet()){
				System.out.println(afpI);
			}
			
			CECalculator calculator = ceA.getCECalculator();
			double[][] m = calculator.getMatMatrix();
			Matrix mat = new Matrix(m);

			DotPlotPanel dotplot = new DotPlotPanel(afpChain, mat);			
			
			//Create JFrame
			
			String title = String.format("Dot plot of %s vs. %s", name1,name2);

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
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
}

