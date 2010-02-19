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
 * Created on Aug 3, 2007
 * 
 */

package org.biojava.bio.structure.gui;

import java.awt.Dimension;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.StrucAligParameters;
import org.biojava.bio.structure.align.StructurePairAligner;
import org.biojava.bio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.bio.structure.align.pairwise.FragmentPair;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava.bio.structure.gui.JMatrixPanel;
import org.biojava.bio.structure.gui.ScaleableMatrixPanel;


/** A JPanel that can display the underlying distance matrix 
 * data of the protein structure alignment algorithm. It adds a 
 * JSlider to a JMatrixPanel.
 * 
 * see also JMatrixPanel.
 * 
 */
public class ScaleableMatrixPanel 
extends JPanel 
implements ChangeListener{

	/**
	 * 
	 */
	private static final long serialVersionUID = -8082261434322968652L;

	JMatrixPanel mPanel;
	JSlider slider;
	JScrollPane scroll;
	
	public static void main(String[] args){

		PDBFileReader pdbr = new PDBFileReader();  
		pdbr.setAutoFetch(true);
		pdbr.setPath("/Users/blivens/pdb/");


		//String pdb1 = "1crl";
		//String pdb2 = "1ede";

		String pdb1 = "1buz";
		String pdb2 = "1ali";            

		//String pdb1 = "5pti";
		//String pdb2 = "5pti";

		// NO NEED TO DO CHANGE ANYTHING BELOW HERE...

		StructurePairAligner sc = new StructurePairAligner();
		StrucAligParameters params = new StrucAligParameters();
		params.setMaxIter(1);
		sc.setParams(params);

		// step1 : read molecules
		try {
			Structure s1 = pdbr.getStructureById(pdb1);
			Structure s2 = pdbr.getStructureById(pdb2);      

			System.out.println("aligning " + pdb1 + " vs. " + pdb2);
			System.out.println(s1);
			System.out.println();
			System.out.println(s2);
			// step 2 : do the calculations
			sc.align(s1,s2);


			ScaleableMatrixPanel smp = new ScaleableMatrixPanel();
			JFrame frame = new JFrame();
			frame.addWindowListener(new WindowAdapter(){
				public void windowClosing(WindowEvent e){
					JFrame f = (JFrame) e.getSource();
					f.setVisible(false);
					f.dispose();
				}

						
				
			});
						
			smp.setMatrix(sc.getDistMat());
			smp.setFragmentPairs(sc.getFragmentPairs());
			smp.setAlternativeAligs(sc.getAlignments());
			
			for (int i = 0; i < sc.getAlignments().length; i++) {
				AlternativeAlignment aa =sc.getAlignments()[i];
				System.out.println(aa);
				
			}
			
			frame.getContentPane().add(smp);

			frame.pack();
			frame.setVisible(true);
			
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public ScaleableMatrixPanel(){

		mPanel   = new JMatrixPanel();
		Box vBox = Box.createVerticalBox();
		
		int RES_MIN  = 1;
		int RES_MAX  = 8;
		int RES_INIT = 1;

		slider = new JSlider(JSlider.HORIZONTAL, RES_MIN,RES_MAX,RES_INIT);
		slider.setInverted(false);
		slider.setPaintTicks(false);
		slider.setPaintLabels(false);
		slider.addChangeListener(this);
		slider.setPreferredSize(new Dimension(100,15));

		vBox.add(slider);

		scroll = new JScrollPane(mPanel);
		scroll.getHorizontalScrollBar().setUnitIncrement(60);
		scroll.getVerticalScrollBar().setUnitIncrement(60);
		scroll.getHorizontalScrollBar().setBlockIncrement(60);
		scroll.getVerticalScrollBar().setBlockIncrement(60);
		vBox.add(scroll);
		this.setPreferredSize(new Dimension(400,400));
		this.add(vBox);
		
		
		mPanel.setLayout(new BoxLayout(mPanel,BoxLayout.Y_AXIS));
		 this.setLayout(new BoxLayout(this,BoxLayout.Y_AXIS));


	}



	public void stateChanged(ChangeEvent e) {
		
		JSlider source = (JSlider)e.getSource();
		
		if ( source.getValueIsAdjusting())
			return;
		
		mPanel.setScale((float)source.getValue());
		
		scroll.repaint();
		scroll.updateUI();
	}

	public Matrix getMatrix() {
		return mPanel.getMatrix();
	}

	public void setMatrix(Matrix matrix) {
		mPanel.setMatrix(matrix);
	
		
	}

	public JMatrixPanel getMatrixPanel(){
		return mPanel;
	}
	
	public FragmentPair[] getFragmentPairs(){
		return mPanel.getFragmentPairs();
	}
	public void setFragmentPairs(FragmentPair[] pairs){
		mPanel.setFragmentPairs(pairs);
	}

	public AlternativeAlignment[] getAlternativeAligs() {
		return mPanel.getAlternativeAligs();
	}



	public void setAlternativeAligs(AlternativeAlignment[] aligs) {
		mPanel.setAlternativeAligs(aligs);
	}
	
	public int getSelectedAlignmentPos() {
		return mPanel.getSelectedAlignmentPos();
	}

	public void setSelectedAlignmentPos(int selectedAlignmentPos) {
		mPanel.setSelectedAlignmentPos(selectedAlignmentPos);
	}


}
