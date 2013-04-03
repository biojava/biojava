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
 * created at May 26, 2008
 */
package org.biojava.bio.structure.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.StructurePairAligner;
import org.biojava.bio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.bio.structure.gui.events.AlignmentPositionListener;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.gui.SequenceDisplay;
import org.biojava.bio.structure.gui.util.AlignedPosition;
import org.biojava.bio.structure.gui.util.SequenceMouseListener;
import org.biojava.bio.structure.gui.util.SequenceScalePanel;

/** A sequence display that can show the results of a protein structure alignment.
 * 
 * @author Andreas Prlic
 * @since 1.7
 */ 
public class SequenceDisplay extends JPanel implements ChangeListener {

	/**
	 * 
	 */
	private static final long serialVersionUID = -1829252532712454236L;

	Structure structure1;
	Structure structure2;
	AlternativeAlignment alig;
	StructurePairAligner structurePairAligner;
	SequenceScalePanel panel1;
	SequenceScalePanel panel2;

	JSlider  residueSizeSlider;
	JLabel  percentageDisplay;	

	int[] idx1;
	int[] idx2;

	/** the maximum value that the scale can get
	 * 
	 */
	public static final int MAX_SCALE               = 10;

	Logger logger = Logger.getLogger("org.biojava");

	List<AlignedPosition> apos;

	float scale;
	SequenceMouseListener mouseListener1;
	SequenceMouseListener mouseListener2;


	JLabel label1;
	JLabel label2;

	public static void main(String[] args){

		try {
			PDBFileReader pdbr = new PDBFileReader();          
			pdbr.setPath("/Users/andreas/WORK/PDB/");


			//String pdb1 = "1crl";
			//String pdb2 = "1ede";

			String pdb1 = "1buz";
			String pdb2 = "5pti";            


			// NO NEED TO DO CHANGE ANYTHING BELOW HERE...

			StructurePairAligner sc = new StructurePairAligner();


			// step1 : read molecules

			System.out.println("aligning " + pdb1 + " vs. " + pdb2);

			Structure s1 = pdbr.getStructureById(pdb1);
			Structure s2 = pdbr.getStructureById(pdb2);                       

			// step 2 : do the calculations
			sc.align(s1,s2);


			AlternativeAlignment[] aligs = sc.getAlignments();
			SequenceDisplay displ = new SequenceDisplay(sc);
			displ.setStructure1(s1);
			displ.setStructure2(s2);
			displ.setAlternativeAlignment(aligs[0]);

			displ.updateDisplay();

			JFrame frame = new JFrame("Sequences for AlternativeAlignment ["+0+"]");

			frame.getContentPane().add(displ);

			frame.pack();
			frame.setVisible(true);
			frame.addWindowListener(new WindowAdapter(){
				public void windowClosing(WindowEvent e){
					JFrame f = (JFrame) e.getSource();
					f.setVisible(false);
					f.dispose();
				}



			});

		} catch (Exception e){
			e.printStackTrace();
		}

	}

	public SequenceDisplay(StructurePairAligner structurePairAligner){
		super();



		structure1 = null;
		structure2 = null;
		alig = null;
		this.structurePairAligner = structurePairAligner;
		panel1 = new SequenceScalePanel(1);
		panel2 = new SequenceScalePanel(2);

		mouseListener1 = new SequenceMouseListener(this);


		panel1.addMouseListener(mouseListener1);
		panel1.addMouseMotionListener(mouseListener1);


		mouseListener2 = new SequenceMouseListener(this);
		panel2.addMouseListener(mouseListener2);
		panel2.addMouseMotionListener(mouseListener2);

		//SequenceMouseListener ml = new SequenceMouseListener(this);
		//this.addMouseListener(ml);
		//this.addMouseMotionListener(ml);

		Box vBox = Box.createVerticalBox();


		Box hBox1 = Box.createHorizontalBox();
		Box hBox2 = Box.createHorizontalBox();

		label1 = new JLabel();
		hBox1.add(label1);

		label2 = new JLabel();
		hBox2.add(label2);


		hBox1.add(panel1);
		hBox2.add(panel2);

		vBox.add(hBox1);
		vBox.add(hBox2);


		int RES_MIN  = 1;
		int RES_MAX  = 100;
		int RES_INIT = 100;
		residueSizeSlider = new JSlider(JSlider.HORIZONTAL,
				RES_MIN, RES_MAX, RES_INIT);
		residueSizeSlider.setInverted(true);
		//residueSizeSlider.setMajorTickSpacing(5);
		//residueSizeSlider.setMinorTickSpacing(2);
		residueSizeSlider.setPaintTicks(false);
		residueSizeSlider.setPaintLabels(false);
		residueSizeSlider.addChangeListener(this);
		//residueSizeSlider.setPreferredSize(new Dimension(100,15));

		percentageDisplay =  new JLabel("100 %");

		Box hBox = Box.createHorizontalBox();
		hBox.setBackground(Color.white);
		hBox.add(Box.createHorizontalGlue());
		hBox.add(residueSizeSlider);
		hBox.add(percentageDisplay);
		hBox.add(Box.createHorizontalGlue());

		//vBox.add(hBox);

		JScrollPane scroll = new JScrollPane(vBox);

		//scroll.setPreferredSize(new Dimension(500,160));

		Box vBox2 = Box.createVerticalBox();
		vBox2.add(scroll);
		vBox2.add(hBox);

		//vBox2.setPreferredSize(new Dimension(500,160));
		//vBox2.setSize(new Dimension(500,160));
		//vBox2.setMinimumSize(new Dimension(500,160));
		//vBox2.setMaximumSize(new Dimension(500,160));
		this.setPreferredSize(new Dimension(500,100));

		this.add(vBox2);

		this.setLayout(new BoxLayout(this,BoxLayout.Y_AXIS));

		apos = new ArrayList<AlignedPosition>();
	}

	public void clearListeners(){

		mouseListener1.clearListeners();
		mouseListener2.clearListeners();

	}
	public void addAlignmentPositionListener(AlignmentPositionListener li){
		mouseListener1.addAlignmentPositionListener(li);
		mouseListener2.addAlignmentPositionListener(li);
	}

	public StructurePairAligner getStructurePairAligner() {
		return structurePairAligner;
	}

	public void setStructurePairAligner(StructurePairAligner structurePairAligner) {
		this.structurePairAligner = structurePairAligner;
	}
	/** get the identical position in the alignment
	 * 
	 * @return identical positions for structure1
	 */
	public int[] getIdx1() {
		return idx1;
	}

	/** set the identical positions in the alignment
	 * 
	 * @param idx identical positions for structure1
	 */
	private void setIdx1(int[] idx) {
		this.idx1 = idx;

	}
	/** get the identical position in the alignment
	 * 
	 * @return identical positions for structure2
	 */
	public int[] getIdx2() {
		return idx2;
	}

	/** set the identical positions in the alignment
	 * 
	 * @param idx identical positions for structure2
	 */
	private void setIdx2(int[] idx) {
		this.idx2 = idx;

	}



	private void buildAligMap(){
		apos.clear();

		int gap   = 0;
		int gpos1 = 0;
		int gpos2 = 0;

		for (int pos = 0 ; pos < idx1.length ; pos ++){

			int p1 = idx1[pos];
			int p2 = idx2[pos];

			int end = Math.max(p1,p2);

			//System.out.println("p1: " + p1 + " p2: " + p2 );
			// fill up  gaps...
			for (;gap<end;gap++){

				AlignedPosition m = new AlignedPosition();
				if ( gpos1 < p1){
					m.setPos1(gpos1);

					gpos1++;
				}
				if ( gpos2 < p2){
					m.setPos2(gpos2);

					gpos2++;
				}
				m.setEquivalent(AlignedPosition.NOT_ALIGNED);

				//System.out.println(m + " => " + end);
				apos.add(m);				
			}

			// add this aligned position
			AlignedPosition m = new AlignedPosition();

			m.setPos1(p1);
			m.setPos2(p2);
			m.setEquivalent(AlignedPosition.EQUIVALENT);

			//System.out.println(m);
			apos.add(m);
			gpos1++;
			gpos2++;
			gap++;

		}

		//System.out.println(apos);

	}


	private void setAtoms(Structure s, SequenceScalePanel panel){
		if ( structurePairAligner == null){
			System.err.println("StructurePairAligner has not been set");
			return;
		}
		Atom[] ca1 = structurePairAligner.getAlignmentAtoms(s);
		Chain c = new ChainImpl();
		c.setChainID("1");
		for (Atom atom : ca1) {
			
			Group g = atom.getGroup();
			
			Chain parentChain = g.getChain();
			
			c.addGroup(g);
			// hack for Jmol?			
			g.setChain(parentChain);
		}
		panel.setChain(c);

	}

	//TODO: add a method to allow the display if the structure alignment
	// has been called with setting the atoms directly



	public void setStructure1(Structure structure){
		this.structure1 = structure;
		if ( structure != null) {
			setAtoms(structure1,panel1);
			label1.setText(structure.getPDBCode());
			label1.repaint();
		}


	}
	public void setStructure2(Structure structure){
		this.structure2 = structure;
		if ( structure != null){
			setAtoms(structure2,panel2);
			label2.setText(structure.getPDBCode());
			label2.repaint();
		}
	}

	public void setAlternativeAlignment(AlternativeAlignment alig){
		this.alig  = alig;
		this.setIdx1(alig.getIdx1());
		this.setIdx2(alig.getIdx2());

		buildAligMap();

		panel1.setAligMap(apos);
		panel2.setAligMap(apos);

		updateDisplay();


	}
	public List<AlignedPosition> getAligMap(){
		return apos;
	}

	public void stateChanged(ChangeEvent e) {

		JSlider source = (JSlider)e.getSource();
		//if (!source.getValueIsAdjusting()) {

		int residueSize = (int)source.getValue();
		calcScale(residueSize);

		updatePercentageDisplay();



		this.repaint();
		this.revalidate();

		//this.updateUI();
		//int width = getTotalWidth();
		//int height = getTotalHeight();
		//Dimension d = new Dimension(width,height);
		//logger.info("setting preferred size" + width + " " + height);
		//this.setPreferredSize(d);
		//this.setSize(d);
		// }


	}
	public void updateDisplay(){

		int residueSize = (int)residueSizeSlider.getValue();
		calcScale(residueSize);
		updatePercentageDisplay();
		this.repaint();
		this.revalidate();

	}
	private void updatePercentageDisplay(){
		int perc = residueSizeSlider.getValue();
		percentageDisplay.setText(perc+ " %");
	}

	private int getMaxSequenceLength(){
		int l1 = panel1.getChain().getAtomGroups("amino").size();
		int l2 = panel2.getChain().getAtomGroups("amino").size();
		if ( l1 > l2)
			return l1;
		else return l2;
	}



	/** calculate the float that is used for display.
	 * 1 * scale = size of 1 amino acid (in pixel).
	 * maximum @see MAX_SCALE
	 * @param zoomFactor
	 * @return a float that is the display "scale" - an internal value required for paintin.
	 * user should only interact with the zoomfactor ...
	 */
	private float getScaleForZoom(int zoomFactor){

		if ( zoomFactor > 100)
			zoomFactor = 100;
		if ( zoomFactor < 1)
			zoomFactor = 1;


		int DEFAULT_X_START 	   
		= SequenceScalePanel.DEFAULT_X_START;
		int DEFAULT_X_RIGHT_BORDER = SequenceScalePanel.DEFAULT_X_RIGHT_BORDER;

		int seqLength = getMaxSequenceLength();
		// the maximum width depends on the size of the parent Component

		int width=getWidth();


		float s = width / (float) ( seqLength + DEFAULT_X_START + DEFAULT_X_RIGHT_BORDER );
		//logger.info("scale for 100% " + s + " " + seqLength + " " + zoomFactor);

		s = (100) * s / ( zoomFactor * 1.0f) ;

		if ( s > MAX_SCALE)
			s = MAX_SCALE;

		//logger.info("but changed to " + s);
		return s;
	}

	/** a value of 100 means that the whole sequence should be displayed in the current visible window
	 * a factor of 1 means that one amino acid shoud be drawn as big as possible   
	 *  
	 * @param zoomFactor - a value between 1 and 100
	 *
	 *  
	 */
	public void calcScale(int zoomFactor){

		float s = getScaleForZoom(zoomFactor);
		scale = s;
		//logger.info("calc scale zoom:"+zoomFactor+ " s: " + s);
		panel1.setScale(s);
		panel2.setScale(s);
		panel1.repaint();
		panel2.repaint();


		//return scale;

	}

	public float getScale(){
		return scale;
	}

}

