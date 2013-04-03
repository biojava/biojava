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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.color.ColorSpace;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.ListCellRenderer;
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
import org.biojava.bio.structure.gui.util.color.*;
import org.biojava.bio.structure.gui.util.color.LinearColorInterpolator.InterpolationDirection;


/** A JPanel that can display the underlying distance matrix 
 * data of the protein structure alignment algorithm. It adds a 
 * JSlider to a JMatrixPanel.
 * 
 * see also JMatrixPanel.
 * 
 */
public class ScaleableMatrixPanel 
extends JPanel 
implements ChangeListener, ActionListener {

	/**
	 * 
	 */
	private static final long serialVersionUID = -8082261434322968652L;

	protected JMatrixPanel mPanel;
	protected JSlider slider;
	protected JScrollPane scroll;
	protected JComboBox coloring;
	
	protected Map<String,ContinuousColorMapper> gradients;
	
	protected static final int SLIDER_STEPS = 8; // Number of minor ticks per unit scaled
	
	
	public static void main(String[] args){

		PDBFileReader pdbr = new PDBFileReader();  
		pdbr.setAutoFetch(true);
		pdbr.setPath("/tmp/");


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
		
		Box gradientBox = Box.createHorizontalBox();
		vBox.add(gradientBox);
		gradientBox.add(new JLabel("Coloring:"));
		gradients = createGradients(); //sets gradients
		//Set first gradient
		this.setCellColor(gradients.values().iterator().next());
		
		coloring = new JComboBox(gradients.keySet().toArray(new String[] {}));
		coloring.setRenderer(new GradientRenderer());
		coloring.addActionListener(this);
		coloring.setMaximumSize(new Dimension(1000,30));
		gradientBox.add(coloring);
		
		int RES_MIN  = 0*SLIDER_STEPS;
		int RES_MAX  = 8*SLIDER_STEPS;
		int RES_INIT = 1*SLIDER_STEPS;

		slider = new JSlider(JSlider.HORIZONTAL, RES_MIN,RES_MAX,RES_INIT);
		slider.setInverted(false);
		slider.setPaintTicks(true);
		//slider.setMinorTickSpacing(1);
		slider.setMajorTickSpacing(SLIDER_STEPS);
		slider.setSnapToTicks(false);
		slider.setPaintLabels(false);
		slider.setPaintTrack(true);
		//slider.setLabelTable(slider.createStandardLabels(SLIDER_STEPS));
		slider.addChangeListener(this);
		slider.setPreferredSize(new Dimension(100,slider.getPreferredSize().height));

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




	protected static Map<String,ContinuousColorMapper> createGradients() {
		SortedMap<String,ContinuousColorMapper> gradients = new TreeMap<String,ContinuousColorMapper>();
		
		int i = 0; //prepend number, since sorted alphabetically
		ColorSpace hsv = HSVColorSpace.getHSVColorSpace();
		LinearColorInterpolator interp;
		GradientMapper gradient;
		
		/*
		DefaultMatrixMapper defaultMapper = new DefaultMatrixMapper(10., .9f);
		gradients.put((++i)+". Default", defaultMapper);
		defaultMapper = new DefaultMatrixMapper(5., .9f);
		gradients.put((++i)+". Sensitive", defaultMapper);
		
		
		gradients.put((++i)+". Rainbow", GradientMapper.getGradientMapper(GradientMapper.RAINBOW_INTENSITY_GRADIENT, 0, 10));
		gradients.put((++i)+". Rainbow", GradientMapper.getGradientMapper(GradientMapper.RAINBOW_INTENSITY_GRADIENT, 10, 0));
		*/
		

		interp = new LinearColorInterpolator(hsv);
		interp.setInterpolationDirection(0, InterpolationDirection.INNER);
		gradient = new GradientMapper(Color.green, Color.black, hsv);
		gradient.put( -50., new Color(hsv,new float[] {0f, .9f, 0f},1f));
		gradient.put( 10., new Color(hsv,new float[] {1f, .9f, 1f},1f));
		gradient.setInterpolator(interp);
		
		gradients.put((++i)+". -50 to 10", gradient);	
		
		
		
		// Mimic DefaultMapper
		interp = new LinearColorInterpolator(hsv);
		interp.setInterpolationDirection(0, InterpolationDirection.INNER);
		gradient = new GradientMapper(Color.green, Color.black, hsv);
		gradient.put( 0., new Color(hsv,new float[] {1f, .9f, 1f},1f));
		gradient.put(10., new Color(hsv,new float[] {0f, .9f, 0f},1f));
		gradient.setInterpolator(interp);
		
		gradients.put((++i)+". Default", gradient);
		
		// Sensitive DefaultMapper
		interp = new LinearColorInterpolator(hsv);
		interp.setInterpolationDirection(0, InterpolationDirection.INNER);
		gradient = new GradientMapper(Color.green, Color.black, hsv);
		gradient.put( 0., new Color(hsv,new float[] {1f, .9f, 1f},1f));
		gradient.put( 5., new Color(hsv,new float[] {0f, .9f, 0f},1f));
		gradient.setInterpolator(interp);
		
		gradients.put((++i)+". Sensitive", gradient);	
		
		// color [0,1]import java.awt.Color;

		interp = new LinearColorInterpolator(hsv);
		interp.setInterpolationDirection(0, InterpolationDirection.INNER);
		gradient = new GradientMapper(Color.green, Color.black, hsv);
		gradient.put( 0., new Color(hsv,new float[] {1f, .9f, 1f},1f));
		gradient.put( 1., new Color(hsv,new float[] {.2f, .9f, 1f},1f));
		gradient.put( 1+1e-6, Color.white);
		gradient.put(5., Color.black);
		gradient.setInterpolator(interp);
		
		gradients.put((++i)+". Emphasize low", gradient);
		
		
		interp = new LinearColorInterpolator(hsv);
		interp.setInterpolationDirection(0, InterpolationDirection.INNER);
		gradient = new GradientMapper(Color.green, Color.black, hsv);
		gradient.put( 0., new Color(hsv,new float[] {0f, .9f, 0f},1f));
		gradient.put( 100., new Color(hsv,new float[] {1f, .9f, 1f},1f));
		gradient.setInterpolator(interp);
		
		gradients.put((++i)+". 0 to 100", gradient);	
		
		// log color
		interp = new LinearColorInterpolator(hsv);
		interp.setInterpolationDirection(0, InterpolationDirection.INNER);
		gradient = new GradientMapper(Color.red, Color.black, hsv);
		gradient.put( 0., new Color(hsv,new float[] {1f, .9f, 1f},1f));
		gradient.put( 10., new Color(hsv,new float[] {0f, .9f, 0f},1f));
		gradient.setInterpolator(interp);
		
		ContinuousColorMapper logGradient = new LogColorMapper(gradient,2);
		gradients.put((++i)+". Logorithmic", logGradient);	
		
		// sqrt color
		interp = new LinearColorInterpolator(hsv);
		interp.setInterpolationDirection(0, InterpolationDirection.INNER);
		gradient = new GradientMapper(Color.red, Color.black, hsv);
		gradient.put( 0., new Color(hsv,new float[] {1f, .9f, 1f},1f));
		gradient.put( 4., new Color(hsv,new float[] {0f, .9f, 0f},1f));
		gradient.setInterpolator(interp);
		
		ContinuousColorMapper sqrtGradient = new SqrtColorMapper(gradient);
		gradients.put((++i)+". Square Root", sqrtGradient);	
		
		return gradients;
	}

	public void stateChanged(ChangeEvent e) {
		
		JSlider source = (JSlider)e.getSource();
		
		if ( source.getValueIsAdjusting()) {
			//return;
		}
		
		//System.out.println("Changed scale to "+source.getValue());
		mPanel.setScale((float)source.getValue()/SLIDER_STEPS);
		
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

	/**
	 * @return the color mapping of the JMatrixPanel
	 */
	public ContinuousColorMapper getCellColor() {
		return mPanel.getCellColor();
	}

	/**
	 * @param cellColor the color mapping of the JMatrixPanel
	 */
	public void setCellColor(ContinuousColorMapper cellColor) {
		mPanel.setCellColor(cellColor);
	}

	/**
	 * A renderer for the the gradient dropdown menu at the top of scaleableMatrixPanes.
	 * Knows how to draw a gradient and nicely label it.
	 * 
	 * @author Spencer Bliven
	 *
	 */
	protected class GradientRenderer extends JPanel
	implements ListCellRenderer {

		private static final long serialVersionUID = -2000575579184232365L;
		private int min,max;
		JLabel title;
		JPanel gradientContainer;
		
		public GradientRenderer() {
			this.setPreferredSize(new Dimension(100,25));
			this.setLayout(new BorderLayout());
			this.min = -1;
			this.max = 10;
			
			JPanel gradientBounds = new JPanel();
			gradientBounds.setLayout(new BorderLayout());
			//gradientBounds.setBorder(BorderFactory.createLineBorder(Color.GRAY));
			gradientBounds.add(new JLabel(Integer.toString(min)),BorderLayout.WEST);
			
			gradientContainer = new JPanel();
			gradientContainer.setLayout(new BorderLayout());
			gradientContainer.setOpaque(false);
			gradientContainer.add(new JLabel("<No gradient>"),BorderLayout.CENTER);
			//gradientContainer.setMinimumSize(new Dimension(50,20));
			//gradientContainer.setBorder(BorderFactory.createLineBorder(Color.DARK_GRAY));
			//gradientBounds.setPreferredSize(new Dimension(200,25));
			gradientBounds.add(gradientContainer,BorderLayout.CENTER);
			gradientBounds.setOpaque(false);
			
			gradientBounds.add(new JLabel(Integer.toString(max)),BorderLayout.EAST);
			
			this.add(gradientBounds, BorderLayout.EAST);
			
			this.title = new JLabel("No gradient");
			//this.title.setOpaque(true);
			this.title.setHorizontalAlignment(JLabel.CENTER);
			this.title.setVerticalAlignment(JLabel.CENTER);
			this.add(this.title,BorderLayout.CENTER);
		}

		/*
		 * This method finds the image and text corresponding
		 * to the selected value and returns the label, set up
		 * to display the text and image.
		 */
		public Component getListCellRendererComponent(
				JList list,
				Object value,
				int index,
				boolean isSelected,
				boolean cellHasFocus) {
			//Get the selected index. (The index param isn't
			//always valid, so just use the value.)
			String gradientLabel = (String)value;
			
			if (isSelected) {
				setBackground(list.getSelectionBackground());
				setForeground(list.getSelectionForeground());
			} else {
				setBackground(list.getBackground());
				setForeground(list.getForeground());
			}
			
			//Set the icon and text.  If icon was null, say so.
			GradientPanel gradPanel = new GradientPanel(gradients.get(gradientLabel),min,max);
			gradPanel.setPreferredSize(new Dimension(100,20));
			//gradPanel.setBorder(BorderFactory.createLineBorder(Color.cyan));
			gradientContainer.removeAll();
			gradientContainer.add(gradPanel,BorderLayout.CENTER);

			title.setText(gradientLabel);

			this.validate();
			
			return this;
		}
	}

	/**
	 * @param e
	 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
	public void actionPerformed(ActionEvent e) {
        JComboBox cb = (JComboBox)e.getSource(); // == coloring
        String gradientName = (String)cb.getSelectedItem();
        ContinuousColorMapper gradient = gradients.get(gradientName);
        assert(gradient != null);
        this.setCellColor(gradient);
        this.repaint();
	}

	
}
