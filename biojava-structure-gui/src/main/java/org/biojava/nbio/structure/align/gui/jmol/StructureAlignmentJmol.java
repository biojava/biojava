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
 * Created on 24.05.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.nbio.structure.align.gui.jmol;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.gui.AlignmentGui;
import org.biojava.nbio.structure.align.gui.DisplayAFP;
import org.biojava.nbio.structure.align.gui.MenuCreator;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.AfpChainWriter;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.ResourceManager;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.align.webstart.AligUIManager;
import org.biojava.nbio.structure.gui.util.color.ColorUtils;
import org.biojava.nbio.structure.jama.Matrix;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import java.awt.*;
import java.awt.event.*;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;

/**
 * A class that provides a simple GUI for Jmol
 *
 * @author Andreas Prlic
 * @since 1.6
 *
 */
public class StructureAlignmentJmol extends AbstractAlignmentJmol implements ChangeListener {

	private Atom[] ca1;
	private Atom[] ca2;
	private AFPChain afpChain;

	private static final String LIGAND_DISPLAY_SCRIPT = ResourceManager.getResourceManager("ce")
			.getString("default.ligand.jmol.script");

	public static void main(String[] args) {

		try {

			UserConfiguration config = new UserConfiguration();
			AtomCache cache = new AtomCache(config);

			Structure struc = cache.getStructure("5pti");

			StructureAlignmentJmol jmolPanel = new StructureAlignmentJmol(null, null, null);

			jmolPanel.setStructure(struc);

			// send some RASMOL style commands to Jmol
			jmolPanel.evalString("select * ; color chain;");
			jmolPanel.evalString("select *; spacefill off; wireframe off; backbone 0.4;  ");

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public StructureAlignmentJmol() {
		// don;t have an afpChain, set it to null...
		this(null, null, null);
	}

	public StructureAlignmentJmol(AFPChain afpChain, Atom[] ca1, Atom[] ca2) {

		AligUIManager.setLookAndFeel();

		nrOpenWindows++;
		jmolPanel = new JmolPanel();

		frame = new JFrame();

		JMenuBar menu = MenuCreator.initJmolMenu(frame,this, afpChain, null);

		frame.setJMenuBar(menu);
		//frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		this.afpChain = afpChain;
		this.ca1 = ca1;
		this.ca2 = ca2;

		frame.addWindowListener( new WindowAdapter()
		{

			@Override
			public void windowClosing(WindowEvent e) {

				nrOpenWindows--;

				destroy();

				if ( nrOpenWindows > 0) {

					frame.dispose();
				}
				else  {
					// check if AlignmentGUI is visible..

					AlignmentGui gui = AlignmentGui.getInstanceNoVisibilityChange();
					if ( gui.isVisible()) {
						frame.dispose();
						gui.requestFocus();
					} else {
						System.exit(0);
					}
				}
			}
		});

		Container contentPane = frame.getContentPane();

		Box vBox = Box.createVerticalBox();

		//try {



		jmolPanel.addMouseMotionListener(this);
		jmolPanel.addMouseListener(this);

		//		} catch (ClassNotFoundException e){
		//			e.printStackTrace();
		//			System.err.println("Could not find Jmol in classpath, please install first. http://www.jmol.org");
		//			return;
		//		}
		jmolPanel.setPreferredSize(new Dimension(DEFAULT_WIDTH,DEFAULT_HEIGHT));
		vBox.add(jmolPanel);


		/// USER SCRIPTING COMMAND
		JTextField field = new JTextField();

		field.setMaximumSize(new Dimension(Short.MAX_VALUE,30));
		field.setText(COMMAND_LINE_HELP);
		RasmolCommandListener listener = new RasmolCommandListener(jmolPanel,field) ;

		field.addActionListener(listener);
		field.addMouseListener(listener);
		field.addKeyListener(listener);
		vBox.add(field);


		/// COMBO BOXES
		Box hBox1 = Box.createHorizontalBox();
		hBox1.add(Box.createGlue());

		String[] styles = new String[] { "Cartoon", "Backbone", "CPK", "Ball and Stick", "Ligands","Ligands and Pocket"};
		JComboBox style = new JComboBox(styles);

		hBox1.setMaximumSize(new Dimension(Short.MAX_VALUE,30));

		hBox1.add(new JLabel("Style"));
		hBox1.add(style);
		vBox.add(hBox1);
		contentPane.add(vBox);

		style.addActionListener(jmolPanel);

		String[] colorModes = new String[] { "Secondary Structure", "By Chain", "Rainbow", "By Element", "By Amino Acid", "Hydrophobicity" ,"Suggest Domains" , "Show SCOP Domains"};
		JComboBox colors = new JComboBox(colorModes);
		colors.addActionListener(jmolPanel);
		hBox1.add(Box.createGlue());
		hBox1.add(new JLabel("Color"));
		hBox1.add(colors);


		// CHeck boxes
		Box hBox2 = Box.createHorizontalBox();
		hBox2.setMaximumSize(new Dimension(Short.MAX_VALUE,30));

		JButton resetDisplay = new JButton("Reset Display");

		resetDisplay.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				System.out.println("reset!!");
				jmolPanel.executeCmd("restore STATE state_1");

			}
		});
			
		hBox2.add(resetDisplay); 
		hBox2.add(Box.createGlue());


		JCheckBox toggleSelection = new JCheckBox("Show Selection");
		toggleSelection.addItemListener(
				new ItemListener() {

					@Override
					public void itemStateChanged(ItemEvent e) {
						boolean showSelection = (e.getStateChange() == ItemEvent.SELECTED);

						if (showSelection){
							jmolPanel.executeCmd("set display selected");
						} else {
							jmolPanel.executeCmd("set display off");
						}

					}
				}
				);



		hBox2.add(toggleSelection);

		hBox2.add(Box.createGlue());
		vBox.add(hBox2);
		
		
		// ZOOM SLIDER
		Box hBox3 = Box.createHorizontalBox();
		hBox3.setMaximumSize(new Dimension(Short.MAX_VALUE,30));
		
		JLabel sliderLabel = new JLabel("Zoom");
		
		hBox3.add(Box.createGlue()); 
		hBox3.add(sliderLabel);
        
		JSlider zoomSlider = new JSlider(JSlider.HORIZONTAL,0,500,100);
		
		zoomSlider.addChangeListener(this);
		
		zoomSlider.setMajorTickSpacing(100);
		zoomSlider.setPaintTicks(true);
		
		Hashtable<Integer, JLabel> labelTable = new Hashtable<Integer, JLabel>();
		labelTable.put(new Integer(0),new JLabel("0%"));
		labelTable.put(new Integer(100),new JLabel("100%"));
		labelTable.put(new Integer(200),new JLabel("200%"));
		labelTable.put(new Integer(300),new JLabel("300%"));
		labelTable.put(new Integer(400),new JLabel("400%"));
		labelTable.put(new Integer(500),new JLabel("500%"));
		
		zoomSlider.setLabelTable(labelTable);
		zoomSlider.setPaintLabels(true);
		
		hBox3.add(zoomSlider); 
		hBox3.add(Box.createGlue());
		
		// SPIN CHECKBOX
		JCheckBox toggleSpin = new JCheckBox("Spin");
		toggleSpin.addItemListener(
				new ItemListener() {

					@Override
					public void itemStateChanged(ItemEvent e) {
						boolean spinOn = (e.getStateChange() == ItemEvent.SELECTED);

						if (spinOn){
							jmolPanel.executeCmd("spin ON");
						} else {
							jmolPanel.executeCmd("spin OFF");
						}

					}
				}
				);
		
		
		hBox3.add(toggleSpin); 
		hBox3.add(Box.createGlue());
		
		vBox.add(hBox3);

		
		// STATUS DISPLAY

		Box hBox = Box.createHorizontalBox();

		status = new JTextField();
		status.setBackground(Color.white);
		status.setEditable(false);
		status.setMaximumSize(new Dimension(Short.MAX_VALUE,30));
		status.setPreferredSize(new Dimension(DEFAULT_WIDTH / 2,30));
		status.setMinimumSize(new Dimension(DEFAULT_WIDTH / 2,30));
		hBox.add(status);
		text = new JTextField();
		text.setBackground(Color.white);
		text.setMaximumSize(new Dimension(Short.MAX_VALUE,30));
		text.setPreferredSize(new Dimension(DEFAULT_WIDTH / 2,30));
		text.setMinimumSize(new Dimension(DEFAULT_WIDTH / 2,30));
		text.setText("Display of Atom info");
		text.setEditable(false);
		hBox.add(text);

		vBox.add(hBox);



		contentPane.add(vBox);
		MyJmolStatusListener li = (MyJmolStatusListener) jmolPanel.getStatusListener();
		li.setTextField(status);
		frame.pack();
		frame.setVisible(true);


		// init coordinates

		initCoords();

		resetDisplay();

	}

	@Override
	protected void initCoords() {
		try {
			if (ca1 == null || ca2 == null) {
				if (structure != null)
					setStructure(structure);
				else {
					// System.err.println("could not find anything to
					// display!");
					return;
				}
			}
			Structure artificial = AlignmentTools.getAlignedStructure(ca1, ca2);
			PDBHeader header = new PDBHeader();
			String title = afpChain.getAlgorithmName() + " V." + afpChain.getVersion() + " : " + afpChain.getName1()
					+ " vs. " + afpChain.getName2();
			header.setTitle(title);
			artificial.setPDBHeader(header);
			setStructure(artificial);
		} catch (StructureException e) {
			e.printStackTrace();
		}
	}

	@Override
	public void destroy() {
		super.destroy();
		afpChain = null;
		ca1 = null;
		ca2 = null;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		String cmd = e.getActionCommand();
		if (cmd.equals(MenuCreator.TEXT_ONLY)) {
			if (afpChain == null) {
				System.err.println("Currently not viewing an alignment!");
				return;
			}
			// Clone the AFPChain to not override the FatCat numbers in alnsymb
			AFPChain textAFP = (AFPChain) afpChain.clone();
			String result = AfpChainWriter.toWebSiteDisplay(textAFP, ca1, ca2);

			DisplayAFP.showAlignmentImage(afpChain, result);

		} else if (cmd.equals(MenuCreator.PAIRS_ONLY)) {
			if (afpChain == null) {
				System.err.println("Currently not viewing an alignment!");
				return;
			}
			String result = AfpChainWriter.toAlignedPairs(afpChain, ca1, ca2);

			DisplayAFP.showAlignmentImage(afpChain, result);

		} else if (cmd.equals(MenuCreator.ALIGNMENT_PANEL)) {
			if (afpChain == null) {
				System.err.println("Currently not viewing an alignment!");
				return;
			}
			try {
				DisplayAFP.showAlignmentPanel(afpChain, ca1, ca2, this);
			} catch (Exception e1) {
				e1.printStackTrace();
				return;
			}

		} else if (cmd.equals(MenuCreator.FATCAT_TEXT)) {
			if (afpChain == null) {
				System.err.println("Currently not viewing an alignment!");
				return;
			}
			String result = afpChain.toFatcat(ca1, ca2);
			result += AFPChain.newline;
			result += afpChain.toRotMat();
			DisplayAFP.showAlignmentImage(afpChain, result);
		}
	}

	public static String getJmolString(AFPChain afpChain, Atom[] ca1, Atom[] ca2) {

		if (afpChain.getBlockNum() > 1) {
			return getMultiBlockJmolScript(afpChain, ca1, ca2);
		}

		StringBuffer j = new StringBuffer();
		j.append(DEFAULT_SCRIPT);

		// now color the equivalent residues ...
		StringBuffer sel = new StringBuffer();
		List<String> pdb1 = DisplayAFP.getPDBresnum(0, afpChain, ca1);
		sel.append("select ");
		int pos = 0;
		for (String res : pdb1) {
			if (pos > 0)
				sel.append(",");
			pos++;

			sel.append(res);
			sel.append("/1");
		}
		if (pos == 0)
			sel.append("none");
		sel.append(";");
		sel.append("backbone 0.6 ;   color orange;");
		sel.append("select */2; color lightgrey; model 2; ");
		// jmol.evalString("select */2; color lightgrey; model 2; ");
		List<String> pdb2 = DisplayAFP.getPDBresnum(1, afpChain, ca2);
		sel.append("select ");
		pos = 0;
		for (String res : pdb2) {
			if (pos > 0)
				sel.append(",");
			pos++;

			sel.append(res);
			sel.append("/2");
		}
		if (pos == 0)
			sel.append("none");
		sel.append("; backbone 0.6 ;   color cyan;");
		// System.out.println(sel);
		j.append(sel);
		// now show both models again.
		j.append("model 0;  ");
		j.append(LIGAND_DISPLAY_SCRIPT);
		// color [object] cpk , set defaultColors Jmol , set defaultColors
		// Rasmol

		// and now select the aligned residues...
		StringBuffer buf = new StringBuffer("select ");
		int count = 0;
		for (String res : pdb1) {
			if (count > 0)
				buf.append(",");
			buf.append(res);
			buf.append("/1");
			count++;
		}

		for (String res : pdb2) {
			buf.append(",");
			buf.append(res);
			buf.append("/2");
		}
		// buf.append("; set display selected;");

		j.append(buf);

		return j.toString();
	}

	public static String getJmolScript4Block(AFPChain afpChain, Atom[] ca1, Atom[] ca2, int blockNr) {
		int blockNum = afpChain.getBlockNum();

		if (blockNr >= blockNum)
			return DEFAULT_SCRIPT;

		int[] optLen = afpChain.getOptLen();
		int[][][] optAln = afpChain.getOptAln();

		if (optLen == null)
			return DEFAULT_SCRIPT;

		StringWriter jmol = new StringWriter();
		jmol.append(DEFAULT_SCRIPT);

		jmol.append("select */2; color lightgrey; model 2; ");

		printJmolScript4Block(ca1, ca2, blockNum, optLen, optAln, jmol, blockNr);

		jmol.append("model 0;  ");
		jmol.append(LIGAND_DISPLAY_SCRIPT);
		// System.out.println(jmol);
		return jmol.toString();

	}

	private static String getMultiBlockJmolScript(AFPChain afpChain, Atom[] ca1, Atom[] ca2) {

		int blockNum = afpChain.getBlockNum();
		int[] optLen = afpChain.getOptLen();
		int[][][] optAln = afpChain.getOptAln();

		if (optLen == null)
			return DEFAULT_SCRIPT;

		StringWriter jmol = new StringWriter();
		jmol.append(DEFAULT_SCRIPT);

		jmol.append("select */2; color lightgrey; model 2; ");

		for (int bk = 0; bk < blockNum; bk++) {

			printJmolScript4Block(ca1, ca2, blockNum, optLen, optAln, jmol, bk);
		}
		jmol.append("model 0;  ");

		jmol.append(LIGAND_DISPLAY_SCRIPT);
		// System.out.println(jmol);
		return jmol.toString();

	}

	private static void printJmolScript4Block(Atom[] ca1, Atom[] ca2, int blockNum, int[] optLen, int[][][] optAln,
			StringWriter jmol, int bk) {
		// the block nr determines the color...
		int colorPos = bk;

		Color c1;
		Color c2;
		// If the colors for the block are specified in AFPChain use them,
		// otherwise the default ones are calculated

		if (colorPos > ColorUtils.colorWheel.length) {
			colorPos = ColorUtils.colorWheel.length % colorPos;
		}

		Color end1 = ColorUtils.rotateHue(ColorUtils.orange, (1.0f / 24.0f) * blockNum);
		Color end2 = ColorUtils.rotateHue(ColorUtils.cyan, (1.0f / 24.0f) * (blockNum + 1));

		c1 = ColorUtils.getIntermediate(ColorUtils.orange, end1, blockNum, bk);
		c2 = ColorUtils.getIntermediate(ColorUtils.cyan, end2, blockNum, bk);

		List<String> pdb1 = new ArrayList<String>();
		List<String> pdb2 = new ArrayList<String>();
		for (int i = 0; i < optLen[bk]; i++) {
			///
			int pos1 = optAln[bk][0][i];
			pdb1.add(JmolTools.getPdbInfo(ca1[pos1]));
			int pos2 = optAln[bk][1][i];
			pdb2.add(JmolTools.getPdbInfo(ca2[pos2]));
		}

		// and now select the aligned residues...
		StringBuffer buf = new StringBuffer("select ");
		int count = 0;
		for (String res : pdb1) {
			if (count > 0)
				buf.append(",");
			buf.append(res);
			buf.append("/1");
			count++;
		}

		buf.append("; backbone 0.6 ; color [" + c1.getRed() + "," + c1.getGreen() + "," + c1.getBlue() + "]; select ");

		count = 0;
		for (String res : pdb2) {
			if (count > 0)
				buf.append(",");

			buf.append(res);
			buf.append("/2");
			count++;
		}
		// buf.append("; set display selected;");

		buf.append("; backbone 0.6 ; color [" + c2.getRed() + "," + c2.getGreen() + "," + c2.getBlue() + "];");

		// now color this block:
		jmol.append(buf);
	}

	@Override
	public void resetDisplay() {

		if (afpChain != null && ca1 != null && ca2 != null) {
			String script = getJmolString(afpChain, ca1, ca2);
			// System.out.println(script);
			evalString(script);
			jmolPanel.evalString("save STATE state_1");
		}
	}

	@Override
	public List<Matrix> getDistanceMatrices() {
		if (afpChain == null)
			return null;
		else
			return Arrays.asList(afpChain.getDisTable1(), afpChain.getDisTable2());
	}

	@Override
	public void stateChanged(ChangeEvent e) {
		JSlider source = (JSlider) e.getSource();
		if (!source.getValueIsAdjusting()) {
			int zoomValue = (int) source.getValue();
			jmolPanel.executeCmd("zoom " + zoomValue);
		}	
	}

	
}
