package org.biojava.nbio.structure.align.gui.jmol;

import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuBar;
import javax.swing.JTextField;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.PDBHeader;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.DisplayAFP;
import org.biojava.nbio.structure.align.gui.MenuCreator;
import org.biojava.nbio.structure.align.gui.MultipleAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.MultipleAlignmentGUI;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentWriter;
import org.biojava.nbio.structure.align.webstart.AligUIManager;
import org.biojava.nbio.structure.jama.Matrix;
import org.jcolorbrewer.ColorBrewer;

/** 
 * A class that provides 3D visualization in Jmol for 
 * {@link MultipleAlignment}s.
 * 
 * @author Aleix Lafita
 * @since 4.1.0
 *
 */
public class MultipleAlignmentJmol extends AbstractAlignmentJmol {

	private MultipleAlignment multAln;
	private List<Atom[]> transformedAtoms;
	private boolean colorByBlocks = false;
	private List<JCheckBox> selectedStructures;

	/**
	 * Default constructor creates an empty JmolPanel window, 
	 * from where alignments can be made through the align menu.
	 */
	public MultipleAlignmentJmol() {
		this(null, null);
	}

	/**
	 * The constructor displays the Mutltiple Alignment in a JmolPanel window.
	 * 
	 * @param multAln: contains the aligned residues.
	 * @param rotatedAtoms: contains the atom coordinates to display, 
	 * so already rotated. 
	 */
	public MultipleAlignmentJmol(MultipleAlignment msa, 
			List<Atom[]> rotatedAtoms) {

		AligUIManager.setLookAndFeel();
		nrOpenWindows++;
		jmolPanel = new JmolPanel();
		frame = new JFrame();
		JMenuBar menu = MenuCreator.initJmolMenu(frame,this,null, msa);

		frame.setJMenuBar(menu);
		this.multAln = msa;
		this.transformedAtoms = rotatedAtoms;
		this.selectedStructures = new ArrayList<JCheckBox>();

		frame.addWindowListener( new WindowAdapter() {

			@Override
			public void windowClosing(WindowEvent e) {

				nrOpenWindows--;
				destroy();

				if ( nrOpenWindows > 0) frame.dispose();
				else {
					MultipleAlignmentGUI gui = MultipleAlignmentGUI.
							getInstanceNoVisibilityChange();
					if (gui.isVisible()) {
						frame.dispose();
						gui.requestFocus();
					} else System.exit(0);
				}
			}
		});

		Container contentPane = frame.getContentPane();
		Box vBox = Box.createVerticalBox();

		jmolPanel.addMouseMotionListener(this);
		jmolPanel.addMouseListener(this);
		jmolPanel.setPreferredSize(
				new Dimension(DEFAULT_WIDTH,DEFAULT_HEIGHT));
		
		vBox.add(jmolPanel);

		/// USER SCRIPTING COMMAND
		JTextField field = new JTextField();

		field.setMaximumSize(new Dimension(Short.MAX_VALUE,30));   
		field.setText(COMMAND_LINE_HELP);
		RasmolCommandListener listener = 
				new RasmolCommandListener(jmolPanel,field);

		field.addActionListener(listener);
		field.addMouseListener(listener);
		field.addKeyListener(listener);
		vBox.add(field);

		/// STRUCTURE SELECTION
		Box hBox0 = Box.createHorizontalBox();
		hBox0.setMaximumSize(new Dimension(Short.MAX_VALUE,30));

		JButton show = new JButton("Show Only: ");
		show.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				jmolPanel.evalString("save selection;");
				String cmd = getJmolString();
				cmd += "; restrict ";
				for (int st=0; st<multAln.size(); st++){
					if (selectedStructures.get(st).isSelected()){
						cmd += "*/"+(st+1)+", ";
					}
				}
				cmd += "none;";
				jmolPanel.executeCmd(cmd+" restore selection;");
			}
		});

		hBox0.add(show);
		hBox0.add(Box.createGlue());

		for (int str=0; str<multAln.size(); str++){
			JCheckBox structureSelection = new JCheckBox(
					multAln.getEnsemble().getStructureNames().get(str));
			hBox0.add(structureSelection);
			hBox0.add(Box.createGlue());
			structureSelection.setSelected(true);
			selectedStructures.add(structureSelection);
		}

		vBox.add(hBox0);	

		/// COMBO BOXES 
		Box hBox1 = Box.createHorizontalBox();
		hBox1.add(Box.createGlue());

		String[] styles = new String[] {"Cartoon", 
				"Backbone", "CPK", "Ball and Stick", 
				"Ligands","Ligands and Pocket"};
		JComboBox style = new JComboBox(styles);

		hBox1.setMaximumSize(new Dimension(Short.MAX_VALUE,30));

		hBox1.add(new JLabel("Style"));
		hBox1.add(style);
		vBox.add(hBox1);
		contentPane.add(vBox);

		style.addActionListener(jmolPanel);

		String[] colorModes = new String[] {"Secondary Structure", 
				"By Chain", "Rainbow", "By Element", "By Amino Acid", 
				"Hydrophobicity" ,"Suggest Domains" , "Show SCOP Domains"};
		JComboBox jcolors = new JComboBox(colorModes);
		jcolors.addActionListener(jmolPanel);
		
		hBox1.add(Box.createGlue());
		hBox1.add(new JLabel("Color"));
		hBox1.add(jcolors);

		String[] cPalette = new String[] {"Spectral", "Set1", "Set2", "Pastel"};
		JComboBox palette = new JComboBox(cPalette);

		palette.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {

				JComboBox source = (JComboBox) e.getSource();
				String value = source.getSelectedItem().toString();
				evalString("save selection; select *; color grey; "
						+ "select ligand; color CPK;");
				
				if (value=="Set1"){
					colorPalette = ColorBrewer.Set1;
				} else if (value=="Set2"){
					colorPalette = ColorBrewer.Set2;
				} else if (value=="Spectral"){
					colorPalette = ColorBrewer.Spectral;
				} else if (value=="Pastel"){
					colorPalette = ColorBrewer.Pastel1;
				}
				String script = getJmolString();
				evalString(script+"; restore selection; ");
			}
		});

		hBox1.add(Box.createGlue());
		hBox1.add(new JLabel("Palette"));
		hBox1.add(palette);


		/// CHECK BOXES
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
		toggleSelection.addItemListener(new ItemListener() {

			@Override
			public void itemStateChanged(ItemEvent e) {
				boolean showSelection = 
						(e.getStateChange() == ItemEvent.SELECTED);

				if (showSelection) jmolPanel.executeCmd("set display selected");
				else jmolPanel.executeCmd("set display off");
			}
		});
		hBox2.add(toggleSelection);
		hBox2.add(Box.createGlue());

		JCheckBox blockColor = new JCheckBox("Color By Block");
		blockColor.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent e) {
				colorByBlocks = (e.getStateChange() == ItemEvent.SELECTED);
				evalString("save selection; "+getJmolString()+
						"; restore selection;");
			}
		});
		
		hBox2.add(blockColor);
		hBox2.add(Box.createGlue());

		vBox.add(hBox2);	

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
		MyJmolStatusListener li = 
				(MyJmolStatusListener) jmolPanel.getStatusListener();
		
		li.setTextField(status);
		frame.pack();
		frame.setVisible(true); 

		initCoords();
		resetDisplay();
	}

	protected void initCoords(){
		try {
			if ( multAln == null ){
				if ( structure != null)
					setStructure(structure);
				else  {
					System.err.println("could not find anything to display!");
					return;
				}
			}
			Structure artificial = 
					DisplayAFP.getAlignedStructure(transformedAtoms);
			
			PDBHeader header = new PDBHeader();
			String title =  multAln.getEnsemble().getAlgorithmName() + 
					" V." +multAln.getEnsemble().getVersion() + " : ";
			
			for (String name:multAln.getEnsemble().getStructureNames()){
				title +=  name + " ";
			}
			
			System.out.println(title);
			header.setTitle(title);
			artificial.setPDBHeader(header);
			setStructure(artificial);
			
		} catch (StructureException e){
			e.printStackTrace();
		}
	}

	@Override
	public void destroy(){
		super.destroy();
		multAln =null;
		transformedAtoms = null;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		String cmd = e.getActionCommand();
		if ( cmd.equals(MenuCreator.TEXT_ONLY)) {
			if ( multAln == null) {
				System.err.println("Currently not viewing an alignment!");
				return;
			}
			System.out.println("Option currently not available for Multiple Alignments");
			//String result = AfpChainWriter.toWebSiteDisplay(textAFP, ca1, ca2) ;
			//DisplayAFP.showAlignmentImage(afpChain, result);

		} else if ( cmd.equals(MenuCreator.PAIRS_ONLY)) {
			if ( multAln == null) {
				System.err.println("Currently not viewing an alignment!");
				return;
			}
			String result = MultipleAlignmentWriter.toAlignedResidues(multAln);
			MultipleAlignmentDisplay.showAlignmentImage(multAln, result);

		} else if (cmd.equals(MenuCreator.ALIGNMENT_PANEL)){
			if ( multAln == null) {
				System.err.println("Currently not viewing an alignment!");
				return;
			}
			try {
				MultipleAlignmentDisplay.showMultipleAligmentPanel(multAln, this, colorPalette);
			} catch (StructureException e1) {
				e1.printStackTrace();
			}

		} else if (cmd.equals(MenuCreator.FATCAT_TEXT)){
			if ( multAln == null) {
				System.err.println("Currently not viewing an alignment!");
				return;
			}
			String result = MultipleAlignmentWriter.toFatCat(multAln);
			result += "\n"+MultipleAlignmentWriter.toTransformMatrices(multAln);
			MultipleAlignmentDisplay.showAlignmentImage(multAln, result);
		}
	}

	/**
	 * Generate a Jmol command String that colors the aligned residues of every structure.
	 */
	private String getJmolString() {

		//Color by blocks if there are flexible alignments (>1 BlockSets) or CPs (>1 Blocks)
		if (colorByBlocks) return getMultiBlockJmolString();
		Color[] colors = colorPalette.getColorPalette(multAln.size());

		StringBuffer j = new StringBuffer();
		j.append(DEFAULT_SCRIPT);

		//Color the equivalent residues of every structure
		StringBuffer sel = new StringBuffer();
		sel.append("select *; color lightgrey; backbone 0.1; ");
		List<List<String>> allPDB = new ArrayList<List<String>>();
		//Loop through all the structures and get the aligned residues
		for (int i=0; i<multAln.size(); i++){

			List<String> pdb = MultipleAlignmentDisplay.getPDBresnum(i,multAln,transformedAtoms.get(i));
			allPDB.add(pdb);
			sel.append("select ");
			int pos = 0;
			for (String res :pdb){
				if (pos > 0)
					sel.append(",");
				pos++;

				sel.append(res);    
				sel.append("/"+(i+1));
			}
			if ( pos == 0)
				sel.append("none");
			sel.append("; backbone 0.4 ; color ["+ colors[i].getRed() +","+ colors[i].getGreen() +","+ colors[i].getBlue() +"]; ");

		}

		j.append(sel);
		j.append("model 0;  ");
		j.append(LIGAND_DISPLAY_SCRIPT);

		//Now select the aligned residues
		StringBuffer buf = new StringBuffer("select ");
		//Loop through all the structures and get the aligned residues
		for (int i=0; i<multAln.size(); i++){
			int count = 0;
			for (String res : allPDB.get(i) ){
				if ( count > 0) buf.append(",");
				buf.append(res);
				buf.append("/"+(i+1));
				count++;
			}
			if (i!=multAln.size()-1) buf.append(",");
		}

		j.append(buf);

		return j.toString();
	}   

	/**
	 * Colors every block of the structures with a different color tonality. It colors each Block differently,
	 * no matter if it is from the same or different BlockSet.
	 */
	private String getMultiBlockJmolString() {

		StringWriter jmol = new StringWriter();
		jmol.append(DEFAULT_SCRIPT);
		jmol.append("select *; color lightgrey; backbone 0.1; ");

		int blockNum = multAln.getBlocks().size();
		Color[] colors = colorPalette.getColorPalette(blockNum);

		//For every structure color all the blocks with the printBlock method
		for (int str=0; str<transformedAtoms.size(); str++){
			jmol.append("select */"+(str+1)+"; color lightgrey; model "+(str+1)+"; ");
			int index = 0;

			for (BlockSet bs : multAln.getBlockSets()) {
				for (Block b : bs.getBlocks() ) {

					List<List<Integer>> alignRes = b.getAlignRes();
					printJmolScript4Block(transformedAtoms.get(str), alignRes, colors[index], jmol, str, index, blockNum);
					index++;
				}
			}
		}

		jmol.append("model 0;  ");
		jmol.append(LIGAND_DISPLAY_SCRIPT);

		return jmol.toString();
	}

	private static void printJmolScript4Block(Atom[] atoms, List<List<Integer>> alignRes, Color blockColor, StringWriter jmol, int str, int colorPos, int blockNum) {

		//Obtain the residues aligned in this block of the structure
		List<String> pdb = new ArrayList<String>();
		for (int i=0;i< alignRes.get(str).size(); i++) {
			//Handle gaps - only color if it is not null
			if (alignRes.get(str).get(i) != null){
				int pos = alignRes.get(str).get(i);
				pdb.add(JmolTools.getPdbInfo(atoms[pos]));
			}
		}

		//Select the aligned residues
		StringBuffer buf = new StringBuffer("select ");
		int count = 0;
		for (String res : pdb){
			if ( count > 0)
				buf.append(",");
			buf.append(res);
			buf.append("/"+(str+1));
			count++;
		}

		buf.append("; backbone 0.4 ; color [" + blockColor.getRed() +"," + blockColor.getGreen() +"," +blockColor.getBlue()+"]; ");

		//Append the string to the global buffer
		jmol.append(buf);
	}

	public void resetDisplay() {

		if (multAln != null) {
			String script = getJmolString();
			//System.out.println(script);
			evalString(script);
			jmolPanel.evalString("hide ligand");
			jmolPanel.evalString("save STATE state_1");
		}
	}

	@Override
	public List<Matrix> getDistanceMatrices() {
		if (multAln == null) return null;
		else return multAln.getEnsemble().getDistanceMatrix();
	}
	
	public void setColorByBlocks(boolean colorByBlocks){
		this.colorByBlocks = colorByBlocks;
		resetDisplay();
	}
}