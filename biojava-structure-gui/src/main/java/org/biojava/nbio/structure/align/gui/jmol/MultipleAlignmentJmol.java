package org.biojava.nbio.structure.align.gui.jmol;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.gui.AlignmentGui;
import org.biojava.nbio.structure.align.gui.DisplayAFP;
import org.biojava.nbio.structure.align.gui.MenuCreator;
import org.biojava.nbio.structure.align.model.BlockSet;
import org.biojava.nbio.structure.align.model.MultipleAlignment;
import org.biojava.nbio.structure.align.model.StructureAlignmentException;
import org.biojava.nbio.structure.align.webstart.AligUIManager;
import org.biojava.nbio.structure.gui.util.color.ColorUtils;
import org.jcolorbrewer.ColorBrewer;

import javax.swing.*;

import java.awt.*;
import java.awt.event.*;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

/** A class that provides GUI for multiple alignments. Code adapted from the StructuralAlignmentJmol class.
 * 
 * @author Aleix Lafita
 *
 */
public class MultipleAlignmentJmol extends AbstractAlignmentJmol {

	private Color[] colors;
	private MultipleAlignment multAln;
	private List<Atom[]> atomArrays;    //already rotated atom arrays of every structure

   /**
    * Default constructor creates an empty window, from where alignments can be made through the align menu.
    * @throws StructureAlignmentException 
    */
   public MultipleAlignmentJmol() throws StructureAlignmentException{
      //don't have a multiple alignment, set it to null...
      this(null, null);
   }
   
   public MultipleAlignmentJmol(MultipleAlignment multAln, List<Atom[]> rotatedAtoms) throws StructureAlignmentException {
	   //Default colors: color set
	   this(multAln, rotatedAtoms, ColorBrewer.Set1.getColorPalette(multAln.size()));
   }

   /**
    * The constructor displays the Mutltiple Alignment in a Jmol panel.
    * @param multAln: contains the aligned residues.
    * @param atomArrays: contains the atom coordinates to display, already rotated.
    * @throws StructureAlignmentException 
    */
   public MultipleAlignmentJmol(MultipleAlignment msa, List<Atom[]> rotatedAtoms, Color[] col) throws StructureAlignmentException {

      AligUIManager.setLookAndFeel();

      nrOpenWindows++;
      jmolPanel = new JmolPanel();

      frame = new JFrame();

      JMenuBar menu = MenuCreator.initMenu(frame,this,null);

      frame.setJMenuBar(menu);
      //frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
      this.multAln = msa;
      this.atomArrays = rotatedAtoms;
      this.colors = col;
      
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
		JComboBox jcolors = new JComboBox(colorModes);
		jcolors.addActionListener(jmolPanel);
		hBox1.add(Box.createGlue());
		hBox1.add(new JLabel("Color"));
		hBox1.add(jcolors);
		
		String[] colorPattelete = new String[] {"Set1", "Set2", "Spectral", "Pastel"};
		JComboBox pattelete = new JComboBox(colorPattelete);
		
		pattelete.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				int size = atomArrays.size();
				JComboBox source = (JComboBox) e.getSource();
				String value = source.getSelectedItem().toString();
				evalString("save selection; select *; color grey; select ligand; color CPK;");
				if (value=="Set1"){
					colors = ColorBrewer.Set1.getColorPalette(size);
				} else if (value=="Set2"){
					colors = ColorBrewer.Set2.getColorPalette(size);
				} else if (value=="Spectral"){
					colors = ColorBrewer.Spectral.getColorPalette(size);
				} else if (value=="Pastel"){
					colors = ColorBrewer.Pastel1.getColorPalette(size);
				}
				String script="";
				try {
					script = getJmolString(multAln, atomArrays, colors);
				} catch (StructureAlignmentException e1) {
					e1.printStackTrace();
				}
				evalString(script+"restore selection; ");
			}
		});

		hBox1.add(Box.createGlue());
		hBox1.add(new JLabel("Pattelete"));
		hBox1.add(pattelete);
		

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

      initCoords();

      resetDisplay();
      
   }
   protected void initCoords(){
      try {
         if ( multAln == null ){
            if ( structure != null)
               setStructure(structure);
            else  {
               //System.err.println("could not find anything to display!");
               return;
            }
         }
         Structure artificial = DisplayAFP.getAlignedStructure(atomArrays);
         PDBHeader header = new PDBHeader();
         String title =  multAln.getAlgorithmName() + " V." +multAln.getVersion() + " : ";
         for (String name:multAln.getStructureNames()) title +=  name + " ";
         System.out.println(title);
         header.setTitle(title);
         artificial.setPDBHeader(header);
         setStructure(artificial);
      } catch (StructureException e){
         e.printStackTrace();
      }
   }

   public void destroy(){
	  super.destroy();
      multAln =null;
      atomArrays = null;
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
         System.out.println("Option currently not available for Multiple Alignments");
         //String result = AfpChainWriter.toAlignedPairs(afpChain, ca1, ca2) ;
         //DisplayAFP.showAlignmentImage(afpChain, result);
         
      } else if (cmd.equals(MenuCreator.ALIGNMENT_PANEL)){
         if ( multAln == null) {
            System.err.println("Currently not viewing an alignment!");
            return;
         }
         System.out.println("Option currently not available for Multiple Alignments");
         //DisplayAFP.showAlignmentImage(afpChain, ca1, ca2, this);

      } else if (cmd.equals(MenuCreator.FATCAT_TEXT)){
         if ( multAln == null) {
            System.err.println("Currently not viewing an alignment!");
            return;
         }
         System.out.println("Option currently not available for Multiple Alignments");
         //String result = afpChain.toFatcat(ca1, ca2) ;
         //result += AFPChain.newline;
         //result += afpChain.toRotMat();
         //DisplayAFP.showAlignmentImage(afpChain, result);
      }
   }

   /**
    * Generate a Jmol command String that colors the aligned residues of every structure.
    * 
    * @param multAln
    * @param atomArrays
    * @throws StructureAlignmentException 
    */
   private static String getJmolString(MultipleAlignment multAln, List<Atom[]> atomArrays, Color[] colors) throws StructureAlignmentException{
     
	  //Color by blocks if there are flexible alignments (>1 BlockSets) or CPs (>1 Blocks)
      if (multAln.getBlockSetNum() > 1) return getMultiBlockJmolString(multAln, atomArrays,colors);
      else if (multAln.getBlockSets().get(0).getBlockNum() > 1) return getMultiBlockJmolString(multAln, atomArrays,colors);

      StringBuffer j = new StringBuffer();
      j.append(DEFAULT_SCRIPT);
      
      //Color the equivalent residues of every structure
      StringBuffer sel = new StringBuffer();
      sel.append("select *; color lightgrey; backbone 0.1; ");
      List<List<String>> allPDB = new ArrayList<List<String>>();
      //Loop through all the structures and get the aligned residues
      for (int i=0; i<multAln.size(); i++){
    	  
    	  List<String> pdb = DisplayAFP.getPDBresnum(i,multAln,atomArrays.get(i));
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
   private static String getMultiBlockJmolString(MultipleAlignment multAln, List<Atom[]> atomArrays, Color[] colors) throws StructureAlignmentException {

	  StringWriter jmol = new StringWriter();
	  jmol.append(DEFAULT_SCRIPT);
	  
	  int blockNum = 0;
	  for (int bs = 0; bs < multAln.getBlockSetNum(); bs ++) {
		  for (int b = 0; b < multAln.getBlockSets().get(bs).getBlockNum(); b ++) {
			  blockNum++;
		  }
	  }
	  
	  //For every structure color all the blocks with the printBlock method
	  for (int str=0; str<atomArrays.size(); str++){
		  jmol.append("select */"+(str+1)+"; color lightgrey; model "+(str+1)+"; ");
		  int index = 0;
		  
		  for (int bs = 0; bs < multAln.getBlockSetNum(); bs ++) {
			  for (int b = 0; b < multAln.getBlockSets().get(bs).getBlockNum(); b ++) {
				  
				  List<List<Integer>> alignRes = multAln.getBlockSets().get(bs).getBlocks().get(b).getAlignRes();
				  printJmolScript4Block(atomArrays.get(str), alignRes, colors, jmol, str, index, blockNum);
				  index++;
			  }
		  }
	  }
	  
	  jmol.append("model 0;  ");
	  jmol.append(LIGAND_DISPLAY_SCRIPT);
	  
	  return jmol.toString();
   }
   
   private static void printJmolScript4Block(Atom[] atoms, List<List<Integer>> alignRes, Color[] colors, StringWriter jmol, int str, int colorPos, int blockNum) {
	 
	 Color start = ColorUtils.lighter(colors[str%colors.length], 1);
	 double fraction = (colorPos*1.0)/blockNum;
	 Color color = ColorUtils.darker(start, fraction);
	 
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

	 buf.append("; backbone 0.6 ; color [" + color.getRed() +"," + color.getGreen() +"," +color.getBlue()+"]; ");
	 
	 //Append the string to the global buffer
	 jmol.append(buf);
}

   public void resetDisplay() throws StructureAlignmentException{
	   
      if (multAln != null) {
         String script = getJmolString(multAln, atomArrays,colors);
         System.out.println(script);
         evalString(script);
         jmolPanel.evalString("save STATE state_1");
         jmolPanel.evalString("hide ligand");
      }
   }

}