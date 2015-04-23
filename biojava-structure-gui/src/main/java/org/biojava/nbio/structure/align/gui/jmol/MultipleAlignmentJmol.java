package org.biojava.nbio.structure.align.gui.jmol;


import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.gui.AlignmentGui;
import org.biojava.nbio.structure.align.gui.DisplayAFP;
import org.biojava.nbio.structure.align.gui.MenuCreator;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.AfpChainWriter;
import org.biojava.nbio.structure.align.model.BlockSet;
import org.biojava.nbio.structure.align.model.MultipleAlignment;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.UserConfiguration;
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
public class MultipleAlignmentJmol extends AlignmentJmol {

   MultipleAlignment multAln;
   List<Atom[]> atomArrays;    //rotated atom arrays of every structure

   public static void main(String[] args){
      try {


         UserConfiguration config = new UserConfiguration();
         AtomCache cache = new AtomCache(config);

         Structure struc = cache.getStructure("5pti");

         MultipleAlignmentJmol jmolPanel = new MultipleAlignmentJmol(null, null);

         jmolPanel.setStructure(struc);

         // send some RASMOL style commands to Jmol
         jmolPanel.evalString("select * ; color chain;");
         jmolPanel.evalString("select *; spacefill off; wireframe off; backbone 0.4;  ");

      } catch (Exception e){
         e.printStackTrace();
      }
   }

   /**
    * Default constructor creates an empty window, from where alignments can be made through the align menu.
    */
   public MultipleAlignmentJmol(){
      //don't have a multiple alignment, set it to null...
      this(null, null);
   }

   /**
    * The constructor displays the Mutltiple Alignment in a Jmol panel.
    * @param multAln: contains the aligned residues.
    * @param atomArrays: contains the atom coordinates to display, already rotated.
    */
   public MultipleAlignmentJmol(MultipleAlignment multAln, List<Atom[]> atomArrays) {

      AligUIManager.setLookAndFeel();

      nrOpenWindows++;
      jmolPanel = new JmolPanel();

      frame = new JFrame();

      JMenuBar menu = MenuCreator.initMenu(frame,this,null);

      frame.setJMenuBar(menu);
      //frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
      this.multAln = multAln;
      this.atomArrays = atomArrays;
      
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
    */
   public String getJmolString(MultipleAlignment multAln, List<Atom[]> atomArrays){
     
	  /* TODO Color by block is disabled right now
      if (multAln.getBlockSets().get(0).getBlockNum() > 1){
         return getMultiBlockJmolString(multAln, atomArrays);
      }*/

      StringBuffer j = new StringBuffer();
      j.append(DEFAULT_SCRIPT);
      
      //Set one color for every structure in the multiple alignment
      Color[] colors = ColorBrewer.Set1.getColorPalette(multAln.getSize());
      
      //Color the equivalent residues of every structure
      StringBuffer sel = new StringBuffer();
      sel.append("select *; color lightgrey; ");
      List<List<String>> allPDB = new ArrayList<List<String>>();
      //Loop through all the structures and get the aligned residues
      for (int i=0; i<multAln.getSize(); i++){
    	  
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
          sel.append("; backbone 0.6 ; color ["+ colors[i].getRed() +","+ colors[i].getGreen() +","+ colors[i].getBlue() +"]; ");
          
      }
      
      j.append(sel);
      j.append("model 0;  ");
      j.append(LIGAND_DISPLAY_SCRIPT);
      
      //Now select the aligned residues
      StringBuffer buf = new StringBuffer("select ");
      //Loop through all the structures and get the aligned residues
      for (int i=0; i<multAln.getSize(); i++){
    	  int count = 0;
	      for (String res : allPDB.get(i) ){
	         if ( count > 0) buf.append(",");
	         buf.append(res);
	         buf.append("/"+(i+1));
	         count++;
	      }
	      if (i!=multAln.getSize()-1) buf.append(",");
      }

      j.append(buf);

      return j.toString();
   }
   
   /**
    * Only select and color one aligned block, indicated by the blockNr.
    */
   public static String getJmolScript4Block(MultipleAlignment multAln, List<Atom[]> atomArrays, int blockNr){
	   
	   BlockSet optAln = multAln.getBlockSets().get(0);
	   int blockNum = optAln.getBlockNum();
	   
	   if ( blockNr >= blockNum)
		   return DEFAULT_SCRIPT;
	   
	   if ( blockNum == 0)
		   return DEFAULT_SCRIPT;

	   StringWriter jmol = new StringWriter();
	   jmol.append(DEFAULT_SCRIPT);

	   jmol.append("select *; color lightgrey; ");
	      
	   printJmolScript4Block(atomArrays, blockNum, optAln, jmol, blockNr);
	   
	   jmol.append("model 0;  ");
	   jmol.append(LIGAND_DISPLAY_SCRIPT);
	   //System.out.println(jmol);
	   return jmol.toString();

	   
   }
   

   private static String getMultiBlockJmolString(MultipleAlignment multAln, List<Atom[]> atomArrays)
   {

	  BlockSet optAln = multAln.getBlockSets().get(0);
	  int blockNum = optAln.getBlockNum();      
	
	  if ( blockNum == 0 )
	     return DEFAULT_SCRIPT;
	
	  StringWriter jmol = new StringWriter();
	  jmol.append(DEFAULT_SCRIPT);
	
	  jmol.append("select */2; color lightgrey; model 2; ");
	  
	  //for(int bk = 0; bk < blockNum; bk ++) printJmolScript4Block(atomArrays, blockNum, optAln, jmol, bk);
	  //TODO Print by block is disabled right now
	  
	  jmol.append("model 0;  ");
	  jmol.append(LIGAND_DISPLAY_SCRIPT);
	  
	  return jmol.toString();


   }

private static void printJmolScript4Block(List<Atom[]> atomArrays, int blockNum,
		BlockSet optAln, StringWriter jmol, int bk) {
	
		//TODO disabled the coloring by block
	 /*//the block nr determines the color...
	 int colorPos = bk;
	 
	 Color c1;
	 Color c2;
	 //If the colors for the block are specified in AFPChain use them, otherwise the default ones are calculated
		 
	 if ( colorPos > ColorUtils.colorWheel.length){
	    colorPos = ColorUtils.colorWheel.length % colorPos ;
	 }
	 
	 Color end1 = ColorUtils.rotateHue(ColorUtils.orange,  (1.0f  / 24.0f) * blockNum  );
	 Color end2 = ColorUtils.rotateHue(ColorUtils.cyan,    (1.0f  / 24.0f) * (blockNum +1)  ) ;
	 	 
	 c1   = ColorUtils.getIntermediate(ColorUtils.orange, end1, blockNum, bk);
	 c2   = ColorUtils.getIntermediate(ColorUtils.cyan, end2, blockNum, bk);
	 
	 List<String> pdb1 = new ArrayList<String>();
	 List<String> pdb2 = new ArrayList<String>();
	 for ( int i=0;i< optLen[bk];i++) {
	    ///
	    int pos1 = optAln[bk][0][i];
	    pdb1.add(JmolTools.getPdbInfo(ca1[pos1]));
	    int pos2 = optAln[bk][1][i];
	    pdb2.add(JmolTools.getPdbInfo(ca2[pos2]));
	 }

	 // and now select the aligned residues...
	 StringBuffer buf = new StringBuffer("select ");
	 int count = 0;
	 for (String res : pdb1 ){
	    if ( count > 0)
	       buf.append(",");
	    buf.append(res);
	    buf.append("/1");
	    count++;
	 }

	 buf.append("; backbone 0.6 ; color [" + c1.getRed() +"," + c1.getGreen() +"," +c1.getBlue()+"]; select ");
	 
	 count = 0;
	 for (String res :pdb2 ){
	    if ( count > 0)
	       buf.append(",");
	 
	    buf.append(res);
	    buf.append("/2");
	    count++;
	 }
	 //buf.append("; set display selected;");

	 buf.append("; backbone 0.6 ; color [" + c2.getRed() +"," + c2.getGreen() +"," +c2.getBlue()+"];");

	 // now color this block:
	 jmol.append(buf);*/
}

   public void resetDisplay(){
	   
      if (multAln != null) {
         String script = getJmolString(multAln, atomArrays);
         //System.out.println(script);
         evalString(script);
         jmolPanel.evalString("save STATE state_1");
      }
   }

}