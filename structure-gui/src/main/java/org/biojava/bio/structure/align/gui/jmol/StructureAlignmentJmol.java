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

package org.biojava.bio.structure.align.gui.jmol;


import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;

import java.awt.event.MouseMotionListener;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import javax.swing.Box;
import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JTextField;


import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.PDBHeader;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.align.gui.AlignmentGui;
import org.biojava.bio.structure.align.gui.DisplayAFP;
import org.biojava.bio.structure.align.gui.MenuCreator;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.model.AfpChainWriter;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.ResourceManager;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.align.webstart.AligUIManager;
import org.biojava.bio.structure.align.gui.aligpanel.AligPanel;
import org.biojava.bio.structure.align.gui.jmol.AtomInfo;
import org.biojava.bio.structure.align.gui.jmol.AtomInfoParser;
import org.biojava.bio.structure.align.gui.jmol.JmolPanel;
import org.biojava.bio.structure.align.gui.jmol.MyJmolStatusListener;
import org.biojava.bio.structure.align.gui.jmol.RasmolCommandListener;


import org.jmol.api.JmolViewer;


/** A class that provides a simple GUI for Jmol
 * 
 * @author Andreas Prlic
 * @since 1.6
 *
 *
 *
 */
public class StructureAlignmentJmol implements MouseMotionListener, MouseListener, WindowListener,ActionListener {

   Structure structure; 

   JmolPanel jmolPanel;
   JFrame frame ;
   JTextField text ;
   JTextField status;

   protected static final String COMMAND_LINE_HELP = "enter Jmol scripting command...";
   Atom[] ca1;
   Atom[] ca2;
   AFPChain afpChain;

   private static final int DEFAULT_HEIGHT = 500;

   private static final int DEFAULT_WIDTH = 500;

   public static final String DEFAULT_SCRIPT = ResourceManager.getResourceManager("ce").getString("default.alignment.jmol.script");

   private static final String LIGAND_DISPLAY_SCRIPT = ResourceManager.getResourceManager("ce").getString("default.ligand.jmol.script");

   static int nrOpenWindows = 0;

   public static void main(String[] args){
      try {


         UserConfiguration config = new UserConfiguration();
         config.setSplit(true);
         config.setAutoFetch(true);
         AtomCache cache = new AtomCache(config);

         Structure struc = cache.getStructure("5pti");

         StructureAlignmentJmol jmolPanel = new StructureAlignmentJmol(null,null,null);

         jmolPanel.setStructure(struc);

         // send some RASMOL style commands to Jmol
         jmolPanel.evalString("select * ; color chain;");
         jmolPanel.evalString("select *; spacefill off; wireframe off; backbone 0.4;  ");

      } catch (Exception e){
         e.printStackTrace();
      }
   }

   public StructureAlignmentJmol(){
      // don;t have an afpChain, set it to null...
      this(null, null, null);
   }


   public StructureAlignmentJmol(AFPChain afpChain, Atom[] ca1, Atom[] ca2) {

      AligUIManager.setLookAndFeel();

      nrOpenWindows++;
      jmolPanel = new JmolPanel();

      frame = new JFrame();

      JMenuBar menu = MenuCreator.initMenu(frame,this, afpChain);

      frame.setJMenuBar(menu);
      //frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
      this.afpChain = afpChain;
      this.ca1 = ca1;
      this.ca2 = ca2;

      frame.addWindowListener( new WindowAdapter()
      {

         public void windowClosing(WindowEvent e) {

            nrOpenWindows--;

            if ( nrOpenWindows > 0)
               frame.dispose();
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


      JTextField field = new JTextField();

      field.setMaximumSize(new Dimension(Short.MAX_VALUE,30));   
      field.setText(COMMAND_LINE_HELP);
      RasmolCommandListener listener = new RasmolCommandListener(jmolPanel,field) ;

      field.addActionListener(listener);
      field.addMouseListener(listener);
      field.addKeyListener(listener);



      vBox.add(field);

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
   private void initCoords(){
      try {
         if ( ca1 == null || ca2 == null ){
            if ( structure != null)
               setStructure(structure);
            else  {
               //System.err.println("could not find anything to display!");
               return;
            }
         }
         Structure artificial = DisplayAFP.getAlignedStructure(ca1,ca2);
         PDBHeader header = new PDBHeader();
         String title =  afpChain.getAlgorithmName() + " V." +afpChain.getVersion() + " : " + afpChain.getName1() + " vs. " + afpChain.getName2();
         header.setTitle(title);
         artificial.setPDBHeader(header);
         setStructure(artificial);
      } catch (StructureException e){
         e.printStackTrace();
      }
   }

   public void destroy(){

      jmolPanel.removeMouseListener(this);
      jmolPanel.removeMouseMotionListener(this);
      afpChain =null;
      ca1 = null;
      ca2 = null;
   }

   public void setAtoms(Atom[] atoms){
      Structure s = new StructureImpl();
      Chain c = new ChainImpl();
      c.setName("A");
      for (Atom a: atoms){
         c.addGroup(a.getParent());
      }
      s.addChain(c);
      setStructure(s);

   }


   public JmolPanel getJmolPanel() {
      return jmolPanel;
   }

   public void setJmolPanel(JmolPanel jmolPanel) {
      this.jmolPanel = jmolPanel;
   }

   public void evalString(String rasmolScript){
      if ( jmolPanel == null ){
         System.err.println("please install Jmol first");
         return;
      }
      jmolPanel.evalString(rasmolScript);
   }

   public void setStructure(Structure s) {

      if ( jmolPanel == null ){
         System.err.println("please install Jmol first");
         return;
      }


      setTitle(s.getPDBCode());

      // actually this is very simple
      // just convert the structure to a PDB file

      String pdb = s.toPDB();	
      //System.out.println(s.isNmr());

      //System.out.println(pdb);
      // Jmol could also read the file directly from your file system
      //viewer.openFile("/Path/To/PDB/1tim.pdb");

      //System.out.println(pdb);
      jmolPanel.openStringInline(pdb);

      // send the PDB file to Jmol.
      // there are also other ways to interact with Jmol, e.g make it directly
      // access the biojava structure object, but they require more
      // code. See the SPICE code repository for how to do this.


      structure = s;

   }

   public Structure getStructure(){
      return structure;
   }

   public void setTitle(String label){
      frame.setTitle(label);
      frame.repaint();
   }
   public String getTitle(){
      return frame.getTitle();
   }

   public void mouseDragged(MouseEvent e) {
      // TODO Auto-generated method stub

   }

   public void mouseMoved(MouseEvent e) {

      JmolViewer viewer = jmolPanel.getViewer();


      // needs latest jmol :-/
      int pos = viewer.findNearestAtomIndex( e.getX(), e.getY() );
      if ( pos == -1 ) { return ; }

      String atomInfo = viewer.getAtomInfo(pos);
      text.setText(atomInfo);


   }

   public void mouseClicked(MouseEvent e) {
      // TODO Auto-generated method stub

   }

   public void mouseEntered(MouseEvent e) {
      // TODO Auto-generated method stub

   }

   public void mouseExited(MouseEvent e) {
      // TODO Auto-generated method stub

   }

   public void mousePressed(MouseEvent e) {
      // TODO Auto-generated method stub

   }

   public void mouseReleased(MouseEvent e) {
      JmolViewer viewer = jmolPanel.getViewer();


      int pos = viewer.findNearestAtomIndex( e.getX(), e.getY() );
      if ( pos == -1 ) { return ; }

      String atomInfo = viewer.getAtomInfo(pos);
      status.setText("clicked: " + atomInfo);

      AtomInfo ai = AtomInfoParser.parse(atomInfo);

      String cmd = "select " + ai.getResidueNumber()+":" +ai.getChainId()+"/"+ai.getModelNumber() + "; set display selected;";

      evalString(cmd);

   }

   public void windowActivated(WindowEvent e) {
      // TODO Auto-generated method stub

   }

   public void windowClosed(WindowEvent e) {
      // TODO Auto-generated method stub

   }

   public void windowClosing(WindowEvent e) {
      destroy();

   }

   public void windowDeactivated(WindowEvent e) {
      // TODO Auto-generated method stub

   }

   public void windowDeiconified(WindowEvent e) {
      // TODO Auto-generated method stub

   }

   public void windowIconified(WindowEvent e) {
      // TODO Auto-generated method stub

   }

   public void windowOpened(WindowEvent e) {
      // TODO Auto-generated method stub

   }

   public void actionPerformed(ActionEvent e) {
      String cmd = e.getActionCommand();
      if ( cmd.equals(MenuCreator.TEXT_ONLY)) {
         if ( afpChain == null) {
            System.err.println("Currently not viewing an alignment!");
            return;
         }
         String result = AfpChainWriter.toWebSiteDisplay(afpChain, ca1, ca2) ;

         DisplayAFP.showAlignmentImage(afpChain, result);
      } else if (cmd.equals(MenuCreator.ALIGNMENT_PANEL)){
         if ( afpChain == null) {
            System.err.println("Currently not viewing an alignment!");
            return;
         }
         DisplayAFP.showAlignmentImage(afpChain, ca1, ca2, this);

      } else if (cmd.equals(MenuCreator.FATCAT_TEXT)){
         if ( afpChain == null) {
            System.err.println("Currently not viewing an alignment!");
            return;
         }
         String result = afpChain.toFatcat(ca1, ca2) ;
         result += AFPChain.newline;
         result += afpChain.toRotMat();
         DisplayAFP.showAlignmentImage(afpChain, result);
      }

   }

   public static String getJmolString(AFPChain afpChain, Atom[] ca1, Atom[] ca2){
     
      if ( afpChain.getBlockNum() > 1){
         return getMultiBlockJmolScript( afpChain,  ca1,  ca2);
      }

      
      
      StringBuffer j = new StringBuffer();
      j.append(DEFAULT_SCRIPT);

     
      // now color the equivalent residues ...
      StringBuffer sel = new StringBuffer();
      List<String> pdb1 = DisplayAFP.getPDBresnum(0,afpChain,ca1);
      sel.append("select ");
      int pos = 0;
      for (String res :pdb1 ){
         if ( pos > 0)
            sel.append(",");
         pos++;

         sel.append(res);    
         sel.append("/1");
      }
      if ( pos == 0)
         sel.append("none");
      sel.append(";");
      sel.append("backbone 0.6 ;   color orange;");
      sel.append("select */2; color lightgrey; model 2; ");
      //jmol.evalString("select */2; color lightgrey; model 2; ");        
      List<String> pdb2 = DisplayAFP.getPDBresnum(1,afpChain,ca2);
      sel.append("select ");
      pos = 0;
      for (String res :pdb2 ){
         if ( pos > 0)
            sel.append(",");
         pos++;

         sel.append(res);
         sel.append("/2");
      }
      if ( pos == 0)
         sel.append("none");
      sel.append("; backbone 0.6 ;   color cyan;");
      //System.out.println(sel);
      j.append(sel);
      // now show both models again.
      j.append("model 0;  ");
      j.append(LIGAND_DISPLAY_SCRIPT);
      //color [object] cpk , set defaultColors Jmol , set defaultColors Rasmol  

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

      for (String res :pdb2 ){
         buf.append(",");
         buf.append(res);
         buf.append("/2");
      }
      //buf.append("; set display selected;");

      j.append(buf);

      return j.toString();
   }
   
   /**
    * Make a color darker.
    * 
    * @param color     Color to make darker.
    * @param fraction  Darkness fraction.
    * @return          Darker color.
    */
   public static Color darker (Color color, double fraction)
   {
     int red   = (int) Math.round (color.getRed()   * (1.0 - fraction));
     int green = (int) Math.round (color.getGreen() * (1.0 - fraction));
     int blue  = (int) Math.round (color.getBlue()  * (1.0 - fraction));

     if (red   < 0) red   = 0; else if (red   > 255) red   = 255;
     if (green < 0) green = 0; else if (green > 255) green = 255;
     if (blue  < 0) blue  = 0; else if (blue  > 255) blue  = 255;    

     int alpha = color.getAlpha();

     return new Color (red, green, blue, alpha);
   }

   

   /**
    * Make a color lighter.
    * 
    * @param color     Color to make lighter.
    * @param fraction  Darkness fraction.
    * @return          Lighter color.
    */
   public static Color lighter (Color color, double fraction)
   {
     int red   = (int) Math.round (color.getRed()   * (1.0 + fraction));
     int green = (int) Math.round (color.getGreen() * (1.0 + fraction));
     int blue  = (int) Math.round (color.getBlue()  * (1.0 + fraction));

     if (red   < 0) red   = 0; else if (red   > 255) red   = 255;
     if (green < 0) green = 0; else if (green > 255) green = 255;
     if (blue  < 0) blue  = 0; else if (blue  > 255) blue  = 255;    

     int alpha = color.getAlpha();

     return new Color (red, green, blue, alpha);
   }
   

   private static String getMultiBlockJmolScript(AFPChain afpChain, Atom[] ca1, Atom[] ca2)
   {

      int blockNum = afpChain.getBlockNum();      
      int[] optLen = afpChain.getOptLen();
      int[][][] optAln = afpChain.getOptAln();

      if ( optLen == null)
         return DEFAULT_SCRIPT;

      StringWriter jmol = new StringWriter();
      jmol.append(DEFAULT_SCRIPT);

      jmol.append("select */2; color lightgrey; model 2; ");
      
      for(int bk = 0; bk < blockNum; bk ++)       {

         //the block nr determines the color...
         int colorPos = bk;
         if ( colorPos > AligPanel.colorWheel.length){
            colorPos = AligPanel.colorWheel.length % colorPos ;
         }
         
         Color c =  AligPanel.colorWheel[colorPos];
         Color cd = darker( AligPanel.colorWheel[colorPos],0.4);
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

         buf.append("; backbone 0.6 ; color [" + cd.getRed() +"," + cd.getGreen() +"," +cd.getBlue()+"]; select ");
         
         count = 0;
         for (String res :pdb2 ){
            if ( count > 0)
               buf.append(",");
         
            buf.append(res);
            buf.append("/2");
            count++;
         }
         //buf.append("; set display selected;");

         buf.append("; backbone 0.6 ; color [" + c.getRed() +"," + c.getGreen() +"," +c.getBlue()+"];");

         // now color this block:
         jmol.append(buf);
      }
      jmol.append("model 0;  ");
      jmol.append(LIGAND_DISPLAY_SCRIPT);
      System.out.println(jmol);
      return jmol.toString();


   }

   public void resetDisplay(){

      if  ( afpChain != null && ca1 != null && ca2 != null) {
         String script = getJmolString( afpChain,ca1,ca2);
         //System.out.println(j.toString());
         evalString(script);
      }
   }


}