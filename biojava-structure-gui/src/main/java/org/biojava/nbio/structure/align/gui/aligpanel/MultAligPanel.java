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
 */
package org.biojava.nbio.structure.align.gui.aligpanel;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.DisplayAFP;
import org.biojava.nbio.structure.align.gui.JPrintPanel;
import org.biojava.nbio.structure.align.gui.MenuCreator;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.AbstractAlignmentJmol;
import org.biojava.nbio.structure.align.gui.jmol.JmolTools;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.MultipleAlignment;
import org.biojava.nbio.structure.align.model.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.model.StructureAlignmentException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.align.webstart.WebStartMain;
import org.biojava.nbio.structure.align.xml.AFPChainXMLParser;
import org.biojava.nbio.structure.gui.events.AlignmentPositionListener;
import org.biojava.nbio.structure.gui.util.AlignedPosition;
import org.jcolorbrewer.ColorBrewer;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.List;


/** A JPanel that can display an AFPChain or a MultipleAlignment in a nice way and interact with Jmol.
 * 	It has been modified from the version specific for AFPChain to include the new MultipleAlignmentDS.
 * 	The AligPanel is initialized with a constructor rather than with the setters now.
 * 
 * @author Andreas Prlic
 * @author Aleix Lafita
 *
 */
public class MultAligPanel  extends JPrintPanel implements AlignmentPositionListener, WindowListener {

   private static final long serialVersionUID = -6892229111166263764L;

   private MultipleAlignment multAln;
   private int size; //number of structures
   private int length; //number of aligned positions
   private Color[] colors;
   
   private Font seqFont;
   private Font eqFont;
   private AbstractAlignmentJmol jmol;
   private MultAligPanelMouseMotionListener mouseMoLi;

   private MultAligmentCoordManager coordManager;
   private BitSet selection;

   private boolean selectionLocked;
   private boolean colorBySimilarity;
   private boolean colorByAlignmentBlock;

   private static final Color COLOR_EQUAL   = Color.decode("#6A93D4");
   private static final Color COLOR_SIMILAR = Color.decode("#D460CF");
   private static final Color[] DEFAULT_COLORS = ColorBrewer.Set1.getColorPalette(10);
   
   /**
    * Default constructor. Empty AligPanel instance.
    */
   public MultAligPanel(){
	      super();
	      this.setBackground(Color.white);
	      coordManager = new MultAligmentCoordManager();
	      seqFont = new Font("SansSerif",Font.PLAIN,12);
	      eqFont = new Font("SansSerif",Font.BOLD,12);

	      mouseMoLi = new MultAligPanelMouseMotionListener(this);
	      this.addMouseMotionListener(mouseMoLi);
	      this.addMouseListener(mouseMoLi);
	      mouseMoLi.addAligPosListener(this);

	      selection = new BitSet();
	      colorBySimilarity = false;
	      colorByAlignmentBlock = false;
	      
	      colors = DEFAULT_COLORS;
   }
   
   /**
    * Constructor using an afpChain and the atom arrays for pairwise alignments.
    * The AFPChain is converted into a MultipleAlignment.
    * 
    * @param afpChain 
    * @param ca1 
    * @param ca2 
    * @param color 
	* @throws StructureException 
	* @throws StructureAlignmentException 
    */
   public MultAligPanel(AFPChain afpChain, Atom[] ca1, Atom[] ca2, Color[] colors) throws StructureAlignmentException, StructureException{
	   this();
	   this.multAln = new MultipleAlignmentImpl(afpChain, ca1, ca2);
	   this.size = multAln.size();
	   this.multAln.updateAlnSequences();
	   this.length = multAln.getAlnSequences().get(0).length();
	   this.colors = colors;
	   if (colors == null) this.colors = DEFAULT_COLORS;
   }
   public MultAligPanel(AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws StructureAlignmentException, StructureException{
	   this(afpChain, ca1, ca2, DEFAULT_COLORS);
   }
   
   /**
    * Constructor using a MultipleAlignment for any other kind of alignment.
    * @param multAln
    * @param colors
 * @throws StructureAlignmentException 
    */
   public MultAligPanel(MultipleAlignment multAln, Color[] colors) throws StructureAlignmentException{
	   this();
	   this.multAln = multAln;
	   this.size = multAln.size();
	   this.length = multAln.getAlnSequences().get(0).length();
	   this.colors = colors;
	   if (colors == null) this.colors = DEFAULT_COLORS;
   }
   public MultAligPanel(MultipleAlignment multAln) throws StructureAlignmentException{
	   this(multAln,DEFAULT_COLORS);
   }
   
   public MultAligmentCoordManager getCoordManager() {
      return coordManager;
   }

   public void addAlignmentPositionListener(AlignmentPositionListener li){
      mouseMoLi.addAligPosListener(li);
   }

   public void destroy(){

      multAln = null;
      mouseMoLi.destroy();	
      jmol = null;
      selection = null;
   }

   public MultipleAlignment getMultipleAlignment(){
      return multAln;
   }

	@Override
	public void paintComponent(Graphics g){

      super.paintComponent(g);

      Graphics2D g2D = (Graphics2D) g;


      g2D.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,
            RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

      g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
            RenderingHints.VALUE_ANTIALIAS_ON);


      // only draw within the ranges of the Clip
      //Rectangle drawHere = g2D.getClipBounds();        


      //int startpos = coordManager.getSeqPos(0,drawHere.x);       
      //int endpos   = coordManager.getSeqPos(0,drawHere.x+drawHere.width-2);

      int startpos = 0;
      int endpos = multAln.getAlnSequences().get(0).length();

      String summary = multAln.toString();
      g2D.drawString(summary, 20, coordManager.getSummaryPos());

      Color significantCol = Color.red;
      //if (multAln.isSignificantResult()) significantCol = Color.green;

      g2D.setPaint(significantCol);
      //draw a darker background
      Rectangle sig = new Rectangle(10,10,10,10);				
      g2D.fill(sig);
      boolean isFATCAT = false;
      if (multAln.getAlgorithmName().startsWith("jFatCat")) isFATCAT = true;
      
      List<Integer> alignedPos = new ArrayList<Integer>();
		try {
			alignedPos = DisplayAFP.getCoreAlignmentPos(multAln);
		} catch (StructureAlignmentException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
      
      for (int i = startpos; i < endpos; i++){

         // color amino acids by hydrophobicity
         boolean isGapped = false;
         g2D.setFont(seqFont);
         
        if (alignedPos.contains(i)) g2D.setFont(eqFont);
        else isGapped = true;

        //Loop through every structure to get all the points
        List<Point> points = new ArrayList<Point>();
        for (int str=0; str<size; str++){
        
        	char c = multAln.getAlnSequences().get(str).charAt(i);
        	if (str==0) points.add(coordManager.getPanelPos(str,i));
        	Point p1 = points.get(points.size()-1);
        	Point p2 = coordManager.getPanelPos(str,i);
	        points.add(p2);
	         
	         if (! isGapped) {
	            Color bg = colors[str];
	            //Color end1 = ColorUtils.rotateHue(colors[str],  (1.0f  / 24.0f) * blockNum  );
	            //Color end2 = ColorUtils.rotateHue(colors[str],    (1.0f  / 24.0f) * (blockNum  +1)) ;
	            
	            //TODO enable coring blocks by rotation of the colors
	            //int blockNum = afpChain.getBlockNum();
	            /*if ( colorByAlignmentBlock) {
	
	               if (! alignedPos.contains(i)){
	
	                  // unaligned!
	                  bg = Color.white;
	                  bg2 = Color.white;
	               } else  {
	                  
	                  int colorPos = 0;
	                  if (isFATCAT) {
	                     int block = 0;
	                     char s = symb[i];
	                     try {
	                        block = Integer.parseInt(s+"") - 1;
	                        bg  = ColorUtils.getIntermediate(ColorUtils.orange, end1, blockNum, block);
	                        bg2   = ColorUtils.getIntermediate(ColorUtils.cyan, end2, blockNum, block);
	                        //bg = ColorUtils.rotateHue(ColorUtils.orange,  (1.0f  / 24.0f) * block  );
	                        //bg2 = ColorUtils.rotateHue(ColorUtils.cyan,  (1.0f  / 16.0f) * block );
	                     } catch (Exception e){}
	                     
	                     if ( colorPos > ColorUtils.colorWheel.length){
	                        colorPos = ColorUtils.colorWheel.length % colorPos ;
	                     }
	                  } else {
	                	  colorPos = AFPAlignmentDisplay.getBlockNrForAlignPos(afpChain, i);
	            		  bg  = ColorUtils.getIntermediate(ColorUtils.orange, end1, blockNum, colorPos);
	            		  bg2   = ColorUtils.getIntermediate(ColorUtils.cyan, end2, blockNum, colorPos);
	            		  //bg = ColorUtils.rotateHue(ColorUtils.orange,  (1.0f  / 24.0f) * colorPos );
	            		  //bg2 = ColorUtils.rotateHue(ColorUtils.cyan,  (1.0f  / 16.0f) * colorPos);
	                  }
	                  
	               }
	            } else {*/
	
	               bg = Color.LIGHT_GRAY;
	            //}
	            
	            // draw a darker background
	            g2D.setPaint(bg);
	            Rectangle rec = new Rectangle(p1.x-1,p1.y-11, (p2.x-p1.x)+12, (p2.y-p1.y)+1);
	            g2D.fill(rec);
	
	            //g2D.setPaint(Color.black);
	            //g2D.draw(rec);
	         }
	         /*if ( colorBySimilarity){
	            if ( c1 == c2){
	               Color bg = COLOR_EQUAL;
	               g2D.setPaint(bg);
	               Rectangle rec = new Rectangle(p1.x-1,p1.y-11, (p2.x-p1.x)+12, (p2.y-p1.y)+12);				
	               g2D.fill(rec);
	            } else if (AFPAlignmentDisplay.aaScore(c1, c2) > 0) {
	               Color bg = COLOR_SIMILAR;
	               g2D.setPaint(bg);
	               Rectangle rec = new Rectangle(p1.x-1,p1.y-11, (p2.x-p1.x)+12, (p2.y-p1.y)+12);				
	               g2D.fill(rec);
	            }
	         }*/
	
	         //if ( selectionStart != null && selectionEnd != null){
	         //	if ( i >= selectionStart.getPos1() && i <= selectionEnd.getPos1()) {
	
	         if ( isSelected(i)){
	            // draw selection
	            Color bg = Color.YELLOW;
	            g2D.setPaint(bg);
	            // draw a darker backgroun
	            Rectangle rec = new Rectangle(p1.x-1,p1.y-11, (p2.x-p1.x)+12, (p2.y-p1.y)+12);
	            g2D.fill(rec);
	            //	}
	         }
	
	         // draw the AA sequence
	         g2D.setColor(Color.black);
	         g2D.drawString(c+"",p2.x,p2.y);
	         
	         //System.out.println(seq1[i] + " " + xpos1 + " " + ypos1 + " " + seq2[i] + xpos2 + " " + ypos2);
        }
      }

      int nrLines = length / MultAligmentCoordManager.DEFAULT_LINE_LENGTH;


      for (int i = 0 ; i < nrLines ; i++){

         try {
            // draw legend at i
        	for (int str=0; str<size; str++){
	            Point p1 = coordManager.getLegendPosition(i,str);
	            
	            int aligPos = i * MultAligmentCoordManager.DEFAULT_LINE_LENGTH;
	            Atom a1 = DisplayAFP.getAtomForAligPos(multAln, str, aligPos);
	            String label1 = JmolTools.getPdbInfo(a1,false);				
	            g2D.drawString(label1, p1.x,p1.y);
	
	            Point p3 = coordManager.getEndLegendPosition(i,str);
	
	            aligPos = i * MultAligmentCoordManager.DEFAULT_LINE_LENGTH + MultAligmentCoordManager.DEFAULT_LINE_LENGTH -1 ;
	            if (aligPos > length) aligPos = length-1;
	            Atom a3 = DisplayAFP.getAtomForAligPos(multAln, str, aligPos);
	
	            String label3 = JmolTools.getPdbInfo(a3,false);
	
	            g2D.drawString(label3, p3.x,p3.y);

        	}
         } catch (StructureAlignmentException e) {
			e.printStackTrace();
		}
      }
   }

 


   private boolean isSelected(int alignmentPosition) {

      return selection.get(alignmentPosition);

   }


   @Override
public void mouseOverPosition(AlignedPosition p) {
      //System.out.println("AligPanel: mouse over position " + p.getPos1() );

      if ( ! selectionLocked)
         selection.clear();
      selection.set(p.getPos1());

      updateJmolDisplay();

      this.repaint();

   }

   private void updateJmolDisplay() {

      if ( jmol == null) return;

      StringBuffer cmd = new StringBuffer("select ");

      int nrSelected = 0;
      try {

         for (int i = 0 ; i< length ; i++){
            if ( selection.get(i)){
            	
            	List<String> select = new ArrayList<String>();
            	Collections.fill(select, "");
            	
            	for (int str=0; str<size; str++){
	               Atom a1 = DisplayAFP.getAtomForAligPos(multAln,str,i);
	               if ( a1 != null ) select.get(str).concat(JmolTools.getPdbInfo(a1));
	               if ( nrSelected > 0) cmd.append(", ");
	                cmd.append(select);
	                cmd.append("/"+(str+1)+", ");
            	}
               nrSelected++;
        	}
        }

      } catch (StructureAlignmentException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
      }	
      if (nrSelected == 0)
         cmd.append(" none;");
      else
         cmd.append("; set display selected;");
      
      jmol.evalString(cmd.toString());
   }


   @Override
public void positionSelected(AlignedPosition p) {
      mouseOverPosition(p);

   }

   @Override
public void rangeSelected(AlignedPosition start, AlignedPosition end) {
      //System.out.println("AligPanel: range selected " + start.getPos1() + " - " + end.getPos1() + " selectionLockedL " + selectionLocked);
      if ( ! selectionLocked )
         selection.clear();
      selection.set(start.getPos1(), end.getPos1()+1);
      updateJmolDisplay();
      this.repaint();

   }

   @Override
public void selectionLocked() {
      selectionLocked = true;

   }

   @Override
public void selectionUnlocked() {
      selectionLocked = false;
      selection.clear();
      this.repaint();

   }


   @Override
public void toggleSelection(AlignedPosition p) {
      selection.flip(p.getPos1());
      //System.out.println("AligPanel: toggle selection " + p.getPos1() + " " + selection.get(p.getPos1()));
      updateJmolDisplay();
      this.repaint();

   }



   public void setStructureAlignmentJmol(AbstractAlignmentJmol jmol) {
      this.jmol = jmol;

   }


   @Override
public void windowActivated(WindowEvent e) {

      // TODO Auto-generated method stub

   }


   @Override
public void windowClosed(WindowEvent e) {
      // TODO Auto-generated method stub

   }


   @Override
public void windowClosing(WindowEvent e) {
      destroy();

   }


   @Override
public void windowDeactivated(WindowEvent e) {
      // TODO Auto-generated method stub

   }


   @Override
public void windowDeiconified(WindowEvent e) {
      // TODO Auto-generated method stub

   }


   @Override
public void windowIconified(WindowEvent e) {
      // TODO Auto-generated method stub

   }


   @Override
public void windowOpened(WindowEvent e) {
      // TODO Auto-generated method stub

   }

   @Override
public void actionPerformed(ActionEvent e) {
      String cmd = e.getActionCommand();
      // print is handled by superclass
      if ( cmd.equals(MenuCreator.PRINT)) {
         super.actionPerformed(e);
      /*} else if (cmd.equals(MenuCreator.TEXT_ONLY)){
         String result = AfpChainWriter.toWebSiteDisplay(afpChain, ca1, ca2);
         DisplayAFP.showAlignmentImage(afpChain, result);
      } else if ( cmd.equals(MenuCreator.PAIRS_ONLY)) {
         String result = AfpChainWriter.toAlignedPairs(afpChain, ca1, ca2) ;
         DisplayAFP.showAlignmentImage(afpChain, result);
      } else if (cmd.equals(MenuCreator.FATCAT_TEXT)){
         String result = afpChain.toFatcat(ca1, ca2);
         result += AFPChain.newline;
         result += afpChain.toRotMat();
         DisplayAFP.showAlignmentImage(afpChain, result);
      } else if ( cmd.equals(MenuCreator.SELECT_EQR)){
         selectEQR();
      } else if ( cmd.equals(MenuCreator.SIMILARITY_COLOR)){
         colorBySimilarity(true);
      } else if ( cmd.equals(MenuCreator.EQR_COLOR)){
         colorBySimilarity(false);
      } else if ( cmd.equals(MenuCreator.FATCAT_BLOCK)){
         colorByAlignmentBlock();*/
      } else {
         System.err.println("Unknown command:" + cmd);
      }

   }


   private void colorByAlignmentBlock() {
      colorByAlignmentBlock = true;
      colorBySimilarity = false;
      this.repaint();

   }


   private void colorBySimilarity(boolean flag) {
      this.colorBySimilarity = flag;
      colorByAlignmentBlock = false;
      this.repaint();
   }


   private void selectEQR() throws StructureAlignmentException {

      selection.clear();

      List<Integer> pos1 = DisplayAFP.getCoreAlignmentPos(multAln);

      for (int pos : pos1){
         selection.flip(pos);
      }
      mouseMoLi.triggerSelectionLocked(true);
      updateJmolDisplay();
      this.repaint();
   }

   public List<Atom[]> getAtomArrays() throws StructureAlignmentException {
      return multAln.getParent().getAtomArrays();
   }

   public static void main(String[] args){

      String file = "/Users/ap3/tmp/4hhb.ce";

      try {
         BufferedReader in = new BufferedReader(new FileReader(file));
         StringBuffer xml = new StringBuffer();
         String str;
         while ((str = in.readLine()) != null) {
            xml.append(str);
         }
         in.close();

         AFPChain[] afps = AFPChainXMLParser.parseMultiXML(xml.toString());
         AFPChain afpChain = afps[0];

         UserConfiguration config = WebStartMain.getWebStartConfig();
         AtomCache cache = new AtomCache(config.getPdbFilePath(),config.getCacheFilePath());

         Atom[] ca1 = cache.getAtoms(afpChain.getName1());
         Atom[] ca2 = cache.getAtoms(afpChain.getName2());

         AFPChainXMLParser.rebuildAFPChain(afpChain, ca1, ca2);


         //StructureAlignment algorithm = StructureAlignmentFactory.getAlgorithm(afpChain.getAlgorithmName());
         StructureAlignmentJmol jmol= StructureAlignmentDisplay.display(afpChain, ca1, ca2);

         DisplayAFP.showAlignmentImage(afpChain, ca1, ca2, jmol);

      } catch (Exception e){
         e.printStackTrace();
      }
   }

}


