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
 * Created on Jun 30, 2010
 * Author: ap3 
 *
 */

package org.biojava.bio.structure.gui.util;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.util.logging.Logger;

import javax.swing.Box;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;

import org.biojava.bio.structure.align.gui.autosuggest.AutoSuggestProvider;
import org.biojava.bio.structure.align.gui.autosuggest.JAutoSuggest;
import org.biojava.bio.structure.align.gui.autosuggest.SCOPAutoSuggestProvider;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.align.webstart.WebStartMain;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;

public class ScopSelectPanel 
extends JPanel
implements StructurePairSelector
{

   /**
    * 
    */
   private static final long serialVersionUID = 757947454156959178L;
   JAutoSuggest dom1;
   JAutoSuggest dom2;

   
   public static Logger logger =  Logger.getLogger("org.biojava");
   
   public ScopSelectPanel(){
      
     this(true);
   }
   
   public ScopSelectPanel(boolean show2boxes){
	   Box vBox = Box.createVerticalBox();
	      
	      //dom1 = new JTextField(10);
	      //dom2 = new JTextField(10);
	      
	      AutoSuggestProvider autoSuggesP = new SCOPAutoSuggestProvider();
			
	      dom1 = new JAutoSuggest(10);		
	      dom1.setAutoSuggestProvider(autoSuggesP);
	      
	      dom2 = new JAutoSuggest(10);
	      dom2.setAutoSuggestProvider(autoSuggesP);
	            
	      Box b1 = getDomainPanel(1,dom1);
	      Box b2 = getDomainPanel(2,dom2);
	      
	      
	      vBox.add(b1);
	      if ( show2boxes)
	    	  vBox.add(b2);
	      
	      this.add(vBox);
   }
   
   private Box getDomainPanel(int pos ,JTextField f){

      //JPanel panel = new JPanel();
      //panel.setBorder(BorderFactory.createLineBorder(Color.black));

      JLabel l01 = new JLabel("SCOP or domain id:");

      //panel.add(l01);
      Box hBox = Box.createHorizontalBox();
      hBox.add(Box.createGlue());
      hBox.add(l01);

      JLabel l11 = new JLabel(pos + ":");
      f.setMaximumSize(new Dimension(Short.MAX_VALUE,30));
      f.setToolTipText("Provide SCOP ID here. Example: d1zyma1");
      hBox.add(l11);
      hBox.add(Box.createVerticalGlue());
      hBox.add(f, BorderLayout.CENTER);
      hBox.add(Box.createGlue());
      

      //hBox21.add(Box.createGlue());

      //panel.add(hBox21);



      return hBox;
  }
   
   public Structure getStructure1() throws StructureException
   {
	   String scop1 = dom1.getText();
      return getStructure(scop1);
   }

   public Structure getStructure2() throws StructureException
   {
      return getStructure(dom2.getText());
   }
   
   private Structure getStructure(String domainID) throws StructureException{
      //PDBFileReader reader = new PDBFileReader();
      
	   if ( domainID == null || domainID.equals(""))
		   return null;
      
      
      
      UserConfiguration config = WebStartMain.getWebStartConfig();
      //String cacheLocation = config.getPdbFilePath();
      
      AtomCache cache = new AtomCache(config);
      
     Structure s = null;
     try {
     	s =	cache.getStructure(domainID);
     	s.setName(domainID);
     } catch (Exception e){
    	 e.printStackTrace();
     }
     return s;
      
//      boolean isSplit = config.isSplit();
//      
//      AtomCache cache = new AtomCache(cacheLocation,isSplit );
//      
//      ScopDatabase scop = ScopInstallationInstance.getInstance().getSCOP();
//      
//      ScopDomain domain = scop.getDomainByScopID(domainID) ;
//      
//      System.out.println("found scop domain :" + domain);
//      
//      if ( domain == null)
//         return null;
//      
//      
//      Structure s = null;
//      try {
//         s =cache.getStructureForDomain(domain);
//         if ( s.getName() == null || s.getName().equals(""))
//        	 s.setName(domainID);
//         s.setPDBCode(domainID);
//      } catch (Exception e){
//         e.printStackTrace();
//         logger.warning(e.getMessage());
//      }
//      
//      return s;

  }

}
