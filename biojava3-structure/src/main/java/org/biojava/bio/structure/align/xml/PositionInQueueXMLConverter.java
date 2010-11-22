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
 * Created on May 10, 2010
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.xml;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;


import org.biojava3.core.util.PrettyXMLWriter;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;

public class PositionInQueueXMLConverter
{

   public String toXML(int position) throws IOException{
      StringWriter swriter = new StringWriter();

      PrintWriter writer = new PrintWriter(swriter);
      PrettyXMLWriter xml = new PrettyXMLWriter(writer);

      xml.openTag("queue");
      xml.attribute("position", position+"");
      xml.closeTag("queue");
      xml.close();
      return swriter.toString(); 
   }
   
   public int fromXML(String xml){
      int position = Integer.MIN_VALUE;
      
      try
      {
         //Convert string to XML document
         DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
         DocumentBuilder db = factory.newDocumentBuilder();
         InputSource inStream = new InputSource();
         inStream.setCharacterStream(new StringReader(xml));
         Document doc = db.parse(inStream);

         // normalize text representation
         doc.getDocumentElement().normalize();


         //Element rootElement = doc.getDocumentElement();

         NodeList listOfAlignments = doc.getElementsByTagName("queue");
         //int numArrays = listOfAlignments.getLength();    
         //System.out.println("got " + numArrays + " alignment results.");
         // go over the blocks
         
       
         for(int afpPos=0; afpPos<listOfAlignments.getLength() ; afpPos++)
         {
            
            Node rootElement       = listOfAlignments.item(afpPos);

            String pos = getAttribute(rootElement,"position");
          
            try {
               position = Integer.parseInt(pos); 
            } catch (NumberFormatException f){
               f.printStackTrace();
            }

         }
      } 
      catch (SAXParseException err) 
      {
         System.out.println ("** Parsing error" + ", line " 
               + err.getLineNumber () + ", uri " + err.getSystemId ());
         System.out.println(" " + err.getMessage ());
      }
      catch (SAXException e)
      {
         Exception x = e.getException ();
         ((x == null) ? e : x).printStackTrace ();
      }
      catch (Throwable t)
      {
         t.printStackTrace ();
      }

      return position;
   }
   

   private static String getAttribute(Node node, String attr){
      if( ! node.hasAttributes()) 
         return null;

      NamedNodeMap atts = node.getAttributes();

      if ( atts == null)
         return null;

      Node att = atts.getNamedItem(attr);
      if ( att == null)
         return null;

      String value = att.getTextContent();

      return value;

   }
}
