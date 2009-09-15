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
 * Created on Aug 3, 2005
 *
 */
package org.biojava.dasobert.das;


import org.xml.sax.helpers.DefaultHandler;
import org.xml.sax.Attributes            ;
import java.util.*;
import java.util.logging.Logger;
import java.awt.Color;

/** a class to parse the XML response of a DAS - stylesheet request.
 * @author Andreas Prlic
 *
 */
public class DAS_Stylesheet_Handler extends DefaultHandler {
    
    List typeGlyphMaps;
    Map  currentType;
    String chars ;
    boolean threeDstyle;
    List threeDGlyphMaps;
    static Logger logger      = Logger.getLogger("org.biojava.spice");
    
    /**
     * 
     */
    public DAS_Stylesheet_Handler() {
        super();
        typeGlyphMaps = new ArrayList();
        currentType = new HashMap();
        threeDGlyphMaps = new ArrayList();
        threeDstyle = false;
    }
    
    
    public Map[] getTypeStyles(){
        return (Map[]) typeGlyphMaps.toArray(new Map[typeGlyphMaps.size()]);
    }
    
    public Map[] get3DStyles(){
        return (Map[]) threeDGlyphMaps.toArray(new Map[threeDGlyphMaps.size()]);
    }
    
    public void startElement (String uri, String name, String qName, Attributes atts){
        chars = "";
        
        if ( qName.equals("CATEGORY")){
            String id = atts.getValue("id");
            if ( id.equals("3D")){
                // here follow the 3D styles...
                threeDstyle = true;
            } else {
                threeDstyle = false;
                
            }
        }
        
        if ( qName.equals("TYPE")){
            // this glyph matches to features of type >id<.
            String id = atts.getValue("id");
            currentType = new HashMap(); 
            currentType.put("type",id);
        } 
        
        else if ( qName.equals("ARROW")){
            currentType.put("style","arrow");
        } else if ( qName.equals("ANCHORED_ARROW")){
            currentType.put("style","anchored_arrow");
        } else if ( qName.equals("BOX")){
            currentType.put("style","box");
        } else if ( qName.equals("CROSS")){
            currentType.put("style","cross");
        } else if ( qName.equals("EX")){
            currentType.put("style","ex");
        } else if ( qName.equals("HELIX")){
            currentType.put("style","helix");
        } else if ( qName.equals("LINE")){
            currentType.put("style","line");
        }  else if ( qName.equals("SPAN")){
            currentType.put("style","span");
        } else if ( qName.equals("TRIANGLE")){
            currentType.put("style","triangle");
        } else if ( qName.equals("GRADIENT")){
            currentType.put("style","gradient");
        }  else if ( qName.equals("HISTOGRAM")){
            currentType.put("style","histogram");
        }  else if ( qName.equals("LINEPLOT")){
            currentType.put("style","lineplot");
        }   else if ( qName.equals("TILING")){
            currentType.put("style","histogram");
        }
        
    }
    
    /**  convert the color provided by the stylesheet to a java Color 
     * 
     * @param chars
     * @return a Color
     */
    private Color getColorFromString(String chars){
        
        
        if (chars.equals("rotate")) {
            // use the default SPICE colors ...
            return null;
        }
        
        try {
            Color col = Color.decode(chars);
            return col;
        } catch ( Exception e) {
            logger.finest("could not decode color from stylesheet " + e.getMessage());
        }
        
        
        // map the string to a build in color...
        // thanks to the Xpm.java script provided by Phil Brown (under LGPL)
        // AP downloaded it from http://www.bolthole.com/java/Xpm.java
        
        // the DAS spec stylesheet only requires 16 VGA colors to be supported, but here we do more... :-)
        
        int[] rgb = Xpm.NameToRGB3(chars);
        if ( rgb != null) {
            Color col = new Color(rgb[0],rgb[1],rgb[2]);
            return col;
        }
        return null ;
    }
    
    public void endElement(String uri, String name, String qName) {
        if ( qName.equals("HEIGHT")){
            currentType.put("height",chars);
        } else if ( qName.equals("COLOR")){
            //System.out.println("got color " + chars);
            Color col = getColorFromString(chars);
            if ( col != null ){
                currentType.put("color",col);
            } else {
                if ( chars.equals("cpk")){
                    currentType.put("cpkcolor","true");
                }
            }
            
        } else if ( qName.equals("OUTLINECOLOR")){
            currentType.put("outlinecolor",chars);
        } else if ( qName.equals("BACKGROUND")){
            currentType.put("background",chars);
        } else if ( qName.equals("BUMP")){
            if ( chars.equals("no"))
                currentType.put("bump","no");
            else 
                currentType.put("bump","yes");
        
            // 3D stuff
        }  else if ( qName.equals("WIREFRAME")){
            currentType.put("display","wireframe");
        } else if ( qName.equals("SPACEFILL")){
            currentType.put("display","spacefill");
        } else if ( qName.equals("BACKBONE")){
            currentType.put("display","backbone");
        } else if ( qName.equals("CARTOON")){
            currentType.put("display","cartoon");
        } else if ( qName.equals("RIBBONS")){
            currentType.put("display","ribbons");            
        } else if ( qName.equals("WIDTH")){
            currentType.put("width",chars);
        } else if ( qName.equals("COLOR1")){
            currentType.put("color1",chars);
        } else if ( qName.equals("COLOR2")){
            currentType.put("color2",chars);
        } else if ( qName.equals("COLOR3")){
            currentType.put("color3",chars);
        } else if ( qName.equals("MAX")){
            currentType.put("max",chars);
        } else if ( qName.equals("MIN")){
            currentType.put("min",chars);
        } else if ( qName.equals("MIN")){
            currentType.put("min",chars);
        } else if ( qName.equals("STEPS")){
            currentType.put("steps",chars);
        } 
        
        else if ( qName.equals("TYPE")){           
            if ( threeDstyle){
             threeDGlyphMaps.add(currentType);   
            } else {
                typeGlyphMaps.add(currentType);
            }
        }
    }
    
    public void characters (char ch[], int start, int length){
        
     
            for (int i = start; i < start + length; i++) {
               chars += ch[i];
            }
        
        
    }
    
}




