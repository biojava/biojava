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
 * Created on May 20, 2010
 * Author: Andreas Prlic 
 *
 */

package org.biojava.nbio.structure.gui.util.color;

import java.awt.*;

public class ColorUtils
{
   
   public static Color orange =   Color.decode("#FFA500");
   public static Color cyan   =   Color.decode("#00FFFF");
   public static Color gold   =   Color.decode("#FFD700");
  
   
   static final Color c1 = Color.decode("#228B22"); // green
   static final Color c2 = Color.decode("#8F8FFF"); // blue   
   static final Color c3 = gold;
   static final Color c4 = Color.decode("#FF8C00"); // orange
   static final Color c5 = Color.decode("#FF00FF"); // pink
   static final Color c6 = Color.decode("#C71585"); // red
   
   public static final Color[] colorWheel = new Color[] {c1, c2, c3, c4 , c5,c6}; 

   
   public static void main(String[] args){
      int i = -1;
      for ( Color color : colorWheel){
         i++;
         float af[] = Color.RGBtoHSB(color.getRed(),color.getGreen(),color.getBlue(), null);
         System.out.println("position "  + i + "  " + af[0] + " " + af[1] + " " + af[2]);
         //System.out.println(rotateHue(color, 0.1f));
      }
   }
   
   public static Color rotateHue (Color color, float fraction) {
      
      float af[] = Color.RGBtoHSB(color.getRed(),color.getGreen(),color.getBlue(), null);
      
      float hue = af[0];
      float saturation = af[1];
      float brightness = af[2];
      
      //System.out.println(hue + " " + saturation + " " + brightness);
      
      float hueNew = hue  + fraction;
      
      //System.out.println(hue + " " + hueNew);
      
      if ( hueNew > 1 ){
         hueNew = hueNew - 1;
      }
     
      return Color.getHSBColor(hueNew, saturation, brightness);
   }
   
   public static Color getIntermediate(Color start, Color end, int stepSize, int position ){
      
      float af1[] = Color.RGBtoHSB(start.getRed(),start.getGreen(),start.getBlue(), null);
      
      float af2[] = Color.RGBtoHSB(end.getRed(),end.getGreen(),end.getBlue(), null);
      
      float hue1 = af1[0];
      float hue2 = af2[0];
     
      if ( hue2 < hue1) {
        
         hue2 = af1[0];
         hue1 = af2[0];
      }
      
      //float saturation = af1[1] + af2[1] / 2f;
      //float brightness = af1[2] + af2[2] / 2f;
      
      float range = Math.abs(hue2-hue1);
      
      while ( position > stepSize){
         position = position - stepSize ; 
      }
      
      float inc = (range * position / (float) stepSize) ;
      float hueNew = hue1 + inc;
      
      //System.out.println(position + " " + hue1 + " " + hue2 + " new: " + hueNew + " inc " + inc + " range " + range); 
      return Color.getHSBColor(hueNew, af1[1], af1[2]);
           
   }
   
   /**
    * Make a color darker. (RGB color scheme)
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
    * Make a color lighter. (RGB color scheme)
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

}
