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
 * Created on Jun 7, 2010
 * Author: ap3 
 *
 */

package org.biojava.bio.structure.align.seq;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;

public class SmithWaterman3DParameters implements ConfigStrucAligParams
{

   short match ;     // match
   short replace ;      // replace 
   short insert ;     // insert
   short delete ;      // delete
   short gapExtend ;      // gapExtend
   
   public SmithWaterman3DParameters(){
      reset();
   }
   
  
   
   public List<String> getUserConfigHelp()
   {
      List<String> params =new ArrayList<String>();
      params.add("The Gap Extension score");
      params.add("The Insert score");
      params.add("The Delete score");
      params.add("The Match score");
      params.add("The Replace score");
      
      // TODO Auto-generated method stub
      return params;
   }

   public List<String> getUserConfigParameterNames()
   {
      List<String> params =new ArrayList<String>();
      params.add("Gap Extension");
      params.add("Insert");
      params.add("Delete");
      params.add("Match");
      params.add("Replace");
      
      return params;
   }

   public List<String> getUserConfigParameters()
   {
      List<String> params =new ArrayList<String>();
      params.add("GapExtend");      
      params.add("Insert");
      params.add("Delete");
      params.add("Match");
      params.add("Replace");
      return params;
   }

   public List<Class> getUserConfigTypes()
   {
      List<Class> params = new ArrayList<Class>();
      params.add(Short.class);
      params.add(Short.class);
      params.add(Short.class);
      params.add(Short.class);
      params.add(Short.class);
      return params;
   }

   public void reset()
   {
       match = (short) -1;     // match
       replace = (short) 3;      // replace 
       insert = (short) 2;     // insert
       delete = (short) 2;      // delete
       gapExtend = (short) 1;      // gapExtend
      
   }



   public Short getMatch()
   {
      return match;
   }



   public void setMatch(Short match)
   {
      this.match = match;
   }



   public Short getReplace()
   {
      return replace;
   }



   public void setReplace(Short replace)
   {
      this.replace = replace;
   }



   public Short getInsert()
   {
      return insert;
   }



   public void setInsert(Short insert)
   {
      this.insert = insert;
   }



   public Short getDelete()
   {
      return delete;
   }



   public void setDelete(Short delete)
   {
      this.delete = delete;
   }



   public Short getGapExtend()
   {
      return gapExtend;
   }



   public void setGapExtend(Short gapExtend)
   {
      this.gapExtend = gapExtend;
   }
   
   

}
