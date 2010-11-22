package org.biojava.bio.structure.scop;

import java.io.StringWriter;

/** Contains data from
 * dir.des.scop.txt_1.75
 * 
 * e.g 
 * <pre>
 * 26154	px	b.47.1.2	d1nrs.1	1nrs L:,H:
 * </pre>
 * <pre>
 * 125030	px	b.47.1.2	d1zgia1	1zgi A:1A-245
 * </pre>
 * 
 * @author Andreas Prlic
 *
 */
public class ScopDescription {

   int sunID;
   ScopCategory category;
   String classificationId;
   String name;
   String description;


   public String toString(){
      StringWriter buf = new StringWriter();
      
      buf.append(sunID+"");
      buf.append("\t");
      buf.append(category.toString());
      buf.append("\t");
      buf.append(classificationId);
      buf.append("\t");
      buf.append(name);
      buf.append("\t");
      buf.append(description);
    
      return buf.toString();
   }
   

   public int getSunID()
   {
      return sunID;
   }
   public void setSunID(int sunID)
   {
      this.sunID = sunID;
   }
   public ScopCategory getCategory()
   {
      return category;
   }
   public void setCategory(ScopCategory category)
   {
      this.category = category;
   }
   public String getClassificationId()
   {
      return classificationId;
   }
   public void setClassificationId(String classificationId)
   {
      this.classificationId = classificationId;
   }
   public String getName()
   {
      return name;
   }
   public void setName(String name)
   {
      this.name = name;
   }
   public String getDescription()
   {
      return description;
   }
   public void setDescription(String description)
   {
      this.description = description;
   }



}
