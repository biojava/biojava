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
package org.biojava.bio.seq.io.agave;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.Annotation;

/**
 * Basic implementation of AGAVEAnnotFilter
 * @author Hanning Ni    Doubletwist Inc
*/
public class SimpleAnnotFilter implements AGAVEAnnotFilter{

  public static final AGAVEAnnotFilterFactory SIMPLE_ANNOT_FILTER_FACTORY
    = new AGAVEAnnotFilterFactory() {
    public AGAVEAnnotFilter getInstance() {
      return new SimpleAnnotFilter();
    }
  };

  SimpleAnnotFilter() {
  }

   /**
     *
     */
   public String getAccession(Annotation annot)
   {
       return (String) null ;
   }

        /**
         *
         */
    public String getKeyword(Annotation annot)
    {
        return (String) null;
    }
    public String getElementId(Annotation annot)
    {
        return (String) null ;
    }
        /**
         *
         */
    public String getOrganism(Annotation annot)
    {
        return (String) null ;
    }

    public String getLabel(Annotation annot)
    {
        return (String) null;
    }

        /**
         *
         */
    public String getDescription(Annotation annot)
    {
        return (String) null;
    }

    public String getMatchAlign(Annotation annot)
    {
        return (String) null ;
    }

    public AGAVEMatchRegion getMatchRegion(Annotation annot)
    {
        return (AGAVEMatchRegion) null ;
    }

    public AGAVEQueryRegion getQueryRegion(Annotation annot)
    {
        return (AGAVEQueryRegion) null ;
    }

    public String getClassifySystem(Annotation annot)
    {
        return (String) null ;
    }

    public String getClassifyId(Annotation annot)
    {
        return (String) null ;
    }

    public String getClassifyType(Annotation annot)
    {
        return (String) null ;
    }

    public String[] getExonIds(Annotation annot)
    {
        return (String[])null ;
    }

    public String getChromNum(Annotation annot)
    {
        return (String)null ;
    }

    public AGAVEIdAlias[] getIdAlias(Annotation annot)
    {
        return (AGAVEIdAlias[])null ;
    }
        /**
         *
         */
    public String getNote(Annotation annot)
    {
        return (String) null ;
    }

    public AGAVEDbId[] getAltIds(Annotation annot)
    {
       return (AGAVEDbId[])null ;
    }

    public AGAVEMapLocation[]  getMapLocation(Annotation annot)
    {
       return (AGAVEMapLocation[]) null  ;
    }

    public AGAVERelatedAnnot[]  getRelatedAnnot(Annotation annot)
    {
       return (AGAVERelatedAnnot[]) null ;
    }

    public String[]  getElementIds(Annotation annot)
    {
       return (String[]) null ;
    }

     public String getGroupOrder(Annotation annot)
     {
         return (String) null ;
     }

     public String getMatchDesc(Annotation annot)
     {
         return (String) null ;
     }

     public String getFeatureType(Annotation annot)
     {
         return (String) null ;
     }
     public String getResultType(Annotation annot)
     {
         return (String) null ;
     }

     public String getConfidence(Annotation annot)
     {
         return (String) null ;
     }

     public String getAlignLength(Annotation annot)
     {
         return (String) null ;
     }

     public String getAlignUnits(Annotation annot)
     {
         return (String) null ;
     }


    public AGAVEXrefs[] getXrefs(Annotation annot)
    {
       return (AGAVEXrefs[])null ;
    }
        /**
         *
         */
     public String getVersion(Annotation annot)
     {
         return (String) null ;
     }

     public String getSequenceId(Annotation annot)
     {
         return (String) null ;
     }

    public String getTaxonId(Annotation annot)
    {
         return (String) null ;
     }

     public String getCloneId(Annotation annot)
     {
         return (String) null ;
     }

     public String getCloneLibrary(Annotation annot)
     {
         return (String) null ;
     }

     public String getChromosome(Annotation annot)
     {
         return (String) null ;
     }

     public String getMapPosition(Annotation annot)
     {
         return (String) null ;
     }

     public String getEcNumber(Annotation annot)
     {
         return (String) null ;
     }

      public String getCreateDate(Annotation annot)
     {
         return (String) null ;
     }

     public String getUpdateDate(Annotation annot)
     {
         return (String) null ;
     }


        /**
         *
         */
     public String getOS(Annotation annot)
     {
         return (String) null ;
     }

        /**
         *
         */
     public String getMolType(Annotation annot)
     {
           return (String) null;
     }

     public AGAVEDbId getDbId(Annotation annot)
     {
          return (AGAVEDbId) null ;
     }
        /**
         *   ThomasD made this a bit safer...
         */
     public AGAVEProperty[] getProperty(Annotation annot, String type)
     {
	 List set = new ArrayList();
	 for (Iterator i = annot.keys().iterator(); i.hasNext();)
	     {
		 Object key = i.next();
		 Object value = annot.getProperty(key);
		 if (key instanceof String && value instanceof String) {
		     set.add(new AGAVEProperty(type, (String) key, null, (String) value ));
		 }
	     }
	 return (AGAVEProperty[]) set.toArray(new AGAVEProperty[set.size()]);
     }
}
