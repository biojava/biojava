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
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.Annotation;

/**
 * Dumping the data from biojava with source of agave into agave format
 * @author Hanning Ni     Doubletwist Inc
 *
 */
public class Agave2AgaveAnnotFilter extends SimpleAnnotFilter{

  public static final AGAVEAnnotFilterFactory AGAVE_AGAVE_ANNOT_FILTER_FACTORY
    = new AGAVEAnnotFilterFactory() {
    public AGAVEAnnotFilter getInstance() {
      return new Agave2AgaveAnnotFilter();
    }
  };

  Agave2AgaveAnnotFilter() {
  }

    /**     */
   public String getAccession(Annotation annot)
   {
       if( annot == null ) return (String)null ;
       return (String) UtilHelper.getProperty(annot, "accession");
   }

   /**       */
    public String getKeyword(Annotation annot)
    {
        if( annot == null ) return (String)null ;
        return (String) UtilHelper.getProperty(annot,"keyword");
    }

    /**      */
    public String getOrganism(Annotation annot)
    {
        if( annot == null ) return (String)null ;
        return (String) UtilHelper.getProperty(annot,"organism_name");
    }

    public String getElementId(Annotation annot)
    {
        if( annot == null ) return (String)null ;
         return (String) UtilHelper.getProperty(annot,"element_id");
    }

    public String getLabel(Annotation annot)
    {
        if( annot == null ) return (String)null ;
         return (String) UtilHelper.getProperty(annot,"label");
    }

    /**      */
    public String getDescription(Annotation annot)
    {
        if( annot == null ) return (String)null ;
         return (String) UtilHelper.getProperty(annot,"description");
    }

    /**        */
    public String getNote(Annotation annot)
    {
        if( annot == null ) return (String)null ;
         return (String) UtilHelper.getProperty(annot,"note");
    }

    /**         */
     public String getVersion(Annotation annot)
     {
        if( annot == null ) return (String)null ;
         return (String) UtilHelper.getProperty(annot,"version");
     }

     /**        */
     public String getOS(Annotation annot)
     {
         if( annot == null ) return (String)null ;
         return (String) UtilHelper.getProperty(annot,"organism_name");
     }

     /**        */
     public String getMolType(Annotation annot)
     {
          if( annot == null ) return (String)null ;
           return (String) UtilHelper.getProperty(annot,"molecule_type");
     }

     public String getTaxonId(Annotation annot)
     {
         if( annot == null ) return (String)null ;
          return (String) UtilHelper.getProperty(annot,"taxon_id");
     }

     public String getCloneId(Annotation annot)
     {
         if( annot == null ) return (String)null ;
          return (String) UtilHelper.getProperty(annot,"clone_id");
     }

     public String getCloneLibrary(Annotation annot)
     {
         if( annot == null ) return (String)null ;
          return (String) UtilHelper.getProperty(annot,"clone_library");
     }

     public String getChromosome(Annotation annot)
     {
         if( annot == null ) return (String)null ;
          return (String) UtilHelper.getProperty(annot,"chromosome");
     }

     public String getMapPosition(Annotation annot)
     {
         if( annot == null ) return (String)null ;
          return (String) UtilHelper.getProperty(annot,"map_position");
     }

     public String getEcNumber(Annotation annot)
     {
         if( annot == null ) return (String)null ;
          return (String) UtilHelper.getProperty(annot,"ec_number");
     }

     public String getCreateDate(Annotation annot)
     {
         if( annot == null ) return (String)null ;
          return (String) UtilHelper.getProperty(annot,"create_date");
     }

     public String getUpdateDate(Annotation annot)
     {
         if( annot == null ) return (String)null ;
          return (String) UtilHelper.getProperty(annot,"update_date");
     }

     public AGAVEXrefs[] getXrefs(Annotation annot)
    {
        if( annot == null ) return (AGAVEXrefs[])null ;
        Object ob = UtilHelper.getProperty(annot,"xrefs") ;
        if( ob != null && ob instanceof List)
        {
            AGAVEXrefs[] set = new AGAVEXrefs[1];
            return (AGAVEXrefs[])((List)ob).toArray( set ) ;
        }
        return (AGAVEXrefs[]) null ;
    }

    public AGAVERelatedAnnot[]  getRelatedAnnot(Annotation annot)
    {
        if( annot == null ) return (AGAVERelatedAnnot[])null ;
        Object ob = UtilHelper.getProperty(annot,"related_annot");
        if( ob != null && ob instanceof List)
        {
            AGAVERelatedAnnot[] set = new AGAVERelatedAnnot[1];
            return (AGAVERelatedAnnot[])((List)ob).toArray( set ) ;
        }
        return (AGAVERelatedAnnot[]) null ;
    }

    public String getGroupOrder(Annotation annot)
    {
         if( annot == null ) return (String)null ;
         return (String) UtilHelper.getProperty(annot,"group_order");
    }

    public String getFeatureType(Annotation annot)
    {
         if( annot == null ) return (String)null ;
         return (String) UtilHelper.getProperty(annot,"feature_type");
    }
    public String getResultType(Annotation annot)
    {
         if( annot == null ) return (String)null ;
         return (String) UtilHelper.getProperty(annot,"result_type");
    }

    public String getConfidence(Annotation annot)
    {
         if( annot == null ) return (String)null ;
         return (String) UtilHelper.getProperty(annot,"confidence");
    }

    public String getMatchAlign(Annotation annot)
    {
         if( annot == null ) return (String)null ;
         return (String) UtilHelper.getProperty(annot,"match_align");
    }

    public AGAVEMatchRegion getMatchRegion(Annotation annot)
    {
        if( annot == null ) return (AGAVEMatchRegion)null ;
        return (AGAVEMatchRegion) UtilHelper.getProperty(annot,"match_region");
    }

    public AGAVEQueryRegion getQueryRegion(Annotation annot)
    {
        if( annot == null ) return (AGAVEQueryRegion)null ;
       return (AGAVEQueryRegion) UtilHelper.getProperty(annot,"query_region") ;
    }

    public String getAlignUnits(Annotation annot)
    {
        if( annot == null ) return (String)null ;
        return (String) UtilHelper.getProperty(annot,"align_units");
    }

    public String getMatchDesc(Annotation annot)
    {
        if( annot == null ) return (String)null ;
        return (String) UtilHelper.getProperty(annot,"match_desc");
    }

    public String[] getElementIds(Annotation annot)
    {
        if( annot == null ) return (String[])null ;
        Object ob = UtilHelper.getProperty(annot,"element_ids");
        if( ob != null && ob instanceof List)
        {
            String[] set = new String[1];
            return (String[])((List)ob).toArray( set ) ;
        }
        return (String[]) null ;
    }

    public AGAVEMapLocation[]  getMapLocation(Annotation annot)
    {
        if( annot == null ) return (AGAVEMapLocation[])null ;
        Object ob = UtilHelper.getProperty(annot,"map_location");
        if( ob != null && ob instanceof List)
        {
            AGAVEMapLocation[] set = new AGAVEMapLocation[1];
            return (AGAVEMapLocation[])((List)ob).toArray( set ) ;
        }
        return (AGAVEMapLocation[]) null ;
    }

    public AGAVEDbId[] getAltIds(Annotation annot)
    {
        if( annot == null ) return (AGAVEDbId[])null ;
        Object ob = UtilHelper.getProperty(annot,"alt_ids");
        if( ob != null && ob instanceof List)
        {
            AGAVEDbId[] set = new AGAVEDbId[1];
            return (AGAVEDbId[])((List)ob).toArray( set ) ;
        }
        return (AGAVEDbId[]) null ;
    }


    public String getClassifySystem(Annotation annot)
    {
        return (String) UtilHelper.getProperty(annot,"system");
    }

    public String getClassifyId(Annotation annot)
    {
         return (String) UtilHelper.getProperty(annot,"id");
    }

    public String getClassifyType(Annotation annot)
    {
         return (String) UtilHelper.getProperty(annot,"type");
    }

    public AGAVEDbId getDbId(Annotation annot)
    {
        return (AGAVEDbId) UtilHelper.getProperty(annot,"db_id");
    }

    public AGAVEIdAlias[] getIdAlias(Annotation annot)
    {
        Object ob = UtilHelper.getProperty(annot,"id_alias");
        if( ob != null && ob instanceof List)
        {
            AGAVEIdAlias[] set = new AGAVEIdAlias[1];
            return (AGAVEIdAlias[])((List)ob).toArray( set ) ;
        }
        return (AGAVEIdAlias[]) null ;
    }

    public String[] getExonIds(Annotation annot)
    {
        Object ob = UtilHelper.getProperty(annot,"exons");
        if( ob != null && ob instanceof List)
        {
            String[] set = new String[1];
            return (String[])((List)ob).toArray( set ) ;
        }
        return (String[]) null ;
    }

    public String getChromNum(Annotation annot)
    {
         return (String) UtilHelper.getProperty(annot,"chromosome_number");
    }

    /**        */
    public AGAVEProperty[] getProperty(Annotation annot, String type)
    {
        for (Iterator i = annot.keys().iterator(); i.hasNext();)
        {
            String key = (String) i.next();
            Object ob = UtilHelper.getProperty(annot,key)  ;
            if( ob instanceof List && ((List)ob).get(0) instanceof AGAVEProperty)
            {
                AGAVEProperty[] tmp = new AGAVEProperty[ 1 ] ;
                return (AGAVEProperty[])((List)ob).toArray( tmp ) ;
            }
         }
         return (AGAVEProperty[]) null ;
     }
}