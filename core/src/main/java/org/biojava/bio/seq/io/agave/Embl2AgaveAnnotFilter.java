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
import org.biojava.bio.Annotation;

/**
 * Map EMBL data into AGAVE format
 *
 * @author Hanning Ni    Doubletwist Inc
 */
public class Embl2AgaveAnnotFilter extends SimpleAnnotFilter{

  public static final AGAVEAnnotFilterFactory EMBL_AGAVE_ANNOT_FILTER_FACTORY
    = new AGAVEAnnotFilterFactory() {
    public AGAVEAnnotFilter getInstance() {
      return new Embl2AgaveAnnotFilter();
    }
  };

  Embl2AgaveAnnotFilter() {
  }

    /**
     *
     */
   public String getAccession(Annotation annot)
   {
       return (String) UtilHelper.getProperty(annot,"AC");
   }

        /**
         *
         */
    public String getKeyword(Annotation annot)
    {
        return (String) UtilHelper.getProperty(annot,"KW");
    }

        /**
         *
         */
    public String getOrganism(Annotation annot)
    {
        return (String) UtilHelper.getProperty(annot,"OC");
    }

        /**
         *
         */
    public String getDescription(Annotation annot)
    {
        return (String) UtilHelper.getProperty(annot,"DE");
    }

        /**
         *
         */
    public String getNote(Annotation annot)
    {
        return
              "RN:" +(String) UtilHelper.getProperty(annot,"RN") + "  \n"  +
              "RP:" + (String) UtilHelper.getProperty(annot,"RP") + "  \n"  +
              "RA:" + (String) UtilHelper.getProperty(annot,"RA") + "  \n"  +
              "RT:" + (String) UtilHelper.getProperty(annot,"RT") + "  \n"  +
              "RL:" + (String) UtilHelper.getProperty(annot,"RL") + "  \n";
    }

        /**
         *
         */
     public String getVersion(Annotation annot)
     {
         return (String) UtilHelper.getProperty(annot,"SV");
     }

        /**
         *
         */
     public String getOS(Annotation annot)
     {
         return (String) UtilHelper.getProperty(annot,"OS");
     }

        /**
         *
         */
     public String getMolType(Annotation annot)
     {
            String id =  (String) UtilHelper.getProperty(annot,"ID");
            if (id.indexOf("rRNA") != -1)
            {
                return "rRNA";
            }
            else if (id.indexOf("tRNA") != -1)
            {
                return "tRNA";
            }
            else if (id.indexOf("mRNA") != -1)
            {
                return "mRNA";
            }
            else if (id.indexOf("RNA") != -1)
            {
                return "RNA";
            }
            else if (id.indexOf("cDNA") != -1)
            {
                return "cDNA";
            }
            else if (id.indexOf("DNA") != -1)
            {
                return "DNA";
            }
            else if (id.indexOf("PROTEIN") != -1)
            {
                return "AA";
            }
            return null;
     }

     public AGAVEDbId getDbId(Annotation annot)
     {
           String AC = (String)UtilHelper.getProperty(annot,"AC")  ;
           String SV = (String)UtilHelper.getProperty(annot,"SV") ;

           if( AC != null && SV != null )
           {
               AGAVEDbId dbid = new AGAVEDbId() ;
               dbid.setId( AC ) ;
               dbid.setDbCode( SV ) ;
               return dbid ;
           }
           else
              return null ;
     }

}