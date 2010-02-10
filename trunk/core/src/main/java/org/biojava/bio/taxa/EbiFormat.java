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
package org.biojava.bio.taxa;


import java.util.Iterator;
import java.util.List;
import java.util.StringTokenizer;

import org.biojava.bio.Annotation;
import org.biojava.utils.ChangeVetoException;

/**
 * Encapsulate the 'EBI' species format used in Embl, Genbank and Swissprot
 * files. The Taxon objects created by this process have the following annotations: <p>
 *
 * <pre>
 * key=PROPERTY_NCBI_TAXON, value=String representation of the NCBI taxon id
 * key=PROPERTY_TAXON_NAMES, value=Map from name-class to name (see ncbi names.dmp)
 * </pre>
 *
 * @author Matthew Pocock
 * @author Len Trigg
 * @deprecated replaced by classes in {@link org.biojavax.bio.taxa org.biojavax.bio.taxa}
 */
public class EbiFormat implements TaxonParser {
  public static final String PROPERTY_NCBI_TAXON = EbiFormat.class + ":NCBI_TAXON";
  public static final String PROPERTY_TAXON_NAMES = EbiFormat.class + ":TAXON_NAMES";
  private static EbiFormat INSTANCE = new EbiFormat();
  
  public static final EbiFormat getInstance() {
    if(INSTANCE == null) {
      INSTANCE = new EbiFormat();
    }
    
    return INSTANCE;
  }
  
  public Taxon parse(TaxonFactory taxonFactory, String taxonString)
    throws
      ChangeVetoException,
      CircularReferenceException
  {
    String name = taxonString.trim();
    if(name.endsWith(".")) {
      name = name.substring(0, name.length() - 1);
    }
    
    Taxon taxon = taxonFactory.getRoot();
    StringTokenizer sTok = new StringTokenizer(name, ";");
    
    if(sTok.countTokens() == 1) {
      return taxonFactory.addChild(taxon, taxonFactory.createTaxon(name, null));
    }
    
    String tok = null;
    CLIMB_TREE:
    while(sTok.hasMoreTokens()) {
      tok = sTok.nextToken().trim();
      for(Iterator i = taxon.getChildren().iterator(); i.hasNext(); ) {
        Taxon child = (Taxon) i.next();
        if(child.getScientificName().equals(tok)) {
          taxon = child;
          continue CLIMB_TREE; // found child by name - go through loop again
        }
      }
      
      break; // couldn't find a child by than name - stop this and move on
    }
    
    for(; sTok.hasMoreTokens(); tok = sTok.nextToken().trim()) {
      taxon = taxonFactory.addChild(
        taxon,
        taxonFactory.createTaxon(tok, null)
      );
    }
    
    return taxon;
  }
  
  public String serialize(Taxon taxon) {
    String name = null;
    
    do {
      String sci = taxon.getScientificName();
      if(name == null) {
        name = sci + ".";
      } else {
        name = sci + "; " + name;
      }
      taxon = taxon.getParent();
    } while(taxon != null && taxon.getParent() != null);
    
    return name;
  }

    public String serializeSource(Taxon taxon) {
        StringBuffer sb = new StringBuffer(taxon.getScientificName());
        String common = taxon.getCommonName();
        if ((common != null) && (common.length() > 0)) {
            sb.append(" (").append(taxon.getCommonName()).append(")");
        }
        sb.append('.');
        return sb.toString();
    }

    public String serializeXRef(Taxon taxon) {
        Annotation anno = taxon.getAnnotation();
        Object t  = anno.getProperty(EbiFormat.PROPERTY_NCBI_TAXON);
        if (t instanceof List) {
            t = (String) ((List) t).get(0);
        }
        return "NCBI_TaxID=" + t + ";";
    }
}
