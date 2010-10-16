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
package org.biojava.bio.seq.io;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import org.biojava.bio.taxa.CircularReferenceException;
import org.biojava.bio.taxa.EbiFormat;
import org.biojava.bio.taxa.Taxon;
import org.biojava.bio.taxa.TaxonFactory;
import org.biojava.bio.taxa.TaxonParser;
import org.biojava.utils.ChangeVetoException;

/**
 * A parser that is able to generate Taxon entries for sequence
 * builder event streams.
 *
 * @author Matthew Pocock
 * @deprecated Use org.biojavax.bio.taxa framework instead
 */
public class OrganismParser
  extends
    SequenceBuilderFilter
//  implements
//    ParseErrorSource
{
  public static final String PROPERTY_ORGANISM = OrganismParser.class + ":organism";
  
  /**
   * Factory which wraps SequenceBuilders in an OrganismParser.
   *
   * @author Matthew Pocock
   */
  public static class Factory
    implements
      SequenceBuilderFactory,
      Serializable
  {
    private SequenceBuilderFactory delegateFactory;
    private String sciNameKey;
    private String commonNameKey;
    private String ncbiTaxonKey;
    private TaxonFactory taxonFactory;
    private TaxonParser taxonParser;
    
    public Factory(
      SequenceBuilderFactory delegateFactory,
      TaxonFactory taxonFactory,
      TaxonParser taxonParser,
      String sciNameKey,
      String commonNameKey,
      String ncbiTaxonKey
    ) {
      this.delegateFactory = delegateFactory;
      this.taxonFactory = taxonFactory;
      this.taxonParser = taxonParser;
      this.sciNameKey = sciNameKey;
      this.commonNameKey = commonNameKey;
      this.ncbiTaxonKey = ncbiTaxonKey;
    }
    
    public SequenceBuilder makeSequenceBuilder() {
      return new OrganismParser(
        delegateFactory.makeSequenceBuilder(),
        taxonFactory,
        taxonParser,
        sciNameKey,
        commonNameKey,
        ncbiTaxonKey
      );
    }
  }
  
  private final TaxonFactory taxonFactory;
  private final TaxonParser taxonParser;
  private final String sciNameKey;
  private final String commonNameKey;
  private final String ncbiTaxonKey;
  private String fullName;
  private String commonName;
  private String ncbiTaxon;
  
  public OrganismParser(
    SequenceBuilder delegate,
    TaxonFactory taxonFactory,
    TaxonParser taxonParser,
    String sciNameKey,
    String commonNameKey,
    String ncbiTaxonKey
  ) {
    super(delegate);
    this.taxonFactory = taxonFactory;
    this.taxonParser = taxonParser;
    this.sciNameKey = sciNameKey;
    this.commonNameKey = commonNameKey;
    this.ncbiTaxonKey = ncbiTaxonKey;
  }
  
  public void addSequenceProperty(Object sciNameKey, Object value)
    throws
      ParseException
  {
    if(this.sciNameKey.equals(sciNameKey)) {
      if(fullName == null) {
        fullName = value.toString();
      } else {
        fullName = fullName + " " + value;
      }
    } else if(this.commonNameKey.equals(sciNameKey)) {
      commonName = value.toString();
    } else if(this.ncbiTaxonKey.equals(sciNameKey)) {
      String tid = value.toString();
      int eq = tid.indexOf("=");
      if(eq >= 0) {
        tid = tid.substring(eq + 1);
      }
      int sc = tid.indexOf(";");
      if(sc >= 0) {
        tid = tid.substring(0, sc);
      }
      if(this.ncbiTaxon == null) {
        this.ncbiTaxon = tid;
      } else {
        this.ncbiTaxon = this.ncbiTaxon + tid;
      }
    } else {
      getDelegate().addSequenceProperty(sciNameKey, value);
    }
  }
  
  public void endSequence()
    throws
      ParseException
  {
    try {
      Taxon taxon = taxonParser.parse(taxonFactory, fullName);
      if(commonName != null && taxon.getCommonName() == null) {
        try {
          taxon.setCommonName(commonName);
        } catch (ChangeVetoException cve) {
          throw new ParseException(cve, "Failed to build Taxon");
        }
      }
      StringTokenizer stok = new StringTokenizer(ncbiTaxon, ",");
      if(stok.countTokens() == 1) {
        taxon.getAnnotation().setProperty(EbiFormat.PROPERTY_NCBI_TAXON, ncbiTaxon);
      } else {
        List tl = new ArrayList();
        while(stok.hasMoreTokens()) {
          tl.add(stok.nextToken());
        }
        taxon.getAnnotation().setProperty(EbiFormat.PROPERTY_NCBI_TAXON, tl);
      }
      getDelegate().addSequenceProperty(PROPERTY_ORGANISM, taxon);
    } catch (ChangeVetoException cve) {
      throw new ParseException(cve, "Could not parse organism: " + fullName);
    } catch (CircularReferenceException cre) {
      throw new ParseException(cre);
    }
  }
}


