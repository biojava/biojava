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

 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */

package org.biojava.bio.seq;

import java.io.InputStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.MissingResourceException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.seq.impl.SimpleGappedSequence;
import org.biojava.bio.seq.impl.SimpleSequenceFactory;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SimpleSymbolPropertyTable;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolPropertyTable;
import org.biojava.utils.ClassTools;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;

/**
 * The central port-of-call for all information and functionality specific to
 * SymbolLists over the protein alphabet.
 *
 * @author Matthew Pocock
 * @author Greg Cox
 * @author Thomas Down
 * @author MarkSchreiber
 * @author Jonathan Warren
 * @author gwaldon (pyrrolysine, pKs)
 */
public class ProteinTools {
    private static final FiniteAlphabet proteinAlpha;
    private static final FiniteAlphabet proteinTAlpha;

    private static final Map tokenToSymbol = new HashMap();

    private static final Map propertyTableMap = new HashMap();

    static {
        try {
            proteinAlpha = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN");
            proteinTAlpha = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN-TERM");
            SymbolTokenization st = proteinTAlpha.getTokenization("token");
            for (Iterator i = proteinTAlpha.iterator(); i.hasNext(); ) {
              AtomicSymbol s = (AtomicSymbol)i.next();
              tokenToSymbol.put(st.tokenizeSymbol(s), s);
            }

        } catch (Exception e) {
            throw new BioError(" Could not initialize ProteinTools", e);
        }
    }


    static {

        Document doc = null;
     /*   try {
            URL proteaseManagerURL = ProteinTools.class.getClassLoader().getResource(
            "org/biojava/bio/symbol/ResidueProperties.xml"
            );
            //If I try and do this here on compile it says "An exception can't be thrown by an initializer"
            InputSource is = Resolver.createInputSource(proteaseManagerURL, true);
            doc = XmlDocument.createXmlDocument(is, true);*/

      try {
          InputStream tablesStream = ClassTools.getClassLoader(ProteinTools.class).getResourceAsStream(
            "org/biojava/bio/symbol/ResidueProperties.xml"
          );
          if(tablesStream == null ) {
            throw new BioError("Couldn't locate ResidueProperties.xml.");
          }

          InputSource is = new InputSource(tablesStream);
          DocumentBuilder parser = DocumentBuilderFactory.newInstance().newDocumentBuilder();
          doc = parser.parse(is);
        }catch (MissingResourceException mre) {
            System.err.println(mre.getMessage());
        }catch(Exception e){//err
            e.printStackTrace();
        }

        try {
            SimpleSymbolPropertyTable monoMassPropertyTable = new SimpleSymbolPropertyTable(
            getAlphabet(),
            SymbolPropertyTable.MONO_MASS
            );

            SimpleSymbolPropertyTable avgMassPropertyTable = new SimpleSymbolPropertyTable(
            getAlphabet(),
            SymbolPropertyTable.AVG_MASS
            );

            SimpleSymbolPropertyTable pK_NtermPropertyTable = new SimpleSymbolPropertyTable(
            getAlphabet(),
            SymbolPropertyTable.PK_Nterm
            );
            
            SimpleSymbolPropertyTable pKPropertyTable = new SimpleSymbolPropertyTable(
            getAlphabet(),
            SymbolPropertyTable.PK
            );

            SimpleSymbolPropertyTable pK_CtermPropertyTable = new SimpleSymbolPropertyTable(
            getAlphabet(),
            SymbolPropertyTable.PK_Cterm
            );
            
            SimpleSymbolPropertyTable HydropathicityTable = new SimpleSymbolPropertyTable(
            getAlphabet(),
            SymbolPropertyTable.HYDROPATHICITY
            );
            
            SymbolTokenization tokens = getAlphabet().getTokenization("token");

            NodeList children = doc.getDocumentElement().getChildNodes();
            for(int i = 0; i < children.getLength(); i++) {
                Node cnode = (Node) children.item(i);
                if(! (cnode instanceof Element)) {
                    continue;
                }
                Element child = (Element) cnode;
                if(child.getNodeName().equals("residue")) {
                    String token = child.getAttribute("token");
                    Symbol s = tokens.parseToken(token);

                    NodeList properyNodes = child.getChildNodes();
                    for(int j = 0; j < properyNodes.getLength(); j++) {
                        cnode = (Node) properyNodes.item(j);
                        if(! (cnode instanceof Element)) {
                            continue;
                        }
                        Element el = (Element) cnode;
                        String name = el.getAttribute("name");
                        if(name.equals(SymbolPropertyTable.MONO_MASS)) {
                            String value = el.getAttribute("value");
                            monoMassPropertyTable.setDoubleProperty(s, value);
                        } else if (name.equals(SymbolPropertyTable.AVG_MASS)) {
                            String value = el.getAttribute("value");
                            avgMassPropertyTable.setDoubleProperty(s, value);
                        } else if (name.equals(SymbolPropertyTable.PK_Nterm)) {
                            String value = el.getAttribute("value");
                            pK_NtermPropertyTable.setDoubleProperty(s, value);
                        } else if (name.equals(SymbolPropertyTable.PK)) {
                            String value = el.getAttribute("value");
                            pKPropertyTable.setDoubleProperty(s, value);
                        } else if (name.equals(SymbolPropertyTable.PK_Cterm)) {
                            String value = el.getAttribute("value");
                            pK_CtermPropertyTable.setDoubleProperty(s, value);
                        }else if (name.equals(SymbolPropertyTable.HYDROPATHICITY)) {
                            String value = el.getAttribute("value");
                            HydropathicityTable.setDoubleProperty(s, value);
                        }
                    }
                }
            }

            propertyTableMap.put(SymbolPropertyTable.MONO_MASS, (SymbolPropertyTable) monoMassPropertyTable);
            propertyTableMap.put(SymbolPropertyTable.AVG_MASS, (SymbolPropertyTable) avgMassPropertyTable);
            propertyTableMap.put(SymbolPropertyTable.PK_Nterm, (SymbolPropertyTable) pK_NtermPropertyTable);
            propertyTableMap.put(SymbolPropertyTable.PK, (SymbolPropertyTable) pKPropertyTable);
            propertyTableMap.put(SymbolPropertyTable.PK_Cterm, (SymbolPropertyTable) pK_CtermPropertyTable);
            propertyTableMap.put(SymbolPropertyTable.HYDROPATHICITY, (SymbolPropertyTable) HydropathicityTable);
        } catch (Exception e) {
            throw new BioError(" Could not initialize ProteinTools", e);
        }
    }
    
    private ProteinTools() {
    }
    
    /**
     *Gets the protein alphabet
     */
    public static final FiniteAlphabet getAlphabet() {
        return proteinAlpha;
    }

    /**
     *Gets the protein alphabet including the translation termination symbols
     */
    public static final FiniteAlphabet getTAlphabet() {
        return proteinTAlpha;
    }

    public static final SymbolPropertyTable getSymbolPropertyTable(String name)
    {
        return (SymbolPropertyTable)propertyTableMap.get(name);
    }

  /**
   * Return a new Protein <span class="type">SymbolList</span> for <span
   * class="arg">protein</span>.
   *
   * @param theProtein a <span class="type">String</span> to parse into Protein
   * @return a <span class="type">SymbolList</span> created form <span
   *         class="arg">Protein</span>
   * @throws IllegalSymbolException if  <span class="arg">dna</span> contains
   *                                any non-Amino Acid characters.
   */
  public static SymbolList createProtein(String theProtein)
          throws IllegalSymbolException
  {
    SymbolTokenization p = null;
    try {
      p = getTAlphabet().getTokenization("token");
    } catch (BioException e) {
      throw new BioError("Something has gone badly wrong with Protein", e);
    }
    return new SimpleSymbolList(p, theProtein);
  }

    /** Get a new protein as a GappedSequence */
    public static GappedSequence createGappedProteinSequence(String theProtein, String name) throws IllegalSymbolException{
        String theProtein1 = theProtein.replaceAll("-", "");
        Sequence protein = createProteinSequence(theProtein1, name);
        GappedSequence protein1 = new SimpleGappedSequence(protein);
        int pos = theProtein.indexOf('-', 0);
        while(pos!=-1){
            protein1.addGapInView(pos+1);
            pos = theProtein.indexOf('-', pos+1);
        }
        return protein1;
    }

  /**
   * Return a new PROTEIN <span class="type">Sequence</span> for
   * <span class="arg">protein</span>.
   *
   * @param protein a <span class="type">String</span> to parse into PROTEIN
   * @param name a <span class="type">String</span> to use as the name
   * @return a <span class="type">Sequence</span> created form
   *         <span class="arg">protein</span>
   * @throws IllegalSymbolException if <span class="arg">protein</span> contains
   *         any non-PROTEIN characters
   */
  public static Sequence createProteinSequence(String protein, String name)
  throws IllegalSymbolException {
    try {
      return new SimpleSequenceFactory().createSequence(
        createProtein(protein),
        "", name, new SimpleAnnotation()
      );
    } catch (BioException se) {
      throw new BioError("Something has gone badly wrong with ProteinTAlpha", se);
    }
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid Alanine
   * (A)
   */
  public static AtomicSymbol ala() {
    return (AtomicSymbol) tokenToSymbol.get("A");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Alanine
   */
  public static AtomicSymbol a() {
    return ala();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Arginine (R)
   */
  public static AtomicSymbol arg() {
    return (AtomicSymbol) tokenToSymbol.get("R");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Arginine
   */
  public static AtomicSymbol r() {
    return arg();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Asparagine (N)
   */
  public static AtomicSymbol asn() {
    return (AtomicSymbol) tokenToSymbol.get("N");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Asparagine
   */
  public static AtomicSymbol n() {
    return asn();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Aspartic Acid (D)
   */
  public static AtomicSymbol asp() {
    return (AtomicSymbol) tokenToSymbol.get("D");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Aspartic Acid
   */
  public static AtomicSymbol d() {
    return asp();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Cysteine (C)
   */
  public static AtomicSymbol cys() {
    return (AtomicSymbol) tokenToSymbol.get("C");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Cysteine
   */
  public static AtomicSymbol c() {
    return cys();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Glutamine (Q)
   */
  public static AtomicSymbol gln() {
    return (AtomicSymbol) tokenToSymbol.get("Q");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Glutamine
   */
  public static AtomicSymbol q() {
    return gln();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Glutamic Acid (E)
   */
  public static AtomicSymbol glu() {
    return (AtomicSymbol) tokenToSymbol.get("E");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Glutamic Acid
   */
  public static AtomicSymbol e() {
    return glu();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Glycine (G)
   */
  public static AtomicSymbol gly() {
    return (AtomicSymbol) tokenToSymbol.get("G");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Glycine
   */
  public static AtomicSymbol g() {
    return gly();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Histidine (H)
   */
  public static AtomicSymbol his() {
    return (AtomicSymbol) tokenToSymbol.get("H");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Histidine
   */
  public static AtomicSymbol h() {
    return his();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Isoleucine (I)
   */
  public static AtomicSymbol ile() {
    return (AtomicSymbol) tokenToSymbol.get("I");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Isoleucine
   */
  public static AtomicSymbol i() {
    return ile();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Leucine (L)
   */
  public static AtomicSymbol leu() {
    return (AtomicSymbol) tokenToSymbol.get("L");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Leucine
   */
  public static AtomicSymbol l() {
    return leu();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Lysine (K)
   */
  public static AtomicSymbol lys() {
    return (AtomicSymbol) tokenToSymbol.get("K");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Lysine
   */
  public static AtomicSymbol k() {
    return lys();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Methionine (M)
   */
  public static AtomicSymbol met() {
    return (AtomicSymbol) tokenToSymbol.get("M");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Methionine
   */
  public static AtomicSymbol m() {
    return met();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Phenylalanine (F)
   */
  public static AtomicSymbol phe() {
    return (AtomicSymbol) tokenToSymbol.get("F");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Phenylalanine
   */
  public static AtomicSymbol f() {
    return phe();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Proline (P)
   */
  public static AtomicSymbol pro() {
    return (AtomicSymbol) tokenToSymbol.get("P");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Proline
   */
  public static AtomicSymbol p() {
    return pro();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Pyrrolysine (O)
   */
  public static AtomicSymbol pyl() {
    return (AtomicSymbol) tokenToSymbol.get("O");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Pyrrolysine
   */
  public static AtomicSymbol o() {
    return pyl();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Selenocysteine (U)
   */
  public static AtomicSymbol sec() {
    return (AtomicSymbol) tokenToSymbol.get("U");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Selenocysteine
   */
   public static AtomicSymbol u(){
     return sec();
   }
   
  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Serine (S)
   */
  public static AtomicSymbol ser() {
    return (AtomicSymbol) tokenToSymbol.get("S");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Serine
   */
  public static AtomicSymbol s() {
    return ser();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Threonine (T)
   */
  public static AtomicSymbol thr() {
    return (AtomicSymbol) tokenToSymbol.get("T");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Threonine
   */
  public static AtomicSymbol t() {
    return thr();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Tryptophan (W)
   */
  public static AtomicSymbol trp() {
    return (AtomicSymbol) tokenToSymbol.get("W");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Tryptophan
   */
  public static AtomicSymbol w() {
    return trp();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Tyrosine (Y)
   */
  public static AtomicSymbol tyr() {
    return (AtomicSymbol) tokenToSymbol.get("Y");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Tyrosine
   */
  public static AtomicSymbol y() {
    return tyr();
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid Valine (V)
   */
  public static AtomicSymbol val() {
    return (AtomicSymbol) tokenToSymbol.get("V");
  }

  /**
   * Returns the <code>AtomicSymbol</code> for the amino acid
   * Valine
   */
  public static AtomicSymbol v() {
    return val();
  }


   /**
    * Returns the <code>AtomicSymbol</code> for the termination (*)
    * placeholder
    */
   public static AtomicSymbol ter() {
     return (AtomicSymbol) tokenToSymbol.get("*");
   }

}
