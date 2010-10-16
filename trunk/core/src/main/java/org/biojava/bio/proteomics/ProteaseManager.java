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
 */


package org.biojava.bio.proteomics;


import java.io.InputStream;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.MissingResourceException;
import java.util.Set;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ClassTools;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;

/**
 * Registry and utility methods for Proteases.
 * @author Mark Schreiber
 */
public final class ProteaseManager {
  private static Map name2Protease = new HashMap();

  static Document doc = null;

  static {
    try {
      InputStream tablesStream =
         ClassTools.getClassLoader(ProteaseManager.class).getResourceAsStream(
            "org/biojava/bio/proteomics/ProteaseManager.xml"
          );

      if(tablesStream == null ) {
        throw new BioException("Couldn't locate ProteaseManager.xml.");
      }

      InputSource is = new InputSource(tablesStream);
      DocumentBuilder parser = DocumentBuilderFactory.newInstance().newDocumentBuilder();
      doc = parser.parse(is);

      NodeList children = doc.getDocumentElement().getChildNodes();
      for(int i = 0; i < children.getLength(); i++) {
        Node cnode = (Node) children.item(i);
        if(! (cnode instanceof Element)) {
          continue;
        }

        Element child = (Element) cnode;
        if(child.getNodeName().equals("protease")) {

          //Parameters
          SymbolList cleavRes = null;
          SymbolList exceptRes = null;
          boolean endo = false;
          String protName = child.getAttribute("name");
          Protease protease = null;

          NodeList proteaseNodes = child.getChildNodes();
          for(int j = 0; j < proteaseNodes.getLength(); j++){
            Node cnode2 = (Node) proteaseNodes.item(j);
            if(! (cnode2 instanceof Element)) {
              continue;
            }
            Element el = (Element) cnode2;
            String name = el.getNodeName();
            String content = el.getFirstChild().getNodeValue();
            if(name.equals("cleaveRes")) {
              cleavRes = createSymbolList(content.trim());
            }else if(name.equals("exceptRes")) {
              exceptRes = createSymbolList(content.trim());
            }else if(name.equals("endo")) {
              endo = new Boolean(content).booleanValue();
            }


            if(cleavRes == null)
              cleavRes = createSymbolList("");
            if(exceptRes == null){
              exceptRes = createSymbolList("");
            }
            protease = new Protease(cleavRes ,endo, exceptRes, protName);
          }
          registerProtease(protease);
        }
      }
    }catch (MissingResourceException mre) {
      System.err.println(mre.getMessage());
    }catch(Exception e){//err
      e.printStackTrace();
    }
  }

  /**
   * Creates and registers a new Protease. In future the Protease can be recovered
   * using the getProteaseByName() method.
   * @param cleaveRes the cleavege residues
   * @param endoProtease is it an endo protease?
   * @param notCleaveRes the exceptions to the cleavage residues
   * @param name the name of the Protease
   * @return a reference to the new Protease
   * @throws IllegalSymbolException if the cleaveRes or notCleaveRes are not
   * from the PROTEIN alphabet
   * @throws BioException if a Protease with the same name already exists.
   */
  public static synchronized Protease createProtease(
      SymbolList cleaveRes,
      boolean endoProtease,
      SymbolList notCleaveRes,
      String name) throws IllegalSymbolException, BioException{

    Protease p = new Protease(cleaveRes, endoProtease, notCleaveRes, name);
    registerProtease(p);
    return p;
  }

  public static synchronized Protease createProtease(
      SymbolList cleaveRes,
      boolean endoProtease,
      String name) throws IllegalSymbolException, BioException{

    Protease p = new Protease(cleaveRes, endoProtease, SymbolList.EMPTY_LIST, name);
    registerProtease(p);
    return p;
  }

  public static synchronized Protease createProtease(
      String cleaveRes,
      boolean endoProtease,
      String notCleaveRes,
      String name) throws BioException, IllegalSymbolException{

    return createProtease(createSymbolList(cleaveRes),
                          endoProtease,
                          createSymbolList(notCleaveRes),
                          name);
  }

  public static synchronized Protease createProtease(
      String cleaveRes,
      boolean endoProtease,
      String name) throws BioException, IllegalSymbolException{

    return createProtease(createSymbolList(cleaveRes),
                          endoProtease,
                          SymbolList.EMPTY_LIST,
                          name);
  }

  /**
   * Registers a protease and ensures its flyweight status
   * @param prot the Protease to register
   * @throws BioException if a Protease with the same name is already registered.
   */
  public static synchronized void registerProtease(Protease prot)throws BioException{
    if(registered(prot.getName()))
       throw new BioException(
           "A Protease has already been registered with the name "
           +prot.getName()
       );

    name2Protease.put(prot.getName(), prot);
  }

  /**
   * Gets a Protease instance by name.
   * @param proteaseName the name of a registered Protease (case sensistive)
   * @return a fly-weight Protease instance
   * @throws BioException if no protease is registered by that name
   */
public static Protease getProteaseByName(String proteaseName)
                             throws BioException {

    Protease protease = (Protease)name2Protease.get(proteaseName);
    if(protease == null){
      throw new BioException("No protease has been registered by that name");
    }
    return protease;
}

/**
 * @return an unmodifiable Set of all the registered Protease names (Strings).
 */
public static Set getNames(){
  return Collections.unmodifiableSet(name2Protease.keySet());
}

/**
 * @return an unmodifiable set of all the registered Protease objects.
 */
public static Set getAllProteases(){
  return Collections.unmodifiableSet(
      new HashSet(name2Protease.values())
  );
}

/**
 * Has a Protease been registered with that name?
 * @param proteaseName the query
 * @return true if one has, false otherwise
 */
public static boolean registered(String proteaseName){
  return name2Protease.containsKey(proteaseName);
}

/**
 * @return a reference to the singleton instance of the ProteaseManager
 */
public static synchronized ProteaseManager getInstance(){
  if(singletonInstance == null){
    singletonInstance = new ProteaseManager();
  }
  return singletonInstance;
}

static private SymbolList createSymbolList(String seq)

                              throws IllegalSymbolException, BioException {
    if(seq == null || seq.trim().equals("")){
      return SymbolList.EMPTY_LIST;
    }
    SymbolList sList;

    FiniteAlphabet prot

             = (FiniteAlphabet)AlphabetManager.alphabetForName("PROTEIN");



    SymbolTokenization tokenization = prot.getTokenization("token");

    sList = new SimpleSymbolList (tokenization, seq);

    return sList;

}

/**
 * @return a flywieght instance of Trypsin
 */
public static Protease getTrypsin(){
  try {
    return getProteaseByName(TRYPSIN);
  }
  catch (BioException ex) {
    throw new BioError("Cannot retreive Trypsin, AlphabetManager.xml may be corrupted", ex);
  }
}

/**
 * @return a flywieght instance of Lys-C
 */
public static Protease getLys_C(){
  try {
    return getProteaseByName(LYS_C);
  }
  catch (BioException ex) {
    throw new BioError("Cannot retreive Lys-C, AlphabetManager.xml may be corrupted", ex);
  }
}

/**
 * @return a flywieght instance of Arg-C
 */
public static Protease getArg_C(){
  try {
    return getProteaseByName(ARG_C);
  }
  catch (BioException ex) {
    throw new BioError("Cannot retreive Arg-C, AlphabetManager.xml may be corrupted",ex);
  }
}

/**
 * @return a flywieght instance of Asp-N
 */
public static Protease getAsp_N(){
  try {
    return getProteaseByName(ASP_N);
  }
  catch (BioException ex) {
    throw new BioError("Cannot retreive Asp-N, AlphabetManager.xml may be corrupted",ex);
  }
}

/**
 * @return a flywieght instance of Glu_C_bicarbonate
 */
public static Protease getGlu_C_bicarbonate(){
  try {
    return getProteaseByName(GLU_C_BICARB);
  }
  catch (BioException ex) {
    throw new BioError("Cannot retreive Glu_C_bicarbonate, AlphabetManager.xml may be corrupted", ex);
  }
}

/**
 * @return a flywieght instance of Glu_C_phosphate
 */
public static Protease getGlu_C_phosphate(){
  try {
    return getProteaseByName(GLU_C_PHOS);
  }
  catch (BioException ex) {
    throw new BioError("Cannot retreive Glu_C_phosphate, AlphabetManager.xml may be corrupted", ex);
  }
}

/**
 * @return a flywieght instance of Chymotrypsin
 */
public static Protease getChymotrypsin(){
  try {
    return getProteaseByName(CHYMOTRYP);
  }
  catch (BioException ex) {
    throw new BioError("Cannot retreive Chymotrypsin, AlphabetManager.xml may be corrupted", ex);
  }
}

/**
 * @return a flywieght instance of CNBr
 */
public static Protease getCNBr(){
  try {
    return getProteaseByName(CNBr);
  }
  catch (BioException ex) {
    throw new BioError("Cannot retreive CNBr, AlphabetManager.xml may be corrupted", ex);
  }
}

private static ProteaseManager singletonInstance;
public static final String TRYPSIN = "Trypsin";
public static final String LYS_C = "Lys-C";
public static final String ARG_C = "Arg-C";
public static final String ASP_N = "Asp-N";
public static final String GLU_C_BICARB = "Glu-C-bicarbonate";
public static final String GLU_C_PHOS = "Glu-C-phosphate";
public static final String CHYMOTRYP = "Chymotrypsin";
public static final String CNBr = "CNBr";

}